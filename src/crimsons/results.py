from __future__ import annotations

import re
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .chemistry import ATOMIC_WEIGHTS, ZSUN, default_solar_abundances
from .config import RunConfig
from .io.hdf5 import load_result, save_result

_RATIO_KEY_RE = re.compile(r"^\[\s*([A-Za-z]+)\s*/\s*([A-Za-z]+)\s*\]$")
_NUMBER_RATIO_KEY_RE = re.compile(r"^([A-Za-z]+)\s*/\s*([A-Za-z]+)$")

# EnrichmentResult.__getitem__ key for the explosion-energy series (see
# EnrichmentResult.energy). Deliberately snake_case with an underscore,
# unlike every element symbol -- [A-Za-z]+ -based, so it can never match
# _RATIO_KEY_RE / _NUMBER_RATIO_KEY_RE and end up treated as a ratio
# numerator/denominator (e.g. '[explosion_energy/Fe]' or
# 'explosion_energy/Fe' both just fail to parse as a ratio and fall
# through to a plain "unknown element" KeyError, same as any other
# nonsense key).
EXPLOSION_ENERGY_KEY = "explosion_energy"

# Elements excluded from the total-metallicity sum ('Z', see
# EnrichmentResult._metallicity_mass).
_METALLICITY_EXCLUDED_ELEMENTS = {"H", "He"}

@dataclass
class ElementSeries:
    """A single element's (or ``[X/Y]`` abundance ratio's) enrichment
    time series across every realization of an `EnrichmentResult`.

    Returned by `EnrichmentResult.__getitem__` -- you normally get one of
    these from ``result['Fe']`` or ``result['[C/Fe]']`` rather than
    constructing it directly.

    Attributes
    ----------
    values : numpy.ndarray
        Shape ``(n_realizations, n_time)``. For a plain element, this is
        cumulative mass returned, in Msun. For a ``[X/Y]`` ratio, this is
        the dimensionless log-ratio (`nan` at times before both elements
        have any nonzero yield, since the ratio is undefined there).
    time : numpy.ndarray
        Shape ``(n_time,)``, Myr -- shared with the parent
        `EnrichmentResult`, not copied.
    label : str
        The key this was looked up with, e.g. ``"Fe"`` or ``"[C/Fe]"``.
    """

    values: np.ndarray
    time: np.ndarray
    label: str

    def mean(self) -> np.ndarray:
        """Mean across realizations, shape ``(n_time,)``.

        For a ``[X/Y]`` ratio, realizations where either element is still
        zero are excluded from the mean at that time step (via
        `numpy.nanmean`) rather than pulling it towards `nan`. Early time
        steps where *every* realization is still `nan` legitimately have
        no defined mean and come back as `nan` themselves.
        """
        if not np.isnan(self.values).any():
            return self.values.mean(axis=0)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            return np.nanmean(self.values, axis=0)

    def std(self) -> np.ndarray:
        """Std dev across realizations, shape ``(n_time,)`` (see `mean`
        for how `nan` values are handled).
        """
        if not np.isnan(self.values).any():
            return self.values.std(axis=0)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            return np.nanstd(self.values, axis=0)

    def __len__(self) -> int:
        return self.values.shape[0]

    def __repr__(self) -> str:
        n_real, n_time = self.values.shape
        return (
            f"ElementSeries({self.label!r}, "
            f"n_realizations={n_real}, n_time={n_time})"
        )


@dataclass
class EnrichmentResult:
    """Output of a Simulation run: N realizations of sampled stellar
    populations (as log-mass bins, not individual stars -- see
    IMF.sample_binned) and time-resolved cumulative element enrichment.

    Index by element symbol or bracket abundance ratio to get an
    `ElementSeries` scoped to just that quantity -- see `__getitem__`.
    Chemistry and explosion energy are kept separate: `elements` /
    `enrichment` only ever cover the 30 tracked elements (never
    'explosion_energy'), and `results['explosion_energy']` is the only
    way to reach the `energy` array through `__getitem__` -- it doesn't
    participate in `[X/Y]` or plain `X/Y` ratios (see `__getitem__`).
    """

    config: RunConfig
    time: np.ndarray
    elements: list  # chemical element symbols only -- never includes 'explosion_energy'
    enrichment: np.ndarray  # (n_realizations, n_time, n_elements) -- chemistry only
    energy: np.ndarray  # (n_realizations, n_time) -- cumulative explosion energy, same units as the yield table's energy column
    fates: list  # list of (n_bins_i,) object arrays, one per realization
    masses: list  # list of (n_bins_i,) float arrays -- mean mass per populated bin
    counts: list  # list of (n_bins_i,) float arrays -- stars represented by each bin
    events_history: list # raw (name, times, yields) tuples per realization

    def __post_init__(self):
        if self.enrichment.shape[-1] != len(self.elements):
            raise ValueError(
                f"enrichment has {self.enrichment.shape[-1]} element columns but "
                f"elements has {len(self.elements)} entries -- these must match. "
                "(If you're loading an old cached result saved before explosion "
                "energy was split out of `elements`, it will have 31 entries "
                "instead of 30 -- regenerate the cache.)"
            )
        if self.energy.shape != self.enrichment.shape[:2]:
            raise ValueError(
                f"energy has shape {self.energy.shape}, expected "
                f"{self.enrichment.shape[:2]} (n_realizations, n_time) to match enrichment"
            )

    def mean(self) -> np.ndarray:
        """Mean chemical enrichment across realizations, shape (n_time,
        n_elements) -- explosion energy isn't included; use
        `self.energy.mean(axis=0)` for that."""
        return self.enrichment.mean(axis=0)

    def std(self) -> np.ndarray:
        """Std dev of chemical enrichment across realizations, shape
        (n_time, n_elements) -- see `mean`."""
        return self.enrichment.std(axis=0)

    def __getitem__(self, key: str) -> ElementSeries:
        """Look up a single element or ``[X/Y]`` abundance ratio.

        Parameters
        ----------
        key : str
            One of: an element symbol present in `elements` (e.g.
            ``"Fe"``); a bracket abundance ratio of two such symbols
            (e.g. ``"[C/Fe]"``); ``"explosion_energy"``, the cumulative
            explosion energy series; ``"Z"``, total metallicity (summed
            ejecta mass of every tracked element except H and He); or
            ``"[Z/Zsun]"``, log10 of the ejecta's own metal mass
            fraction relative to solar (see `_metallicity_solar_ratio`
            for exactly how that's defined -- it is *not* simply
            ``log10(result['Z'].values / ZSUN)``, since `Z` there is a
            mass in Msun and ZSUN is a mass fraction; see the Notes).

        Returns
        -------
        ElementSeries
            ``result[key].mean()`` and ``.std()`` give ``(n_time,)``
            arrays; ``result[key].values`` gives the raw
            ``(n_realizations, n_time)`` array.

        Raises
        ------
        KeyError
            If an element symbol isn't in `elements`. Only ``"Z"``
            paired with ``"Zsun"`` is wired up as a ratio -- e.g.
            ``"[Z/Fe]"`` or ``"Z/Zsun"`` (no brackets) raise KeyError
            too, since ``"Z"`` isn't itself a tracked element.
        TypeError
            If `key` isn't a string.

        Notes
        -----
        A ``[X/Y]`` ratio is computed from the *cumulative ejecta mass*
        of X and Y (this library's enrichment source term -- see the
        Physics docs), not a gas-phase abundance diluted into an ISM
        reservoir, since `EnrichmentResult` doesn't track total gas mass.
        It's :math:`\\log_{10}(M_X/M_Y) - \\log_{10}(M_X/M_Y)_\\odot`,
        which is algebraically identical to the standard number-abundance
        definition :math:`\\log_{10}(N_X/N_Y) - \\log_{10}(N_X/N_Y)_\\odot`
        regardless of atomic weights, since they cancel in the
        ejecta-to-solar ratio. The solar reference comes from
        `crimsons.chemistry.default_solar_abundances` unless you've
        called `crimsons.chemistry.set_default_solar_abundances`.

        ``"[Z/Zsun]"`` follows the same "no gas reservoir tracked" logic,
        but ZSUN (`crimsons.chemistry.ZSUN` = 0.0142) is a mass
        *fraction*, not a mass -- so unlike every other ``[X/Y]``, this
        one can't just take a ratio of two ejecta masses. Instead it
        computes the ejecta's own metallicity mass fraction (metal mass
        ejected so far / *all* tracked-element mass ejected so far,
        i.e. H+He+metals) and compares that fraction to ZSUN. This means
        ``result['[Z/Zsun]']`` is *not* ``log10(result['Z'].values /
        ZSUN)`` -- ``result['Z']`` is a raw mass (Msun), consistent with
        every other single-key lookup, while the ``Z`` inside
        ``[Z/Zsun]`` is a dimensionless fraction. See
        `_metallicity_solar_ratio`.

        Examples
        --------
        >>> mean_fe = result['Fe'].mean()          # doctest: +SKIP
        >>> mean_c_fe = result['[C/Fe]'].mean()     # doctest: +SKIP
        >>> mean_energy = result['explosion_energy'].mean()  # doctest: +SKIP
        >>> mean_z = result['[Z/Zsun]'].mean()      # doctest: +SKIP
        """
        if not isinstance(key, str):
            raise TypeError(
                "EnrichmentResult indices must be an element symbol "
                f"(e.g. 'Fe') or a bracket ratio (e.g. '[C/Fe]'), not {key!r}"
            )

        if key == EXPLOSION_ENERGY_KEY:
            return self._energy_series()

        if key == "Z":
            return self._metallicity_series(label=key)

        if key == "ejected_mass":
            return self._ejected_mass_series(label=key)

        # Check for [X/Y] bracket notation (scaled to solar)
        match = _RATIO_KEY_RE.match(key)
        if match:
            numerator, denominator = match.group(1), match.group(2)
            if numerator == "Z" and denominator == "Zsun":
                return self._metallicity_solar_ratio(label=key)
            return self._abundance_ratio(numerator, denominator, label=key)

        # Check for X/Y plain notation (number of atoms ratio)
        match_number = _NUMBER_RATIO_KEY_RE.match(key)
        if match_number:
            return self._number_ratio(match_number.group(1), match_number.group(2), label=key)

        return self._element_series(key)
        

    def __contains__(self, key: str) -> bool:
        try:
            self[key]
        except (KeyError, TypeError):
            return False
        return True

    def _element_index(self, symbol: str) -> int:
        try:
            return self.elements.index(symbol)
        except ValueError:
            raise KeyError(
                f"Unknown element {symbol!r}. Available: {', '.join(self.elements)}"
            ) from None

    def _element_series(self, symbol: str) -> ElementSeries:
        idx = self._element_index(symbol)
        return ElementSeries(self.enrichment[:, :, idx], self.time, label=symbol)

    def _energy_series(self) -> ElementSeries:
        return ElementSeries(self.energy, self.time, label=EXPLOSION_ENERGY_KEY)

    def _metallicity_mass(self) -> np.ndarray:
        """Cumulative ejecta mass summed over every tracked element
        except H and He, shape (n_realizations, n_time), in Msun --
        this is `results['Z']`. Same convention as any single-element
        series: a raw mass, never NaN (0 before any metal enrichment,
        same as e.g. `results['Fe']` before the first SN)."""
        idx = [i for i, el in enumerate(self.elements) if el not in _METALLICITY_EXCLUDED_ELEMENTS]
        return self.enrichment[:, :, idx].sum(axis=2)

    def _metallicity_series(self, label: str) -> ElementSeries:
        return ElementSeries(self._metallicity_mass(), self.time, label=label)

    def _ejected_mass(self) -> np.ndarray:
        """Cumulative ejecta mass summed over ALL tracked elements, 
        shape (n_realizations, n_time), in Msun.
        """
        # Sum across the element axis (axis=2) without any filtering
        return self.enrichment.sum(axis=2)

    def _ejected_mass_series(self, label: str) -> ElementSeries:
        return ElementSeries(self._ejected_mass(), self.time, label=label)

    def _metallicity_solar_ratio(self, label: str) -> ElementSeries:
        """log10(Z/Zsun) -- see the __getitem__ Notes for why this needs
        a mass *fraction*, not the raw mass `results['Z']` returns.

        Z here is the ejecta's own metallicity mass fraction: metal mass
        ejected so far, divided by *all* tracked-element mass ejected so
        far (H+He+metals) -- i.e. "what fraction, by mass, of everything
        returned to the ISM so far is metals", not a gas-phase ISM
        abundance (EnrichmentResult doesn't track a diluting gas
        reservoir; same caveat as every `[X/Y]` ratio). NaN before any
        mass has been ejected at all (0/0), same convention as `[X/Y]`.
        """
        z_mass = self._metallicity_mass()
        total_mass = self.enrichment.sum(axis=2)  # every tracked element: H+He+metals
        with np.errstate(divide="ignore", invalid="ignore"):
            z_fraction = np.where(total_mass > 0, z_mass / total_mass, np.nan)
            values = np.where(z_fraction > 0, np.log10(z_fraction) - np.log10(ZSUN), np.nan)
        return ElementSeries(values, self.time, label=label)

    def _abundance_ratio(self, numerator: str, denominator: str, label: str) -> ElementSeries:
        i_num = self._element_index(numerator)
        i_den = self._element_index(denominator)
        
        m_num = self.enrichment[:, :, i_num]
        m_den = self.enrichment[:, :, i_den]

        solar_ratio = default_solar_abundances().ratio(numerator, denominator)

        with np.errstate(divide="ignore", invalid="ignore"):
            ejecta_ratio = np.where((m_num > 0) & (m_den > 0), m_num / m_den, np.nan)
            log_ratio = np.log10(ejecta_ratio)

        solar_ratio = default_solar_abundances().ratio(numerator, denominator)
        values = log_ratio - np.log10(solar_ratio)
        return ElementSeries(values, self.time, label=label)

    def _number_ratio(self, numerator: str, denominator: str, label: str) -> ElementSeries:
        """Computes the linear ratio of the number of atoms for two elements."""
        i_num = self._element_index(numerator)
        i_den = self._element_index(denominator)
        
        m_num = self.enrichment[:, :, i_num]
        m_den = self.enrichment[:, :, i_den]

        # Fetch atomic weights to convert mass to number of atoms
        try:
            m_atomic_num = ATOMIC_WEIGHTS[numerator]
            m_atomic_den = ATOMIC_WEIGHTS[denominator]
        except KeyError as e:
            raise KeyError(f"Atomic mass for {e.args[0]} not found in ATOMIC_MASSES dictionary.") from None

        # N = Mass / Atomic Mass
        n_num = m_num / m_atomic_num
        n_den = m_den / m_atomic_den

        with np.errstate(divide="ignore", invalid="ignore"):
            # Set to NaN where either element hasn't been produced yet to avoid div-by-zero
            ratio = np.where((n_num > 0) & (n_den > 0), np.log10(n_num / n_den), np.nan)

        return ElementSeries(ratio, self.time, label=label)
    
    def save(self, path):
        save_result(self, Path(path))

    @classmethod
    def load(cls, path) -> EnrichmentResult:
        return load_result(cls, Path(path))