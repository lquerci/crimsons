from __future__ import annotations

import re
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .chemistry import default_solar_abundances
from .config import RunConfig
from .io.hdf5 import load_result, save_result

_RATIO_KEY_RE = re.compile(r"^\[\s*([A-Za-z]+)\s*/\s*([A-Za-z]+)\s*\]$")


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
    `ElementSeries` scoped to just that quantity -- see
    `__getitem__`.
    """

    config: RunConfig
    time: np.ndarray
    elements: list
    enrichment: np.ndarray  # (n_realizations, n_time, n_elements)
    fates: list  # list of (n_bins_i,) object arrays, one per realization
    masses: list  # list of (n_bins_i,) float arrays -- mean mass per populated bin
    counts: list  # list of (n_bins_i,) float arrays -- stars represented by each bin
    events_history: list # raw (name, times, yields) tuples per realization

    def mean(self) -> np.ndarray:
        """Mean enrichment across realizations, shape (n_time, n_elements)."""
        return self.enrichment.mean(axis=0)

    def std(self) -> np.ndarray:
        """Std dev across realizations, shape (n_time, n_elements)."""
        return self.enrichment.std(axis=0)

    def __getitem__(self, key: str) -> ElementSeries:
        """Look up a single element or ``[X/Y]`` abundance ratio.

        Parameters
        ----------
        key : str
            Either an element symbol present in `elements` (e.g.
            ``"Fe"``), or a bracket abundance ratio of two such symbols
            (e.g. ``"[C/Fe]"``).

        Returns
        -------
        ElementSeries
            ``result[key].mean()`` and ``.std()`` give ``(n_time,)``
            arrays; ``result[key].values`` gives the raw
            ``(n_realizations, n_time)`` array.

        Raises
        ------
        KeyError
            If an element symbol isn't in `elements`.
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

        Examples
        --------
        >>> mean_fe = result['Fe'].mean()          # doctest: +SKIP
        >>> mean_c_fe = result['[C/Fe]'].mean()     # doctest: +SKIP
        """
        if not isinstance(key, str):
            raise TypeError(
                "EnrichmentResult indices must be an element symbol "
                f"(e.g. 'Fe') or a bracket ratio (e.g. '[C/Fe]'), not {key!r}"
            )
        match = _RATIO_KEY_RE.match(key)
        if match:
            return self._abundance_ratio(match.group(1), match.group(2), label=key)
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

    def _abundance_ratio(self, numerator: str, denominator: str, label: str) -> ElementSeries:
        i_num = self._element_index(numerator)
        i_den = self._element_index(denominator)
        
        m_num = self.enrichment[:, :, i_num]
        m_den = self.enrichment[:, :, i_den]

        solar_ratio = default_solar_abundances().ratio(numerator, denominator)

        # add floors
        num_clamped = np.maximum(1.0e-20, m_num)
        den_clamped = np.maximum(1.0e-10, m_den)

        # Compute abundance ratio
        log_ratio = np.log10(num_clamped / den_clamped) - np.log10(solar_ratio)

        #with np.errstate(divide="ignore", invalid="ignore"):
        #    #ejecta_ratio = np.where((m_num > 0) & (m_den > 0), m_num / m_den, np.nan)
        #    ejecta_ratio = np.where((m_num > 0) & (m_den > 0), m_num / m_den, 1e-15)
        #    log_ratio = np.log10(ejecta_ratio)

        # Compute [X/Y] ratio with log terms
        values = np.where(m_num > 1.0e-20, log_ratio, -15.0)
        #values = log_ratio - np.log10(solar_ratio)

        # Apply -15.0 floor when element is below threshold
        #values = np.where(m_num > 1.0e-20, log_ratio, -15.0)

        # Handle Fe == Fe case explicitly (jj == 26 in Fortran)
        #if numerator == denominator:
        #    values = np.zeros_like(m_num)
#
        return ElementSeries(values, self.time, label=label)


        with np.errstate(divide="ignore", invalid="ignore"):
            ejecta_ratio = np.where((m_num > 0) & (m_den > 0), m_num / m_den, np.nan)
            log_ratio = np.log10(ejecta_ratio)

        solar_ratio = default_solar_abundances().ratio(numerator, denominator)
        values = log_ratio - np.log10(solar_ratio)
        return ElementSeries(values, self.time, label=label)

    def save(self, path):
        save_result(self, Path(path))

    @classmethod
    def load(cls, path) -> EnrichmentResult:
        return load_result(cls, Path(path))