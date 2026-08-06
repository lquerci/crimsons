from __future__ import annotations

import warnings
from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from ..chemistry import POPIII_THRESHOLD


class MetallicityOutOfRangeWarning(UserWarning):
    """A queried metallicity fell within a population regime (Population
    III or Population II/I) but outside that regime's own tabulated
    range -- the nearest available metallicity bin was used instead."""


def _regime_lookup(masses, metallicities, yields):
    """Build a callable(mass, metallicity) -> (n_mass, n_elements) array
    for ONE population regime's (mass, metallicity) grid.

    Mass is always silently clamped to this regime's range. Metallicity
    is also clamped, but *warns* when it has to -- unlike mass, landing
    outside a regime's tabulated Z range despite being in the right
    regime (PopIII vs PopII/I) usually means the model's coverage is
    sparser than you assumed, which is worth knowing about, whereas
    mass being slightly outside a table's range is routine. Falls back
    to 1D mass-only interpolation if this regime has only one
    metallicity point (e.g. a Population III table computed at Z=0
    only).
    """
    if len(metallicities) == 1:
        z_only = metallicities[0]
        mass_only_yields = yields[:, 0, :]

        def call(mass, metallicity):
            mass = np.atleast_1d(np.asarray(mass, dtype=float))
            mass_c = np.clip(mass, masses.min(), masses.max())
            if not np.isclose(metallicity, z_only):
                warnings.warn(
                    f"metallicity {metallicity:.3g} doesn't match this regime's "
                    f"only tabulated value ({z_only:.3g}) -- using it anyway",
                    MetallicityOutOfRangeWarning,
                    stacklevel=4,
                )
            return np.stack(
                [
                    np.interp(mass_c, masses, mass_only_yields[:, ei])
                    for ei in range(mass_only_yields.shape[-1])
                ],
                axis=-1,
            )

        return call

    interp = RegularGridInterpolator(
        (masses, metallicities), yields, bounds_error=False, fill_value=None
    )
    z_min, z_max = metallicities.min(), metallicities.max()

    def call(mass, metallicity):
        mass = np.atleast_1d(np.asarray(mass, dtype=float))
        mass_c = np.clip(mass, masses.min(), masses.max())
        if metallicity < z_min or metallicity > z_max:
            warnings.warn(
                f"metallicity {metallicity:.3g} is outside this regime's tabulated "
                f"range [{z_min:.3g}, {z_max:.3g}] -- clamping to the nearest bin",
                MetallicityOutOfRangeWarning,
                stacklevel=4,
            )
        z_c = np.clip(np.full(mass.shape, metallicity, dtype=float), z_min, z_max)
        return interp(np.stack([mass_c, z_c], axis=-1))

    return call


@dataclass
class YieldTable:
    """A mass-metallicity grid of element yields for one enrichment channel.

    `yields` has shape (n_mass, n_metallicity, n_elements), in solar
    masses of each element ejected by a star of that initial mass and
    metallicity.

    Population III and Population II/I are treated as genuinely
    different kinds of stars, never interpolated across: the
    metallicity grid is split at `popiii_threshold` into two independent
    regimes, and a query is routed to whichever regime it falls in.
    Querying a regime this table has no data for raises a ValueError
    (better than silently extrapolating one population's yields into the
    other); querying within a regime but outside that regime's own
    tabulated range warns and clamps to the nearest bin (see
    MetallicityOutOfRangeWarning).
    """

    masses: np.ndarray
    metallicities: np.ndarray
    elements: list
    yields: np.ndarray
    popiii_threshold: float = POPIII_THRESHOLD

    def __post_init__(self):
        self.masses = np.asarray(self.masses, dtype=float)
        self.metallicities = np.asarray(self.metallicities, dtype=float)
        self.yields = np.asarray(self.yields, dtype=float)
        self.elements = list(self.elements)

        expected_shape = (len(self.masses), len(self.metallicities), len(self.elements))
        if self.yields.shape != expected_shape:
            raise ValueError(
                f"yields shape {self.yields.shape} doesn't match "
                f"(n_mass={len(self.masses)}, n_metallicity={len(self.metallicities)}, "
                f"n_elements={len(self.elements)}) = {expected_shape}"
            )

        is_popiii = self.metallicities < self.popiii_threshold
        self._regimes = {}
        if np.any(is_popiii):
            self._regimes["Population III"] = _regime_lookup(
                self.masses, self.metallicities[is_popiii], self.yields[:, is_popiii, :]
            )
        if np.any(~is_popiii):
            self._regimes["Population II/I"] = _regime_lookup(
                self.masses, self.metallicities[~is_popiii], self.yields[:, ~is_popiii, :]
            )

    def __call__(self, mass, metallicity, rng=None):
        # rng is accepted (and ignored) only so this has the same call
        # signature as StochasticYieldTable -- Channel.yields doesn't need
        # to know which kind of table it's calling
        regime_name = "Population III" if metallicity < self.popiii_threshold else "Population II/I"
        regime = self._regimes.get(regime_name)
        if regime is None:
            raise ValueError(
                f"no yield data for {regime_name} (Z={metallicity:.3g}, "
                f"threshold={self.popiii_threshold:.3g}) in this table -- "
                f"available: {sorted(self._regimes) or 'none'}"
            )
        return regime(mass, metallicity)


def _is_distribution(value) -> bool:
    """True if `value` should be treated as a per-star distribution to
    sample from, rather than a fixed value: anything with an `.rvs`
    method (scipy.stats-style frozen distributions) or anything plain
    callable (f(rng, n) -> array of n draws)."""
    return value is not None and (hasattr(value, "rvs") or callable(value))


def _sample_distribution(value, rng, n: int) -> np.ndarray:
    """Draw `n` samples from `value` -- a scipy.stats-like frozen
    distribution, a plain callable(rng, n) -> array, or (as a fallback)
    a fixed scalar broadcast to n identical draws."""
    if hasattr(value, "rvs"):
        return np.atleast_1d(np.asarray(value.rvs(size=n, random_state=rng), dtype=float))
    if callable(value):
        return np.atleast_1d(np.asarray(value(rng, n), dtype=float))
    return np.full(n, float(value))


class StochasticYieldTable:
    """A yield table with one or more "extra" model-parameter axes (e.g.
    SN II rotation velocity, explosion energy) resolved *per star* from a
    supplied distribution, rather than fixed once for the whole table.

    Holds one population-regime-aware YieldTable per combination of
    nearest-grid-index along each stochastic axis (built once, up front
    -- cheap, since these axes have only a handful of grid values). At
    call time, draws a fresh value per input mass from each axis's
    distribution, snaps it to that axis's nearest tabulated value (these
    are distinct model calculations, not a smoothly interpolable physical
    dimension -- consistent with how a *fixed* extra parameter is already
    resolved by nearest-value, not interpolation), and looks up the
    matching table for each star, grouped by shared index for efficiency.

    Not constructed directly -- see yields.io.load_yield_table_hdf5,
    which builds one of these when any model_params entry is a
    distribution instead of a plain value.
    """

    def __init__(self, tables_by_index: dict, axis_names: list, axis_grids: dict, axis_distributions: dict):
        self._tables_by_index = tables_by_index
        self._axis_names = list(axis_names)
        self._axis_grids = axis_grids
        self._axis_distributions = axis_distributions

    def __call__(self, mass, metallicity, rng=None):
        if rng is None:
            raise ValueError(
                "this yield table has stochastic extra parameters and needs "
                "an rng to draw them -- call it via Channel.yields(mass, "
                "metallicity, rng), not directly"
            )
        mass = np.atleast_1d(np.asarray(mass, dtype=float))
        n = len(mass)

        if not self._axis_names:
            nearest_per_star = [()] * n
        else:
            nearest_idx_per_axis = []
            for name in self._axis_names:
                grid = self._axis_grids[name]
                draws = _sample_distribution(self._axis_distributions[name], rng, n)
                nearest = np.array([int(np.argmin(np.abs(grid - d))) for d in draws])
                nearest_idx_per_axis.append(nearest)
            nearest_per_star = list(zip(*nearest_idx_per_axis))

        n_elements = len(next(iter(self._tables_by_index.values())).elements)
        result = np.empty((n, n_elements))
        nearest_per_star = np.array(nearest_per_star)
        for key in {tuple(row) for row in nearest_per_star}:
            mask = np.all(nearest_per_star == np.asarray(key), axis=1) if key else np.ones(n, dtype=bool)
            table = self._tables_by_index[key]
            result[mask] = table(mass[mask], metallicity)
        return result

    @property
    def elements(self):
        return next(iter(self._tables_by_index.values())).elements


class Channel(ABC):
    """Base class for an enrichment channel (SNII, SNIa, AGB, PISN, or a
    custom one). A channel decides which stars go through it, when their
    enrichment reaches the ISM, and how much of each element they release.
    """

    name: str

    @abstractmethod
    def contributes(self, mass, metallicity, rng) -> np.ndarray:
        """Boolean mask over `mass`: which stars go through this channel."""

    @abstractmethod
    def delay_time(self, mass, metallicity, lifetime, rng) -> np.ndarray:
        """Time (Gyr) after birth at which this channel's enrichment
        reaches the ISM, for the stars selected by `contributes`. For a
        mass-triggered channel this is usually just the stellar lifetime;
        a delay-time-distribution channel (like SN Ia) can differ.
        """

    @abstractmethod
    def yield_table(self):
        """A YieldTable or StochasticYieldTable."""

    def yields(self, mass, metallicity, rng=None):
        return self.yield_table()(mass, metallicity, rng=rng)


class MassRangeChannel(Channel):
    """A channel triggered purely by initial stellar mass (e.g. SNII, AGB,
    PISN), where enrichment happens at the end of the star's life."""

    def __init__(self, name: str, mass_min: float, mass_max: float, yield_table):
        self.name = name
        self.mass_min = mass_min
        self.mass_max = mass_max
        self._yield_table = yield_table

    def contributes(self, mass, metallicity, rng):
        return (mass >= self.mass_min) & (mass <= self.mass_max)

    def delay_time(self, mass, metallicity, lifetime, rng):
        return lifetime.copy()

    def yield_table(self):
        return self._yield_table


class PopulationChannel(Channel):
    """A channel whose events come from the stellar population as a
    whole (e.g. integrated over a delay-time distribution, or a fixed
    fraction of the total mass formed) rather than from selecting
    individual stars by mass -- currently just SN Ia.

    This doesn't fit MassRangeChannel's "which stars, at their own
    lifetime" pattern at all: a DTD's whole point is that explosions are
    spread over cosmic time independent of any single star's lifetime,
    and the number of explosions in a time bin is a population-level
    statistic, not a per-star yes/no. So instead of contributes/
    delay_time, the engine calls `population_events` directly for any
    channel of this kind (dispatched via isinstance, in
    enrichment/engine.py), and skips the per-star mass-selection path
    entirely -- accordingly, PopulationChannel instances never appear in
    a realization's `fates` array (that array describes what happened to
    each mass bin; population events aren't attributed to one).

    contributes/delay_time are given harmless/defensive default
    implementations below (never called by the engine for this channel
    kind) purely so this remains a valid Channel.
    """

    def contributes(self, mass, metallicity, rng):
        return np.zeros(np.asarray(mass).shape, dtype=bool)

    def delay_time(self, mass, metallicity, lifetime, rng):
        raise NotImplementedError(
            f"{type(self).__name__} is a PopulationChannel -- events come from "
            "population_events(), not the per-star contributes/delay_time path"
        )

    @abstractmethod
    def population_events(self, mass_formed, metallicity, lifetime_fn, time_grid, rng):
        """(times, yields) for this realization, or None if nothing
        happened -- `times` a subset of `time_grid`'s points, `yields`
        shape (len(times), n_elements), already scaled by however many
        discrete events occurred at each time."""