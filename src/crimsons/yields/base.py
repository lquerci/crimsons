from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np
from scipy.interpolate import RegularGridInterpolator


@dataclass
class YieldTable:
    """A mass-metallicity grid of element yields for one enrichment channel.

    `yields` has shape (n_mass, n_metallicity, n_elements), in solar
    masses of each element ejected by a star of that initial mass and
    metallicity. Calling the table interpolates (and clamps at the grid
    edges rather than extrapolating wildly) for arbitrary masses/Z.

    If only one metallicity is present (e.g. a Population III / PISN
    table computed at Z=0 only), this falls back to 1D interpolation over
    mass alone, since scipy's RegularGridInterpolator needs >=2 points
    per axis.
    """

    masses: np.ndarray
    metallicities: np.ndarray
    elements: list
    yields: np.ndarray

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

        if len(self.metallicities) == 1:
            self._interp = None
            self._mass_only_yields = self.yields[:, 0, :]
        else:
            self._interp = RegularGridInterpolator(
                (self.masses, self.metallicities),
                self.yields,
                bounds_error=False,
                fill_value=None,  # extrapolate; we clamp inputs below instead
            )

    def __call__(self, mass, metallicity):
        mass = np.atleast_1d(np.asarray(mass, dtype=float))
        mass_c = np.clip(mass, self.masses.min(), self.masses.max())

        if self._interp is None:
            return np.stack(
                [
                    np.interp(mass_c, self.masses, self._mass_only_yields[:, ei])
                    for ei in range(len(self.elements))
                ],
                axis=-1,
            )

        z = np.full(mass.shape, metallicity, dtype=float)
        z_c = np.clip(z, self.metallicities.min(), self.metallicities.max())
        points = np.stack([mass_c, z_c], axis=-1)
        return self._interp(points)


class Channel(ABC):
    """Base class for an enrichment channel (SNII, SNIa, AGB, or a custom
    one). A channel decides which stars go through it, when their
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
    def yield_table(self) -> YieldTable: ...

    def yields(self, mass, metallicity):
        return self.yield_table()(mass, metallicity)


class MassRangeChannel(Channel):
    """A channel triggered purely by initial stellar mass (e.g. SNII, AGB),
    where enrichment happens at the end of the star's life."""

    def __init__(self, name: str, mass_min: float, mass_max: float, yield_table: YieldTable):
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