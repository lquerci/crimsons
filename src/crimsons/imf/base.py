from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np

from .mass_range import STANDARD_MASS_RANGE, default_mass_range


class IMF(ABC):
    """Base class for stellar initial mass functions.

    Subclasses implement `_build` (construct whatever internal state
    depends on m_min/m_max -- normalization, grids, etc.), `pdf`, and
    `inverse_cdf`. `sample` (below) uses `inverse_cdf` to draw stars
    until their summed mass reaches a target.

    Mass range resolution
    ----------------------
    Pass m_min/m_max explicitly to fix the range regardless of
    metallicity -- this is locked in permanently, even if the IMF is
    later handed to a Simulation with a different metallicity.

    Otherwise, the range defaults from metallicity via `mass_range_fn`
    (see imf.mass_range.default_mass_range). This default
    is *implicit* and gets re-resolved whenever a new metallicity becomes
    known -- in particular, a Simulation automatically calls
    `bind_metallicity(self.metallicity)` before running, so you only need
    to specify metallicity once, on the Simulation, not redundantly on
    the IMF too.
    """

    def __init__(self, m_min=None, m_max=None, metallicity=None, mass_range_fn=None):
        self._mass_range_fn = mass_range_fn or default_mass_range
        self._bounds_are_explicit = m_min is not None or m_max is not None

        if self._bounds_are_explicit:
            std_min, std_max = STANDARD_MASS_RANGE
            self.m_min = m_min if m_min is not None else std_min
            self.m_max = m_max if m_max is not None else std_max
        else:
            self.m_min, self.m_max = self._mass_range_fn(metallicity)

        self._build()

    def bind_metallicity(self, metallicity):
        """Re-resolve the *implicit* default mass range for `metallicity`.

        No-op if m_min/m_max were given explicitly at construction --
        explicit bounds always win. Called automatically by Simulation;
        call it yourself if using an IMF standalone and you want its
        metallicity-adaptive default to reflect a metallicity you didn't
        know at construction time.
        """
        if self._bounds_are_explicit:
            return
        self.m_min, self.m_max = self._mass_range_fn(metallicity)
        self._build()

    @abstractmethod
    def _build(self):
        """Construct whatever internal state depends on m_min/m_max
        (normalization constants, sampling grids, ...). Called once at
        construction and again any time bind_metallicity resolves a new
        range."""

    @abstractmethod
    def pdf(self, mass: np.ndarray) -> np.ndarray:
        """dN/dM at the given masses (unnormalized is fine)."""

    @abstractmethod
    def inverse_cdf(self, u: np.ndarray) -> np.ndarray:
        """Map uniform(0, 1) samples to stellar masses. Must be vectorized."""

    def sample(
        self,
        rng: np.random.Generator,
        target_mass: float,
        n_bins : float = 500, 
        chunk_size: int = 1_000_000,
    ) -> np.ndarray:
        """Sample stellar masses until they sum to about `target_mass`.

        Draws in chunks (vectorized per chunk, not star-by-star) and stops
        as soon as the cumulative mass reaches the target, trimming the
        overshoot. The final star may push the total slightly over
        target_mass -- that's standard for this kind of stochastic
        sampling, since you can't include a fractional star.
        """

        edges = np.geomspace(self.m_min, self.m_max, n_bins + 1)
        mass_bin_centers = np.sqrt(edges[:-1] * edges[1:]) # geometric center
        count_sum = np.zeros(n_bins)
        mass_sum = np.zeros(n_bins)
        total = 0.0
 
        while total < target_mass:
            u = rng.random(chunk_size)
            chunk = self.inverse_cdf(u)
            chunk_total = float(chunk.sum())
 
            if total + chunk_total >= target_mass:
                # boundary chunk: trim star-by-star to the exact star that
                # crosses target_mass, same logic as `sample`'s global
                # trim, just applied locally to this chunk
                cum = np.cumsum(chunk)
                cutoff = int(np.searchsorted(cum, target_mass - total)) + 1
                chunk = chunk[:cutoff]
 
            counts, _ = np.histogram(chunk, bins=edges)
            sums, _ = np.histogram(chunk, bins=edges, weights=chunk)
            count_sum += counts
            mass_sum += sums
            total += float(chunk.sum())

        # filtering results in uneven arrays between relization, removed
        #populated = count_sum > 0  
        #bin_counts = count_sum[populated]
        #bin_masses = mass_sum[populated] / bin_counts

        return mass_bin_centers , count_sum