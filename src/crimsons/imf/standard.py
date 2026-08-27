from __future__ import annotations

import numpy as np

from .base import IMF
from .functional import FunctionalIMF


class BrokenPowerLawIMF(IMF):
    """dN/dM proportional to M^-alpha_i on each segment, continuous across
    breakpoints. `slopes` gives alpha_i for each segment; `interior_breakpoints`
    gives the break masses *between* m_min and m_max (so len(interior_breakpoints)
    must be len(slopes) - 1 -- m_min and m_max themselves are resolved by the
    base IMF class, not passed here as breakpoints).

    Both Salpeter (one segment, no interior breaks) and Kroupa (three
    segments, two interior breaks) are special cases of this -- add your
    own broken power law by subclassing with different slopes/breaks. For
    anything that isn't a power law, use FunctionalIMF instead.
    """

    def __init__(
        self,
        slopes,
        interior_breakpoints=(),
        m_min=None,
        m_max=None,
        metallicity=None,
        mass_range_fn=None,
    ):
        slopes = np.asarray(slopes, dtype=float)
        interior_breakpoints = list(interior_breakpoints)
        if len(interior_breakpoints) != len(slopes) - 1:
            raise ValueError(
                f"need exactly len(slopes)-1={len(slopes) - 1} interior "
                f"breakpoints, got {len(interior_breakpoints)}"
            )
        self.slopes = slopes
        self._interior_breakpoints = interior_breakpoints
        super().__init__(
            m_min=m_min, m_max=m_max, metallicity=metallicity, mass_range_fn=mass_range_fn
        )

    def _build(self):
        # consider only breakpoints between m_min and m_max
        valid_interior = [
            b for b in self._interior_breakpoints if self.m_min < b < self.m_max
        ]

        # sort breakpoints
        self.breakpoints = np.array(
            [self.m_min, *valid_interior, self.m_max], dtype=float
        )

        # keep slopes to the remaining active segments
        slopes = []
        full_breaks = [self.m_min, *self._interior_breakpoints, self.m_max]
        for i in range(len(self.slopes)):
            lo, hi = full_breaks[i], full_breaks[i + 1]
            if hi > self.m_min and lo < self.m_max:
                slopes.append(self.slopes[i])

        self.active_slopes = np.array(slopes, dtype=float)

        n_seg = len(self.active_slopes)
        # amplitude of each segment, chosen so dN/dM is continuous across
        # breakpoints (amplitude[0] is an arbitrary overall scale)
        amplitude = np.ones(n_seg)
        for i in range(1, n_seg):
            amplitude[i] = amplitude[i - 1] * self.breakpoints[i] ** (
                self.active_slopes[i] - self.active_slopes[i - 1]
            )

        # integral of dN/dM over each segment, used to pick which segment
        # a given uniform draw falls into
        seg_integral = np.empty(n_seg)
        for i in range(n_seg):
            a = self.active_slopes[i]
            lo, hi = self.breakpoints[i], self.breakpoints[i + 1]
            if np.isclose(a, 1.0):
                seg_integral[i] = amplitude[i] * (np.log(hi) - np.log(lo))
            else:
                n = 1.0 - a
                seg_integral[i] = amplitude[i] / n * (hi**n - lo**n)

        cum_weight = np.concatenate([[0.0], np.cumsum(seg_integral)])
        self._cum_weight = cum_weight / cum_weight[-1]
        self._amplitude = amplitude / cum_weight[-1]

    def pdf(self, mass):
        mass = np.atleast_1d(np.asarray(mass, dtype=float))
        out = np.zeros_like(mass)
        for i in range(len(self.active_slopes)):
            lo, hi = self.breakpoints[i], self.breakpoints[i + 1]
            mask = (mass >= lo) & (mass <= hi)
            out[mask] = self._amplitude[i] * mass[mask] ** (-self.active_slopes[i])
        return out

    def inverse_cdf(self, u):
        u = np.atleast_1d(np.asarray(u, dtype=float))

        seg_idx = np.searchsorted(self._cum_weight, u, side="right") - 1
        seg_idx = np.clip(seg_idx, 0, len(self.active_slopes) - 1)

        lo = self.breakpoints[seg_idx]
        hi = self.breakpoints[seg_idx + 1]
        a = self.active_slopes[seg_idx]
        u_lo = self._cum_weight[seg_idx]
        u_hi = self._cum_weight[seg_idx + 1]
        u_local = (u - u_lo) / (u_hi - u_lo)

        n = 1.0 - a
        flat = np.isclose(n, 0.0)
        n_safe = np.where(flat, 1.0, n)

        pow_result = (u_local * (hi**n_safe - lo**n_safe) + lo**n_safe) ** (
            1.0 / n_safe
        )
        log_result = lo * np.exp(u_local * np.log(hi / lo))
        return np.where(flat, log_result, pow_result)


class Salpeter1955(BrokenPowerLawIMF):
    """Single power law, dN/dM ~ M^-2.35 (Salpeter 1955)."""

    def __init__(
        self, m_min=None, m_max=None, metallicity=None, mass_range_fn=None, slope=2.35
    ):
        super().__init__(
            slopes=[slope],
            interior_breakpoints=(),
            m_min=m_min,
            m_max=m_max,
            metallicity=metallicity,
            mass_range_fn=mass_range_fn,
        )


class Kroupa2001(BrokenPowerLawIMF):
    """Three-segment broken power law (Kroupa 2001).

    Note: m_min/m_max follow the shared metallicity-adaptive default
    (imf.mass_range.default_mass_range) like every other IMF in this
    package, rather than Kroupa's own commonly-cited brown-dwarf cutoff
    of 0.01 Msun. Pass m_min=0.01 explicitly if you want that back.
    """

    def __init__(
        self,
        m_min=None,
        m_max=None,
        metallicity=None,
        mass_range_fn=None,
        m_break1=0.08,
        m_break2=0.5,
        slopes=(0.3, 1.3, 2.3),
    ):
        super().__init__(
            slopes=list(slopes),
            interior_breakpoints=(m_break1, m_break2),
            m_min=m_min,
            m_max=m_max,
            metallicity=metallicity,
            mass_range_fn=mass_range_fn,
        )

def chabrier_lognormal_shape(m, m_ch: float = 0.35, alpha: float = 1.35):
    """dN/dM ~ M^-(1+alpha) * exp(-m_ch/M) -- a Chabrier (2003)-type
    shape: turns over below the characteristic mass `m_ch`, approaches a
    M^-(1+alpha) power law well above it.

    Not the literal piecewise Chabrier (2003) form (a separate
    log-normal below 1 Msun spliced to a power law above) -- this is one
    continuous function with the same qualitative shape, meant to be
    used with `FunctionalIMF` (see `Chabrier2003`).
    """
    m = np.asarray(m, dtype=float)
    return m ** -(1.0 + alpha) * np.exp(-m_ch / m)


class Chabrier2003(FunctionalIMF):
    """FunctionalIMF built from `chabrier_lognormal_shape`.

    `m_ch` is the characteristic turnover mass -- 0.35 Msun (the usual
    solar-neighborhood value) by default. Population III star formation
    is often modeled as much more top-heavy, with a characteristic mass
    of order several to ~10 Msun; see `crimsons.imf.defaults.default_imf`,
    which raises `m_ch` to 10 automatically below
    `crimsons.chemistry.POPIII_THRESHOLD`.
    """

    def __init__(
        self,
        m_min=None,
        m_max=None,
        metallicity=None,
        mass_range_fn=None,
        m_ch: float = 0.35,
        alpha: float = 1.35,
        n_grid: int = 2000,
    ):
        self.m_ch = m_ch
        self.alpha = alpha
        super().__init__(
            chabrier_lognormal_shape,
            m_min=m_min,
            m_max=m_max,
            metallicity=metallicity,
            mass_range_fn=mass_range_fn,
            n_grid=n_grid,
            imf_params={"m_ch": m_ch, "alpha": alpha},
        )


class FlatIMF(BrokenPowerLawIMF):
    """Flat Initial Mass Function: dN/dM = constant on [m_min, m_max]."""


    def __init__(
        self, m_min=None, m_max=None, metallicity=None, mass_range_fn=None, slope=0.0
    ):
        super().__init__(
            slopes=[slope],
            interior_breakpoints=(),
            m_min=m_min,
            m_max=m_max,
            metallicity=metallicity,
            mass_range_fn=mass_range_fn,
        )