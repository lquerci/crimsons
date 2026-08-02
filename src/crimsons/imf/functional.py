from __future__ import annotations

import numpy as np

from .base import IMF


class FunctionalIMF(IMF):
    """An IMF built from an arbitrary shape function.

    Give it any callable `shape(mass) -> dN/dM`, up to an overall
    normalization -- it doesn't need to integrate to 1, or to anything in
    particular, and it doesn't need to be vectorized (FunctionalIMF calls
    it elementwise as a fallback if a vectorized call doesn't work). This
    class normalizes it into a proper pdf and builds a numerical inverse
    CDF on a log-spaced mass grid, so you don't need a closed-form inverse
    to use a custom shape: log-normal-plus-power-law forms (Chabrier-
    style), tapered (Larson-style) forms, top-heavy Population III
    shapes, or anything else from a paper.

    For a plain broken power law, BrokenPowerLawIMF (and Salpeter1955 /
    Kroupa2001) are exact and slightly faster since they skip the
    numerical grid; this is the general-purpose option for everything
    else.

    Parameters
    ----------
    shape : callable
        mass -> dN/dM (unnormalized).
    n_grid : int
        Number of log-spaced points used to build the internal CDF grid.
        Higher is more accurate but a bit slower to construct (sampling
        speed itself is unaffected -- it's just table lookups either way).

    Example
    -------
    >>> # a simple two-part shape: flat below 1 Msun, Salpeter slope above
    >>> def my_shape(m):
    ...     m = np.atleast_1d(m)
    ...     return np.where(m < 1.0, 1.0, m ** -2.35)
    >>> imf = FunctionalIMF(my_shape, m_min=0.1, m_max=100.0)
    >>> masses = imf.sample(np.random.default_rng(0), target_mass=1000.0)
    """

    def __init__(
        self,
        shape,
        m_min=None,
        m_max=None,
        metallicity=None,
        mass_range_fn=None,
        n_grid: int = 2000,
        # argument to pass to the shape function
        imf_params: dict | None = None,  # <-- Added parameter dictionary
        **kwargs,
    ):
        self._shape = shape
        self._n_grid = n_grid
        self.imf_params = dict(imf_params or {})
        self.imf_params.update(kwargs)

        super().__init__(
            m_min=m_min, m_max=m_max, metallicity=metallicity, mass_range_fn=mass_range_fn
        )

    def _evaluate(self, mass):
        """Call the shape function, vectorized if possible, elementwise
        if not."""
        try:
            # Try vectorized call with params
            result = np.asarray(self._shape(mass, **self.imf_params), dtype=float)
            if result.shape == mass.shape:
                return result
        except TypeError:
            # Fallback for shapes that don't take imf_params
            try:
                result = np.asarray(self._shape(mass), dtype=float)
                if result.shape == mass.shape:
                    return result
            except Exception: # noqa: BLE001, S110 -- deliberate fallback for non-vectorized evaluation
                pass
        except Exception:  # noqa: BLE001, S110 -- shape() is arbitrary user
            # code; any failure here just means "not vectorized", so we
            # fall through to the elementwise path below.
            pass
        try:
            return np.array([float(self._shape(m, **self.imf_params)) for m in mass])
        except TypeError:
            return np.array([float(self._shape(m)) for m in mass])

    def _build(self):
        mass_grid = np.geomspace(self.m_min, self.m_max, self._n_grid)
        density = self._evaluate(mass_grid)

        if np.any(density < 0):
            raise ValueError("IMF shape function returned negative values")

        cdf_grid = np.concatenate(
            [[0.0], np.cumsum(0.5 * (density[1:] + density[:-1]) * np.diff(mass_grid))]
        )
        total = cdf_grid[-1]
        if total <= 0:
            raise ValueError(
                "IMF shape function integrates to ~zero over "
                f"[{self.m_min}, {self.m_max}] -- check it's positive somewhere in range"
            )

        self._mass_grid = mass_grid
        self._density_grid = density / total
        self._cdf_grid = cdf_grid / total

    def pdf(self, mass):
        mass = np.atleast_1d(np.asarray(mass, dtype=float))
        return np.interp(mass, self._mass_grid, self._density_grid, left=0.0, right=0.0)

    def inverse_cdf(self, u):
        u = np.clip(np.atleast_1d(np.asarray(u, dtype=float)), 0.0, 1.0)
        return np.interp(u, self._cdf_grid, self._mass_grid)