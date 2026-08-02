from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np

from ..chemistry import POPIII_THRESHOLD, ZSUN


class LifetimeFunction(ABC):
    """Base class for stellar lifetime formulas.

    Subclasses implement a vectorized mapping from (mass, metallicity) to
    lifetime in Gyr. This is where you plug in whatever fit your Fortran
    code used (e.g. Padovani & Matteucci 1993, Raiteri et al. 1996).
    """

    @abstractmethod
    def __call__(self, mass: np.ndarray, metallicity: float) -> np.ndarray: ...

class StellarLifetime(LifetimeFunction):
    """
    Stellar lifetime model combining Raiteri et al. (1996) and Scharer et al. (2002).
    """

    def __init__(self):
        self.zsun = ZSUN

    def _get_raiteri_coeffs(self, abs_z: float) -> tuple[float, float, float]:
        """Calculates a_0, a_1, a_2 polynomial coefficients."""
        a = 7.0e-5
        b = 3.0e-2

        # Clamp absolute metallicity between [a, b] as implemented in Fortran
        z_eff = np.clip(abs_z, a, b)
        log_z = np.log10(z_eff)

        a_0 = 10.13 + 0.07547 * log_z - 0.008084 * (log_z**2.0)
        a_1 = -4.424 - 0.7939 * log_z - 0.1187 * (log_z**2.0)
        a_2 = 1.262 + 0.3385 * log_z + 0.05417 * (log_z**2.0)

        return a_0, a_1, a_2

    def __call__(self, mass: np.ndarray | float, metallicity: float) -> np.ndarray:
        """Computes lifetime in Gyr for an array of masses and absolute metallicity."""
        m = np.atleast_1d(np.asarray(mass, dtype=float))
        t_life = np.zeros_like(m)

        # metallicity is absolute
        abs_z = metallicity

        log_m = np.log10(np.maximum(m, 1e-10))  # Protection against non-positive masses

        if metallicity <= POPIII_THRESHOLD:
            # Pop III stars: Raiteri+1996 for m < 5, Scharer+2002 for m >= 5
            mask_low = m < 5.0
            mask_mid = (m >= 5.0) & (m <= 500.0)
            mask_high = m > 500.0

            # Low mass (< 5 Msun): Raiteri+1996
            if np.any(mask_low):
                a_0, a_1, a_2 = self._get_raiteri_coeffs(abs_z)
                t_log = a_0 + a_1 * log_m[mask_low] + a_2 * (log_m[mask_low] ** 2.0)
                t_life[mask_low] = 10.0**t_log  # in years

            # Medium mass (5 - 500 Msun): Scharer+2002
            if np.any(mask_mid):
                log_m_mid = log_m[mask_mid]
                t_log = (
                    9.785
                    - 3.759 * log_m_mid
                    + 1.413 * (log_m_mid**2.0)
                    - 0.186 * (log_m_mid**3.0)
                )
                t_life[mask_mid] = 10.0**t_log  # in years

            # High mass (> 500 Msun): Clamped to 500 Msun
            if np.any(mask_high):
                log_m_500 = np.log10(500.0)
                t_log = (
                    9.785
                    - 3.759 * log_m_500
                    + 1.413 * (log_m_500**2.0)
                    - 0.186 * (log_m_500**3.0)
                )
                t_life[mask_high] = 10.0**t_log  # in years

        else:
            # Pop II/I stars: Raiteri+1996 for all valid masses (m > 0)
            valid_mask = m > 0.0
            if np.any(valid_mask):
                a_0, a_1, a_2 = self._get_raiteri_coeffs(abs_z)
                t_log = (
                    a_0
                    + a_1 * log_m[valid_mask]
                    + a_2 * (log_m[valid_mask] ** 2.0)
                )
                t_life[valid_mask] = 10.0**t_log  # in years

        # Convert years to Myr (1 Myr = 10^6 years)
        return t_life / 1.0e6