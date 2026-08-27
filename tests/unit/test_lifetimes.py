import numpy as np

from crimsons import ZSUN
from crimsons.stars.lifetimes import StellarLifetime


def test_lifetime_decreases_with_mass():
    lifetime = StellarLifetime()
    masses = np.array([1.0, 5.0, 20.0])
    tau = lifetime(masses, metallicity=0.02)
    assert np.all(np.diff(tau) < 0)


def test_lifetime_is_vectorized_and_positive():
    lifetime = StellarLifetime()
    masses = np.linspace(0.5, 50, 100)
    tau = lifetime(masses, metallicity=0.0)
    assert tau.shape == masses.shape
    assert np.all(tau > 0)

def test_popIII_explode_earlier():
    lifetime = StellarLifetime()
    masses = np.linspace(0.5, 50, 100)
    tau_popIII = lifetime(masses, metallicity=-5)
    tau_popII  = lifetime(masses, metallicity=ZSUN)
    assert np.all(tau_popIII <= tau_popII + 1e-12)
