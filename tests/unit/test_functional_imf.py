import numpy as np
import pytest

from crimsons.imf.functional import FunctionalIMF
from crimsons.imf.standard import Salpeter1955


def test_matches_analytic_salpeter_within_tolerance():
    """The numerical grid-based inverse CDF should closely reproduce the
    exact analytic one for the same underlying shape -- this validates
    the numerics, not just that the code runs."""
    slope = 2.35
    functional = FunctionalIMF(lambda m: m**-slope, m_min=0.1, m_max=100.0, n_grid=4000)
    analytic = Salpeter1955(m_min=0.1, m_max=100.0, slope=slope)

    u = np.linspace(0.01, 0.99, 25)
    m_functional = functional.inverse_cdf(u)
    m_analytic = analytic.inverse_cdf(u)

    np.testing.assert_allclose(m_functional, m_analytic, rtol=0.01)


def test_sampling_reaches_target_mass():
    imf = FunctionalIMF(lambda m: m**-2.35, m_min=0.1, m_max=100.0)
    rng = np.random.default_rng(0)
    bin_masses, bin_counts = imf.sample(rng, target_mass=1000.0)
    assert (bin_masses*bin_counts).sum() >= 1000.0
    assert np.all(bin_masses >= 0.1) and np.all(bin_counts <= 100.0)


def test_accepts_non_vectorized_scalar_only_function():
    def scalar_shape(m):
        # deliberately not numpy-friendly: math.exp on a plain float only
        import math

        return math.exp(-m)

    imf = FunctionalIMF(scalar_shape, m_min=0.5, m_max=5.0, n_grid=500)
    rng = np.random.default_rng(1)
    bin_masses, bin_counts = imf.sample(rng, target_mass=50.0)
    assert (bin_masses*bin_counts).sum() >= 50.0


def test_flat_shape_samples_roughly_uniformly():
    imf = FunctionalIMF(lambda m: np.ones_like(m), m_min=1.0, m_max=3.0, n_grid=1000)
    rng = np.random.default_rng(2)
    u = rng.random(20_000)
    masses = imf.inverse_cdf(u)
    # uniform on [1, 3] has mean 2.0; with 20k draws the sample mean
    # should land well within a generous tolerance
    assert masses.mean() == pytest.approx(2.0, abs=0.05)


def test_negative_density_raises():
    with pytest.raises(ValueError):
        FunctionalIMF(lambda m: -np.ones_like(m), m_min=0.1, m_max=10.0)


def test_metallicity_adaptive_range_also_works_for_functional_imf():
    imf = FunctionalIMF(lambda m: m**-2.35, metallicity=1e-7)
    assert (imf.m_min, imf.m_max) == (0.8, 1000.0)
