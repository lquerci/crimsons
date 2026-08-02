import numpy as np
import pytest

from crimsons.imf.standard import Kroupa2001, Salpeter1955, FlatIMF


@pytest.mark.parametrize("imf_cls", [Salpeter1955, Kroupa2001, FlatIMF])
def test_sample_reaches_target_mass_without_wild_overshoot(imf_cls):
    imf = imf_cls()
    rng = np.random.default_rng(0)
    bin_masses, bin_counts = imf.sample(rng, target_mass=1000.0)

    assert (bin_masses*bin_counts).sum() >= 1000.0
    assert (bin_masses*bin_counts).sum() < 1000.0 + bin_masses.max()
    assert np.all(bin_masses >= imf.m_min)
    assert np.all(bin_masses <= imf.m_max)
    assert np.all(bin_counts >= 0)


@pytest.mark.parametrize("imf_cls", [Salpeter1955, Kroupa2001, FlatIMF])
def test_sampling_is_reproducible_with_same_seed(imf_cls):
    imf = imf_cls()
    m1 = imf.sample(np.random.default_rng(42), target_mass=500.0)
    m2 = imf.sample(np.random.default_rng(42), target_mass=500.0)
    np.testing.assert_array_equal(m1, m2)


@pytest.mark.parametrize("imf_cls", [Salpeter1955, Kroupa2001, FlatIMF])
def test_inverse_cdf_covers_full_mass_range(imf_cls):
    imf = imf_cls()
    lo = imf.inverse_cdf(np.array([1e-6]))[0]
    hi = imf.inverse_cdf(np.array([1 - 1e-6]))[0]
    assert lo == pytest.approx(imf.m_min, rel=0.05)
    assert hi == pytest.approx(imf.m_max, rel=0.05)
