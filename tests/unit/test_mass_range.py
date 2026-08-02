import numpy as np

from crimsons.imf.mass_range import (
    EXTREMELY_METAL_POOR_MASS_RANGE,
    POPIII_THRESHOLD,
    STANDARD_MASS_RANGE,
    default_mass_range,
)
from crimsons.imf.standard import Salpeter1955


def test_default_mass_range_switches_below_threshold():
    assert default_mass_range(0.02) == STANDARD_MASS_RANGE
    assert default_mass_range(POPIII_THRESHOLD * 10) == STANDARD_MASS_RANGE
    assert default_mass_range(POPIII_THRESHOLD / 10) == EXTREMELY_METAL_POOR_MASS_RANGE
    assert default_mass_range(None) == STANDARD_MASS_RANGE


def test_imf_with_no_args_uses_standard_range():
    imf = Salpeter1955()
    assert (imf.m_min, imf.m_max) == STANDARD_MASS_RANGE


def test_imf_constructed_with_low_metallicity_gets_wide_range():
    imf = Salpeter1955(metallicity=1e-7)
    assert (imf.m_min, imf.m_max) == EXTREMELY_METAL_POOR_MASS_RANGE


def test_explicit_bounds_are_never_overridden_by_metallicity():
    imf = Salpeter1955(m_min=0.5, m_max=50.0, metallicity=1e-3)
    assert (imf.m_min, imf.m_max) == (0.5, 50.0)
    imf.bind_metallicity(1e-8)  # should be a no-op
    assert (imf.m_min, imf.m_max) == (0.5, 50.0)


def test_bind_metallicity_updates_implicit_default_and_rebuilds():
    imf = Salpeter1955()  # implicit default, resolves to standard range
    assert (imf.m_min, imf.m_max) == STANDARD_MASS_RANGE

    imf.bind_metallicity(1e-7)
    assert (imf.m_min, imf.m_max) == EXTREMELY_METAL_POOR_MASS_RANGE

    # the rebuilt IMF should actually sample within the new range, not
    # the stale one from before bind_metallicity
    rng = np.random.default_rng(0)
    bin_masses, bin_counts = imf.sample(rng, target_mass=2000.0)
    assert np.all(bin_masses >= EXTREMELY_METAL_POOR_MASS_RANGE[0])
    assert np.all(bin_counts <= EXTREMELY_METAL_POOR_MASS_RANGE[1])


def test_custom_mass_range_fn_is_honored():
    def always_narrow(metallicity):
        return (1.0, 2.0)

    imf = Salpeter1955(mass_range_fn=always_narrow, metallicity=1e-6)
    assert (imf.m_min, imf.m_max) == (1.0, 2.0)
