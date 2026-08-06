from importlib.resources import files

import numpy as np
import pytest

from crimsons.chemistry import ELEMENTS, ZSUN
from crimsons.yields.channels import (
    SNIa,
    )

H5_PATH = files("crimsons.yields") / "data" / "stellar_yields.h5"

# --- SNIa INITIALIZATION TESTS ---

def test_snia_initializes_with_default_model_and_params():
    """Test the fix for the SNIa/Iwamoto default model_params bug."""
    snia = SNIa()
    assert snia.model == "Iwamoto"
    assert snia.name == "SNIa"
    # Ensure it successfully picked W7 as the default model
    assert snia._yield_table.model_params["model"] == "W7"

def test_snia_initializes_with_explicit_model_params():
    """Test that explicit model parameters (like WDD2) are respected."""
    snia = SNIa(model="Iwamoto", model_params={"model": "WDD2"})
    assert snia.model == "Iwamoto"
    assert snia._yield_table.model_params["model"] == "WDD2"

# --- POPULATION CHANNEL ARCHITECTURE TESTS ---

def test_snia_is_population_channel_and_overrides_base_methods():
    """SNIa must bypass per-star evaluations (contributes/delay_time)."""
    snia = SNIa()
    rng = np.random.default_rng(42)
    masses = np.array([1.0, 3.0, 5.0, 8.0, 10.0])
    
    # contributes() should strictly return all False for a PopulationChannel
    mask = snia.contributes(masses, ZSUN, rng)
    assert not np.any(mask)
    assert mask.shape == masses.shape
    assert mask.dtype == bool
    
    # delay_time() should explicitly raise an error to prevent accidental usage
    with pytest.raises(NotImplementedError, match="PopulationChannel"):
        snia.delay_time(masses, ZSUN, np.ones_like(masses), rng)

# --- POPULATION EVENTS EXECUTION TESTS ---

def mock_lifetime(mass, metallicity):
    """Simple mock lifetime function returning Myr (t ~ 10000 * M^-2.5)."""
    return 10000.0 * (np.asarray(mass) / 1.0)**-2.5

def test_snia_population_events_dtd_mode():
    """Test the Maoz DTD event generation and exact return shapes."""
    # mass_min=3.0, mass_max=8.0 is the corrected intermediate mass range
    snia = SNIa(mode="dtd", dtd_shape="maoz", mass_min=3.0, mass_max=8.0)
    time_grid = np.logspace(1, 4, 50)  # 10 Myr to 10,000 Myr
    rng = np.random.default_rng(42)
    mass_formed = 1e6  # Large enough to guarantee multiple events
    
    pop_res = snia.population_events(mass_formed, ZSUN, mock_lifetime, time_grid, rng)
    
    assert pop_res is not None
    times, yields, counts = pop_res  # Expecting the updated 3-tuple
    
    # 1. Times should be a subset of the time_grid
    assert len(times) > 0
    assert np.all(np.isin(times, time_grid))
    
    # 2. Counts should be integers > 0
    assert np.all(counts > 0)
    assert np.all(counts == np.floor(counts))  # checks that e.g. 46.0 == 46.0
    
    # 3. Yields should be appropriately shaped (n_events, n_elements)
    assert yields.shape == (len(times), len(ELEMENTS))
    assert np.all(yields >= 0)
    
    # SNIa produce lots of Iron -- check that Fe is present in the yields
    fe_idx = ELEMENTS.index("Fe")
    assert np.sum(yields[:, fe_idx]) > 0

def test_snia_population_events_single_burst_mode():
    """Test that a single burst localizes all SNIa to one time bin."""
    snia = SNIa(mode="single_burst", burst_delay_myr=300.0, rate_per_msun=1e-3)
    time_grid = np.logspace(1, 4, 50) 
    rng = np.random.default_rng(42)
    mass_formed = 1e5
    
    pop_res = snia.population_events(mass_formed, ZSUN, mock_lifetime, time_grid, rng)
    
    assert pop_res is not None
    times, _yields, counts = pop_res
    
    # Single burst should only occur at exactly one timestep
    assert len(times) == 1
    
    # It should clamp to the closest time grid point <= 300 Myr
    expected_idx = min(int(np.searchsorted(time_grid, 300.0)), len(time_grid) - 1)
    assert times[0] == time_grid[expected_idx]
    
    # Total counts should equal mass_formed * rate_per_msun (integerized)
    expected_counts = int(mass_formed * 1e-3)
    assert counts[0] == pytest.approx(expected_counts, abs=1)

def test_population_events_returns_none_if_mass_too_low():
    """If the population mass is tiny, no integer events should occur."""
    snia = SNIa()
    time_grid = np.logspace(1, 4, 50)
    rng = np.random.default_rng(42)
    
    # Forming 1 solar mass should yield 0 SNIa events
    pop_res = snia.population_events(1.0, ZSUN, mock_lifetime, time_grid, rng)
    assert pop_res is None