from importlib.resources import files

import numpy as np
import pytest

from crimsons.chemistry import HDF_COLUMNS, ZSUN
from crimsons.yields.channels import (
    SNIa,
    )
from crimsons.yields.io import describe_model

H5_PATH = files("crimsons.yields") / "data" / "stellar_yields.h5"

# --- SNIa INITIALIZATION TESTS ---

def test_snia_initializes_with_default_model_and_params():
    """Test the fix for the SNIa/Iwamoto default model_params bug."""
    snia = SNIa()
    assert snia.model == "Iwamoto99"
    assert snia.name == "SNIa"
    # Ensure it successfully picked W7 as the default model
    assert snia._yield_table.model_params["model"] == "W7"

def test_snia_initializes_with_explicit_model_params():
    """Test that explicit model parameters (like WDD2) are respected."""
    snia = SNIa(model="Iwamoto99", model_params={"model": "WDD2"})
    assert snia.model == "Iwamoto99"
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
    assert yields.shape == (len(times), len(HDF_COLUMNS))
    assert np.all(yields >= 0)
    
    # SNIa produce lots of Iron -- check that Fe is present in the yields
    fe_idx = HDF_COLUMNS.index("Fe")
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

def test_describe_iwamoto_shows_mapped_metallicity():
    """Test that describe_model exposes 'metallicity' even though it's 
    not in the 'axes' attribute for SNIa/Iwamoto.
    """
    info = describe_model(H5_PATH, "SNIa", "Iwamoto99")
    
    # It should not be an independent grid axis
    assert "metallicity" not in info["axes"]
    assert info["axes"] == ["model", "mass", "elements"]
    
    # But it must be present in the metadata payload
    assert "metallicity" in info
    # The metallicity array should perfectly map 1-to-1 with the models
    assert len(info["metallicity"]) == len(info["model"])


def test_snia_routes_to_w70_for_pop_iii():
    """Test that initializing SNIa with Pop III metallicity automatically 
    selects the W70 model and isolates the table to the Pop III regime.
    """
    snia = SNIa(metallicity=1e-7)
    
    # Did it pick the right fallback model?
    assert snia._yield_table.model_params["model"] == "W70"
    
    # Did it correctly build the 1D metallicity grid from the mapped dataset?
    table_z = snia._yield_table.metallicities
    assert len(table_z) == 1
    assert table_z[0] == 1e-7
    
    # Did YieldTable file this under Population III?
    assert "Population III" in snia._yield_table._regimes
    assert "Population II/I" not in snia._yield_table._regimes


def test_snia_routes_to_w7_for_pop_ii():
    """Test that initializing SNIa with Pop II/I metallicity automatically 
    selects the W7 model and isolates the table to the Pop II/I regime.
    """
    snia = SNIa(metallicity=ZSUN)
    
    assert snia._yield_table.model_params["model"] == "W7"
    
    # YieldTable should file this under Population II/I
    assert "Population II/I" in snia._yield_table._regimes
    assert "Population III" not in snia._yield_table._regimes


def test_snia_yield_table_raises_wrong_regime():
    """Test that querying a loaded SNIa table with a metallicity belonging 
    to a different population regime raises a ValueError.
    """
    # Initialize with Pop III (loads W70)
    snia = SNIa(metallicity=1e-7) 
    rng = np.random.default_rng(42)
    
    # Attempting to query the W70 Pop III table at Z=ZSUN (Pop II) must crash
    with pytest.raises(ValueError, match="no yield data for Population II/I"):
        snia._yield_table(mass=1.2, metallicity=ZSUN, rng=rng)


def test_snia_mapped_metallicity_shape_and_content():
    """Test that the loaded 3D matrix collapses correctly and contains
    the specific explosion energy expected from the append script.
    """
    snia = SNIa(model="Iwamoto99", model_params={"model": "W7"})
    table = snia.yield_table()
    
    # Ensure dimensions collapsed appropriately: (1 mass, 1 metallicity, 31 elements)
    assert table.yields.shape == (1, 1, 31)
    
    # Energy is at index 0 and should be 1.2e51 (from SNIA_EXPLOSION_ENERGY)
    energy_idx = HDF_COLUMNS.index("Energy")
    assert table.yields[0, 0, energy_idx] == pytest.approx(1.2e51)