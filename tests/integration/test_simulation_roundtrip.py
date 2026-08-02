import numpy as np

from crimsons import (
    PISN,
    EnrichmentResult,
    Salpeter1955,
    Simulation,
    StellarLifetime,
    default_channels,
)


def _make_simulation(**overrides):
    kwargs = {
        "imf": Salpeter1955(),
        "lifetime_fn": StellarLifetime(),
        "channels": default_channels(),
        "mass_formed": 2000.0,
        "metallicity": 0.02,
        "n_realizations": 3,
        "seed": 123,
        "time_grid": np.linspace(0, 13.8, 30),
    }
    kwargs.update(overrides)
    return Simulation(**kwargs)


def test_simulation_runs_and_aggregates():
    result = _make_simulation().run()
    assert result.enrichment.shape == (3, 30, len(result.elements))
    assert result.mean().shape == (30, len(result.elements))
    assert result.std().shape == (30, len(result.elements))


def test_save_and_load_roundtrip(tmp_path):
    result = _make_simulation(seed=7, n_realizations=2).run()
    path = tmp_path / "result.h5"
    result.save(path)

    loaded = EnrichmentResult.load(path)

    np.testing.assert_allclose(loaded.enrichment, result.enrichment)
    assert loaded.elements == result.elements
    assert loaded.config.run_id() == result.config.run_id()
    assert len(loaded.masses) == len(result.masses)
    assert len(loaded.counts) == len(result.counts)
    for i in range(len(result.masses)):
        np.testing.assert_allclose(loaded.masses[i], result.masses[i])
        np.testing.assert_allclose(loaded.counts[i], result.counts[i])
        assert np.array_equal(loaded.fates[i], result.fates[i])


def test_simulation_resolves_metallicity_adaptive_imf_range():
    # imf is constructed with no explicit bounds and no metallicity of its
    # own -- Simulation's metallicity should drive the default range
    sim = _make_simulation(metallicity=1e-7, imf=Salpeter1955())
    result = sim.run()
    assert (sim.imf.m_min, sim.imf.m_max) == (0.8, 1000.0)
    assert result.config.imf_mass_min == 0.8
    assert result.config.imf_mass_max == 1000.0


def test_simulation_with_pisn_runs_end_to_end():
    """PISN alongside SNII/AGB/SNIa, with a large enough mass_formed and
    wide enough IMF range that some bins should actually land in PISN's
    140-260 Msun window."""
    sim = _make_simulation(
        channels=[*default_channels(), PISN()],
        imf=Salpeter1955(m_min=0.8, m_max=1000.0),
        metallicity=1e-6,
        mass_formed=500_000.0,
        n_realizations=2,
    )
    result = sim.run()
    all_fates = np.concatenate(result.fates)
    assert set(np.unique(all_fates)) <= {"none", "SNII", "AGB", "SNIa", "PISN"}
    assert "PISN" in all_fates  # confirms the channel actually fired, not just wired up
    assert np.all(result.enrichment >= 0)


def test_cache_skips_recomputation_and_returns_same_result(tmp_path):
    sim = _make_simulation(seed=99, mass_formed=800.0, n_realizations=2)
    cache_dir = tmp_path / "cache"

    result1 = sim.run(cache_dir=cache_dir)
    result2 = sim.run(cache_dir=cache_dir)  # should load from cache

    np.testing.assert_allclose(result1.enrichment, result2.enrichment)
    assert list(cache_dir.glob("*.h5"))  # cache file was actually written


def test_different_n_bins_are_different_cache_entries(tmp_path):
    cache_dir = tmp_path / "cache"
    sim_coarse = _make_simulation(seed=5, n_bins=50)
    sim_fine = _make_simulation(seed=5, n_bins=200)

    sim_coarse.run(cache_dir=cache_dir)
    sim_fine.run(cache_dir=cache_dir)

    assert sim_coarse.config().run_id() != sim_fine.config().run_id()
    assert len(list(cache_dir.glob("*.h5"))) == 2