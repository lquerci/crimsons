import numpy as np

from crimsons.chemistry import HDF_COLUMNS, ZSUN
from crimsons.enrichment.engine import (
    bin_enrichment,
    run_and_bin_realization,
    run_realization,
)
from crimsons.imf.standard import Salpeter1955
from crimsons.stars.lifetimes import StellarLifetime
from crimsons.yields.channels import default_channels


def test_run_realization_smoke():
    imf = Salpeter1955()
    lifetime_fn = StellarLifetime()
    channels = default_channels(ZSUN)
    rng = np.random.default_rng(0)

    bin_masses, bin_counts, fates, events = run_realization(
        imf,
        mass_formed=5000.0,
        metallicity=ZSUN,
        channels=channels,
        lifetime_fn=lifetime_fn,
        time_grid=np.logspace(0,4,50),
        rng=rng,
    )

    assert (bin_masses * bin_counts).sum() >= 5000.0
    assert bin_counts.shape == bin_masses.shape
    assert fates.shape == bin_masses.shape
    assert set(np.unique(fates)) <= {"none", "SNII", "AGB", "SNIa"}

    time_grid = np.linspace(0, 13.8, 50)
    enrichment = bin_enrichment(events, time_grid, len(HDF_COLUMNS))

    assert enrichment.shape == (50, len(HDF_COLUMNS))
    assert np.all(enrichment >= 0)
    # cumulative enrichment must be non-decreasing in time
    assert np.all(np.diff(enrichment, axis=0) >= -1e-12)


def test_event_yields_are_scaled_by_bin_count():
    """A bin representing many stars should contribute proportionally
    more enrichment than one representing few -- this is the actual
    mechanical point of threading counts through the engine."""
    imf = Salpeter1955()
    lifetime_fn = StellarLifetime()
    channels = default_channels(metallicity=ZSUN)

    # a small population: some bins will have small counts
    _, small_counts, _, small_events = run_realization(
        imf, mass_formed=2000.0, metallicity=ZSUN, channels=channels,
        lifetime_fn=lifetime_fn, rng=np.random.default_rng(3), 
        time_grid=np.logspace(0,4,50),
    )
    # a much larger population with the same seed-derived shape: bins
    # should carry larger counts and proportionally larger event yields
    _, large_counts, _, large_events = run_realization(
        imf, mass_formed=200_000.0, metallicity=ZSUN, channels=channels,
        lifetime_fn=lifetime_fn, rng=np.random.default_rng(3),
        time_grid=np.logspace(0,4,50),
    )

    assert large_counts.sum() > small_counts.sum()

    def total_yield(events):
        return sum(float(np.sum(yields)) for _, _, yields, _ in events)

    assert total_yield(large_events) > total_yield(small_events)


def test_run_and_bin_realization_matches_manual_run_then_bin():
    """run_and_bin_realization (the unit of work Simulation.run hands
    out, serially or to a worker process) must be exactly equivalent to
    calling run_realization then bin_enrichment by hand and stripping
    yields from events -- it's a pure refactor, not new behavior."""
    imf = Salpeter1955()
    lifetime_fn = StellarLifetime()
    channels = default_channels(ZSUN)
    time_grid = np.logspace(0, 4, 50)
    seed = 12345

    rng = np.random.default_rng(seed)
    exp_bin_masses, exp_bin_counts, exp_fates, exp_events = run_realization(
        imf, mass_formed=5000.0, metallicity=ZSUN, channels=channels,
        lifetime_fn=lifetime_fn, time_grid=time_grid, rng=rng,
    )
    exp_enrichment_row = bin_enrichment(exp_events, time_grid, len(HDF_COLUMNS))
    exp_memory_safe_events = [(name, times, counts) for name, times, yields, counts in exp_events]

    bin_masses, bin_counts, fates, memory_safe_events, enrichment_row = run_and_bin_realization(
        imf, mass_formed=5000.0, metallicity=ZSUN, channels=channels,
        lifetime_fn=lifetime_fn, time_grid=time_grid, seed=seed,
        n_bins=1000, chunk_size=1_000_000, n_columns=len(HDF_COLUMNS),
    )

    np.testing.assert_allclose(bin_masses, exp_bin_masses)
    np.testing.assert_allclose(bin_counts, exp_bin_counts)
    assert np.array_equal(fates, exp_fates)
    np.testing.assert_allclose(enrichment_row, exp_enrichment_row)
    assert len(memory_safe_events) == len(exp_memory_safe_events)
    for (name, times, counts), (exp_name, exp_times, exp_counts) in zip(
        memory_safe_events, exp_memory_safe_events
    ):
        assert name == exp_name
        np.testing.assert_allclose(times, exp_times)
        np.testing.assert_allclose(counts, exp_counts)
        # yields must actually be stripped -- only 3 fields, not 4
        assert len(memory_safe_events[0]) == 3


def test_run_and_bin_realization_accepts_a_seed_not_just_a_generator():
    """This is the whole point of taking `seed` instead of a ready-made
    Generator: a freshly-spawned worker process needs to build its own
    rng from a seed (or SeedSequence), not reuse parent-process state."""
    imf = Salpeter1955()
    lifetime_fn = StellarLifetime()
    channels = default_channels(ZSUN)
    time_grid = np.logspace(0, 4, 20)

    seed_seq = np.random.SeedSequence(0).spawn(1)[0]
    result_int = run_and_bin_realization(
        imf, mass_formed=1000.0, metallicity=ZSUN, channels=channels,
        lifetime_fn=lifetime_fn, time_grid=time_grid, seed=seed_seq,
        n_bins=200, chunk_size=1_000_000, n_columns=len(HDF_COLUMNS),
    )
    # same SeedSequence given twice -> identical results (determinism is
    # what makes n_jobs safe to change without affecting physics)
    result_int_2 = run_and_bin_realization(
        imf, mass_formed=1000.0, metallicity=ZSUN, channels=channels,
        lifetime_fn=lifetime_fn, time_grid=time_grid, seed=seed_seq,
        n_bins=200, chunk_size=1_000_000, n_columns=len(HDF_COLUMNS),
    )
    np.testing.assert_allclose(result_int[-1], result_int_2[-1])