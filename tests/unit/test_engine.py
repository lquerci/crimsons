import numpy as np

from crimsons.chemistry import HDF_COLUMNS, ZSUN
from crimsons.enrichment.engine import bin_enrichment, run_realization
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