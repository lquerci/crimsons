from __future__ import annotations

import numpy as np

from ..chemistry import check_metallicity
from ..yields.base import PopulationChannel, StochasticYieldTable


def run_realization(
    imf,
    mass_formed,
    metallicity,
    channels,
    lifetime_fn,
    time_grid,
    rng,
    n_bins: int = 1000,
    chunk_size: int = 1_000_000,
):
    """Sample one stellar population -- as a log-mass-binned histogram,
    not one entry per star (see IMF.sample_binned) -- and compute its
    enrichment events.

    Each "star" here is really a bin: `bin_masses[i]` is the mean mass of
    `bin_counts[i]` stars. Channel selection (`contributes`) and the
    per-selected-star calculations (`delay_time`, `yields`) run once per
    bin exactly as they used to run once per star; the only mechanical
    difference is that each bin's yield is scaled by its count before
    being recorded as an enrichment event, so a bin of 50,000 stars
    contributes 50,000x what one of them would.

    Caveat -- this is exact for a deterministic mass-range channel (SNII,
    AGB): a bin is unambiguously either in range or not, so scaling by
    count reproduces sampling every star individually. It is NOT exact
    for a *stochastic* channel like the current SNIa, whose `contributes`
    /`delay_time` draw once per input element: applied to a bin, that
    means the bin's entire count either all becomes SN Ia progenitors
    exploding at one shared delay time, or none do, rather than each
    underlying star getting its own independent binary/DTD draw. That's a
    known, temporary approximation -- SN Ia's channel implementation is
    being reworked separately -- and doesn't affect SNII/AGB.

    The same "once per bin, not once per star" caveat applies to yield
    tables with a *stochastic* extra parameter (e.g. SN II rotation
    velocity drawn from a distribution rather than fixed): every star in
    a bin shares that bin's one draw. With the default n_bins=1000 this
    is a fine approximation (a bin at the top of a steep IMF's mass range
    might represent thousands of stars, but even there each bin spans a
    narrow slice in log-mass); it stops being fine if n_bins is small
    relative to how fast the extra parameter's effect on yields varies.

    Returns
    -------
    bin_masses, bin_counts : as returned by IMF.sample_binned
    fates : (n_populated_bins,) object array naming which channel each
        bin went through ("none" if it doesn't contribute to any tracked
        channel)
    events : list of (channel_name, delay_times, yields) per channel that
        contributed -- yields are already scaled by each bin's count.
    """

    check_metallicity(metallicity=metallicity)


    bin_masses, bin_counts = imf.sample(
        rng, mass_formed, n_bins=n_bins, chunk_size=chunk_size
    )
    lifetimes = lifetime_fn(bin_masses, metallicity)
    fates = np.full(bin_masses.shape, "none", dtype=object)
    events = []

    for channel in channels:
        # BRANCH 1: PopulationChannel (e.g. SNIa DTD) ---
        if isinstance(channel, PopulationChannel):
            pop_res = channel.population_events(
                mass_formed=mass_formed,
                metallicity=metallicity,
                lifetime_fn=lifetime_fn,
                time_grid=time_grid,
                rng=rng,
            )
            if pop_res is not None:
                t_pop, y_pop, n_events_pop = pop_res
                events.append((channel.name, np.asarray(t_pop), np.asarray(y_pop), np.asarray(n_events_pop)))
        # BRANCH 2: MassRangeChannel (AGB, SNII, PISN)
        else:
            mask = channel.contributes(bin_masses, metallicity, rng)
            if not np.any(mask):
                continue
            m_sub = bin_masses[mask]
            c_sub = bin_counts[mask]
            t_sub = channel.delay_time(m_sub, metallicity, lifetimes[mask], rng)

            is_stochastic = isinstance(channel._yield_table, StochasticYieldTable)            

            # Decide route: Use slow (individual star) sampling only if it's stochastic 
            # AND the total number of stars in this subset is small enough to avoid memory spikes.
            total_stars_in_subset = np.sum(c_sub)
            use_individual_sampling = is_stochastic and (total_stars_in_subset < 50_000)

            if not use_individual_sampling: 
                # FAST : sample the yields extra dimension once per mass bin
                y_sub = channel.yields(m_sub, metallicity, rng) * c_sub[:, None]
            else:                    
                # SLOW: sample the yields extra dimension once per star

                # ensure counts are integers
                counts = np.round(c_sub).astype(int)

                # unroll bin into a flat array
                unrolled_m = np.repeat(m_sub, counts)

                # get the yield for each individual stellar particle
                unrolled_y = channel.yields(unrolled_m, metallicity, rng)

                # group sum the yields
                bin_indices = np.repeat(np.arange(len(m_sub)), counts)
                y_sub = np.zeros((len(m_sub), unrolled_y.shape[1]))
                np.add.at(y_sub, bin_indices, unrolled_y)

            
            events.append((channel.name, np.asarray(t_sub), np.asarray(y_sub), np.asarray(c_sub)))
            fates[mask] = channel.name

    return bin_masses, bin_counts, fates, events


def bin_enrichment(events, time_grid, n_elements):
    """Bin per-star-group enrichment events onto a time grid.

    Returns the *cumulative* mass of each element returned to the ISM by
    each point on time_grid, shape (n_time, n_elements). This is the raw
    enrichment injected by the stellar population -- combining it with a
    star formation history, gas inflow/outflow, etc. into a full one-zone
    chemical evolution model is left to the caller; this library computes
    the enrichment *source term*, not a full GCE model.

    Unchanged by binning: `events`' yields are already scaled by each
    bin's star count (see run_realization), so this just sums them onto
    the time grid exactly as it always has.
    """
    time_grid = np.asarray(time_grid)
    increments = np.zeros((len(time_grid), n_elements))
    for _name, times, yields, n_events in events:
        # Find the first index where time_grid[idx] >= times
        idx = np.searchsorted(time_grid, times, side="left")
        idx = np.clip(idx, 0, len(time_grid) - 1)
        np.add.at(increments, idx, yields)
    return np.cumsum(increments, axis=0)


def run_and_bin_realization(
    imf,
    mass_formed,
    metallicity,
    channels,
    lifetime_fn,
    time_grid,
    seed,
    n_bins: int,
    chunk_size: int,
    n_columns: int,
):
    """Run one full realization (`run_realization`) and bin it onto
    `time_grid` (`bin_enrichment`) in a single call -- this is the unit
    of work `Simulation.run` hands out, whether it's running serially
    (a plain for loop calling this once per realization) or in parallel
    (a worker process calling this once per submitted task -- see
    `crimsons.simulation.ensemble`).

    Doing both steps here, rather than in the caller, keeps a parallel
    worker's return value small regardless of population size: yields
    are already reduced to a per-time-bin enrichment array and events
    have already had their (potentially large) per-star-group yields
    arrays stripped (keeping only name/times/counts) before anything
    needs to cross a process boundary -- the same memory-saving trick
    `Simulation.run` already did in its serial loop.

    `seed` is anything `numpy.random.default_rng` accepts (an int, a
    `numpy.random.SeedSequence`, ...) rather than a ready-made
    `Generator`, since a `Generator`'s state doesn't need to (and for a
    freshly-spawned worker process, can't meaningfully) come from the
    parent process -- each call gets its own independent stream.

    Returns
    -------
    bin_masses, bin_counts, fates : as returned by `run_realization`
    memory_safe_events : list of (channel_name, times, counts) -- same
        as `run_realization`'s `events`, minus the yields array
    enrichment_row : (len(time_grid), n_columns) array, as returned by
        `bin_enrichment`
    """
    rng = np.random.default_rng(seed)
    bin_masses, bin_counts, fates, events = run_realization(
        imf,
        mass_formed,
        metallicity,
        channels,
        lifetime_fn,
        time_grid,
        rng,
        n_bins=n_bins,
        chunk_size=chunk_size,
    )
    enrichment_row = bin_enrichment(events, time_grid, n_columns)
    memory_safe_events = [(name, times, counts) for name, times, yields, counts in events]
    return bin_masses, bin_counts, fates, memory_safe_events, enrichment_row