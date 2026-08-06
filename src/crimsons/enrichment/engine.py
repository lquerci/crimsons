from __future__ import annotations

import numpy as np

from ..yields.base import PopulationChannel


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
            y_sub = channel.yields(m_sub, metallicity, rng) * c_sub[:, None]
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
        idx = np.searchsorted(time_grid, times, side="right") - 1
        idx = np.clip(idx, 0, len(time_grid) - 1)
        np.add.at(increments, idx, yields)
    return np.cumsum(increments, axis=0)