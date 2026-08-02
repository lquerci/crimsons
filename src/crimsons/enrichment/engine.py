from __future__ import annotations

import numpy as np


def run_realization(
        imf,
        mass_formed,
        metallicity,
        channels,
        lifetime_fn,
        rng,
        n_bins: int = 1000,
        chunk_size: int = 1_000_000,        
    ):
    """Sample one stellar population and compute its enrichment events.

    Everything here operates on whole arrays of stars at once (sample all
    masses -> bin them ->  vectorized lifetimes -> vectorized fate masks -> vectorized
    yield interpolation) rather than looping star by star -- that's what
    keeps a pure-Python/numpy implementation fast.

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

    # add a filter for empty bins?

    lifetimes = lifetime_fn(bin_masses, metallicity)
    fates = np.full(bin_masses.shape, "none", dtype=object)
    events = []
    for channel in channels:
        mask = channel.contributes(bin_masses, metallicity, rng)
        if not np.any(mask):
            continue
        m_sub = bin_masses[mask]
        c_sub = bin_counts[mask]
        t_sub = channel.delay_time(m_sub, metallicity, lifetimes[mask], rng)
        y_sub = channel.yields(m_sub, metallicity) * c_sub[:, None]
        events.append((channel.name, np.asarray(t_sub), np.asarray(y_sub)))
        fates[mask] = channel.name

    return bin_masses, bin_counts, fates, events


def bin_enrichment(events, time_grid, n_elements):
    """Bin per-star enrichment events onto a time grid.

    Returns the *cumulative* mass of each element returned to the ISM by
    each point on time_grid, shape (n_time, n_elements). This is the raw
    enrichment injected by the stellar population -- combining it with a
    star formation history, gas inflow/outflow, etc. into a full one-zone
    chemical evolution model is left to the caller; this library computes
    the enrichment *source term*, not a full GCE model.
    """
    time_grid = np.asarray(time_grid)
    increments = np.zeros((len(time_grid), n_elements))
    for _name, times, yields in events:
        idx = np.searchsorted(time_grid, times, side="right") - 1
        idx = np.clip(idx, 0, len(time_grid) - 1)
        np.add.at(increments, idx, yields)
    return np.cumsum(increments, axis=0)