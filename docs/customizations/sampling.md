# Stochastic Sampling

CRIMSONS samples stellar populations **stochastically but in bins**, not
star-by-star -- this is what keeps a $10^6\,M_\odot$ (or larger) realization
tractable.

## Drawing to a target mass

`IMF.sample(rng, target_mass, n_bins=1000, chunk_size=1_000_000)` draws
stars in vectorized chunks from `inverse_cdf(uniform(0, 1))`, accumulating
their summed mass, until the cumulative total reaches `target_mass`. The
chunk that crosses the target is trimmed star-by-star to the exact draw
that pushes the total over -- the final star may push the sum slightly
past `target_mass` (you can't include a fractional star), which is
standard for this kind of sampling.

## Binning, not enumerating

Rather than returning every individual sampled mass, `IMF.sample` bins
them onto `n_bins` **log-spaced** mass edges spanning $[m_{\min}, m_{\max}]$
and returns two arrays:

- `bin_masses` -- each bin's geometric center, $\sqrt{m_{\rm lo}\, m_{\rm hi}}$
- `bin_counts` -- how many sampled stars landed in that bin

Every downstream calculation -- `contributes`, `delay_time`, `yields` --
runs once **per bin** rather than once per star, and each bin's yield is
scaled by its star count before being recorded as an enrichment event.
Unpopulated bins are kept (with `count = 0`) rather than filtered out, so
every realization's arrays share the same fixed length and shape.

## Where the bin approximation is exact, and where it isn't

For a **deterministic mass-range channel** (`SNII`, `AGB`, `PISN`): a bin
is unambiguously either inside or outside `[mass_min, mass_max]`, so
scaling its yield by its star count reproduces sampling every star
individually -- no approximation.

For a **yield table with a stochastic extra model parameter** (e.g. an SN
II rotation velocity drawn from a distribution rather than fixed via
`model_params`, see [Selecting Yield Models](../examples/yield-models.md)):
every star sharing a bin also shares that bin's *one* draw, rather than
each getting an independent draw. With the default `n_bins=1000` this is a
good approximation (even a bin at the top of a steep IMF's mass range,
representing many thousands of stars, still spans a narrow slice in
log-mass) -- but it degrades if `n_bins` is small relative to how quickly
the extra parameter's effect on yields varies with mass.

`SNIa`, by contrast, doesn't go through this per-bin path at all -- it's a
`PopulationChannel`, handled by `population_events` directly on the
population as a whole (see [Enrichment Channels](channels.md)).

## Reproducibility across an ensemble

`Simulation(..., seed=...)` uses `numpy.random.SeedSequence(seed).spawn(n_realizations)`
to derive one independent child seed per realization. The same master
seed always reproduces the same set of realizations, both individually
and as an ensemble; `seed=None` draws fresh entropy each time.

## Tuning `n_bins`

`Simulation(..., n_bins=...)` (default `1000`) trades accuracy for speed:
more bins track the IMF's shape more finely (and better resolve any
per-bin stochastic yield draws) at the cost of more work per realization.
`chunk_size` (`Simulation(..., sample_chunk_size=...)`, default
`1_000_000`) controls how many candidate stars are drawn per vectorized
batch while filling the target mass, and rarely needs tuning.
