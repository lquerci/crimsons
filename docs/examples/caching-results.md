# Caching & Persistence

Continuing from [Output Interaction](output-interaction.md): beyond
reading a result's elements and ratios, you can also save it, reload it,
and skip re-running a configuration you've already computed.

## Saving and loading a result directly

```python
result.save("my_run.h5")

from crimsons import EnrichmentResult
loaded = EnrichmentResult.load("my_run.h5")
```

This round-trips everything: the time grid, per-element enrichment array,
the `RunConfig` used to produce it, and per-realization detail (sampled
masses, counts, fates, and event history).

## Automatic caching by configuration

`Simulation.run(cache_dir=...)` skips re-running a simulation whose exact
configuration has already been computed:

```python
sim = Simulation(imf=Kroupa2001(), mass_formed=1e6, metallicity=-1.0, seed=42)

result = sim.run(cache_dir="./cache")   # runs, then saves to ./cache/<run_id>.h5
result = sim.run(cache_dir="./cache")   # instant -- loaded from cache, not re-run
```

Force a fresh run (and overwrite the cache entry) with `overwrite=True`:

```python
result = sim.run(cache_dir="./cache", overwrite=True)
```

## How the cache key is computed

Every `Simulation` derives a [`RunConfig`][crimsons.config.RunConfig] --
IMF type and mass bounds, mass formed, metallicity, realization count,
seed, bin count, and channel names -- and hashes it to a 16-character
`run_id` via `config.run_id()`. Two `Simulation`s with matching configs
(down to the seed) hash to the same ID and therefore the same cache file,
`<cache_dir>/<run_id>.h5`.

Because `RunConfig` doesn't capture the *contents* of custom channels or
yield tables (only `channel_names`), swapping in a different yield model
or a hand-built `Channel` while keeping the same channel names will **not**
change the cache key -- use separate `cache_dir`s (or `overwrite=True`) if
you're iterating on channel internals rather than the parameters
`RunConfig` tracks. The same is true for the time grid: the properties of the time grid, such as the extent or the number of sterps, are not stored in the cache key and thus any change in the time grid requires rerunning the setup.  

## What gets stored in the HDF5 file

`EnrichmentResult.save` writes:

- `time`, `elements`, `enrichment` (the ensemble array) at the top level
- `config_json`, an attribute holding the serialized `RunConfig`
- one `realizations/<i>/` group per realization, with that realization's
  sampled `masses`, `counts`, `fates`, and an `events/` subgroup recording
  each contributing channel's event times and star counts (per-event
  *yields* are intentionally **not** persisted here -- only counts -- to
  keep cache files a reasonable size; the cumulative `enrichment` array
  already captures the yield-weighted result)

Next: back to [Basic Simulation](basic-simulation.md), or dig into
[Customize](../customize/index.md) to change what these runs actually
compute.