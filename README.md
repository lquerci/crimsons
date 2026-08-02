# crimsons

Stochastic chemical enrichment of stellar populations from a sampled IMF,
tracking SN II, SN Ia, and AGB contributions across 30 elements.

> **Placeholder name.** `crimsons` is very likely taken on PyPI — pick a
> real name and rename the `src/crimsons` directory plus every `crimsons`
> import before you publish. Same goes for the placeholder lifetimes and
> yield tables described below: this repo runs end to end out of the box,
> but the *science* in it is illustrative, not validated. See "What's a
> placeholder" below before using this for real results.

## Install (development)

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
pytest
```

## Quickstart

```python
from crimsons import Simulation, Salpeter1955, PowerLawLifetime, default_channels

sim = Simulation(
    imf=Salpeter1955(),
    lifetime_fn=PowerLawLifetime(),
    channels=default_channels(),       # SNII + AGB + SNIa, bundled placeholder yields
    mass_formed=1e5,                   # Msun formed in this stellar population
    metallicity=0.02,
    n_realizations=20,
    seed=42,
)

result = sim.run(cache_dir="./cache")  # reruns with the same config load from cache instead of recomputing

mean_enrichment = result.mean()        # (n_time, n_elements)
result.save("my_run.h5")

# later, in another script:
from crimsons import EnrichmentResult
result = EnrichmentResult.load("my_run.h5")
```

See `examples/quickstart.py` for a runnable version.

## What this library computes (and what it doesn't)

Given an IMF, a target stellar mass, and a metallicity, this library
samples a stochastic stellar population, evolves each star through its
lifetime, and reports the *cumulative mass of each element returned to
the interstellar medium as a function of time* by SN II, SN Ia, and AGB
channels. That's the enrichment **source term**.

It does **not** run a full one-zone galactic chemical evolution model
(star formation history, gas inflow/outflow, dilution, multiple stellar
generations). Combining this source term with those ingredients is left
to the caller — this library's job is to get the stochastic
IMF-sampling-and-yields part right and reusable.

### Populations are stored as log-mass bins, not individual stars

For large `mass_formed`, the actual number of stars can reach into the
hundreds of millions to billions — holding that many individual masses
in memory isn't viable. `IMF.sample_binned` (used internally by
`Simulation`) draws stars exactly as before, chunk by chunk, but
immediately folds each chunk into a histogram over `n_bins` (default
1000) log-spaced mass bins and discards the individual values, so peak
memory stays roughly constant regardless of `mass_formed`. Each
`EnrichmentResult.masses[i]`/`.counts[i]` pair gives the mean mass and
star count of every populated bin in realization `i`, rather than one
entry per star.

This is exact for a mass-range channel (SN II, AGB): a bin is
unambiguously either in the channel's range or not, so scaling its yield
by its star count reproduces sampling every star individually. It is
**not** currently exact for SN Ia, whose `contributes`/`delay_time` draw
once per input element — applied to a bin, that means the bin's entire
star count either all becomes SN Ia progenitors exploding at one shared
delay time, or none do, rather than each star getting an independent
binary/DTD draw. This is a known, temporary approximation while the SN
Ia channel's approach gets reworked; it doesn't affect SN II/AGB.

## Architecture

```
src/crimsons/
├── imf/            pluggable IMFs (Salpeter1955, Kroupa2001, FunctionalIMF, ...)
├── stars/          pluggable stellar lifetime functions
├── yields/         pluggable enrichment channels + yield table interpolation
│   └── data/       stellar_yields.h5 (SNII/AGB/PISN, multi-model) +
│                   snia_placeholder.csv + manifest.yaml
├── enrichment/      the vectorized per-realization engine
├── simulation/      Simulation: loops over N realizations, seeds, aggregates
├── io/               HDF5 save/load, config-hash-based caching
├── results.py       EnrichmentResult container
├── config.py         RunConfig: the hashable, JSON-serializable run spec
└── elements.py       the 30 tracked element names
```

The guiding principle: **physics plugins never touch the engine, and the
engine never touches I/O.** If you're adding a new IMF, lifetime formula,
or yield set, you're editing exactly one small file and nothing else.

## The yields module: SNII, AGB, PISN

These three are all `MassRangeChannel`s (a star's initial mass alone
decides whether it goes through the channel), and all load from a single
bundled `stellar_yields.h5`. Each channel group in that file holds one or
more named **models** (literature sources, e.g. `LC`, `NK`, `HW`), and
each model is an N-dimensional grid: always `metallicity` × `mass` ×
`elements`, plus whatever extra physical parameters that source has (SN
II's `LC` has `rotation`; `HW` has `energy` and `mixing`). The schema is
documented in `yields/io.py`.

```python
from crimsons import SNII, AGB, PISN
from crimsons.yields.channels import list_available_models, describe_available_model

list_available_models("SNII")                    # ['HW', 'LC', 'NK', 'NK_HN']
describe_available_model("SNII", "LC")            # axes, grid values, shape

SNII()                                             # default model (NK) -- no extra params needed
SNII(model="LC", model_params={"rotation": 300})   # LC needs rotation specified: 3 values, no safe default
AGB(model="MM")
PISN()                                             # not in default_channels() -- add explicitly for
                                                    # Population III / extremely metal-poor runs
```

Extra parameter axes are resolved to the *nearest* available grid value
at load time, not interpolated across — they're discrete model variants
(different literature calculations), not a smooth dimension the engine
has any way to vary per star. An axis with only one grid value doesn't
need `model_params` at all; one with several raises a clear error
telling you what to pick if you don't specify it. This is also why
`SNII()`'s default model is `NK` rather than `LC`: `LC` needs a rotation
choice with no principled silent default, so it stays opt-in.

`SNIa` is unaffected by any of this — it's still the CSV-based
placeholder in `snia_placeholder.csv`, since its whole approach is a
known TODO to rework separately (see the binned-sampling caveat above).

## What's a placeholder (replace before doing real science)

- **`PowerLawLifetime`** (`stars/lifetimes.py`) — a toy `tau = t0 * M^-slope`
  fit, metallicity-independent. Replace with your Fortran code's actual
  formula by subclassing `LifetimeFunction`.
- **`stellar_yields.h5`** (`yields/data/`) — I don't have your real data,
  so this is a synthetic file generated by
  `scripts/generate_stellar_yields_h5.py`, matching the schema in
  `yields/io.py` but with made-up numbers. Regenerate it with your real
  tables converted into the same schema (one HDF5 group per model, an
  `axes` attribute naming the `yields` dataset's dimensions, one 1-D
  dataset per axis) — `generate_stellar_yields_h5.py`'s `write_model`
  helper is a template for the writing side. If your real data's layout
  differs from this schema, the loader in `yields/io.py` is the one
  place to adjust.
- **`snia_placeholder.csv`** — synthetic, from `generate_placeholder_data.py`.
- **`SNIa`'s progenitor mass range, binary fraction, and DTD shape**
  (`yields/channels.py`) — illustrative numbers. If your Fortran code
  handles SN Ia differently (e.g. a fixed delay rather than a DTD, or a
  different progenitor scenario), adjust `contributes`/`delay_time`
  accordingly; the `Channel` interface doesn't assume any particular
  scheme.

## Extending

- **New IMF from an arbitrary function**: `FunctionalIMF(shape, m_min=..., m_max=...)`
  takes any `shape(mass) -> dN/dM` (vectorized or not) and builds a
  normalized pdf and a numerical inverse CDF automatically — use this for
  anything that isn't a pure power law (log-normal + power-law tails,
  tapered forms, top-heavy Population III shapes, ...).
- **New broken-power-law IMF**: subclass `crimsons.imf.base.IMF` directly,
  or `BrokenPowerLawIMF` for a piecewise power law with a closed-form
  (exact, slightly faster) inverse CDF — see `Salpeter1955`/`Kroupa2001`
  for examples.
- **Metallicity-adaptive default mass range**: every IMF resolves its
  default `m_min`/`m_max` from metallicity via `imf.mass_range.default_mass_range`
  (standard `0.1–100 Msun`, widened to `0.8–1000 Msun` below `Z = 3e-5`)
  *unless* you pass `m_min`/`m_max` explicitly, which always wins. A
  `Simulation` resolves this automatically from its own `metallicity`, so
  you only specify metallicity once. Pass `mass_range_fn=` to any IMF to
  use your own policy instead of the built-in threshold.
- **New lifetime formula**: subclass `crimsons.stars.lifetimes.LifetimeFunction`.
- **New channel**: subclass `crimsons.yields.base.Channel` (or
  `MassRangeChannel` for anything triggered purely by initial mass) and
  give it a `YieldTable`.
- **New yield model**: add a group to `stellar_yields.h5` following the
  schema above — `SNII(model="your_new_model")` picks it up with no code
  changes. Or bypass the HDF5 file entirely: pass `yield_table=` to
  `SNII`/`AGB`/`PISN`/`SNIa` with your own `YieldTable` built from arrays
  or loaded via `yields.io.load_yield_table` (CSV).

## Testing against the original Fortran code

The most valuable test you can add is a regression test: run the Fortran
code with a fixed seed and a small stellar population, save its output,
and add a test in `tests/integration/` asserting the Python port
reproduces it within tolerance. `tests/integration/test_simulation_roundtrip.py`
is a template for the plumbing (running a `Simulation`, saving/loading);
swap in the real reference numbers once you have the physics ported.

## Roadmap ideas

- Parallelize realizations across cores (`joblib`) once the physics is
  validated — realizations are embarrassingly parallel by construction.
- Rework `SNIa` to draw individual binary/DTD outcomes per bin count
  rather than once per bin (see the caveat above) — this is the next
  planned change, not yet implemented.
- Optional `[plot]` extra with a few standard plots (enrichment vs time,
  IMF sampling histogram, fate-by-mass-bin).