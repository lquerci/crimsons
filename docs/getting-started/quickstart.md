# Quickstart

This walks through the three objects you touch in almost every CRIMSONS
run: an **IMF**, a **`Simulation`**, and the **`EnrichmentResult`** it
returns.

## 1. Pick an IMF

```python
from crimsons import Kroupa2001

imf = Kroupa2001()
```

`Kroupa2001()`, `Salpeter1955()`, and `FunctionalIMF(...)` (for an
arbitrary custom shape -- see
[Custom IMF Shapes](../examples/custom-imf.md)) are the built-in options.
You don't need to set a mass range here: if you don't pass `m_min`/`m_max`
explicitly, the IMF resolves a metallicity-appropriate default range once
it knows the run's metallicity (see step 2, and
[Metallicity & Population III](../physics/metallicity.md)).

## 2. Configure a `Simulation`

```python
from crimsons import Simulation

sim = Simulation(
    imf=imf,
    mass_formed=1e6,       # Msun of stars formed, per realization
    metallicity=-1.0,      # see convention below
    n_realizations=20,     # independent Monte Carlo realizations
    seed=42,               # reproducible: each realization gets its own child seed
)
```

**Metallicity convention:** pass `metallicity <= 0` for `log10(Z / Zsun)`,
or `metallicity > 0` for absolute `Z`. `metallicity=-1.0` above means
`Z = 0.1 * Zsun`; `metallicity=0.02` would mean absolute `Z = 0.02`. Both
forms are always available on the result via `result.config.metallicity`
(stored as absolute `Z`).

If you don't pass `lifetime_fn` or `channels`, `Simulation` defaults to
`StellarLifetime()` and `default_channels()` (SNII + AGB + SNIa -- add
`PISN()` yourself for Population III / extremely metal-poor runs).

## 3. Run it

```python
result = sim.run()
```

Each of the `n_realizations` realizations samples its own stellar
population from the IMF, evolves it through every configured channel, and
accumulates enrichment onto a shared time grid (`Myr`, log-spaced from 1
to 10,000 by default -- pass `time_grid=` to `Simulation` to override).

Want a progress bar for long runs? Pass `verbose=True` (requires
`pip install tqdm`).

## 4. Read the result

`result` is an [`EnrichmentResult`](../api/core.md):

```python
result.elements          # ["H", "He", "Li", ..., "Zn"] -- 30 tracked elements
result.time               # (n_time,) array, Myr
result.enrichment.shape   # (n_realizations, n_time, n_elements) -- cumulative mass returned, Msun

mean = result.mean()      # (n_time, n_elements) -- mean across realizations
std = result.std()        # (n_time, n_elements) -- scatter across realizations
```

Per-realization detail is also available: `result.masses`, `result.counts`
and `result.fates` give the sampled mass bins, how many stars each bin
represents, and which channel (or `"none"`) each bin's stars went through.

## 5. Save and reload

```python
result.save("my_run.h5")

from crimsons import EnrichmentResult
loaded = EnrichmentResult.load("my_run.h5")
```

Or let `Simulation` manage this for you via `run(cache_dir=...)` -- see
[Caching & Persistence](../examples/caching-results.md).

## Full example

```python
import numpy as np
from crimsons import Kroupa2001, Simulation

imf = Kroupa2001()
sim = Simulation(
    imf=imf,
    mass_formed=1e6,
    metallicity=-1.0,
    n_realizations=20,
    seed=42,
)
result = sim.run()

mean = result.mean()
i_fe = result.elements.index("Fe")
print("Final mean Fe returned to the ISM:", mean[-1, i_fe], "Msun")
```

Next: [worked examples](../examples/basic-simulation.md), or the
[physics](../physics/index.md) behind each piece.
