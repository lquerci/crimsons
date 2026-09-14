# Quickstart

This walks through the three objects you touch in almost every CRIMSONS
run: an **IMF**, a **`Simulation`**, and the **`EnrichmentResult`** it
returns.

## 1. Pick an IMF

```python
from crimsons import Kroupa2001

imf = Kroupa2001()
```

`Kroupa2001()`, `Salpeter1955()`, and `Chabrier2003()` are the built-in options. For the full 
list of IMFs see [Initial Mass Function](../physics/imf.md). It is also possible to pass an 
arbitrary custom IMF shape with `FunctionalIMF(...)`,  see
[Custom IMF Shapes](../examples/custom-imf.md). 

If the IMF is not  specified, CRIMSONS defaults to `Chabrier2003()`, with `m_ch = 0.35` or  `10` based on metallicity.
Similarly, if no `m_min`/`m_max` is passed explicitly, the IMF resolves a metallicity-appropriate default range once
it knows the run's metallicity (see step 2, and
[Metallicity & Population III](../physics/metallicity.md)).

## 2. Configure a `Simulation`

```python
from crimsons import Simulation

sim = Simulation(
    imf=imf,
    mass_formed=1e6,       # Msun of stars formed, per realization
    metallicity=0.00142,   # absolute metallicity 
    n_realizations=20,     # independent Monte Carlo realizations
    seed=42,               # reproducible: each realization gets its own child seed
)
```

**Metallicity convention:** the passed metallicity, $Z$, is expected to be absolute. The assumed solar metallicity in CRIMSONS is $Z_\odot = 0.0142$ (Asplund et al., 2009), therefore the value in the example corresponds to $\log(Z/Z_\odot) = -1$. You can convert $\log(Z/Z_\odot) = -1$ to $Z$ using the `z_from_logz()` function. 

```python
from crimsons import z_from_logz

absolute_metallicity = z_from_logz(-1)

```

The metallicity distinguishes between Pop III and Pop II single stellar populations with the threshold value set to $Z_{crit} = 10^{-4.5} Z_\odot$. Once the metallicity is specified, all parameters defaults to the fiducial parameters of Rossi et al., 2026. Additionally, if you don't pass `lifetime_fn`, CRIMSONS default is `StellarLifetime()` which uses Raitieri+97 and Schaerer 2002 for stellar lifetimes of PopII and PopIII, respectively.  Similarly, you can specify the chemical enrichment channels passing `channels` as arguments, otherwise CRIMSONS uses `default_channels()` which are: Supernovae type II (SNII), type Ia (SNIa), and Asymptiotic Giant Branch (AGB), for Pop II and PopIII, with the addition of Pair-instability SNe (PISN) for the latter.

## 3. Run it

```python
result = sim.run()
```

Each of the `n_realizations` realizations samples its own stellar
population from the IMF, evolves it through every configured channel, and
accumulates enrichment onto a shared time grid (`Myr`, log-spaced from 1
to 10,000 by default -- pass `time_grid=` to `Simulation` to override).

Want a progress bar for long runs? Pass `verbose=True` as `Simulation` argument.

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

Or let `Setup` manage this for you via `run(cache_dir=...)` -- see
[Caching & Persistence](../examples/caching-results.md).

## Full example

The following is an example of computin the mean iron enrichemnt from a Pop II single stellar population considering the enrichment from SNII and AGB only.  

```python
import numpy as np
from crimsons import Kroupa2001, Setup, SNII, AGB

imf = Kroupa2001()
sim = Setup(
    imf=imf,
    mass_formed=1e4,
    metallicity=0.0142,
    n_realizations=20,
    seed=42,
    channels=[SNII(), AGB()]
)
result = sim.run()

mean = result.mean()
i_fe = result.elements.index("Fe")
print("Final mean Fe returned to the ISM:", mean[-1, i_fe], "Msun")
```

Next: [worked examples](../examples/basic-setup.md), or the
[physics](../physics/index.md) behind each piece.
