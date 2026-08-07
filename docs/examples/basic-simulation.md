# Basic Simulation

A complete run: sample a Kroupa IMF, evolve it through the default
channels, and plot the mean enrichment history with its realization-to-
realization scatter.

```python
import matplotlib.pyplot as plt
import numpy as np

from crimsons import Kroupa2001, Simulation

sim = Simulation(
    imf=Kroupa2001(),
    mass_formed=1e6,       # Msun
    metallicity=-1.0,      # log10(Z / Zsun)
    n_realizations=30,
    seed=0,
    verbose=True,          # progress bar (needs `pip install tqdm`)
)
result = sim.run()

mean = result.mean()       # (n_time, n_elements)
std = result.std()         # (n_time, n_elements)

for name in ["O", "Mg", "Fe"]:
    i = result.elements.index(name)
    plt.plot(result.time, mean[:, i], label=name)
    plt.fill_between(
        result.time, mean[:, i] - std[:, i], mean[:, i] + std[:, i], alpha=0.2
    )

plt.xscale("log")
plt.xlabel("Time [Myr]")
plt.ylabel("Cumulative mass returned [Msun]")
plt.legend()
plt.tight_layout()
plt.savefig("enrichment.png", dpi=150)
```

`matplotlib` isn't a CRIMSONS dependency -- install it separately
(`pip install matplotlib`) if you want to reproduce this.

## Inspecting individual realizations

`result.enrichment[i]` is realization `i`'s own `(n_time, n_elements)`
array. Combined with `result.fates[i]` (which channel each sampled mass
bin went through) and `result.masses[i]` / `result.counts[i]` (the bin
centers and star counts behind that realization), you can drill into a
single draw rather than only the ensemble mean:

```python
i = 0  # first realization
fates, masses, counts = result.fates[i], result.masses[i], result.counts[i]

for channel_name in set(fates):
    n_stars = counts[fates == channel_name].sum()
    print(f"{channel_name:>6}: {n_stars:,.0f} stars")
```

`"none"` is a valid fate: it covers mass bins that didn't fall into any
configured channel (e.g. stars below the AGB mass window, or above PISN's
if you haven't added it).

## Changing the time grid

By default, enrichment is reported on `np.logspace(0, 4, 200)` Myr (1 Myr
to 10 Gyr). Pass your own via `Simulation(..., time_grid=...)`:

```python
sim = Simulation(
    imf=Kroupa2001(),
    mass_formed=1e6,
    metallicity=-1.0,
    time_grid=np.linspace(0, 1000, 500),  # 0-1 Gyr, linear, Myr
)
```

Next: [custom IMF shapes](custom-imf.md), or
[selecting yield models](yield-models.md).
