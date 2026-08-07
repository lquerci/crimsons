---
title: CRIMSONS
description: Chemical Evolution Yield & IMF Stochastic Engine
---

# CRIMSONS

Chemical Evolution Yield & IMF Stochastic Engine

CRIMSONS is a Python engine for stochastic chemical-enrichment simulations of
stellar populations: sample an initial mass function, evolve it through
core-collapse supernovae, AGB winds, Type Ia supernovae and (optionally)
pair-instability supernovae, and get back the time-resolved abundances of 30
elements released to the interstellar medium.

!!! warning "Bundled yield tables are placeholders"
    The `stellar_yields.h5` file shipped with this package is
    **synthetically generated**, not real published nucleosynthesis data --
    see [Selecting Yield Models](examples/yield-models.md) for how to point
    CRIMSONS at your own tables before using it for science.

## What it does

- **Sample an IMF stochastically.** Draw a stellar population up to a target
  formed mass, binned in log-mass rather than star-by-star, so a single
  realization stays fast even at $10^6\,M_\odot$ and above.
- **Evolve it through multiple enrichment channels.** SNII and AGB are
  triggered by initial mass; SNIa is driven by a delay-time distribution (or
  a single-burst approximation) over the population as a whole; PISN covers
  the (typically Population III) pair-instability window.
- **Track metallicity- and population-dependent physics.** Stellar
  lifetimes, IMF mass ranges, and yield-table lookups all branch between
  Population III and Population II/I regimes at a configurable metallicity
  threshold.
- **Run an ensemble, not a single draw.** `Simulation` runs $N$ independent,
  reproducibly-seeded realizations and reports the mean and scatter across
  them -- the "stochastic" in the name.
- **Cache and persist results.** Every run's configuration hashes to a
  stable ID; results round-trip to HDF5 so repeated runs of the same
  configuration are free.

## A 30-second example

```python
from crimsons import Kroupa2001, Simulation

imf = Kroupa2001()

sim = Simulation(
    imf=imf,
    mass_formed=1e6,      # Msun of stars formed in this population
    metallicity=-1.0,     # log10(Z / Zsun) = -1
    n_realizations=20,
    seed=42,
)

result = sim.run()
mean_enrichment = result.mean()  # shape (n_time, n_elements)
```

Head to [Getting Started](getting-started/installation.md) for the full
walkthrough, or straight to [Examples](examples/basic-simulation.md) for more
worked scenarios.

## Where to go next

<div class="grid cards" markdown>

- :material-download:{ .lg .middle } **Getting Started**

    ---

    Install CRIMSONS and run your first simulation in a few minutes.

    [:octicons-arrow-right-24: Installation](getting-started/installation.md)

- :material-flask-outline:{ .lg .middle } **Examples**

    ---

    Worked examples: custom IMF shapes, yield-model selection, caching.

    [:octicons-arrow-right-24: Basic simulation](examples/basic-simulation.md)

- :material-atom:{ .lg .middle } **Physics**

    ---

    The IMF, lifetime, channel, and metallicity models behind the code.

    [:octicons-arrow-right-24: Physics overview](physics/index.md)

- :material-api:{ .lg .middle } **API Reference**

    ---

    Full, auto-generated reference for every public class and function.

    [:octicons-arrow-right-24: API overview](api/index.md)

</div>

## Citing

If CRIMSONS contributes to a publication, please cite the repository (add
your Zenodo DOI / paper reference here once available) and cite the
literature sources for whichever yield tables and IMF/lifetime
prescriptions you actually used -- see [Physics](physics/index.md) and
[Selecting Yield Models](examples/yield-models.md) for the references baked
into each default.
