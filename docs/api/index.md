# API Reference

This section is generated directly from docstrings in `src/crimsons` by
[mkdocstrings](https://mkdocstrings.github.io/) -- every page's "Source
code" toggle shows the real implementation alongside its documentation,
so it never drifts out of sync with what's actually shipped.

## Package layout

| Page | Module(s) | Covers |
|---|---|---|
| [Simulation](simulation.md) | `crimsons.simulation.ensemble` | The `Simulation` class -- ties everything below together into an ensemble run. |
| [Initial Mass Functions](imf.md) | `crimsons.imf.*` | `IMF` base class, `Salpeter1955`/`Kroupa2001`/`FlatIMF`/`BrokenPowerLawIMF`, `FunctionalIMF`, and the metallicity-adaptive mass-range policy. |
| [Stellar Lifetimes](stars.md) | `crimsons.stars.lifetimes` | `LifetimeFunction` base class and the bundled `StellarLifetime` model. |
| [Yields & Channels](yields.md) | `crimsons.yields.*` | `Channel`/`MassRangeChannel`/`PopulationChannel`, `YieldTable`/`StochasticYieldTable`, the bundled `SNII`/`AGB`/`PISN`/`SNIa` channels, and yield-table I/O. |
| [Enrichment Engine](enrichment.md) | `crimsons.enrichment.engine` | The per-realization sampling/enrichment loop `Simulation` drives. |
| [Results & Config](core.md) | `crimsons.results`, `crimsons.config`, `crimsons.chemistry` | `EnrichmentResult`, `RunConfig`, and the tracked-element/metallicity constants. |
| [I/O & Caching](io.md) | `crimsons.io.*` | HDF5 read/write for `EnrichmentResult`, and the config-hash-based run cache. |

## Top-level imports

Everything most users need is re-exported from the `crimsons` package
itself:

```python
from crimsons import (
    AGB, PISN, SNII, SNIa,               # channels
    ELEMENTS, ZSUN,                       # constants
    EnrichmentResult,
    FunctionalIMF, Kroupa2001, Salpeter1955,  # IMFs
    Simulation,
    StellarLifetime,
    default_channels, describe_available_model, list_available_models,
)
```

`FlatIMF`, `BrokenPowerLawIMF`, and the individual `Channel`/`YieldTable`
base classes aren't re-exported at the top level -- import them from their
own submodule (e.g. `from crimsons.imf.standard import BrokenPowerLawIMF`)
when subclassing or introspecting directly.
