# Selecting Yield Models

!!! warning "Synthetic bundled data"
    `stellar_yields.h5` (SNII/AGB/PISN) and the SNIa placeholder table
    shipped with CRIMSONS are **synthetically generated** for development
    and testing -- not real nucleosynthesis calculations. The model names
    below (`NK`, `VAN`, `HW`, `LC`, `Iwamoto`, ...) exist in the bundled
    file, but swap in your own `h5_path=` before drawing scientific
    conclusions.

Each channel (`SNII`, `AGB`, `PISN`, `SNIa`) picks a `model` from the
bundled `stellar_yields.h5`, defaulting to:

| Channel | Default model |
|---|---|
| `SNII` | `NK` |
| `AGB` | `VAN` |
| `PISN` | `HW` |
| `SNIa` | `Iwamoto` (sub-model `W7`) |

`SNII` defaults to `NK` rather than `LC` deliberately: `LC` needs a
rotation-velocity choice with no principled silent default, so it's opt-in.

## Discovering what's available

```python
from crimsons import list_available_models, describe_available_model

list_available_models("SNII")
# ['HW', 'LC', 'NK', 'NK_HN']

describe_available_model("SNII", "LC")
# {'axes': ['metallicity', 'rotation', 'mass', 'elements'],
#  'yields_shape': (4, 3, 9, 30),
#  'metallicity': [1.42e-05, 1.42e-04, 1.42e-03, 1.42e-02],
#  'rotation': [0.0, 150.0, 300.0],
#  'mass': [13.0, 15.0, 20.0, 25.0, 30.0, 40.0, ...],
#  'elements': ['H', 'He', ..., 'Zn']}
```

Every model has `metallicity`, `mass`, and `elements` axes. Any *other*
axis (like `LC`'s `rotation`, or `HW`'s `energy`/`mixing`) is a discrete
model parameter you resolve via `model_params`.

## Picking a model with extra axes

An axis with only one grid value is picked automatically. An axis with
several values needs an explicit choice in `model_params`, snapped to the
*nearest* tabulated grid point (these are separate model calculations, not
a smoothly interpolable physical dimension):

```python
from crimsons import SNII

# LC needs a rotation choice -- nearest match to 300 among [0, 150, 300]
snii = SNII(model="LC", model_params={"rotation": 300})
```

Omitting a required `model_params` entry raises a `ValueError` that names
the axis and its available values, so you don't have to memorize
`describe_available_model`'s output.

## Stochastic (per-star) model parameters

Instead of a fixed value, `model_params` can hold a **distribution** to
draw from independently per star -- either a `scipy.stats`-style frozen
distribution (anything with `.rvs`) or a plain `callable(rng, n) -> array`:

```python
import scipy.stats as st
from crimsons import SNII

snii = SNII(model="LC", model_params={"rotation": st.norm(loc=200, scale=50)})
```

This returns a table that draws a fresh rotation velocity per star at
sampling time (still snapped to the nearest of `[0, 150, 300]`), rather
than fixing one rotation for the whole run.

## Custom mass windows

Each channel's mass range is also a constructor argument, independent of
the yield model's own tabulated mass grid (queried masses outside the
table's grid are simply clamped to the nearest tabulated mass):

```python
from crimsons import SNII

snii = SNII(mass_min=10.0, mass_max=25.0)
```

## SN Ia: delay-time distribution vs. single burst

`SNIa` doesn't select by mass at all -- it's driven by the population as a
whole:

```python
from crimsons import SNIa

# Default: Maoz+12 power-law DTD, literature rate_per_msun=0.0013
snia = SNIa()

# Mannucci+06 / Matteucci+06 double-Gaussian DTD instead
snia = SNIa(dtd_shape="mannucci")  # rate_per_msun defaults to 0.0025

# Delta-function limit: every eligible SN Ia goes off at one fixed delay
snia = SNIa(mode="single_burst", burst_delay_myr=1000.0, rate_per_msun=0.001)
```

See [Enrichment Channels](../physics/channels.md) for the DTD formulas and
[`SNIa`][crimsons.yields.channels.SNIa] for the full parameter list.

## Using your own yield tables

!!! warning "CHECK THE load CSV function"
    we might have to check the presene of `load_csv_function` and in case describe briefely the expected shape for the yields 

Bypass the bundled HDF5 file entirely by passing `yield_table=` (a
[`YieldTable`][crimsons.yields.base.YieldTable] or
[`StochasticYieldTable`][crimsons.yields.base.StochasticYieldTable]) or
point at a different file with `h5_path=`:

```python
from crimsons.yields.io import load_yield_table_hdf5
from crimsons import AGB

table = load_yield_table_hdf5("/path/to/your_yields.h5", "AGB", "YOUR_MODEL")
agb = AGB(yield_table=table)

# or, keep using the model-name interface against your own file:
agb = AGB(model="YOUR_MODEL", h5_path="/path/to/your_yields.h5")
```

A simple CSV format (`mass, metallicity, <one column per element>`) is also
supported via
[`load_yield_table`][crimsons.yields.io.load_yield_table] -- see
`yields/data/snia_placeholder.csv` in the source for the expected shape.

Next: [Caching & Persistence](caching-results.md).
