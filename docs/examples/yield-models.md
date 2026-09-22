# Selecting Yield Models

This page walks through picking and swapping yield models in practice;
see [Customize: Enrichment Channels](../customize/channels.md) for the
physics each channel's defaults represent, and how to change a channel's
mass window or SN Ia's DTD shape instead of its yield model.

Each channel (`SNII`, `AGB`, `PISN`, `SNIa`) picks a `model` from the
bundled `stellar_yields.h5` and defaults to the metallicity-related set of yields:

| Channel | Pop III default |  Pop II default |
|---|---|---|
| `SNII` | `Heager10` | `Limongi18` rotation 0 km/s |
| `AGB` | `Meynet02` rotation 300 km/s| `vandenHoek97` |
| `PISN` | `Heager02` | `Heager02` |
| `SNIa` | `Iwamoto` (model `W70`) | `Iwamoto` (model `W7`) |

The default population is Pop II, if no metallicity is passed.  


## Discovering what's available

```python
from crimsons import list_available_models, describe_available_model

list_available_models("SNII")
# ['Heager10', 'Limongi18', 'Nomoto13', 'Nomoto13_Hypernovae']

describe_available_model("SNII", "Limongi18")
# {'axes': ['metallicity', 'rotation', 'mass', 'elements'],
#  'yields_shape': (4, 3, 9, 30),
#  'metallicity': [1.42e-05, 1.42e-04, 1.42e-03, 1.42e-02],
#  'rotation': [0.0, 150.0, 300.0],
#  'mass': [13.0, 15.0, 20.0, 25.0, 30.0, 40.0, ...],
#  'elements': ['H', 'He', ..., 'Zn']}
```

Every model has `metallicity`, `mass`, and `elements` axes. Any *other*
axis (like `Limongi18`'s `rotation`, or `Heager10`'s `energy`/`mixing`) is a discrete
model parameter you resolve via `model_params` (see [Extra axes](yield-models.md#picking-a-model-with-extra-axes)) 


## Extra axis

Any extra property of the model, i.e. not `metallicity`, `mass`, and `elements`, can be selected with two different approaches. 

### Picking a model with extra axes

You can pick a specific value making an explicity choice in `model_params`. For example you can pass 300 to "rotation" to select the 
rotational velocity of 300 km/S in the `Limongi18` model:

```python
from crimsons import SNII

# Limongi18 needs a rotation choice -- nearest match to 300 among [0, 150, 300]
snii = SNII(model="Limongi18", model_params={"rotation": 300})
```

The choice is snapped to the *nearest* value. In the `Limongi18` example, it means that passing a value of 300, 250 or 500 for the rotational velocity
results in the same value of 300 km/s.

Omitting a required `model_params` entry raises a `ValueError` that names the axis and its available values, so you don't have to memorize
`describe_available_model`'s output.


### Stochastic (per-star) model parameters

Instead of a fixed value, `model_params` can hold a **distribution** to
draw from independently per star -- either a `scipy.stats`-style frozen
distribution (anything with `.rvs`) or a plain `callable(rng, n) -> array`:

```python
import scipy.stats as st
from crimsons import SNII

snii = SNII(model="Limongi18", model_params={"rotation": st.norm(loc=200, scale=50)})
```

This returns a SNII table that draws a fresh rotation velocity per star from a gaussian distribution centered at 200 with
standard deviation of 50. The sampling is still snapped to the nearest of `[0, 150, 300]` but this option enables to 
explore different distributions rather than fixing the rotation. 


```python
import scipy.stats as st
from crimsons import SNII


def heager10_energy_distribution(rng, n):
    # tabulated explosion energies in 10^51 erg
    energies = [3.0, 6.0, 12.0, 15.0, 18.0, 30.0, 50.0, 100.0]
    return rng.choice(energies, size=n)

def heager10_mixing_distribution(rng, n):
    # tabulated internal mixing, f_mix, multiplied by 100 
    mixing = [39.8, 63.1, 100.0, 158.5]
    return rng.choice(mixing, size=n)


snii = SNII(model="Heager10", 
        model_params={"energy": heager10_energy_distribution,
                      "mixing": heager10_mixing_distribution})
```

Similarly to the previous example, this returns a SNII table that draws a random explosion energy and internal mixing per star from an
uniform probability distribution.    


## Custom mass windows

`SNII`, `AGB`, and `PISN` are [Mass Range Channel](../api/yields.md#channels--yield-table-data-structures), which means that each 
star contributes as to the chemical channle if its mass falls within the channel's mass range.
The metallicity-dependent default mass ranges are:

| Channel | Pop III mass range \[$M_\odot$\] |  Pop II mass range \[$M_\odot$\] |
|---|---|---|
| `AGB` | 2 - 8| 2 - 8 |
| `SNII` | 10 - 100 | 13 - 120 |
| `PISN` | 140 - 260 | 140 - 260 |

Each channel's mass range can be overwritten independently of
the yield model's own tabulated mass grid (masses outside the table's grid are simply scaled to the nearest tabulated mass):W

```python
from crimsons import SNII

snii = SNII(mass_min=10.0, mass_max=25.0)
```

Here we define the mass range of the default SNII channel between 10 and 25 $M_\odot$. 

## SN Ia: delay-time distribution vs. single burst

`SNIa` chemical channel is a [population channel](../api/yields.md#channels--yield-table-data-structures) which is driven by the population as a whole. There are two possible enrichment paths for SNIa: delay time distribution "DTD" or single burst. In the former SNIa events follows a distribution that can be specified with the dtd_shape argument, currently "Maoz12" and "Mannucci06" are the two possible choices. In the latter the SNIa enrichment happens after a specific delay time and uses.   

```python
from crimsons import SNIa

# Default: Maoz+12 power-law DTD, literature rate_per_msun=0.0013
snia = SNIa()

# Mannucci+06 / Matteucci+06 double-Gaussian DTD instead
snia = SNIa(dtd_shape="mannucci")  # rate_per_msun defaults to 0.0025

```

This is the example of the DTD case for both the distribution implemented. 

```python
from crimsons import SNIa

# Delta-function limit: every eligible SN Ia goes off at one fixed delay
snia = SNIa(mode="single_burst", burst_delay_myr=1000.0, rate_per_msun=0.001)
```

Here instead we select the single burst mode with a delay of 1 Gyr and a rate of 0.001 $N_{SNIa}/M_\odot$.

See [Customize: Enrichment Channels](../customize/channels.md) for the DTD
formulas and mass-window/mode customization, and
[`SNIa`][crimsons.yields.channels.SNIa] for the full parameter list.

## Using your own yield tables

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

Next: [Output Interaction](output-interaction.md), or see
[Customize: Enrichment Channels](../customize/channels.md) for the
physics behind each channel's defaults.