# Metallicity & Population III

A single metallicity value, and a single threshold, reshape three
independent parts of a `Simulation`: the IMF's default mass range, the
stellar lifetime formula, and which yield-table data a channel is allowed
to use.

## What it is

### Specifying metallicity

The metallicity convention is that metallicity, $Z$, is positive and is interpreted as **absolute** metallicity. This is true for all metallicity-dependent functions and call. In case there is the necessity to convert the absolute metallicity into $\log_{10}(Z/Z_\odot)$, the solar metallicity is $Z_\odot = 0.0142$ (Asplund et al, 2009) throughout (`crimsons.chemistry.ZSUN`). 
### The Population III threshold

```python
POPIII_THRESHOLD = 3e-7  # absolute Z, roughly log10(Z/Zsun) < -4.5
```

Any absolute metallicity at or below this value is treated as Population
III; above it, Population II/I. This one constant
(`crimsons.chemistry.POPIII_THRESHOLD`) drives three separate branches:

**1. IMF default mass range** (`imf.mass_range.default_mass_range`, only
when `m_min`/`m_max` weren't given explicitly):

| Regime | Default range |
|---|---|
| $Z > 3\times10^{-7}$ | $0.1$–$100\,M_\odot$ |
| $Z \le 3\times10^{-7}$ | $0.8$–$1000\,M_\odot$ |

**2. Stellar lifetime formula** (see [Stellar Lifetimes](lifetimes.md)):
Population II/I uses the Raiteri et al. (1996) fit for all masses;
Population III switches to a Raiteri+Schaerer combination, with a
different fit above $5\,M_\odot$.

**3. Yield table regime routing** (`YieldTable`, see
[Selecting Yield Models](../examples/yield-models.md)): a table's
tabulated metallicities are split at this same threshold into two
**independent** regimes that are never interpolated across each other.
Population III and Population II/I are treated as genuinely different
kinds of stars in this library, not points on a continuum -- querying a
regime a table has no data for raises a `ValueError` rather than silently
extrapolating one population's yields into the other.

### Two different kinds of "outside the tabulated range"

Within the *correct* regime, a query landing outside that regime's own
tabulated metallicity range doesn't error -- it **warns** and clamps to
the nearest tabulated bin (`MetallicityOutOfRangeWarning`). This is
deliberately different from a query landing in the *wrong* regime
entirely (a hard `ValueError`): a table's Z-coverage being sparser than
you assumed is worth a warning; asking a Population II/I table for
Population III yields is a modeling error worth stopping on. Mass, by
contrast, is silently clamped either way -- a queried mass slightly
outside a table's grid is considered routine.

## Customizing

### Overriding the Population III threshold per table

`YieldTable.popiii_threshold` overrides the package-wide
`POPIII_THRESHOLD` for a single table, if you need a different Pop
III/II boundary for a specific yield set without changing it everywhere
else:

```python
from crimsons.yields.io import load_yield_table_hdf5

table = load_yield_table_hdf5("/path/to/your_yields.h5", "SNII", "YOUR_MODEL")
table.popiii_threshold = 1e-6
```

### Bypassing the metallicity-adaptive IMF range

Pass `m_min`/`m_max` explicitly to any IMF constructor to opt out of the
metallicity-adaptive default range described above entirely -- see
[Initial Mass Function: changing the mass range](imf.md#changing-the-mass-range).

### Binding metallicity to an IMF standalone

`Simulation.__init__` automatically calls `imf.bind_metallicity(self.metallicity)`
before running, so in normal use you only specify metallicity once, on the
`Simulation`. `bind_metallicity` is a no-op if the IMF was given explicit
`m_min`/`m_max` at construction (explicit bounds always win); otherwise it
re-resolves the default range for the new metallicity. Call it yourself if
you're using an IMF standalone (outside a `Simulation`) and want its
metallicity-adaptive default to reflect a metallicity you didn't know at
construction time:

```python
imf = Kroupa2001()          # mass range not yet resolved
imf.bind_metallicity(0.001) # now resolved for this metallicity
```
