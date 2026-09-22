# Output Interaction

Everything a `Simulation.run()` call returns lives on `EnrichmentResult`.
Beyond `.mean()`/`.std()` (the ensemble aggregates, shape
`(n_time, n_elements)`), the result also supports indexing by element
symbol or abundance-ratio, so you don't have to look up an element's
index by hand.

## Element access

```python
result = sim.run()

# Element index
idx = result.elements.index("Fe")
mean_fe = result.mean()[:, idx]

# same thing
mean_fe = result["Fe"].mean()
```

`result["Fe"]` returns an `ElementSeries`: a thin wrapper around that
one element's `(n_realizations, n_time)` slice of `result.enrichment`,
with its own `.mean()` and `.std()`:

```python
series = result["C"]
series.values   # (n_realizations, n_time) -- raw cumulative mass, Msun
series.mean()   # (n_time,)
series.std()    # (n_time,)
len(series)     # n_realizations
```

An unknown symbol raises `KeyError` listing what's available:

```python
result["Xx"]
# KeyError: "Unknown element 'Xx'. Available: H, He, Li, Be, B, C, N, O, F, Ne, ..."
```

## Abundance ratios: `[X/Y]`

Bracket syntax computes a standard `[X/Y]` abundance ratio instead of a
raw mass:

```python
mean_cfe = result["[C/Fe]"].mean()
```

This works for **any** two tracked
elements, including `[X/O]`:

```python
result["[C/O]"]    # carbon relative to oxygen
result["[Fe/O]"]    # iron relative to oxygen
result["[Mg/H]"]    # magnesium relative to hydrogen
result["[O/Fe]"]    # the usual alpha-element-vs-iron diagnostic
```

`[X/Y]` is computed as

$$
[\mathrm{X/Y}] = \log_{10}\!\left(\frac{M_X}{M_Y}\right)_{\!\rm ejecta} - \log_{10}\!\left(\frac{M_X}{M_Y}\right)_{\!\odot}
$$

which is algebraically identical to the number-abundance
definition $\log_{10}(N_X/N_Y) - \log_{10}(N_X/N_Y)_\odot$. This apporach removes the atomic-weight table which is still available to convert a log-abundance solar reference into mass units when the solar table is loaded.

!!! note "This is the ejecta ratio, not a diluted gas-phase abundance"
    `[X/Y]` here describes the *cumulative nucleosynthetic yield* -- this
    library's enrichment source term (see [Customize](../customize/index.md#scope))
    -- not an abundance diluted into an ISM gas reservoir, since
    `EnrichmentResult` doesn't track total gas mass. Combining this with a
    gas-mass/inflow model to get a gas-phase `[X/Y]` is left to the caller.

Before any of an element pair has been produced, that time step's ratio is
`nan` rather than `-inf`/`inf` -- `.mean()`/`.std()` already handle this
via `numpy.nanmean`/`numpy.nanstd` internally, so realizations that
haven't produced a given element yet are excluded from the mean at that
time step rather than pulling it towards `nan`:

```python
result["[C/Fe]"].values[:, 0]   # likely all nan -- nothing has exploded yet at t=0
result["[C/Fe]"].mean()[0]      # nan only if *every* realization is still nan there
```



## Number ratios: `X/Y`

Ratios syntax computes a standard `X/Y` number ratio instead of a
raw mass:

```python
mean_cfe = result["C/Fe"].mean()
```

This works for **any** two tracked
elements:

```python
result["C/O"]    # carbon relative to oxygen
result["Fe/O"]    # iron relative to oxygen
result["Mg/H"]    # magnesium relative to hydrogen
result["O/Fe"]    # oxigen relative to iron
```

`X/Y` is computed as

$$
\mathrm{X/Y} = \log_{10}\!\left(\frac{M_X A_Y}{M_Y A_X}\right)
$$

where $M_X$ is the cumulative mass of the chemical element $X$ and $A_X$ is its atomic-weight 

!!! note "This is the ejecta ratio, not a diluted gas-phase abundance"
    `X/Y` here describes the *cumulative nucleosynthetic yield* -- this
    library's enrichment source term (see [Customize](../customize/index.md#scope))
    -- not an abundance diluted into an ISM gas reservoir, since
    `EnrichmentResult` doesn't track total gas mass. Combining this with a
    gas-mass/inflow model to get a gas-phase `X/Y` is left to the caller.

Similarly to `[X/Y]`, before any of an element pair has been produced, that time step's ratio is
`nan` rather than `-inf`/`inf` -- `.mean()`/`.std()` already handle this
via `numpy.nanmean`/`numpy.nanstd` internally, so realizations that
haven't produced a given element yet are excluded from the mean at that
time step rather than pulling it towards `nan`:

```python
result["C/O"].values[:, 0]   # likely all nan -- nothing has exploded yet at t=0
result["C/O"].mean()[0]      # nan only if *every* realization is still nan there
```

## The solar reference

`[X/Y]` uses a bundled default solar abundance table (Asplund et al. 2009,
Table 1 photospheric composition). Swap in your own for the rest of the
process with `set_default_solar_abundances`:

```python
from crimsons import set_default_solar_abundances

# a CSV with an `element` column, and either a log number-abundance
# column (A(X) = log10(N_X/N_H) + 12, auto-detected by name) or a
# mass-fraction column (any column name containing "mass")
set_default_solar_abundances("my_solar_abundances.csv")

result["[C/Fe]"].mean()  # now uses your table
```

or build one in memory without a file, via `SolarAbundances`:

```python
from crimsons import SolarAbundances, set_default_solar_abundances

table = SolarAbundances.from_log_eps({"H": 12.00, "He": 10.93, "C": 8.43, ...})
set_default_solar_abundances(table)
```

Both `load_solar_abundances(path)` and `SolarAbundances` are also
importable directly if you want a table without changing the process-wide
default (e.g. to compare two solar references side by side).

## Checking what's available

```python
"Fe" in result        # True
"[C/O]" in result      # True -- both elements exist
"[C/Xx]" in result     # False
```

## Next

- [Caching & Persistence](caching-results.md) -- `result.save()`/`.load()`,
  and `Simulation.run(cache_dir=...)` to skip re-running an identical
  configuration.
- [Customize](../customize/index.md) -- change what goes *into* a run,
  rather than how you read what comes out.
