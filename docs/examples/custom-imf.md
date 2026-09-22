# Custom IMF Shapes

The built-in options for the IMF are `Salpeter1955`, `Kroupa2001`,
`FlatIMF`, and `Chabrier2003`. `Salpeter1955`, `Kroupa2001`, and `FlatIMF`
are exact, fast implementations of
[broken power laws](../customize/imf.md#broken-power-laws); `Chabrier2003`
-- log-normal below $1\,M_\odot$, a power law above -- isn't a broken
power law, so it's implemented on top of
[`FunctionalIMF`][crimsons.imf.functional.FunctionalIMF] instead. For any
other shape -- tapered (Larson-style) forms, top-heavy Population III
shapes, or something lifted straight from a paper -- use `FunctionalIMF`
directly, the same way, as shown below.

Give it any callable `shape(mass) -> dN/dM`, up to an overall
normalization: it doesn't need to integrate to anything in particular, and
it doesn't need to be vectorized (`FunctionalIMF` falls back to an
elementwise call if a vectorized call doesn't work).

```python
import numpy as np
from crimsons import FunctionalIMF


def my_shape(m):
    """Flat below 1 Msun, Salpeter slope above."""
    m = np.atleast_1d(m)
    return np.where(m < 1.0, 1.0, m ** -2.35)


imf = FunctionalIMF(shape=my_shape, m_min=0.1, m_max=100.0, n_grid=2000)

masses, counts = imf.sample(np.random.default_rng(0), target_mass=1000.0)
```

`m_min` is the minimum stellar mass (in $M_\odot$) while `m_max` is the maximum stellar mass (in $M_\odot$). `n_grid` (default `2000`) controls how many log-spaced points are used to
build the internal CDF grid that `inverse_cdf` interpolates on -- higher is
more accurate but slower to *construct* (sampling speed itself is
unaffected, since it's table lookups either way).

Any IMF has a `sample` method that can be used to display the sampled distribution

```python
# performs the IMF sampling  
masses, counts = imf.sample(np.random.default_rng(42), target_mass=1000.0)

# total numeber of stars
total_stars = np.sum(counts)

# mass bins for visualization 
macro_edges = np.geomspace(imf.m_min, imf.m_max, 36)
macro_centers = np.sqrt(macro_edges[:-1] * macro_edges[1:])
macro_dm = np.diff(macro_edges)

# sampled distribution
counts, _ = np.histogram(masses, bins=macro_edges, weights=counts)
pdf_macro = counts / (total_stars * macro_dm)
plt.step(macro_centers, pdf_macro, where="mid", lw=2.0)

# assumed IMF
mass_grid = np.logspace(np.log10(imf.m_min), np.log10(imf.m_max), 500)
pdf_values = imf.pdf(mass_grid)
plt.plot(mass_grid, pdf_values, linewidth=2)

plt.xscale("log")
plt.yscale("log")
plt.show()
```

## Passing extra parameters to your shape

If your IMF function requires extra parameters, you can pass them either as arguments or as part of the  `imf_params` dictionary. Both are forwarded to `shape` on every call, which is convenient for a parametrized family of shapes:

```python
def chabrier_like(m, mc=0.2, sigma=0.55, alpha=2.3):
    m = np.atleast_1d(m)
    lognormal = np.exp(-((np.log10(m) - np.log10(mc)) ** 2) / (2 * sigma**2))
    power_law = m ** -alpha
    return np.where(m < 1.0, lognormal, power_law * lognormal[np.searchsorted(m, 1.0) - 1])

# extra parameters passed as function arguments
imf = FunctionalIMF(shape=chabrier_like, m_min=0.05, m_max=120.0, mc=0.25, sigma=0.6)

# extra parameters passed in imf_params
imf = FunctionalIMF(shape=chabrier_like, m_min=0.05, m_max=120.0, imf_params={'mc':0.25, 'sigma':0.6})

```

## Using it in a `Simulation`

`FunctionalIMF` (like every IMF) works as imf agrument:

```python
sim = Simulation(imf=imf, mass_formed=1e5, metallicity=0.02, n_realizations=10)
result = sim.run()
```

## Metallicity-adaptive mass ranges

If you *don't* pass `m_min`/`m_max`, every IMF (including `FunctionalIMF`)
falls back to a metallicity-dependent default range instead -- see
[Metallicity & Population III](../customize/metallicity.md) for exactly
what that range is and how to opt out of it.

## Building your own exact broken power law

If your custom shape actually *is* a broken power law, subclassing
`BrokenPowerLawIMF` (as `Salpeter1955` and `Kroupa2001` do) is exact and a
little faster than `FunctionalIMF`, since it skips the numerical CDF grid:

```python
from crimsons.imf.standard import BrokenPowerLawIMF

class MyThreeSlopeIMF(BrokenPowerLawIMF):
    def __init__(self, m_min=None, m_max=None, metallicity=None, mass_range_fn=None):
        super().__init__(
            slopes=[0.3, 1.3, 2.7],
            interior_breakpoints=(0.1, 0.6),
            m_min=m_min, m_max=m_max,
            metallicity=metallicity, mass_range_fn=mass_range_fn,
        )

imf = MyThreeSlopeIMF()
```

In this example the broken power law has three *negative* slopes of 0.3, 1.3, and 2.7 and the intervals for the three slopes are [m_min,0.1), [0.1,0.6), and [0.6,M_max), respectively.  

Next: [Selecting Yield Models](yield-models.md), or see
[Customize: Initial Mass Function](../customize/imf.md) for the
normalization and mass-range details behind all of the above.