# Custom IMF Shapes

The built-in options for the IMF are `Salpeter1955`, `Kroupa2001`, `Chabrier2003`, `Flat`.  
`Salpeter1955`, `Kroupa2001` and `Flat` are exact, fast
implementations of [broken power laws](../physics/imf.md#broken-power-laws). `Chabrier2003` is implemented using the `FunctionalIMF`[`FunctionalIMF`][crimsons.imf.functional.FunctionalIMF], which computes the .   For anything else -- log-normal
forms (Chabrier-style), tapered (Larson-style) forms, top-heavy
Population III shapes, or a shape lifted straight from a paper --
use [`FunctionalIMF`][crimsons.imf.functional.FunctionalIMF].

Give it any callable `shape(mass) -> dN/dM`, up to an overall
normalization: it doesn't need to integrate to anything in particular, and
it doesn't need to be vectorized (`FunctionalIMF` falls back to an
elementwise call if a vectorized call doesn't work).

```python
import numpy as np
from crimsons import FunctionalIMF, Simulation


def my_shape(m):
    """Flat below 1 Msun, Salpeter slope above."""
    m = np.atleast_1d(m)
    return np.where(m < 1.0, 1.0, m ** -2.35)


imf = FunctionalIMF(my_shape, m_min=0.1, m_max=100.0)

masses, counts = imf.sample(np.random.default_rng(0), target_mass=1000.0)
```

`n_grid` (default `2000`) controls how many log-spaced points are used to
build the internal CDF grid that `inverse_cdf` interpolates on -- higher is
more accurate but slower to *construct* (sampling speed itself is
unaffected, since it's table lookups either way).

## Passing extra parameters to your shape

`imf_params` (or arbitrary keyword arguments) are forwarded to `shape` on
every call, which is convenient for a parametrized family of shapes:

```python
def chabrier_like(m, mc=0.2, sigma=0.55, alpha=2.3):
    m = np.atleast_1d(m)
    lognormal = np.exp(-((np.log10(m) - np.log10(mc)) ** 2) / (2 * sigma**2))
    power_law = m ** -alpha
    return np.where(m < 1.0, lognormal, power_law * lognormal[np.searchsorted(m, 1.0) - 1])


imf = FunctionalIMF(chabrier_like, m_min=0.05, m_max=120.0, mc=0.25, sigma=0.6)
```

## Using it in a `Simulation`

`FunctionalIMF` (like every IMF) works as a drop-in for `Kroupa2001` /
`Salpeter1955`:

```python
sim = Simulation(imf=imf, mass_formed=1e5, metallicity=0.02, n_realizations=10)
result = sim.run()
```

## Metallicity-adaptive mass ranges

If you *don't* pass `m_min`/`m_max`, every IMF (including `FunctionalIMF`)
falls back to a metallicity-dependent default range the first time its
metallicity is known -- see
[Metallicity & Population III](../physics/metallicity.md). Passing
`m_min`/`m_max` explicitly always overrides that default, permanently, even
if the IMF is later handed to a `Simulation` with a different metallicity.

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
```

Next: [Selecting Yield Models](yield-models.md).
