# Initial Mass Function (IMF)

The Initial Mass Function $\xi(m) = dN/dm$ describes the relative number of
stars formed per unit mass interval in a single stellar population. Every
IMF in CRIMSONS exposes a probability density `pdf(mass)` and an
`inverse_cdf(u)` used to draw stellar masses (see
[Stochastic Sampling](sampling.md)).

## What it is

### Normalization

CRIMSONS normalizes $\xi(m)$ as a proper *number*-density probability
distribution over the active mass range $[m_{\min}, m_{\max}]$:

$$
\int_{m_{\min}}^{m_{\max}} \xi(m)\, dm = 1
$$

This is a **number** normalization, not a mass normalization -- $\xi(m)$
integrates to one star drawn, not to one solar mass formed. (Sampling
still targets a total *mass*, `mass_formed`, by drawing stars from this
normalized distribution until their summed mass reaches the target -- see
[Stochastic Sampling](sampling.md).)

### Broken power laws

`Salpeter1955`, `Kroupa2001`, and `FlatIMF` are  exact special cases of `BrokenPowerLawIMF`: $\xi(m) \propto m^{-\alpha_k}$ on each segment
$[m_{k-1}, m_k)$, with amplitudes chosen so $\xi(m)$ is continuous across
breakpoints.

- **Salpeter (1955):** a single segment, $\alpha = 2.35$.
- **Kroupa (2001):** three segments (default breakpoints at $0.08$ and
  $0.5\,M_\odot$), with slopes $\alpha = (0.3,\ 1.3,\ 2.3)$.
- **Flat (`FlatIMF`):** a single segment with $\alpha = 0$, i.e.
  $\xi(m) = \text{const}$.

## Customizing

### Choosing a different shape

The public, importable IMFs are `Salpeter1955`, `Kroupa2001`, `FlatIMF`, `Chabrier2003` and
`FunctionalIMF`. `FlatIMF` and the `BrokenPowerLawIMF` base class are also
available from `crimsons.imf.standard` for building your own broken power
law by subclassing with different slopes/breakpoints.

For anything that isn't a power law -- log-normal-plus-power-law forms
(Chabrier-style), tapered (Larson-style) forms, top-heavy Population III
shapes, or a shape lifted straight from a paper -- use `FunctionalIMF`,
which builds a numerical inverse CDF from any shape function you supply.
See [Custom IMF Shapes](../examples/custom-imf.md) for full worked
examples of both routes (a custom `FunctionalIMF` shape, and subclassing
`BrokenPowerLawIMF` directly).

### Changing the mass range

$[m_{\min}, m_{\max}]$ is either:

- **explicit** -- pass `m_min` and/or `m_max` to any IMF constructor; this
  is locked in permanently, even if the IMF is later handed to a
  `Simulation` with a different metallicity, or
- **implicit and metallicity-adaptive** -- if you don't pass them, the
  range is resolved from metallicity via
  `imf.mass_range.default_mass_range` (see
  [Metallicity & Population III](metallicity.md)) every time metallicity
  becomes known or changes.

!!! note "Kroupa's own mass range"
    `Kroupa2001`'s default range follows the same metallicity-adaptive
    policy as every other IMF here, **not** Kroupa's commonly-cited
    brown-dwarf cutoff of $0.01\,M_\odot$. Pass `m_min=0.01` explicitly if
    you want that back.

See the [API reference](../api/imf.md) for the full parameter list of
every IMF class.
