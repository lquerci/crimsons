# Enrichment Channels

A **channel** decides which stars enrich the ISM, when their enrichment
arrives, and how much of each element they release. CRIMSONS has two kinds:

- **`MassRangeChannel`** -- triggered purely by a star's initial mass;
  enrichment arrives at the end of that star's life (`SNII`, `AGB`,
  `PISN`).
- **`PopulationChannel`** -- driven by the population as a whole rather
  than any single star's mass (`SNIa`, via a delay-time distribution).

`default_channels()` returns `SNII() + AGB() + SNIa()`. Add `PISN()`
yourself for Population III / extremely metal-poor runs.

## What it is

### Mass-triggered channels

| Channel | Default mass window | Physical picture |
|---|---|---|
| `AGB` | $2$–$8\,M_\odot$ | Low/intermediate-mass stars enriching via stellar winds; primary $s$-process and light elements (C, N, F). |
| `SNII` | $8$–$40\,M_\odot$ | Core-collapse supernovae; hydrostatic burning up to iron-group elements plus explosive nucleosynthesis. |
| `PISN` | $140$–$260\,M_\odot$ | Pair-instability supernovae: complete disruption, no compact remnant, of very massive, metal-free/extremely metal-poor stars. Effectively Population III only. |

A star contributes to a channel if its mass falls in `[mass_min,
mass_max]`; its delay time is simply its stellar lifetime (see
[Stellar Lifetimes](lifetimes.md)); and its yield comes from that
channel's yield table (see [Selecting Yield Models](../examples/yield-models.md)),
evaluated at its mass and the population's metallicity.

### SN Ia: a delay-time distribution over the population

SN Ia doesn't select individual progenitor stars by mass -- it treats the
whole population's eventual SN Ia rate as a population-level statistic.
Explosions follow a delay-time distribution (DTD) which is normalized to
integrate to 1 over the run's `time_grid` and is scaled by `rate_per_msun`
(SN Ia per $M_\odot$ of stars formed) to give an absolute expected number
of explosions per time bin:

$$
\dot N_{\rm Ia}(t) \propto \Psi(\tau)
$$

Two shapes are available for $\Psi(\tau)$ ($\tau$ in Myr):

=== "Maoz+12 (default)"

    A power law, with support bounded below by `min_time`:

    $$
    \Psi(\tau) \propto \left(\frac{\tau}{\tau_{\min}}\right)^{-1.12}, \quad \tau > \tau_{\min}
    $$

    Default literature rate: `rate_per_msun = 0.0013`.

=== "Mannucci+06 / Matteucci+06"

    A double-Gaussian-in-log-time shape, split at $t_0 = 10^{7.93}\,\mathrm{yr}$
    into a "prompt" and a "delayed" component:

    $$
    \log_{10}\Psi(\tau) =
    \begin{cases}
    1.4 - 50\,(\log_{10}t + 7.7)^2 & t \le t_0 \\
    -0.8 - 0.9\,(\log_{10}t - 8.7)^2 & t > t_0
    \end{cases}
    \quad (t = \tau\ \text{in yr})
    $$

    Default literature rate: `rate_per_msun = 0.0025`.

Either shape's support is bounded by the progenitor mass range's
(default $3$–$8\,M_\odot$) lifetimes: the most massive progenitor sets
$\tau_{\min}$, the least massive sets $\tau_{\max}$, both evaluated through
whichever `lifetime_fn` the `Simulation` uses. The (generally fractional)
expected explosion count per time bin is converted to an integer count via
**deterministic remainder carry-over** -- each bin's leftover fraction
rolls into the next bin's expected count -- rather than a random (e.g.
Poisson) draw.

!!! note "SN Ia is deterministic across realizations, by default"
    Because neither the DTD shape, its normalization, nor the
    integerization step involves the realization's random number
    generator, every realization in an ensemble gets an **identical**
    SN Ia contribution (same `mass_formed`, `metallicity`, and
    `time_grid` in, same result out) -- unless you configure a yield
    table with a stochastic model parameter (see
    [Selecting Yield Models](../examples/yield-models.md)). The
    realization-to-realization *scatter* you see in an
    `EnrichmentResult` comes from the stochastic IMF sampling behind
    `SNII`/`AGB`/`PISN`, not from `SNIa`.

## Customizing

### Changing a channel's mass window

The mass windows in the table above are constructor defaults, not
physical constants -- pass `mass_min=`/`mass_max=` to any of the three
mass-triggered channels to change the window (e.g. for a different
progenitor model or literature convention):

```python
from crimsons import SNII

snii = SNII(mass_min=10.0, mass_max=25.0)
```

### Choosing SN Ia's shape and mode

`SNIa` supports `dtd_shape="maoz"` (default) or `dtd_shape="mannucci"`
for the two formulas above, plus a `mode="single_burst"` alternative --
the discretized limit of a delta-function DTD, where every eligible SN Ia
explodes at one fixed delay (`burst_delay_myr`) after formation instead of
following a distribution:

```python
from crimsons import SNIa

snia = SNIa(dtd_shape="mannucci")
snia = SNIa(mode="single_burst", burst_delay_myr=1000.0, rate_per_msun=0.001)
```

See [Selecting Yield Models](../examples/yield-models.md#sn-ia-delay-time-distribution-vs-single-burst)
for the additional parameter walkthrough.

### Choosing a yield model

Each channel also takes `model=`/`model_params=` to select *which*
tabulated nucleosynthesis calculation it draws yields from (independent
of the mass window above) -- see
[Selecting Yield Models](../examples/yield-models.md) for the full
walkthrough, including stochastic (per-star) model parameters and
supplying your own yield tables entirely.

### Writing your own channel

Subclass `MassRangeChannel` for another mass-triggered channel, or
`Channel`/`PopulationChannel` directly for something that doesn't fit
either existing pattern -- see
[`Channel`][crimsons.yields.base.Channel] in the API reference for the
methods a channel must implement.
