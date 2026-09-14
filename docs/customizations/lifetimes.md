# Stellar Lifetimes

A star's main-sequence(-ish) lifetime sets *when* its channel's enrichment
reaches the ISM (for mass-triggered channels -- see
[Enrichment Channels](channels.md)) and bounds the delay-time distribution
for SN Ia. `StellarLifetime` (the default `lifetime_fn` for `Simulation`)
implements a two-source, mass-and-metallicity-dependent fit, ported from a
Fortran chemical-evolution code.

## Population II/I: Raiteri et al. (1996)

For metallicity above the Population III threshold (see
[Metallicity & Population III](metallicity.md)), lifetime is a
metallicity-dependent quadratic in $\log_{10} m$:

$$
\log_{10}\!\left(\frac{t}{\mathrm{yr}}\right) = a_0 + a_1 \log_{10} m + a_2 \left(\log_{10} m\right)^2
$$

with coefficients that are themselves quadratic in $\log_{10} Z$:

$$
\begin{aligned}
a_0 &= 10.13 + 0.07547\,\log_{10}Z - 0.008084\,(\log_{10}Z)^2 \\
a_1 &= -4.424 - 0.7939\,\log_{10}Z - 0.1187\,(\log_{10}Z)^2 \\
a_2 &= 1.262 + 0.3385\,\log_{10}Z + 0.05417\,(\log_{10}Z)^2
\end{aligned}
$$

$Z$ here is the **absolute** metallicity, clamped to $[7\times10^{-5},\ 3\times10^{-2}]$
before taking its log (matching the range the original fit was calibrated
over) -- values outside that range use the coefficients at the nearest
boundary rather than extrapolating the quadratic.

## Population III: Raiteri (1996) + Schaerer (2002)

Below the Population III threshold, the single Raiteri+1996 fit is
replaced by a piecewise combination, split by mass:

$$
\log_{10}\!\left(\frac{t}{\mathrm{yr}}\right) =
\begin{cases}
\text{Raiteri+1996 (above), evaluated at this star's mass} & m < 5\,M_\odot \\[4pt]
9.785 - 3.759\log_{10}m + 1.413(\log_{10}m)^2 - 0.186(\log_{10}m)^3 & 5 \le m \le 500\,M_\odot \\[4pt]
\text{the }m=500\,M_\odot\text{ value of the line above} & m > 500\,M_\odot
\end{cases}
$$

The middle branch is the Schaerer (2002) fit for very massive,
metal-free-ish stars; above $500\,M_\odot$ it's held flat at the
$500\,M_\odot$ value rather than extrapolated.

## Units

Internally the fits are evaluated in years and then converted to **Myr**
(the unit `StellarLifetime.__call__` returns, and the unit CRIMSONS'
`time_grid` and delay times are expressed in throughout).

## Plugging in your own fit

`lifetime_fn` just needs to be a callable
`(mass: np.ndarray, metallicity: float) -> np.ndarray` (in Myr, vectorized
over mass) -- subclass
[`LifetimeFunction`][crimsons.stars.lifetimes.LifetimeFunction] to plug in
whatever fit you need (e.g. Padovani & Matteucci 1993, or your own
isochrone-based table), and pass it as `Simulation(..., lifetime_fn=...)`.
