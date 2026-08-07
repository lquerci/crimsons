# Enrichment Engine (`crimsons.enrichment`)

The per-realization loop `Simulation.run()` calls once per realization:
sample the IMF into log-mass bins, route each bin (or the population as a
whole, for `SNIa`) through the configured channels, and bin the resulting
events onto a shared time grid. See
[Stochastic Sampling](../physics/sampling.md) for the binned Monte Carlo
approach this implements.

::: crimsons.enrichment.engine
    options:
      members_order: source
