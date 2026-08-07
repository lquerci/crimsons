# Results & Config (`crimsons.results`, `crimsons.config`, `crimsons.chemistry`)

The data returned by a run (`EnrichmentResult`), the configuration that
produced it (`RunConfig`, also the basis for cache keys -- see
[Caching & Persistence](../examples/caching-results.md)), and the
tracked-element / solar-metallicity constants shared across the package.

## Results

::: crimsons.results
    options:
      members_order: source

## Run configuration

::: crimsons.config
    options:
      members_order: source

## Chemistry constants

::: crimsons.chemistry
    options:
      members_order: source
      show_source: true
