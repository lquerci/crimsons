# Changelog

All notable changes to CRIMSONS are documented here, following the spirit
of [Keep a Changelog](https://keepachangelog.com/en/1.1.0/){: target="_blank" } and
[Semantic Versioning](https://semver.org/){: target="_blank" }.

## [0.1.0]

Initial development release.

- Core simulation loop (`Simulation`, `EnrichmentResult`, `RunConfig`).
- IMFs: `Salpeter1955`, `Kroupa2001`, `FlatIMF`, `FunctionalIMF`, `Chabrier2003`.
- Channels: `SNII`, `AGB`, `PISN`, `SNIa` (DTD and single-burst modes).
- `StellarLifetime` (Raiteri et al. 1996 / Schaerer et al. 2002).
- HDF5-backed yield-table loading with fixed and stochastic model
  parameters.
- Result caching keyed on a hashed `RunConfig`.

[Unreleased]: https://github.com/your-org/crimsons/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/your-org/crimsons/releases/tag/v0.1.0
