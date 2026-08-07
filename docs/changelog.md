# Changelog

All notable changes to CRIMSONS are documented here, following the spirit
of [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and
[Semantic Versioning](https://semver.org/).

!!! info "Template"
    This page is a starting point, not a historical record -- back-fill it
    with your actual release history, and update it as part of each
    release going forward. Consider automating it (e.g.
    [towncrier](https://towncrier.readthedocs.io/) or
    [git-cliff](https://git-cliff.org/)) once release cadence picks up.

## [Unreleased]

### Added

- Initial public documentation site.

## [0.1.0]

Initial development release.

- Core simulation loop (`Simulation`, `EnrichmentResult`, `RunConfig`).
- IMFs: `Salpeter1955`, `Kroupa2001`, `FlatIMF`, `FunctionalIMF`.
- Channels: `SNII`, `AGB`, `PISN`, `SNIa` (DTD and single-burst modes).
- `StellarLifetime` (Raiteri et al. 1996 / Schaerer et al. 2002).
- HDF5-backed yield-table loading with fixed and stochastic model
  parameters.
- Result caching keyed on a hashed `RunConfig`.

[Unreleased]: https://github.com/your-org/crimsons/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/your-org/crimsons/releases/tag/v0.1.0
