# Changelog

All notable changes to DockMate-VS are documented here. The project
uses semantic versioning after the initial public release.

## [0.1.0] - Unreleased

### Added

- Spreadsheet-driven Protocol Development and Screening workflows.
- Sweeps over docking engine, box geometry, water handling, exhaustiveness,
  seed, and supported Smina score-only rescoring functions.
- Matched-control, labelled-assay, and unlabelled-screening input models.
- Per-structure and per-target enrichment, early-enrichment, pose-recovery,
  runtime, and failure reporting.
- Symmetry-aware heavy-atom RMSD and separate native, best-score, and
  best-RMSD pose inspection.
- Resumable campaigns with manifests, machine-readable results, rendered
  summaries, and provenance for code, dependencies, and external engines.
- Offline docking after preflight retrieval of required public structures.
- Automated unit and integration tests for the screening and reporting paths.
- Optional linux/amd64 Docker backend with bundled headless docking tools.
- Active-environment installer for external docking, preparation,
  pocket-detection, and optional visualization programs, with automatic fast
  solver selection and a guard against prolonged classic-Conda solves.
- Architecture, output-contract, third-party licensing, and release-validation
  documentation.

### Changed

- Docking executables are discovered from the active `PATH` or explicit settings
  instead of machine-specific default paths.
- Ligand preparation now retains one normalized charge state rather than
  redocking chemically identical copies labelled as different protonation states.
  Version 0.1 explicitly does not provide pH-aware ionization-state enumeration.

### Fixed

- The binding-site self-docking helper now calculates RMSD in the receptor
  coordinate frame instead of fitting a displaced pose onto the crystal ligand.

[0.1.0]: https://github.com/gkoorsen/DockMate-VS/releases/tag/v0.1.0
