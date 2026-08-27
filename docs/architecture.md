# DockMate-VS architecture

This document describes the release architecture, campaign data flow, output
contracts, and intended extension points. User-facing operation is documented in
[`USER_GUIDE.md`](../USER_GUIDE.md).

## Entry points

- `dockmate-vs` starts the Tkinter GUI when no subcommand is supplied.
- `dockmate-vs protocol --config FILE` runs a protocol-development campaign.
- `dockmate-vs screen --config FILE` runs a frozen screening protocol.
- `dockmate-vs report --run PATH` regenerates a report from saved raw results.
- `dockmate-vs doctor [--strict]` reports available Python and external tools.

`dockmate_vs/cli.py` owns command parsing. `dockmate_vs/headless.py` validates
YAML/JSON campaign files and adapts the same campaign implementation used by the
GUI; headless execution is not a separate docking pipeline.

## Components

| Component | Responsibility |
| --- | --- |
| `dockmate_vs/gui/app.py` | Campaign orchestration, spreadsheet expansion, restart state, reporting, and desktop views |
| `dockmate_vs/headless.py` | Configuration loading, non-interactive execution, report regeneration, and installation diagnostics |
| `dockmate_vs/preparation/` | Protein repair/water policies and ligand variant/conformer preparation |
| `dockmate_vs/binding_site/` | Co-crystal and configured binding-site definitions |
| `dockmate_vs/docking/` | Docking-engine interfaces and Smina/Vina execution |
| `dockmate_vs/adaptive_docking.py` | Ordered adaptive redocking attempts and attempt-level provenance |
| `dockmate_vs/analysis/` | Enrichment and bootstrap analyses |
| `dockmate_vs/utils/` | RMSD, structure alignment, result analysis, and shared file utilities |
| `dockmate_vs/config/` | Typed configuration models used by lower-level components |

The platform delegates pose generation and primary scoring to external docking
engines. DockMate-VS prepares inputs, invokes those engines, preserves outputs,
and calculates workflow-level pose-recovery and enrichment statistics. It does
not implement a new docking score.

## Campaign flow

1. The spreadsheet loader normalizes receptor, compound, control, and optional
   label columns.
2. Matched-control rows expand to separate active and decoy cases. Explicit
   assay labels remain attached to their compounds; unlabelled compounds remain
   screening samples.
3. Preflight downloads all required public structures before docking starts.
4. Protein and ligand preparation create the receptor and enumerated ligand
   variants for each case.
5. Protocol Development expands the requested engine, box, water, search, seed,
   and rescoring combinations. Screening applies one frozen condition.
6. Each completed case is written to restart state before the next case runs.
7. Raw records are converted into structured JSON, tabular CSV, Markdown, and
   rendered GUI summaries. Reports can be regenerated without redocking.

Network access is required only for unresolved public structures during
preflight. Locally available inputs, preparation, docking, analysis, and report
generation do not require continuous connectivity.

## Output contracts

### Screening

- `run_manifest.json` records creation time, software and dependency versions,
  external binary paths/versions, normalized configuration, and case identities.
- `redock_progress.json` contains incrementally completed case records used by a
  compatible restart.
- `redock_results.json` contains a `results` array of case-level records.
- `redock_results.csv` provides the same case-level fields in tabular form.
- `redock_summary.json` is the machine-readable aggregate report.
- `redock_summary.md` is generated from the same aggregate report.

Case records distinguish run status, scores, control labels, RMSD metrics,
rescoring metrics, preparation descriptors, output paths, and error messages.
Fields that do not apply, such as RMSD for an unknown screening compound, remain
null rather than being assigned a numerical sentinel.

### Protocol Development

The `protocol_development/` directory contains:

- `protocol_development_manifest.json` with the sweep and software provenance;
- `protocol_development_results.csv` with one row per expanded condition;
- the generated Markdown summary, recommendations, and charts.

Baseline columns describe the docking engine's original pose ordering. Rescored
columns describe the selected score-only method on those same saved poses.

Consumers should prefer named fields over column position and preserve unknown
fields when transforming records. Additive fields may be introduced in a minor
release; incompatible schema changes must be documented in `CHANGELOG.md`.

## Restart and provenance

A campaign resumes only when its saved manifest is compatible with the current
input and settings. Completed cases with their expected output files are skipped;
failed or incomplete cases are retried. Changing a preparation, docking, or
scoring condition invalidates reuse for the affected campaign.

For publication, retain the input checksum, manifest, raw results, summaries,
external-tool versions, git revision, seed, CPU count, failures, and data-source
terms. A favourable summary alone is not sufficient provenance.

## Extension points

### Docking engines

New engines should implement the interface in `dockmate_vs/docking/base.py`,
return explicit pose/score records, expose version information, and add tests for
command construction, failure handling, and ranking direction. Add the engine to
configuration validation and both GUI/headless selection paths.

### Preparation and binding sites

Preparation policies belong in `dockmate_vs/preparation/`; site-definition
strategies belong in `dockmate_vs/binding_site/`. New options must be represented
in the run manifest so a resumed or published campaign is unambiguous.

### Analyses and reports

New metrics should be calculated from raw case records, define their valid
denominator and ranking direction, and return null when required data are absent.
Rendered text and charts should consume the same structured summary to avoid
divergent values.

## Verification

The focused test suite is under `dockmate_vs/tests/`. It covers spreadsheet
classification, campaign configuration, control expansion, restart behavior,
RMSD handling, rescoring, enrichment, reporting, result loading, and external
viewer integration. External binaries are mocked in automated integration tests;
release verification should additionally run the documented example with the
actual engine build named in its manifest.

Run the suite from the repository root:

```bash
python -m pip install -e ".[test]"
python -m pytest dockmate_vs/tests -q
```
