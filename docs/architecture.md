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
4. Protein preparation creates one signature-validated receptor per unique
   source/site/chain/water condition and materializes identical case-local
   copies. Ligand preparation creates or safely reuses tool- and
   configuration-versioned enumerated variants.
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
- `prepared_receptors/` stores campaign-level receptor preparations and the
  signatures used to validate safe reuse; each case still receives local
  receptor files for inspection and pose viewing.

Case records distinguish run status, scores, control labels, RMSD metrics,
rescoring metrics, preparation descriptors, output paths, and error messages.
Fields that do not apply, such as RMSD for an unknown screening compound, remain
null rather than being assigned a numerical sentinel.

`utils/result_records.py` normalizes all exported native-pose RMSD fields,
including adaptive attempt records and restart checkpoints. `rmsd_status` is
`ok`, `not_computable`, or `not_applicable`; `error_message` explains applicable
failures. CSV uses empty cells and JSON uses null. Legacy sentinels remain
readable. A completed docking with failed RMSD can still resume, and the failed
comparison remains in the applicable recovery denominator. Do not infer
experimental activity from successful pose mapping.

The molecular coordinate formats and downstream simulation boundary are
documented in [Outputs and MD handoff](md_handoff.md).

### Protocol Development

The `protocol_development/` directory contains:

- `protocol_development_manifest.json` with the sweep and software provenance;
- `protocol_development_results.csv` with one row per expanded condition;
- `protocol_development_candidates.json` with near-equivalent pose-recovery
  candidates for subsequent labelled enrichment runs;
- the generated Markdown summary, recommendations, and charts.

Baseline columns describe the docking engine's original pose ordering. Rescored
columns describe the selected score-only method on those same saved poses.

Consumers should prefer named fields over column position and preserve unknown
fields when transforming records. Additive fields may be introduced in a minor
release; incompatible schema changes must be documented in `CHANGELOG.md`.

## Restart and provenance

A campaign resumes only when its saved manifest and restart-schema version are
compatible with the current input and settings. Completed cases with their
expected output files are skipped; failed or incomplete cases are retried.
Changing a preparation, docking, or scoring condition invalidates reuse for the
affected campaign. Receptor and ligand caches additionally include source,
policy, schema, and preparation-tool identities. Changing a non-swept Protocol
Development setting also invalidates reuse. Protocol campaigns reject result
files containing crystal complexes from another workbook rather than merging
unrelated targets.

For publication, retain the input checksum, manifest, raw results, summaries,
external-tool versions, git revision, seed, CPU count, failures, and data-source
terms. A favourable summary alone is not sufficient provenance.

## Extension points

### Adding a docking engine

Engine extensions currently require source changes. Installing another binary
or adding a name to a YAML file does not register an engine. There is no public
plugin registry or entry-point discovery API in this release.

1. Implement `DockingEngine` in `dockmate_vs/docking/base.py`: `dock()` accepts
   receptor/ligand paths, box center and size, and an output path;
   `parse_output()` reads the saved poses; `get_version()` reports the binary
   version. Return the `DockingPose`/`DockingResult` types from that same module.
   `docking/smina.py` is the existing adapter example. The similarly named
   dataclass in `adaptive_docking.py` represents an adaptive protocol attempt,
   not this adapter interface.
2. Preserve pose order and the score for each pose. The current adapter contract
   and campaign baseline ranking assume **lower is better**. Do not pass a
   higher-is-better score through unchanged. Document any normalization, retain
   raw values in metadata, and avoid calling unitless scores kcal/mol. Save
   poses in a format supported by `_load_poses_and_scores()` in `gui/app.py`, or
   extend that reader and the viewer/topology conversion tests.
3. Wire the adapter into `DockMateVSApp._run_single_case()`, used by both GUI and
   `HeadlessDockMateRunner` campaigns. Handle its actual receptor/ligand format,
   binding-site preparation, sampling controls, timeout, and output path.
   Headless campaigns reuse this method without constructing Tk. Separately
   update `AdaptiveDockingPipeline._run_single_protocol()` if the engine will
   participate in the adaptive protocol cascade.
4. Add selection and validation in `headless.load_campaign_config()` and the
   GUI's engine selectors, `_collect_single_config()`, and protocol-sweep
   expansion. Update the `DockingEngine` enum in `config/schema.py` for the
   standalone configuration API. Exhaustiveness does not apply to all engines;
   rDock additionally needs cavity preparation, runs, and cavity provenance.
   Implement only the controls the new engine supports. Score-only rescoring
   remains a separate Smina operation and needs compatible pose files.
5. Record version, preparation, sampling, and score semantics in the manifest.
   Review `_protocol_resume_signature()` and condition identities in `gui/app.py`
   so incompatible results cannot be reused. Update score-source labels and
   units in `gui/unknown_charts.py`, dependency reporting in `headless.doctor()`,
   installation instructions, and third-party notices where applicable.

At minimum, add two kinds of tests:

- An adapter test with mocked subprocess calls covering command arguments,
  timeout/nonzero exit, parsed pose order/scores, and version reporting.
- A campaign test using `HeadlessDockMateRunner` and a tiny local workbook,
  covering Protocol Development expansion, exported results, and resume. Follow
  `tests/test_headless_cli.py` and `tests/test_screening_integration.py`; prove
  the engine is reached through the campaign path without initializing Tk.

Also perform a real-binary smoke run and a GUI selection check before claiming
support. A dummy adapter test establishes integration, not scientific validation
of the external program.

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
