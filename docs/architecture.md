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

Engine registration is explicit in `dockmate_vs/docking/__init__.py`. Installing
a binary alone does not register it; automatic entry-point discovery is not
implemented. GUI selectors, headless validation, protocol sweeps, adaptive
dispatch, score labels, provenance, and `doctor` use the shared registry.

1. Implement the `DockingEngine` ABC in `docking/base.py`: `dock()`,
   `parse_output()`, and `get_version()`. Return its `DockingResult` and
   `DockingPose` dataclasses, not the adaptive-attempt dataclass.
2. Override `from_config(config)` when the adapter needs binaries or custom
   options. It receives common campaign settings plus the selected engine's
   options in `config["engine_options"]`; its `binaries` mapping contains
   resolved paths. The base implementation forwards common sampling settings
   and custom options to the constructor. Adapters must validate their options,
   impose subprocess timeouts, and report failures.
3. Declare capabilities and register the adapter once:

   ```python
   # In docking/__init__.py, after importing the adapter:
   from .my_engine import MyDockEngine
   register_engine(
       "mydock", MyDockEngine,
       EngineCapabilities(
           binary_names=("mydock",),
           supports_exhaustiveness=False,
           requires_cavity_preparation=True,
           site_definition_label="mydock_cavity_v1",
           supports_score_only_rescoring=False,
           receptor_format="pdb", ligand_format="sdf", output_suffix=".sdf",
           docking_score_label="MyDock", score_units="",
       ),
   )
   ```

| Capability | Contract |
| --- | --- |
| `supports_exhaustiveness` | Enables the GUI control and exhaustiveness sweep; otherwise each seed/box/water condition runs once. |
| `requires_cavity_preparation` | The adapter owns cavity preparation; exposes cavity settings and marks legacy site provenance. |
| `site_definition_label` | Versioned site-definition identity used for resume; change it when cavity semantics change. The CSV field retains the legacy name `rdock_cavity_version`. |
| `supports_score_only_rescoring` | Output is compatible with the workflow's separate Smina rescoring step, not a claim that the adapter has native score-only mode. Unsupported rescore sweeps are rejected. |
| `binary_names` | Executables shown by `doctor` and captured in run provenance. |
| `receptor_format`, `ligand_format` | Prepared receptor PDB/PDBQT and ligand SDF/PDBQT, respectively. |
| `output_suffix` | `.pdbqt`, `.sdf`, or `.sd`. The adapter must write the requested output path. |
| `scoring_functions` | Optional selectable scoring names; one entry fixes the choice, zero means no scoring selector. |
| `score_units`, `docking_score_label`, `color` | Score-axis labeling and engine comparison styling; do not label unitless scores kcal/mol. |

Common legacy settings (`smina_binary`, `vina_binary`, `rdock_root`,
`rdock_runs`, and `rdock_radius`) still work. Custom options and binary overrides
go under `single.engine_options` in YAML:

```yaml
single:
  engine: mydock
  engine_options:
    mydock:
      binaries:
        mydock: /absolute/path/to/mydock
      custom_sampling: 20
protocol_sweep:
  engines: [mydock]
```

The GUI's **Engine options (JSON)** field accepts the corresponding mapping,
for example `{"mydock":{"binaries":{"mydock":"/absolute/path/to/mydock"},"custom_sampling":20}}`.
Options are saved in the manifest and included in protocol-resume compatibility.
The adapter defines which custom options it accepts. The cavity panel retains
rDock's legacy settings; other cavity engines use their own engine-options keys.
`doctor` reports registered binaries using PATH and RBT_ROOT; it does not load a
campaign's custom binary overrides. Optional engines do not make `doctor
--strict` fail unless registered with `required_by_doctor=True`.

Preserve output pose order, one-based mode numbers, and scores. The campaign
contract is **lower is better**. Normalize higher-is-better scores in both
the saved output and the returned result; document that normalization and retain
raw scores in metadata. PDBQT output needs Vina-style score remarks (or Smina
minimizedAffinity); SDF/SD needs a `SCORE` or `SCORE.INTER` property. These are the
formats understood by the shared results reader and pose viewer.

`dock()` also receives workflow context in keyword arguments:
`ligand_smiles`, `ligand_name`, `prepared_receptor_pdb`, and
`reference_ligand_pdb`. Consume or ignore these explicitly, never forward them
as arbitrary command-line flags. Preserve PDBQT atom serials, names, types and
the prepared input geometry before docking for strict pose export. SDF output
must preserve the prepared ligand's chemistry.

The adaptive cascade's default composition remains Smina/Vina and optional
rDock. A programmatic caller can supply `DockingProtocol(engine="mydock", ...)`
entries through `pipeline.protocols`; execution uses the same registry.
Registration does not automatically add an engine to that scientific cascade.

Add at least two tests:

- Adapter tests with mocked subprocess calls: arguments, timeout/nonzero exit,
  output order/scores and version handling.
- A headless protocol campaign over a tiny local workbook: real dispatch,
  condition expansion, result export, resume, and refusal of incompatible
  options. `tests/test_engine_registry.py` demonstrates a registered dummy
  adapter without constructing Tk.

Also perform a real-binary smoke run and GUI selection check before claiming
scientific support for a new engine. Registration supplies integration, not
validation of an external docking program.

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
