# DockMate-VS

[![Tests](https://github.com/gkoorsen/DockMate-VS/actions/workflows/tests.yml/badge.svg)](https://github.com/gkoorsen/DockMate-VS/actions/workflows/tests.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-0B6E69.svg)](LICENSE)

DockMate-VS is a spreadsheet-driven desktop application for developing,
recording, and applying molecular-docking protocols. It keeps native-pose
recovery, active/inactive enrichment, and score-only screening as separate
questions while producing resumable, machine-readable campaigns.

![DockMate-VS workflow](docs/images/dockmate_vs_workflow.png)

## Why use it?

- Sweep docking engine, box definition, water handling, exhaustiveness, random
  seed, and rescoring conditions before committing to a screen.
- Compare baseline and Smina score-only rankings on the same saved poses.
- Calculate symmetry-aware heavy-atom RMSD without superimposing a displaced
  docking pose onto its native ligand.
- Evaluate labelled controls per receptor structure and target while keeping
  unlabelled compounds out of enrichment statistics.
- Triage top screening hits with PoseBusters plausibility checks and PLIP
  native-contact similarity.
- Describe a complete batch campaign in one reviewable spreadsheet, including
  explicit compound-receptor pairings for multi-receptor docking.
- Resume compatible campaigns and retry missing or failed outputs.
- Inspect native, best-score, and best-RMSD poses and export CSV, JSON, Markdown,
  charts, and viewer files.

## Installation

### 1. Install conda (if you don't have it)

We recommend **Miniforge**, which defaults to the conda-forge channel that DockMate-VS builds against. The one-liner below picks the right installer for your OS and architecture:

```bash
cd ~
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash "Miniforge3-$(uname)-$(uname -m).sh"
```

Accept the licence, keep the default location, and answer `yes` when it offers to run `conda init`. Then open a new terminal, or reload:

```bash
source ~/.bashrc      # macOS: source ~/.zshrc
conda --version
```

### 2. Install DockMate-VS

From Github:
```bash
git clone https://github.com/gkoorsen/DockMate-VS.git
cd DockMate-VS
```
Alternatively, download the Zenodo distribution at https://zenodo.org/records/22283782

Set up the environment:
```bash
conda env create -f environment.yml
conda activate dockmate-vs
python -m pip install -e .
```

Pose plausibility and native-contact similarity reports use PoseBusters and
PLIP. They are installed by `environment.yml`. If you are updating an existing
environment instead of recreating it, run:

```bash
python -m pip install -e ".[pose-quality]"
```

New runs save each prepared ligand variant as both PDBQT and SDF. PoseBusters
uses the prepared SDF topology with the selected docked coordinates when
available, avoiding misleading chemistry inferred from PDBQT alone.

Your prompt should now begin with `(dockmate-vs)`. Verify:

```bash
python -c "import dockmate_vs; print(dockmate_vs.__version__)"
```

### Already using Anaconda or Miniconda?

The environment will still build, but conda may refuse with `CondaToSNonInteractiveError` if Anaconda's default channels are in your
configuration. Either remove them:

```bash
conda config --add channels conda-forge
conda config --remove channels defaults
conda config --set channel_priority strict
```

or accept the channel terms as prompted by the error message.

### Slow solve?

`conda env create` can take several minutes. Miniforge includes `mamba`, which is much faster:

```bash
mamba env create -f environment.yml
```

### 3. Install the external tools

*Optional*: Download LigPlot+
Download LigPlot+ from https://www.ebi.ac.uk/thornton-srv/software/LigPlus/. This will allow invoking LigPlot+ to view 2D representations in the DockMate-VS Pose Viewer.

Install the external tools required:
```bash
scripts/install_external_tools.sh --with-pymol --ligplus-archive PATH_TO_LIGPLOT_PLUS_ARCHIVE
```
If LigPlot+ was downloaded, replace `PATH_TO_LIGPLOT_PLUS_ARCHIVE` with its path
(quote paths containing spaces). Otherwise, omit `--ligplus-archive` and its
argument. LigPlot+ installation is not required.

Docking engines are installed into separate conda environments (their
native dependencies conflict) and linked into `dockmate-vs`.

### 4. Reactivate the environment

Required — the installer writes conda activation hooks that set `RBT_ROOT`
for rDock, and these only take effect on a fresh activation:

```bash
conda deactivate && conda activate dockmate-vs
```

Confirm:

```bash
echo "$RBT_ROOT"    # should print a path ending in /envs/dockmate-rdock
vina --version
```

### Optional core container

The core Docker image provides a reproducible headless environment containing
DockMate-VS, Vina, Smina, rDock, Open Babel, fpocket, OpenMM, and PDBFixer:

```bash
scripts/dockmate-docker build
scripts/dockmate-docker doctor
scripts/dockmate-docker protocol examples/campaign.protocol.yml
```

The native GUI can select **Docker** as its execution backend while keeping
results and the PyMOL/LigPlot+ launchers on the host. PyMOL and LigPlot+ are not
redistributed in the image; see [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md).
The current image targets `linux/amd64` because the available Smina package is
not published for every architecture.

## 4. Launch

```bash
dockmate-vs
```

or, from a source checkout:

```bash
python scripts/launch_dockmate_vs.py
```

Headless campaigns use YAML or JSON configuration files:

```bash
dockmate-vs protocol --config examples/campaign.protocol.yml
dockmate-vs screen --config examples/campaign.screen.yml
dockmate-vs report --run results/dude_aces_screen
dockmate-vs doctor
```

## First reproducible example

1. Open the **Protocol Development** tab.
2. Select `examples/dude_aces_protocol_development_1E66.xlsx`.
3. Use the settings in `examples/campaign.protocol.yml`, or run that campaign
   directly with the headless command shown above.
4. Select a clean output directory and start the run.
5. Review best-pose, Top-1/5/10, best-pose-rank, rescoring, factor effects, and
   the near-equivalent pose-recovery candidate set.

The protocol workbook contains the public DUD-E ACES 1E66/HUX native-pose
control. A successful run writes a protocol manifest, condition-level CSV,
Markdown summary, machine-readable candidate protocols, and interactive charts
below the selected output directory. Compare qualified candidates in separate
labelled enrichment runs and freeze one protocol before screening unknowns.

For enrichment, select
`examples/dude_aces_screening_subset_1E66_seed42.xlsx` in the **Screening** tab
or run `dockmate-vs screen --config examples/campaign.screen.yml`. It contains a
fixed, score-independent subset of 20 DUD-E clustered actives and 200 DUD-E
property-matched decoys at one receptor. Matched decoys are presumed
non-binders, not experimentally confirmed inactives. This compact campaign
demonstrates the complete workflow rather than establishing a new scoring-
function benchmark. The frozen example uses Vina scoring, a 4 A margin,
selective water retention, exhaustiveness 8, 20 requested modes, seed 42, and
thorough score-based ligand-variant selection. Exact scores and poses can vary
with docking-engine and preparation-tool builds, hardware, and random seed. See
`examples/README.md` for reference results, source checksums, reconstruction,
and interpretation guidance.

## Input models

The spreadsheet is a campaign definition rather than a one-compound form. Each
row identifies a receptor structure, binding-site reference, compound, and
optional experimental class. This makes large batches auditable and allows the
same library to be applied to several receptor structures by listing each
compound-receptor pairing explicitly. Reports retain the receptor and target
grouping instead of assuming that scores from different structures are directly
comparable.

Only the **first worksheet** is read, regardless of its name. Put the campaign
table there, with column headings in the first row. Additional worksheets, such
as metadata or validated-structure lists, are not loaded as docking jobs.

Common column headings are:

| Column | Meaning |
| --- | --- |
| `Protein` | Target name used to group results. |
| `PDB_ID` | Receptor structure identifier. |
| `Ligand` | Co-crystal ligand residue name used as the binding-site reference. |
| `Target_Ligand` | Name of the compound to dock, not its activity class. |
| `SMILES` | Structure of the compound to dock. Supply this for unknown compounds. |
| `decoy compound` | Comma-separated decoy names, in the same order as their SMILES. |
| `decoy SMILES` | Comma-separated decoy structures for a matched-control row. |
| `label` or `control_label` | Explicit class: `1` active, `0` decoy/inactive, blank unknown. |

Column headings are detected automatically, ignoring case, spaces, underscores,
and punctuation. The application supports two ways to specify controls:

- **Matched controls:** omit explicit label columns from the worksheet. A row
  with nonblank `decoy SMILES` expands to one active plus one docking case for
  each listed decoy. For example, 30 decoy SMILES produce one active and 30
  decoys. The active is named by `Target_Ligand` (or `Ligand` when no compound
  name is supplied). Missing decoy names receive names such as `decoy_1`.
  A row with blank decoy SMILES is an unknown screening compound, even if it
  contains decoy names. Active-control and unknown rows can share this worksheet.
- **Explicit labels:** put one compound per row and use `label` or
  `control_label` to mark its class. `active` and `inactive`/`decoy` are also
  accepted values. Multiple actives and negatives can share a receptor
  structure. Leave a compound's label blank to treat it as unknown.

**Label-column precedence:** if a recognized label column exists anywhere in the
table, its values determine classification for every row and the decoy columns
are not expanded. This applies even if the entire label column is blank. Remove
the label column completely when using the matched-control format; alternatively,
list every active, decoy, and unknown as a separate explicitly labelled or blank-
labelled row. Other recognized label headings include `activity_label`,
`active_decoy`, `class`, `activity`, `is_active`, `active`, and `actives`.

Unknown screening compounds are ranked separately and never enter ROC AUC or
enrichment calculations. A worksheet with neither decoy SMILES nor explicit
labels therefore defines an unknown-compound screen. Compound names, workbook
names, and output-folder names do not establish activity labels.

The **Filters** tab can exclude known
additives/cofactors and optionally sample unlabelled screening compounds; all
labelled or matched controls are always retained.

## Outputs

Each screening run writes:

- `run_manifest.json`: settings, cases, code revision, dependency versions, and
  external docking-tool versions.
- `redock_progress.json`: incremental state used for compatible restart.
- `redock_results.csv` and `redock_results.json`: case-level scores, status,
  preparation descriptors, pose metrics where applicable, and provenance.
- `redock_summary.json` and `redock_summary.md`: structured and human-readable
  summaries regenerated from raw result records.

Protocol-development runs produce an analogous manifest, condition-level CSV,
summary, recommendations, and plots.

For runs containing unknown compounds, **Results > Charts > Unknown docking
scores** plots every available selected ranking score by receptor structure and
highlights the lowest-scoring compound in each group. It uses the same score
selection as the top-ranked-compound summary: GNINA score, then rescoring score,
then the original docking score. Hover over a point for its compound name,
target, structure, engine, score, and score source. Use the structure selector
to focus on one group and the plot toolbar to zoom or export an image. Active
controls and decoys are excluded; scored/total counts show missing or failed
cases.

The component boundaries, campaign flow, output contracts, and extension points
are described in [`docs/architecture.md`](docs/architecture.md).

## Testing

```bash
python -m pip install -e ".[test]"
python -m pytest dockmate_vs/tests -q
```

GitHub Actions runs the focused suite on Python 3.9, 3.11, and 3.12 under Linux.
External docking binaries are mocked in automated integration tests. Before a
release, maintainers also verify the documented example with an installed
docking engine and record the engine version in the generated manifest.

## Interpretation and limitations

The platform orchestrates established docking engines; it does not introduce a
new scoring function. A recoverable native-like pose does not imply that the
engine ranks it first, and enrichment does not establish biochemical activity.
Raw scores should normally be compared within one receptor structure. Standard
Vina/Smina workflows are not valid for covalently linked complexes, which are
detected and skipped.

See [`USER_GUIDE.md`](USER_GUIDE.md) for complete workflow instructions,
troubleshooting, and results interpretation.

## ESR1 case-study materials

The reconstruction-based [ESR1 supplement](supplementary/README.md) provides
the exact subset identifiers, historical campaign settings, structure-free
results, and scripts supporting the SoftwareX case study.
Download the [supplement ZIP](supplementary/packages/DockMate-VS_ESR1_Supplement_S1.zip)
and its [SHA-256 checksum](supplementary/packages/DockMate-VS_ESR1_Supplement_S1.zip.sha256).
LIT-PCBA structures must be obtained separately; they are not redistributed.
This supplement is hosted on GitHub separately from the Zenodo software archive.

## Citation and support

Citation metadata are provided in [`CITATION.cff`](CITATION.cff). Please cite
the version used in your work. The current software version is **0.1.1**:

> Koorsen, G. (2026). DockMate-VS (v0.1.1). Zenodo.
> https://doi.org/10.5281/zenodo.22685196

When the SoftwareX article is available, cite both the article and the archived
software release.

- Problems and feature requests: [GitHub Issues](https://github.com/gkoorsen/DockMate-VS/issues)
- Support contact: [gkoorsen@uj.ac.za](mailto:gkoorsen@uj.ac.za)
- Contribution guide: [`CONTRIBUTING.md`](CONTRIBUTING.md)
- Release history: [`CHANGELOG.md`](CHANGELOG.md)

## License

DockMate-VS is distributed under the [MIT License](LICENSE). External
programs and benchmark datasets retain their own licenses and terms.
