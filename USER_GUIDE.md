# DockMate-VS User Guide

## 1. Scope

DockMate-VS supports two stages:

1. **Protocol Development** compares docking and ranking conditions against
   native co-crystal poses.
2. **Screening** applies one frozen protocol to labelled controls, assay
   benchmarks, and/or unlabelled compounds.

Use Protocol Development before Screening when changing receptor preparation,
water treatment, box geometry, docking engine, search effort, or scoring.

## 2. Requirements

### Python environment

Python 3.9-3.12 is supported. From the repository root:

```bash
conda env create -f environment.yml
conda activate dockmate-vs
python -m pip install -e .
```

The GUI uses Tk. Linux users may need their distribution's `python3-tk`
package. On macOS, use a Python build that includes Tk support.

### External programs

| Program | Role | Requirement | Installation instructions |
| --- | --- | --- | --- |
| AutoDock Vina or Smina | Docking | At least one is required | [Vina](https://autodock-vina.readthedocs.io/en/latest/installation.html) / [Smina](https://github.com/mwojcikowski/smina) |
| Smina | Score-only rescoring | Required only when rescoring is enabled | [Smina](https://github.com/mwojcikowski/smina) |
| Open Babel (`obabel`) | Rigid receptor and ligand PDBQT conversion | Required by current Vina/Smina workflows | [Open Babel](https://openbabel.org/docs/Installation/install.html) |
| OpenMM and PDBFixer | Missing-atom and missing-loop repair | Optional but recommended | [OpenMM](https://docs.openmm.org/latest/userguide/application/01_getting_started.html#installing-openmm) / [PDBFixer](https://github.com/openmm/pdbfixer#installation) |
| Reduce and PROPKA | Optional H-bond optimization and pKa diagnostics | Optional pipeline components; see molecular preparation below | [Reduce](https://github.com/rlabduke/reduce#building) / [PROPKA](https://propka.readthedocs.io/en/latest/installation.html) |
| rDock | Alternative docking engine | Optional | [rDock](https://rdock.github.io/installation/) |
| Meeko (`mk_prepare_receptor.py`) | Flexible-receptor PDBQT preparation | Optional | [Meeko](https://meeko.readthedocs.io/en/develop/installation.html) |
| fpocket | Apo binding-site prediction | Optional | [fpocket](https://github.com/Discngine/fpocket#installing) |
| PyMOL and LigPlot+ | Pose and interaction visualization | Optional | [PyMOL](https://pymol.org/) / [LigPlot+](https://www.ebi.ac.uk/thornton-srv/software/LigPlus/) |

After activating `dockmate-vs`, the bundled installer can install all
package-manageable tools and keep their commands on the conda environment's
`PATH`:

```bash
scripts/install_external_tools.sh
```

The installer supports Micromamba, Mamba, or Conda configured with the
libmamba solver. It automatically prefers `micromamba`, then `mamba`, before
using `conda`. This avoids prolonged dependency resolution with Conda's classic
solver. macOS users can install Micromamba with:

```bash
brew install micromamba
```

See the [Micromamba installation
instructions](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html)
for other platforms. Use `--package-manager /path/to/micromamba` to select a
specific executable. The installer stops before invoking Conda's classic
solver by default; `--allow-classic-conda` overrides this guard when the delay
is acceptable.

Use `--with-pymol` to add open-source PyMOL. LigPlot+ must first be obtained
under its own licence; register an extracted installation without copying it
into the repository:

```bash
scripts/install_external_tools.sh --with-pymol \
  --smina-bin /path/to/smina \
  --ligplus-root /path/to/LigPlus
```

The script accepts `--vina-bin`, `--smina-bin`, `--rdock-root`, `--pymol-bin`,
and `--ligplus-root` for existing installations. It validates these paths
before solving the environment, links explicitly supplied binaries into the
active environment, and writes conda activation hooks for `RBT_ROOT` and
LigPlot+ variables. It does not edit `.zshrc`, `.bashrc`, or other shell startup
files. The installer supports macOS and Linux. The current bioconda rDock
package is Linux-only; on macOS, provide an existing rDock root or use the
Docker backend. Windows users should use Docker, WSL, or the linked manual
installation instructions.

If an older copy of the installer repeatedly displays `Solving environment`
and mentions frozen or flexible solves, cancel it with `Ctrl+C`. Install
Micromamba, update the repository, and rerun the installer. Do not run Conda and
Micromamba transactions against the same environment concurrently.

Install tools according to their own documentation and licenses. DockMate-VS
discovers `vina`, `smina`, and `rbdock` from the active executable `PATH` when
they are available. The GUI can override Vina/Smina executable paths and the
rDock root (`RBT_ROOT`). Commands without dedicated GUI fields, including
`obabel`, `reduce`, and `mk_prepare_receptor.py`, must be on `PATH`. Example
Vina/Smina overrides are:

```text
/usr/local/bin/vina
/home/user/apps/smina/smina
```

Test the selected binary before a campaign:

```bash
/path/to/vina --version
/path/to/smina --version
```

LigPlot+ is discovered through `LIGPLOT_BIN`, the executable `PATH`, or a
LigPlus installation named by `LIGPLUS_ROOT`/`LIGPLOT_HOME`. Set
`HET_GROUP_DICTIONARY` when `components.cif` is stored outside that installation.

### Optional Docker backend

Docker can provide the computation environment when installing all docking
tools locally is inconvenient. From the repository root:

```bash
scripts/dockmate-docker build
scripts/dockmate-docker doctor
```

The core image includes Vina, Smina, rDock, Open Babel, fpocket, OpenMM, and
PDBFixer. It deliberately excludes PyMOL and LigPlot+: those remain optional
host applications launched by the native GUI. LigPlot+ cannot be redistributed
in the public image under its standard academic licence.

The image currently targets `linux/amd64`. Apple Silicon systems can run it
through Docker Desktop emulation by retaining the default platform setting.
Docker Desktop must be running before a Docker-backed campaign starts.

## 3. Launching

```bash
dockmate-vs
```

or:

```bash
python scripts/launch_dockmate_vs.py
```

Select an input workbook and a dedicated output directory. Avoid reusing an
unrelated run directory.

The main notebook has four tabs: **Protocol Development**, **Screening**,
**Filters**, and **Pose Viewer**. The first two select the workflow. Filters are
shared input policies, while Pose Viewer inspects completed runs without changing
the selected workflow.

Choose **Local** under Execution backend to use binaries configured in the GUI.
Choose **Docker** to run the campaign in the core image. The GUI mounts the input
workbook read-only, mounts the selected output directory read/write, streams the
container log into the normal progress window, and stores the generated campaign
YAML with the run. Result paths remain valid for the native Pose Viewer.

### Headless commands

The GUI and command line call the same campaign implementation:

```bash
dockmate-vs protocol --config examples/campaign.protocol.yml
dockmate-vs screen --config examples/campaign.screen.yml
dockmate-vs report --run /path/to/completed/run
dockmate-vs doctor
```

Configuration paths are resolved relative to the YAML file. The supplied files
`examples/campaign.protocol.yml` and `examples/campaign.screen.yml` document all
common settings. The protocol campaign uses
`examples/dude_aces_protocol_development_1E66.xlsx`; the screening campaign
uses `examples/dude_aces_screening_subset_1E66_seed42.xlsx`. The latter contains
20 DUD-E clustered ACES actives and 200 property-matched decoys selected with a
fixed seed before docking. DUD-E matched decoys are presumed non-binders rather
than experimentally confirmed inactive compounds. The compact subset exercises
assay-benchmark reporting and must not be presented as a new scoring-function
validation. For container execution, keep the campaign file, workbook, and
output path below the current working directory and run:

```bash
scripts/dockmate-docker protocol examples/campaign.protocol.yml
scripts/dockmate-docker screen examples/campaign.screen.yml
```

See `docs/architecture.md` for component ownership, campaign data flow, output
contracts, and the supported extension points for engines and analyses.

## 4. Spreadsheet schema

Column names are matched case-insensitively where supported. Recommended names
are shown below.

### Why a spreadsheet?

The workbook is an explicit campaign manifest, not merely a convenient way to
enter one ligand. Each row associates a compound with a receptor structure,
binding-site definition, optional experimental label, and provenance notes.
This supports reproducible batch experiments without requiring users to write
Python code:

- hundreds or thousands of docking cases can be queued in one reviewed file;
- one compound library can be docked against several conformations, structures,
  sites, or targets;
- matched controls and assay labels remain attached to the receptor against
  which they are valid; and
- the original workbook can be archived with the manifest and results as an
  independently inspectable description of the campaign.

Multi-receptor docking is explicit rather than implicit: repeat a compound on
one row for every `PDB_ID`/site-ligand combination to be evaluated. This avoids
an accidental Cartesian product and makes the intended jobs visible before the
run. DockMate-VS expands those rows into cases, reuses compatible repeated
native controls, and reports results by receptor structure and target. Raw
scores from different receptor structures are not assumed to share a common
scale.

### Shared receptor columns

| Column | Meaning |
| --- | --- |
| `Protein` | Target name used in grouped reports |
| `Target_Ligand` | Human-readable compound or case name |
| `PDB_ID` | Four-character PDB identifier |
| `Ligand` | Co-crystal/site ligand residue name |
| `Chain` | Ligand/receptor chain where needed |
| `Notes` | Free-text provenance notes |

### Matched-control rows

| Column | Meaning |
| --- | --- |
| `decoy compound` | Matched-decoy identifier |
| `decoy SMILES` | Matched-decoy structure |

A populated decoy name and SMILES create two cases: the native active and the
decoy. Repeated PDB/native-ligand pairs are docked once per compatible condition.
A row with empty decoy fields is **not** automatically a control.

### Assay-benchmark rows

| Column | Meaning |
| --- | --- |
| `SMILES` | Compound structure |
| `label` | `1` for assay active; `0` for assay inactive |

This format supports multiple actives and inactives at one receptor structure.
Assay inactives are not tested as if they were property-matched decoys to one
crystal ligand.

### Unlabelled screening rows

A row with a compound SMILES but no label and no decoy is an unknown screening
sample. It receives scores and ranks but is excluded from ROC AUC, average
precision, enrichment factors, and RMSD success rates.

DockMate-VS detects and uses populated SMILES columns automatically. There is no
separate SMILES toggle because disabling SMILES would silently change a compound
screen into a co-crystal redocking task.

### Filters and random sampling

The **Filters** tab can exclude configured additive and cofactor residue names.
Controls are never removed by these filters. Random sampling applies only in the
Screening workflow and its sample size counts unlabelled compounds only. Every
labelled or matched active/decoy control is retained, while sampled unknowns are
distributed across receptor structures using the supplied seed.

## 5. Molecular preparation

This section documents the version 0.1 implementation, not a generic docking
preparation recipe. The implementation is in
[`dockmate_vs/preparation/protein.py`](dockmate_vs/preparation/protein.py) and
[`dockmate_vs/preparation/ligand.py`](dockmate_vs/preparation/ligand.py).
User-selected preparation settings and external-tool versions are captured in
the run manifest; its git commit identifies fixed implementation defaults.
Prepared files are retained inside each case directory for inspection.

### Receptor preparation

For a co-crystal campaign, preflight obtains the requested PDB structure before
docking begins. The receptor pipeline then performs the following operations:

1. BioPython parses the PDB file. PDBFixer attempts to identify missing
   residues, add missing heavy atoms, and model missing loops when supported.
   If this optional repair fails, the event is logged and the original parsed
   structure is retained rather than silently aborting the whole campaign.
2. The selected crystallographic-water policy is applied to `HOH` and `WAT`
   residues. **Remove all** deletes them; **retain all** preserves them;
   **selective** retains waters satisfying protein-contact, occupancy, and
   optional B-factor criteria. Current campaign defaults require a water oxygen
   to be within 3.0 A of protein atoms, make at least two such contacts, and
   have occupancy at least 0.5. The lower-level API can additionally require a
   water to be within 3.0 A of a supplied site centre; the standard GUI/headless
   co-crystal path does not currently pass that optional centre to receptor
   preparation.
3. The co-crystal ligand named in the workbook's `Ligand` column is removed from
   the receptor before docking. The original ligand is extracted separately to
   define the box and, in Protocol Development, provide the RMSD reference.
4. The pipeline attempts PROPKA analysis at pH 7.4, PDBFixer hydrogen addition,
   and Reduce H-bond optimization when those optional tools are available.
   Missing optional tools and non-fatal failures are written to the run log.
5. Open Babel converts the cleaned receptor to a rigid PDBQT. Only the first
   model and blank/`A` alternate locations are retained for this conversion.
   The current rigid conversion deliberately uses the heavy-atom-cleaned
   structure to avoid Open Babel valence failures; PROPKA output is therefore
   diagnostic and must not be interpreted as an explicitly curated PDBQT
   protonation-state assignment.

The selected site ligand is stripped automatically, but unrelated hetero groups
and additional protein chains are not universally deleted. The **Filters** tab
controls which workbook cases are admitted; it is not a general receptor-cleaning
editor. Inspect `receptor_prepared.pdbqt`, especially for metals, cofactors,
alternate sites, unusual residues, or systems requiring a specific protonation
state. Such choices are target-specific protocol variables and should be
validated rather than assumed.

### Ligand preparation

DockMate-VS prepares each populated SMILES with RDKit and Open Babel:

1. RDKit parses the SMILES, keeps the largest disconnected fragment, and applies
   its charge normalizer. This can remove salts and neutralize removable formal
   charges.
2. RDKit enumerates tautomers and retains at most **Max tautomers** after a
   deterministic heuristic ranking that favours aromatic, carbonyl, and amide
   forms and penalizes charge separation and selected enol/iminol patterns.
3. Explicit hydrogens and three-dimensional coordinates are generated with
   ETKDGv3 using preparation seed 42. Rigid molecules produce one conformer;
   flexible molecules are oversampled, MMFF94s-minimized, and reduced to at most
   **Max conformers** by greedy 0.5 A RMSD-diversity selection beginning with
   the lowest-energy conformer.
4. Each retained tautomer/conformer is minimized with MMFF94s when parameters
   are available and converted from SDF to PDBQT by Open Babel. The normal
   partial-charge order is Gasteiger, MMFF94, then an explicitly warned
   charge-model-free fallback. Metal-containing ligands skip tautomer
   enumeration and use the Open Babel `qeq`, `qtpie`, and `eem` charge fallbacks.
5. The selected ligand-variant policy determines how many prepared variants are
   docked. Adaptive selection uses rotatable-bond and heavy-atom counts; the
   fast option uses the lowest MMFF-energy variant; thorough/all modes dock a
   larger subset. Screening always chooses the final retained result by docking
   score, never by native-pose RMSD.

Version 0.1 does **not** enumerate ligand ionization states as a function of pH.
The configured pH range is reserved for future use, and the normalized SMILES
charge state is the starting state. Curate input charge states externally when
ionization is important, and inspect the saved `ligand_variants/*.pdbqt` files.
The ligand cache key includes the molecule identifier, input SMILES, and relevant
preparation limits, so repeated compatible compounds can avoid regeneration.

### Preparation records

Useful case-level files include:

| Path | Meaning |
| --- | --- |
| `receptor_prepared.pdbqt` | Rigid receptor supplied to Vina/Smina |
| `receptor_prepared.pdb` | PDB representation used by viewers and other engines |
| `crystal_ligand.pdb` | Extracted site ligand and RMSD reference when applicable |
| `ligand_variants/*.pdbqt` | Prepared tautomer/conformer candidates |
| `variants/<variant>/docked.pdbqt` | Engine poses for a docked variant |
| `run_manifest.json` | Frozen settings, versions, paths, and case identities |

Always inspect representative prepared receptors and ligands before a large
campaign. A completed conversion proves file compatibility, not chemical
correctness.

## 6. Protocol Development

### Select conditions

The protocol sweep can vary:

- Vina, Smina, and installed rDock engines;
- co-crystal box margins or fixed box dimensions;
- remove-all, retain-all, and selective crystallographic-water handling;
- exhaustiveness and random seed;
- baseline score ordering and supported Smina score-only methods.

Rescoring does not redock the ligand. It evaluates the selected scoring function
on every saved pose and compares baseline and rescored rankings.

### Ligand variant limits

Set explicit maximum tautomer and conformer counts. Large values multiply the
number of dockings. Begin with two tautomers and two conformers for exploration,
then justify any increase. Screening must select variants by score, not RMSD.

The detailed preparation and variant-selection rules are described in
**Molecular preparation** above. Inspect prepared forms for charge-sensitive
compounds and document any externally curated charge-state policy used to create
the input library.

### Interpret pose recovery

- **Best RMSD:** lowest RMSD generated anywhere in the saved poses.
- **Top-1 RMSD:** RMSD of the top-ranked pose.
- **Top-5/10 success:** whether a pose below the selected RMSD threshold occurs
  within that score-ranked prefix.
- **Best-pose rank:** score rank of the lowest-RMSD pose.

A low best RMSD with a poor Top-1 RMSD means sampling found a native-like pose
but scoring did not prioritize it. Do not report this as Top-1 validation.

Protocol Development reports **pose-recovery candidate protocols**, not an
artificial single winner. Conditions with near-best success and RMSD values
within 0.25 A of the best result are retained as an equivalence group. The
report builds a shortlist of up to eight protocols while preserving eligible
alternatives in scoring, water handling, engine, box, and exhaustiveness. It
keeps original and rescored scoring as separate candidates because enrichment
rankings can differ even when pose rankings do not. The report records selection
confidence; one crystal complex and one seed provide low-confidence evidence
even when RMSD is excellent. The machine-readable
`protocol_development_candidates.json` contains the same shortlist.

Run near-equivalent candidates as separate campaigns against the same labelled
active/inactive or active/decoy enrichment set. Keep the compound set fixed and
compare ROC AUC, average precision, EF1/5/10, failures, score pathologies, and
runtime. Do not pool scores from different candidate protocols into one AUC.

## 7. Screening

### Freeze the protocol

Labelled enrichment is the second protocol-selection stage. Freeze one protocol
only after comparing the pose-recovery candidates on the labelled set. Record
the selected receptor, binding site, water treatment, engine, scoring, box,
exhaustiveness, random seed, number of modes, energy range, ligand-variant
limits, and CPU count before screening unknown compounds. Do not select a
favourable receptor or seed after examining unknown-compound outcomes.

### Run and resume

Remote structures are prefetched before docking. Once preflight succeeds, the
remaining preparation, docking, scoring, and analysis can proceed without
continuous network connectivity.

Progress is saved incrementally. Restart with the same input and compatible
settings to skip completed cases. Changed docking or scoring settings invalidate
reuse. Failed cases and missing output files are retried. Protocol Development
also verifies workbook identity, crystal complexes, and non-swept settings
before resuming; an incompatible output directory is refused rather than mixed
with the new campaign.

### Interpret enrichment

- **ROC AUC:** pairwise discrimination across all score thresholds.
- **Average precision:** precision-recall summary; compare it with active
  prevalence.
- **EF1/5/10:** early enrichment. Tied scores at a cutoff are averaged and a
  possible min-max range is reported.
- **Per-structure AUC:** preferred when receptor structures have different score
  scales.
- **Target macro AUC:** mean of valid per-structure AUC values.
- **Target-pooled AUC:** diagnostic when multiple receptor structures are pooled.

Unknown compounds have no experimental pose and therefore no RMSD. Assay
benchmark labels support enrichment but do not create reference poses.

## 8. Results and inspection

The Results window renders summary sections and charts from the raw result
records. Reopening a run regenerates the summary with current reporting logic.

Use **Load Run Folder...** and select one completed campaign directory directly,
not the parent `output` directory that contains several campaigns. For a
Screening campaign, the selected folder contains `redock_results.json` or
`redock_results.csv`. For Protocol Development, select either the campaign
folder containing the `protocol_development` subfolder or that subfolder itself.
Use **Load Results File...** only when selecting the exact
`redock_results.json`, `redock_results.csv`, or
`protocol_development_results.csv` file.

The **Pose Viewer** tab loads a completed run folder or exact results file and can
display:

- native co-crystal ligand;
- best-scoring pose;
- lowest-RMSD pose when an experimental reference exists.

LigPlot+ interaction diagrams use the pose identified by the selected metric.
For screening compounds, inspect the selected best-score or rescored pose rather
than a nonexistent best-RMSD pose.

## 9. Reproducibility checklist

Before publishing a campaign, retain:

- exact input workbook and its checksum;
- run manifest and software git commit;
- docking-engine and dependency versions;
- result CSV/JSON and generated summary;
- failures and excluded compounds;
- hardware, CPU count, campaign wall time, and random seed;
- data-source DOI, cleaning script, and redistribution terms.

## 10. Troubleshooting

### Docking executable not found

Use the complete path and verify it from the same shell that launches the GUI.
Conda-installed binaries can be located with `which vina` or `which smina` after
activating the environment.

### RDKit is reported as unavailable on Apple Silicon

If startup reports an RDKit import failure containing `code signature invalid`,
replace the affected conda-forge build with the tested release and retry:

```bash
conda install -n dockmate-vs -c conda-forge "rdkit=2025.03.3"
conda activate dockmate-vs
python -c "import rdkit; print(rdkit.__version__)"
dockmate-vs
```

The supplied `environment.yml` excludes the affected RDKit 2026 series for new
environments. A message saying RDKit cannot be imported does not necessarily
mean that the package is absent; DockMate-VS also reports the underlying native
library error to distinguish these cases.

### Multi-model rescoring error

Current rescoring splits a multi-model PDBQT and scores each pose independently.
If an older run reports `Unexpected multi-MODEL input`, reopen or rerun it with a
current release.

### No RMSD in screening

This is expected for non-native compounds. RMSD is only meaningful when the
docked molecule matches an experimental reference ligand.

### Covalent complex skipped

Vina/Smina conventional noncovalent docking is not a valid substitute for a
dedicated covalent-docking method. Select a noncovalent validation structure or
use an appropriate external covalent workflow.

### Run is slow

Reduce tautomer/conformer limits, requested modes, or exhaustiveness during
exploration. Use the CPU field to control engine cores. Preparation and variant
enumeration can dominate wall time even when the recorded docking timer is
short.

### Results differ across computers

Confirm engine version, random seed, CPU count, receptor preparation, input
order, dependency versions, and architecture. Small score or pose differences
can arise from external engines and numerical libraries; report compatible,
rather than bit-identical, replication where appropriate.
