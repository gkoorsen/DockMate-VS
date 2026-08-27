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
| OpenMM and PDBFixer | Protein repair | Optional but recommended | [OpenMM](https://docs.openmm.org/latest/userguide/application/01_getting_started.html#installing-openmm) / [PDBFixer](https://github.com/openmm/pdbfixer#installation) |
| Reduce and PROPKA | Hydrogen/protonation preparation | Optional pipeline components | [Reduce](https://github.com/rlabduke/reduce#building) / [PROPKA](https://propka.readthedocs.io/en/latest/installation.html) |
| rDock | Alternative docking engine | Optional | [rDock](https://rdock.github.io/installation/) |
| Meeko (`mk_prepare_receptor.py`) | Flexible-receptor PDBQT preparation | Optional | [Meeko](https://meeko.readthedocs.io/en/develop/installation.html) |
| PyMOL and LigPlot+ | Pose and interaction visualization | Optional | [PyMOL](https://pymol.org/) / [LigPlot+](https://www.ebi.ac.uk/thornton-srv/software/LigPlus/) |

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
common settings. For container execution, keep the campaign file, workbook, and
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

## 5. Protocol Development

### Select conditions

The protocol sweep can vary:

- Vina, Smina, and installed rDock engines;
- co-crystal box margins or fixed box dimensions;
- remove-all, retain-all, and selective crystallographic-water handling;
- exhaustiveness and random seed;
- baseline score ordering and supported Smina score-only methods.

Rescoring does not redock the ligand. It evaluates the selected scoring function
on every saved pose and compares baseline and rescored rankings.

### Ligand variants

Set explicit maximum tautomer and conformer counts. Large values multiply the
number of dockings. Begin with two tautomers and two conformers for exploration,
then justify any increase. Screening must select variants by score, not RMSD.

Version 0.1 applies RDKit largest-fragment selection and charge normalization,
then enumerates tautomers and conformers. It does not perform pH-aware ligand
ionization-state enumeration. Inspect prepared forms for charge-sensitive
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

## 6. Screening

### Freeze the protocol

Record the selected receptor, binding site, water treatment, engine, scoring,
box, exhaustiveness, random seed, number of modes, energy range, ligand-variant
limits, and CPU count before screening. Do not select a favourable receptor or
seed after examining screening labels.

### Run and resume

Remote structures are prefetched before docking. Once preflight succeeds, the
remaining preparation, docking, scoring, and analysis can proceed without
continuous network connectivity.

Progress is saved incrementally. Restart with the same input and compatible
settings to skip completed cases. Changed docking or scoring settings invalidate
reuse. Failed cases and missing output files are retried.

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

## 7. Results and inspection

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

## 8. Reproducibility checklist

Before publishing a campaign, retain:

- exact input workbook and its checksum;
- run manifest and software git commit;
- docking-engine and dependency versions;
- result CSV/JSON and generated summary;
- failures and excluded compounds;
- hardware, CPU count, campaign wall time, and random seed;
- data-source DOI, cleaning script, and redistribution terms.

## 9. Troubleshooting

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
