# Supplementary Material S1: ESR1 reconstruction and analysis

This supplement supports the DockMate-VS SoftwareX case study for ESR1
1XP1/AIH: the deterministic 80-active/800-inactive screen and the protocol sweep.
It distributes source identifiers, calculated results, settings and analysis
code, **not** the LIT-PCBA compound structures or original input workbooks.
No SMILES, molecular coordinate files, docking poses or third-party binaries
are included. Source data must be obtained separately from their providers.

## Contents

- `selection.csv`: the 880 source identifiers in the exact screening order.
  These are LIT-PCBA source identifiers (PubChem substance IDs), not PubChem CIDs.
- `source.json`: upstream location, archive SHA-256, cleaning counts and expected
  checksum of the reconstructed compound CSV.
- `campaigns/`: historical screening and protocol-sweep settings; host-specific
  paths are replaced with portable references. These are documentation of the
  original runs, not configurations guaranteed to run in the current release.
- `results/redock_results.csv`: per-compound scores, computed descriptors,
  activity labels and execution statistics. Original machine paths and free-text
  fields are excluded. These are result records, not raw engine logs or poses.
- `protocol_development/`: structure-free sweep and seed-check result tables.
- `publication/tables/`: reference statistical outputs for comparison.
- `provenance/`: original preparation/analysis scripts, historical environment
  information and the checksums of the frozen source evidence.
- `reconstruct.py`: reconstructs the cleaned subset locally from an archive you
  downloaded separately, verifying both selection order and exact CSV checksum.
- `verify.py`: checks package integrity and compares reconstructed/reanalysed
  outputs with the archived evidence.

## Environment

The original analyses used Python 3.9.15. In a separate Python 3.9 environment,
install `analysis-requirements.txt` for analysis. Reconstruction also requires
RDKit 2022.03.5 (for example, from conda-forge). Full historical environment
information is in `provenance/provenance.json`; this is not a complete docking
environment lockfile. Different library versions may change canonicalization,
sampling, model fitting or figures. The validation commands detect differences
in the recorded numerical outputs rather than assuming equivalence.

Run the following commands from the extracted supplement directory.

### 1. Verify the archive contents

```sh
python verify.py
```

### 2. Reconstruct the screening inputs

Obtain `LIT-PCBA_full.tar.gz` from the original provider:
https://lab.drugdesign.unistra.fr/datasets/lit-pcba/

The exact download URL and SHA-256 are recorded in `source.json`. A newer or
otherwise different archive is rejected, even if it has the same filename.
The script does not download anything or extract the full archive. It reads
only the two ESR1 active/inactive source files into a temporary directory.

```sh
python reconstruct.py --archive /path/to/LIT-PCBA_full.tar.gz --output reconstructed
python verify.py --reconstructed reconstructed
```

This recreates the same canonicalization, within-class deduplication, conflict
exclusion and seed-42 sampling used in the paper. The output CSV must match the
original CSV byte-for-byte. The screening workbook has the same job rows, but
its metadata and ZIP timestamps differ; workbook byte identity is not claimed.
A one-row native-redocking workbook is also generated. Reconstructed structures
and workbooks are local outputs, not part of the distributed supplement.

The receptor and reference ligand can be obtained separately using PDB accession
1XP1, chain A, ligand AIH. Checksums of the original downloaded PDB and ideal SDF
are recorded in `source.json`. No full docking rerun is performed by this script.

### 3. Reproduce the statistics without downloading structures

```sh
python provenance/analyze_esr1_seed42.py --frozen . --manuscript-root reanalysis
python verify.py --reanalysis reanalysis
```

This reruns the 10,000-replicate stratified bootstrap, enrichment calculations,
descriptor diagnostics and property-only model on the distributed result table.
It writes tables and illustrative analytical figures under `reanalysis/`.
Those figures are not the final GUI screenshots used in the manuscript.
The expected metrics include ROC AUC 0.695, average precision 0.174 and mean
property-only cross-validated ROC AUC 0.743 (rounded values).

DockMate-VS 0.1.1 can load `results/redock_results.csv` for the screening charts
and `protocol_development/protocol_development_results.csv` for the protocol
dashboard. Pose viewing is unavailable because pose files are not supplied.
Use Top-5 and a 2 A threshold for the protocol figure; the 162 saved rows
represent 126 conditions after collapsing repeated rDock exhaustiveness values.

## Statistical methods

Screening ranks use negative docking score so larger values rank first.
Confidence intervals are the 2.5th and 97.5th percentiles of 10,000 bootstrap
replicates, sampling actives and inactives separately with replacement, retaining
their original class sizes and using NumPy `default_rng(42)`. Expected enrichment
accounts for compounds tied at each cutoff; the associated minimum/maximum
columns are tie bounds, not bootstrap confidence intervals.

The property-only classifier uses molecular weight, logP, topological polar
surface area, rotatable-bond count and ligand charge as saved by the historical
preparation workflow. StandardScaler and logistic regression are fitted together
inside each training fold. Logistic regression uses L2 regularization, C=1,
the lbfgs solver, max_iter=2000, random_state=42 and otherwise scikit-learn 1.1.1
defaults. Evaluation uses 10 stratified, shuffled folds with random_state=42.
The reported AUC is the unweighted mean of the ten fold AUCs; the reported SD is
the sample standard deviation. This diagnostic does not establish prospective
or receptor-specific prediction performance.

## Historical docking versus current software

The docking runs predate version 0.1.1 and used the former `adaptive` ligand
variant-selection mode and legacy charge handling. Version 0.1.1 removed that
mode and changed the default charge policy. Do not substitute `all` or claim
that rerunning 0.1.1 recreates the historical scores. The recorded commits and
environment identify the historical development period, but do not establish
a bitwise-exact snapshot of every tool or intervening working-tree change.
The analysis is reproducible directly from the archived result records without
rerunning docking; a new docking campaign requires its own provenance.

Absolute paths in the original manifests are not portable and have been removed
from this supplement. Unmodified calculation values and selected numerical
tables are retained; packaging provenance documents these transformations.
The complete original LIT-PCBA archive and per-case docking outputs are excluded.

## Attribution and publication

Cite LIT-PCBA: https://doi.org/10.1021/acs.jcim.0c00155 and the Protein Data Bank:
https://doi.org/10.1093/nar/28.1.235. DockMate-VS 0.1.1 is identified by
https://doi.org/10.5281/zenodo.22685196.

`LICENSE` applies to the accompanying DockMate-VS code, not to third-party data.
This reconstruction-based distribution does not assert or grant redistribution
rights over LIT-PCBA. Activity labels and source identifiers are attributed to
LIT-PCBA; calculated scores and diagnostics are outputs of the case study.

This ZIP is distributed through the DockMate-VS GitHub repository under
`supplementary/packages/DockMate-VS_ESR1_Supplement_S1.zip`, separately from the
Zenodo software archive. Cite a commit-specific GitHub permalink to identify
the exact package used. The S1 filename is a package identifier, not a claim
that the journal has assigned a supplement number.
