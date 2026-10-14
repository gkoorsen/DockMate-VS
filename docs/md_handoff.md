# Outputs and MD handoff

DockMate-VS writes docking coordinates, scores, preparation records, and analysis
tables. These provide starting structures for downstream molecular dynamics
(MD) and binding free-energy work. Force-field parameters, solvent/ions,
equilibration protocols, restraints, and simulation input decks must be prepared
with the downstream tools.

## File inventory

Paths below are relative to a run directory or an individual case directory.
Use each result's `output_file` to identify its selected ligand variant.

| File | Contents and interpretation |
| --- | --- |
| `redock_results.csv` | One case per row; identifiers, control labels, baseline/rescored scores, RMSDs, errors, and paths. Empty cells represent missing values. |
| `redock_results.json` | The same case records under `results`; unavailable values are null and adaptive attempts are nested records. |
| `protocol_development/protocol_development_results.csv` | One receptor/condition/ranking-method record per row; engine, box, water handling, sampling, seed, completion, and recovery metrics. |
| `redock_summary.json` / `.md` | Derived analyses; the JSON includes selected unknown-hit pose-quality records and their generated file paths. |
| `run_manifest.json` / `protocol_development_manifest.json` | Configuration, identifiers, software/tool versions, and provenance. |
| `redock_progress.json` | Restart checkpoints, using the case-result schema. |
| `receptor_prepared.pdbqt` | Rigid Vina/Smina receptor coordinates, docking atom types, and partial charges. Not an MD topology. |
| `receptor_prepared.pdb` | Prepared receptor coordinates in PDB format. Inspect the protein, waters, cofactors, metals, and chain selection for the intended MD model. |
| `crystal_ligand.pdb` | Extracted crystallographic reference; its coordinates are not the newly docked pose. |
| `ligand_variants/<variant>.sdf` | Prepared ligand molecular graph, bond orders, formal charges, and pre-docking coordinates. |
| `ligand_variants/<variant>.pdbqt` | Prepared ligand with docking atom types, partial charges, and torsion tree. |
| `variants/<variant>/docked.pdbqt` | Multi-model Vina/Smina pose coordinates and score remarks. |
| `variants/<variant>/docked.sd` | rDock molecular records and score properties. |
| `pose_quality/<id>/selected_pose_prepared_topology.sdf` | Selected docked coordinates using prepared SDF chemistry, when topology transfer succeeds. |
| `pose_quality/<id>/selected_pose.sdf` | Fallback reconstructed chemistry; inspect `posebusters_topology_source` and any topology error before reuse. |
| `pose_quality/<id>/selected_pose.pdb` | Selected ligand coordinates with residue name `ZQX`. PDB alone does not preserve ligand bond-order information. |
| `pose_quality/<id>/selected_pose_complex.pdb` | Receptor and selected ligand in the same coordinate frame, generated for contact analysis. |

Coordinates and RMSDs are in angstroms. Vina/Smina affinity scores are reported
in kcal/mol; rDock and CNN score fields have their own semantics. Baseline and
rescore values describe different scoring methods and may rank the same poses
differently. Record both the score source and the actual selected pose number.
`rmsd_status` distinguishes `ok`, `not_computable`, and `not_applicable`; a failed
RMSD is not evidence that docking itself failed.

## Worked coordinate handoff

1. Load a screening run in Results, or run `dockmate-vs report --run /path/to/run`.
   Pose-quality analysis considers the top five unknown compounds per structure;
   it does not generate these files for every saved pose. In
   `redock_summary.json`, locate the relevant `screening_pose_quality` record and
   record its `case_id`, `output_file`, one-based `pose_index`, and score source.
2. Use the path in `posebusters_ligand_sdf`. Prefer a record whose
   `posebusters_topology_source` is `prepared_sdf` or `prepared_sdf_mcs`, with no
   topology error. This file contains docked coordinates. The original
   `ligand_variants/*.sdf` supplies chemistry but contains pre-docking coordinates.
3. Inspect this ligand together with `selected_pose_complex.pdb` and the case's
   receptor. Retain the SDF for downstream ligand parameterization. If a PDB
   ligand is needed, convert the selected SDF with Open Babel:

   ```bash
   obabel -isdf "/path/to/selected_pose_prepared_topology.sdf" \
     -opdb -O "/path/to/ligand_for_md.pdb"
   ```

   Format conversion is not parameterization. Check the atom count, elements,
   stereochemistry, bond orders, formal charge, hydrogens, and coordinates after
   conversion. Preserve the receptor-relative pose during any atom renaming or
   coordinate reordering.
4. Use `selected_pose_complex.pdb` as a coordinate starting point when available.
   Its ligand residue is `ZQX`. If assembling another pose manually, combine its
   ligand HETATM records with the receptor using a molecular structure reader/
   writer; place the final `END` after both components and avoid duplicate serial
   or residue identifiers. Do not independently align the ligand or regenerate
   its 3D coordinates. Keep the selected SDF alongside the complex PDB.

## Exporting a saved pose

The `export` command automates the coordinate/topology handoff without docking
again, minimising the pose, or changing saved results:

```bash
dockmate-vs export --run "/path/to/run" --list-cases
dockmate-vs export --run "/path/to/run" --case "CASE_ID" --pose 1 \
  --format sdf --output "/path/to/exports/pose1.sdf"
dockmate-vs export --run "/path/to/run" --case "CASE_ID" --pose 1 \
  --format complex-pdb --output "/path/to/exports/complex.pdb"
```

`--pose` is **one-based docking-file order**, not a rescored rank or the pose with
the lowest RMSD. Consult the results table before choosing it. Older protocol
results without case IDs use their case-directory names. Both screening and
protocol run folders are supported, including relocated copies; `--list-cases`
marks missing outputs. Ambiguous case IDs and existing destination files are
rejected. Without `--output`, new files go into the current working directory.

For PDBQT, the matching prepared variant SDF **and** PDBQT are required. Atom
coordinates from that preparation pair establish the mapping, and preserved
PDBQT atom identities transfer the docked coordinates onto the SDF topology.
Mismatches fail rather than infer new bonds. A regenerated SDF whose coordinates
do not match the original prepared PDBQT is insufficient for this strict route.
For SDF/SD output, export verifies chemical equivalence and maps the saved pose
onto its prepared topology. Historical rDock output uses the actual input SDF
inside the docking directory when available.

Mapped atom coordinates are retained. Ordinary hydrogens omitted from PDBQT are
rebuilt with RDKit without minimisation; unsupported/missing nonstandard atoms
cause an error. A JSON sidecar records source files, canonical SMILES, pose
number and the count of rebuilt hydrogens. Review hydrogen geometry downstream.
Export also checks atom count and canonical SMILES through an SDF round trip
before writing files; inconsistent chemistry or stereochemistry is rejected.

`complex-pdb` also requires `receptor_prepared.pdb`. It adds the ligand as residue
`LIG` on an unused chain with unique atom serials; `--ligand-resname DRG` changes
the name. It writes a companion `complex.ligand.sdf` retaining bond orders,
because PDB alone does not preserve the complete ligand chemistry.
These are **geometry/topology exports, not simulation-ready systems**: no
force-field parameters, solvation, ion placement, MD decks or free-energy
protocols are generated.

For manual inspection of older PDBQT-only output, basic conversion can extract
model 1:

```bash
obabel -ipdbqt "/path/to/docked.pdbqt" -f 1 -l 1 \
  -osdf -O "/path/to/pose1_inferred.sdf"
```

Change both model selectors to the desired pose number. PDBQT lacks a complete
chemical graph, so this route can infer incorrect bond orders or aromaticity.
Validate the inferred molecule against the matching prepared SDF and transfer
docked coordinates with an explicit atom mapping. Do not use raw conversion as
proof of chemical equivalence. Restore the matching preparation files to use
the verified export route above.

## Downstream preparation

| Destination | Additional work and required inputs |
| --- | --- |
| GROMACS | Prepare protein topology and ligand parameters compatible with the selected force field; combine coordinate files and include the ligand `.itp` in the system `.top`. Check atom order, molecule counts, and charge; then solvate, add ions, and prepare `.mdp`/`.tpr` inputs. `pdb2gmx` does not parameterize arbitrary drug-like ligands. |
| NAMD | Build a consistent PSF/PDB pair plus force-field parameter files, including the ligand. Set solvent, periodic boundaries, and simulation configuration. A PDB without the corresponding molecular topology is insufficient. |
| AMBER | Generate ligand parameters, for example with Antechamber/GAFF and `parmchk2`, then assemble with a compatible protein force field using LEaP. Write the system topology (`prmtop`) and coordinates (`inpcrd`/`rst7`); review charge and missing parameters before simulation. |
| BFEE2 | Supply the parameterized complex in a supported topology/coordinate combination, such as PSF/PDB with parameters or AMBER topology/coordinates. Prepare the system and define the receptor/ligand selections and free-energy protocol required by that version. |
| BAT.py | Follow its protein/ligand preparation and pose-PDB naming conventions. Its workflow includes ligand parameterization and simulation preparation; provide the required settings, selections, and external dependencies. DockMate-VS does not generate BAT configuration files. |

Protonation, tautomer state, cofactors, retained waters, and the receptor model
must be reviewed for the downstream conditions. DockMate-VS preserves supplied
formal charges by default, but does not enumerate ionization states at a chosen
pH. Docking partial charges should not be substituted for MD force-field charges.

Official preparation guides:

- [GROMACS pdb2gmx](https://manual.gromacs.org/current/onlinehelp/gmx-pdb2gmx.html)
- [NAMD system configuration](https://www.ks.uiuc.edu/Research/namd/2.14/ug/node12.html)
- [AMBER ligand preparation tutorial](https://ambermd.org/tutorials/basic/tutorial4b/)
- [BFEE project and versioned releases](https://github.com/fhh2626/BFEE2)
- [BAT.py user guide and docking examples](https://github.com/GHeinzelmann/BAT.py)
