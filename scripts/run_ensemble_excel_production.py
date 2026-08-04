#!/usr/bin/env python3
"""
Complete ensemble validation from Excel template - PRODUCTION VERSION.

Uses your existing docking platform infrastructure:
- LigandPreparation from ligand.py
- SminaDockingEngine from smina.py
- Proper PDBQT handling

Usage:
    python run_ensemble_excel_production.py \\
        --excel anathi_template_complete.xlsx \\
        --pdb-dir structures/ \\
        --output results/
"""

import argparse
import sys
from pathlib import Path
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from loguru import logger

# Add to path if needed
sys.path.insert(0, str(Path(__file__).parent.parent))

from docking_platform_gui.experiments.ensemble_validation import EnsembleValidationExperiment
from docking_platform_gui.analysis.enrichment_analysis import EnrichmentAnalysis, BootstrapValidation
from docking_platform_gui.experiments.docking_adapter_final import (
    DockingManagerAdapter,
    CachedDockingAdapter
)


def read_excel_and_validate(excel_file):
    """Read Excel and validate it has all required data."""
    logger.info(f"Reading Excel: {excel_file}")
    
    df = pd.read_excel(excel_file)
    df.columns = df.columns.str.strip()
    
    # Validate required columns
    required = ['Protein', 'Target_Ligand', 'PDB_ID', 'Ligand', 'SMILES', 'decoy SMILES']
    missing = [col for col in required if col not in df.columns]
    if missing:
        logger.error(f"Missing required columns: {missing}")
        logger.error("Please run prepare_excel_with_smiles.py first!")
        return None
    
    # Check for missing SMILES
    missing_smiles = df[df['SMILES'].isna()]
    if len(missing_smiles) > 0:
        logger.error(f"⚠️  {len(missing_smiles)} entries missing SMILES:")
        for idx, row in missing_smiles.iterrows():
            logger.error(f"  Row {idx}: {row['Target_Ligand']}")
        logger.error("\nPlease add SMILES for all compounds!")
        return None
    
    # Every PDB_ID row is pooled into ONE ensemble sharing ONE grid box, so a template
    # spanning several proteins is not a valid input - it would dock each compound
    # against unrelated receptors. Use one workbook per target.
    ens_rows = df[df['PDB_ID'].notna()]
    if 'Protein' in df.columns and len(ens_rows) > 0:
        proteins = sorted(set(ens_rows['Protein'].dropna().astype(str).str.strip()) - {''})
        if len(proteins) > 1:
            logger.error(f"⚠️  Structure rows span {len(proteins)} proteins: {', '.join(proteins)}")
            for prot in proteins:
                ids = ens_rows[ens_rows['Protein'].astype(str).str.strip() == prot]['PDB_ID']
                logger.error(f"    {prot}: {', '.join(str(p) for p in ids)}")
            logger.error("")
            logger.error("All structures with a PDB_ID are pooled into ONE ensemble sharing ONE")
            logger.error("grid box. Mixing proteins would dock every compound against unrelated")
            logger.error("receptors. Split this into one workbook per target and run each separately.")
            return None

    logger.info(f"✅ Excel validation passed")
    logger.info(f"  - {len(df)} total compounds")
    logger.info(f"  - {df['PDB_ID'].notna().sum()} with PDB structures (ensemble)")
    logger.info(f"  - {df['PDB_ID'].isna().sum()} unknowns to screen")
    logger.info(f"  - {df['decoy SMILES'].notna().sum()} with decoys")
    
    return df


def create_sdf_files(df, output_dir):
    """Create actives.sdf and decoys.sdf from Excel data."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    logger.info("\nGenerating 3D structures from SMILES...")
    
    # Actives
    actives_file = output_dir / 'actives.sdf'
    writer = Chem.SDWriter(str(actives_file))
    
    actives_count = 0
    for idx, row in df.iterrows():
        smiles = row['SMILES']
        mol = Chem.MolFromSmiles(smiles)
        
        if mol is None:
            logger.warning(f"  Invalid SMILES for {row['Target_Ligand']}: {smiles}")
            continue
        
        # Generate 3D
        mol = Chem.AddHs(mol)
        result = AllChem.EmbedMolecule(mol, randomSeed=42)
        if result != 0:
            logger.warning(f"  Could not generate 3D for {row['Target_Ligand']}")
            continue
        
        AllChem.UFFOptimizeMolecule(mol)
        
        # Set properties
        compound_id = f"active_{idx}_{row['Target_Ligand'].replace(' ', '_')}"
        mol.SetProp('_Name', compound_id)
        mol.SetProp('compound_name', row['Target_Ligand'])
        mol.SetProp('original_smiles', smiles)
        
        if pd.notna(row['PDB_ID']):
            mol.SetProp('has_structure', 'true')
            mol.SetProp('pdb_id', row['PDB_ID'])
        
        writer.write(mol)
        actives_count += 1
    
    writer.close()
    logger.info(f"  ✅ Created: {actives_file} ({actives_count} compounds)")
    
    # Decoys
    decoys_file = output_dir / 'decoys.sdf'
    writer = Chem.SDWriter(str(decoys_file))
    
    decoys_count = 0
    def _split_list(value):
        if pd.isna(value):
            return []
        if isinstance(value, str):
            return [v.strip() for v in value.split(",") if v.strip()]
        text = str(value).strip()
        return [text] if text else []

    for idx, row in df.iterrows():
        decoy_smiles_list = _split_list(row['decoy SMILES'])
        if not decoy_smiles_list:
            continue

        decoy_name_list = _split_list(row.get('decoy compound', ''))

        for j, smiles in enumerate(decoy_smiles_list, 1):
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                logger.warning(f"  Invalid decoy SMILES for row {idx}: {smiles}")
                continue

            mol = Chem.AddHs(mol)
            result = AllChem.EmbedMolecule(mol, randomSeed=42)
            if result != 0:
                continue

            AllChem.UFFOptimizeMolecule(mol)

            decoy_id = decoy_name_list[j - 1] if j - 1 < len(decoy_name_list) else f"decoy_{idx}_{j}"
            mol.SetProp('_Name', f"decoy_{idx}_{decoy_id}")
            mol.SetProp('paired_with', row['Target_Ligand'])
            mol.SetProp('original_smiles', smiles)

            writer.write(mol)
            decoys_count += 1
    
    writer.close()
    logger.info(f"  ✅ Created: {decoys_file} ({decoys_count} compounds)")
    
    return actives_file, decoys_file


def get_ensemble_structures(df, pdb_dir):
    """Get paths to ensemble PDB structures (must be PDBQT format)."""
    pdb_dir = Path(pdb_dir)
    
    ensemble_df = df[df['PDB_ID'].notna()].copy()
    
    structures = []
    logger.info("\nEnsemble structures:")
    
    for idx, row in ensemble_df.iterrows():
        pdb_id = row['PDB_ID']
        
        # Try to find PDBQT file (required for docking)
        possible_names = [
            f"{pdb_id}_prepared.pdbqt",
            f"{pdb_id}_prep.pdbqt",
            f"{pdb_id}.pdbqt",
            f"{pdb_id}_receptor.pdbqt"
        ]
        
        found = None
        for name in possible_names:
            path = pdb_dir / name
            if path.exists():
                found = str(path)
                break
        
        if found is None:
            logger.error(f"  ✗ NOT FOUND: {pdb_id}.pdbqt")
            logger.error(f"    Tried: {possible_names}")
            logger.error(f"    Please convert PDB to PDBQT format:")
            logger.error(f"      prepare_receptor4.py -r {pdb_id}.pdb -o {pdb_id}_prepared.pdbqt")
            return None
        
        structures.append(found)
        logger.info(f"  ✓ {pdb_id} - {row['Target_Ligand']} ({Path(found).name})")
    
    return structures


def _pdb_fallback_candidates(pdbqt_path: Path):
    """Return likely sibling PDB files for a receptor PDBQT path."""
    name = pdbqt_path.name
    base = None
    for suffix in ("_prepared.pdbqt", "_prep.pdbqt", "_receptor.pdbqt"):
        if name.endswith(suffix):
            base = name[:-len(suffix)]
            break
    if base is None and name.endswith(".pdbqt"):
        base = name[:-6]

    candidates = []
    if base:
        candidates.append(pdbqt_path.with_name(f"{base}.pdb"))
    candidates.append(pdbqt_path.with_suffix(".pdb"))
    return candidates


def _structure_source_for_grid(pdbqt_path):
    """
    Locate the PDB/mmCIF that carries the ligand for a receptor PDBQT.

    Receptor PDBQT files produced by prepare_receptor4.py have the ligand stripped, so
    site geometry has to come from the sibling PDB. Returns None if none is present.
    """
    pdbqt_file = Path(pdbqt_path)
    candidates = list(_pdb_fallback_candidates(pdbqt_file))
    candidates += [pdbqt_file.with_suffix(".cif"), pdbqt_file.with_suffix(".mmcif")]
    base = pdbqt_file.name
    for suffix in ("_prepared.pdbqt", "_prep.pdbqt", "_receptor.pdbqt", ".pdbqt"):
        if base.endswith(suffix):
            stem = base[:-len(suffix)]
            for ext in (".pdb", ".cif", ".mmcif"):
                candidates.append(pdbqt_file.with_name(f"{stem}_aligned{ext}"))
                candidates.append(pdbqt_file.with_name(f"{stem}{ext}"))
            break
    return next((p for p in candidates if p.exists()), None)


def check_ensemble_coframed(structures, grid_center, grid_size, comp_ids=None):
    """
    Verify every ensemble member's binding site falls inside the grid box.

    This is the gate that catches an unaligned ensemble. One grid box is shared by all
    receptors, which is only valid if they occupy a common coordinate frame. PDB entries
    generally do NOT: on a real CK2-alpha ensemble the pairwise ligand-centroid
    distances were 396.5, 387.1 and 88.2 A, and the mean of those centres sat in bulk
    solvent. Docking still completes and returns plausible-looking scores.

    Returns True when the ensemble is usable, False when it is not.
    """
    sources = {}
    for struct in structures:
        src = _structure_source_for_grid(struct)
        if src is not None:
            sources[struct] = src

    if not sources:
        logger.warning(
            "\nCannot verify ensemble coordinate frames: no sibling PDB/CIF found "
            "next to the receptor PDBQT files."
        )
        logger.warning(
            "  Place the corresponding .pdb files alongside them so the grid can be checked."
        )
        return True

    try:
        from docking_platform_gui.utils.structure_alignment import validate_grid
    except ImportError:
        logger.warning("\ngemmi not installed - skipping ensemble frame check.")
        logger.warning("  Install with: pip install gemmi")
        return True

    logger.info("\nVerifying ensemble structures share a coordinate frame...")
    ok, rows = validate_grid(
        list(sources.values()), grid_center, grid_size, comp_ids=comp_ids
    )

    for row in rows:
        label = Path(row["structure"]).name
        if row["status"] == "no_ligand":
            logger.info(f"  {label}: no ligand found (skipped)")
        elif row["status"] == "OUTSIDE":
            logger.error(
                f"  ✗ {label}: site ({row['comp_id']}) is {row['distance']} A from the "
                f"grid centre - OUTSIDE the box"
            )
        elif row["status"] == "near_edge":
            logger.warning(
                f"  ! {label}: site ({row['comp_id']}) is {row['distance']} A from the "
                f"grid centre - close to the box face"
            )
        else:
            logger.info(f"  ✓ {label}: site ({row['comp_id']}) {row['distance']} A from centre")

    if not ok:
        logger.error("")
        logger.error("=" * 70)
        logger.error("ABORTING: ensemble structures are not in a common coordinate frame.")
        logger.error("=" * 70)
        logger.error("One grid box is shared by every receptor, so at least one structure")
        logger.error("would be docked into empty solvent. The run would complete and report")
        logger.error("scores that mean nothing.")
        logger.error("")
        logger.error("Fix: superpose the structures first, then prepare receptors from the")
        logger.error("aligned copies:")
        logger.error("")
        logger.error("  python -m docking_platform_gui.utils.align_structures \\")
        logger.error("      --structures 1ABC.pdb 2DEF.pdb 3GHI.pdb \\")
        logger.error("      --out-dir structures_aligned")
        logger.error("")
        logger.error("Then re-run with --pdb-dir pointing at the aligned receptors.")
        logger.error("To proceed anyway (not recommended), pass --skip-frame-check.")

    return ok


def calculate_consensus_center(structures, comp_ids=None):
    """
    Consensus grid center from the binding-site ligand of each ensemble member.

    Each structure contributes ONE ligand copy, chosen as the largest non-additive
    component (or an explicit comp_id), and for multimers the copy nearest the reference
    site. Averaging every HETATM instead - as this function previously did - lets
    crystallisation additives and glycans dominate: in 6TGU ethylene glycol contributes
    70 atoms against the inhibitor's 66 and pulls the centre ~17 A off-site; in 1MX1 the
    NAG/SIA/NDG glycans pull it off the steroid pocket entirely.
    """
    logger.info("\nCalculating consensus grid center...")

    sources = {}
    for struct in structures:
        src = _structure_source_for_grid(struct)
        if src is not None:
            sources[struct] = src

    if sources:
        try:
            from docking_platform_gui.utils.structure_alignment import derive_grid

            spec = derive_grid(list(sources.values()), comp_ids=comp_ids)
            for site in spec.sites:
                logger.info(
                    f"  {site.path.name}: {site.comp_id} "
                    f"({site.center[0]:.2f}, {site.center[1]:.2f}, {site.center[2]:.2f})"
                    + (f"  [copy 1 of {site.n_copies}]" if site.n_copies > 1 else "")
                )
            logger.info(
                f"\n  Consensus: ({spec.center[0]:.2f}, {spec.center[1]:.2f}, {spec.center[2]:.2f})"
            )
            logger.info(f"  Std dev: {tuple(spec.spread_std)}")
            logger.info(f"  Max deviation from consensus: {spec.max_deviation} A")
            logger.info(
                f"  Suggested box size: {spec.size[0]:.0f} A "
                f"(pass --grid-size to override)"
            )
            if spec.max_deviation > 5.0:
                logger.warning(
                    f"  Sites disagree by up to {spec.max_deviation} A - the structures may "
                    "not be superposed, or one ligand may occupy a different pocket."
                )
            return list(spec.center)
        except ImportError:
            logger.warning("  gemmi not installed - falling back to raw HETATM averaging.")
            logger.warning("  Install gemmi for correct site selection: pip install gemmi")
        except ValueError as exc:
            logger.error(f"  {exc}")
            return None

    centers = []
    exclude_resnames = {"HOH", "WAT", "SOL", "TIP3", "CL", "NA", "K", "CA", "MG", "ZN", "FE", "MN", "BR", "I", "SO4", "PO4", "GOL"}
    
    for pdbqt_path in structures:
        # Extract ligand coordinates from PDBQT
        # PDBQT files have HETATM lines for ligands
        coords = []

        with open(pdbqt_path, 'r') as f:
            for line in f:
                if line.startswith('HETATM'):
                    try:
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        coords.append([x, y, z])
                    except ValueError:
                        continue

        if not coords:
            pdbqt_file = Path(pdbqt_path)
            pdb_path = next((p for p in _pdb_fallback_candidates(pdbqt_file) if p.exists()), None)
            if pdb_path:
                res_counts = {}
                with open(pdb_path, 'r') as f:
                    for line in f:
                        if not line.startswith('HETATM'):
                            continue
                        resname = line[17:20].strip()
                        if resname in exclude_resnames:
                            continue
                        try:
                            x = float(line[30:38].strip())
                            y = float(line[38:46].strip())
                            z = float(line[46:54].strip())
                        except ValueError:
                            continue
                        coords.append([x, y, z])
                        res_counts[resname] = res_counts.get(resname, 0) + 1
                if coords and res_counts:
                    res_list = ", ".join(sorted(res_counts))
                    logger.info(f"  {pdbqt_file.name}: using HETATM from {pdb_path.name} ({res_list})")
        
        if coords:
            center = np.mean(coords, axis=0)
            centers.append(center)
            logger.info(f"  {Path(pdbqt_path).name}: ({center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f})")
        else:
            logger.warning(f"  Could not extract ligand coords from {Path(pdbqt_path).name}")
    
    if len(centers) == 0:
        logger.error("Could not extract any ligand positions.")
        return None
    
    consensus = np.mean(centers, axis=0)
    std = np.std(centers, axis=0)
    
    logger.info(f"\n  Consensus: ({consensus[0]:.2f}, {consensus[1]:.2f}, {consensus[2]:.2f})")
    logger.info(f"  Std dev: ({std[0]:.2f}, {std[1]:.2f}, {std[2]:.2f})")

    return consensus.tolist()


def calculate_protein_centroid(structures):
    """Calculate consensus center from protein ATOM coordinates (apo-friendly fallback)."""
    logger.info("\nCalculating consensus grid center from protein centroid...")

    centers = []
    for structure_path in structures:
        structure_file = Path(structure_path)
        coords = []

        with open(structure_file, 'r') as f:
            for line in f:
                if not line.startswith('ATOM'):
                    continue
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    coords.append([x, y, z])
                except ValueError:
                    continue

        if not coords:
            pdb_path = next((p for p in _pdb_fallback_candidates(structure_file) if p.exists()), None)
            if pdb_path:
                with open(pdb_path, 'r') as f:
                    for line in f:
                        if not line.startswith('ATOM'):
                            continue
                        try:
                            x = float(line[30:38].strip())
                            y = float(line[38:46].strip())
                            z = float(line[46:54].strip())
                            coords.append([x, y, z])
                        except ValueError:
                            continue

        if not coords:
            logger.warning(f"  Could not extract protein ATOM coordinates from {structure_file.name}")
            continue

        center = np.mean(coords, axis=0)
        centers.append(center)
        logger.info(f"  {structure_file.name}: ({center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f})")

    if not centers:
        logger.error("Could not extract protein centroids from any structure.")
        return None

    consensus = np.mean(centers, axis=0)
    std = np.std(centers, axis=0)
    logger.info(f"\n  Consensus protein centroid: ({consensus[0]:.2f}, {consensus[1]:.2f}, {consensus[2]:.2f})")
    logger.info(f"  Std dev: ({std[0]:.2f}, {std[1]:.2f}, {std[2]:.2f})")
    logger.warning("Protein centroid is a fallback. Prefer explicit known pocket center when possible.")
    return consensus.tolist()


def main():
    parser = argparse.ArgumentParser(
        description="Run ensemble validation from completed Excel template (PRODUCTION)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        '--excel', required=True,
        help='Excel template WITH SMILES column filled'
    )
    parser.add_argument(
        '--pdb-dir', required=True,
        help='Directory with prepared PDBQT structures'
    )
    parser.add_argument(
        '--output', default='ensemble_validation_results',
        help='Output directory'
    )
    parser.add_argument(
        '--grid-size', nargs=3, type=float, default=[22, 22, 22],
        metavar=('X', 'Y', 'Z'),
        help='Grid box size (Angstroms)'
    )
    parser.add_argument(
        '--grid-center', nargs=3, type=float, default=None,
        metavar=('X', 'Y', 'Z'),
        help='Explicit grid center. Recommended for apo targets without co-crystal ligand.'
    )
    parser.add_argument(
        '--center-mode',
        choices=['ligand', 'protein'],
        default='ligand',
        help='How to infer center when --grid-center is not provided.'
    )
    parser.add_argument(
        '--bootstrap', type=int, default=1000,
        help='Bootstrap iterations for statistical validation'
    )
    parser.add_argument(
        '--exhaustiveness', type=int, default=16,
        help='Docking exhaustiveness (higher = more thorough)'
    )
    parser.add_argument(
        '--cpu', type=int, default=4,
        help='Number of CPU cores'
    )
    parser.add_argument(
        '--use-cache', action='store_true',
        help='Enable ligand caching (speeds up re-runs)'
    )
    parser.add_argument(
        '--smina-binary', default=None,
        help='Path to smina executable (uses system PATH if not specified)'
    )
    parser.add_argument(
        '--compare-all-singles', action='store_true',
        help='Compare ensemble against each single structure (adds summary CSV)'
    )
    parser.add_argument(
        '--skip-frame-check', action='store_true',
        help='Skip verifying that ensemble structures share a coordinate frame. '
             'Only use this if you have already superposed them by other means - '
             'an unaligned ensemble docks into solvent and returns meaningless scores.'
    )
    
    args = parser.parse_args()
    
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Setup logging
    logger.add(
        output_dir / "experiment.log",
        level="INFO",
        format="{time:YYYY-MM-DD HH:mm:ss} | {level} | {message}"
    )
    
    logger.info("="*70)
    logger.info("ENSEMBLE VALIDATION FROM EXCEL - PRODUCTION VERSION")
    logger.info("="*70)
    logger.info(f"Started: {pd.Timestamp.now()}")
    logger.info(f"Using your existing docking platform infrastructure:")
    logger.info(f"  - LigandPreparation (SMILES → PDBQT)")
    logger.info(f"  - SminaDockingEngine")
    logger.info(f"  - DockingResult/DockingPose")
    
    # 1. Read and validate Excel
    df = read_excel_and_validate(args.excel)
    if df is None:
        return 1
    
    # 2. Create SDF files
    try:
        actives_file, decoys_file = create_sdf_files(df, output_dir / 'compounds')
    except Exception as e:
        logger.exception(f"Failed to create compound files: {e}")
        return 1
    
    # 3. Get ensemble structures
    ensemble_structures = get_ensemble_structures(df, args.pdb_dir)
    if ensemble_structures is None:
        return 1
    
    # 4. Calculate grid center
    if args.grid_center is not None:
        grid_center = [float(v) for v in args.grid_center]
        logger.info("\nUsing explicit grid center from CLI.")
    elif args.center_mode == "protein":
        grid_center = calculate_protein_centroid(ensemble_structures)
    else:
        grid_center = calculate_consensus_center(ensemble_structures)

    if grid_center is None:
        logger.error("Failed to determine grid center.")
        logger.error("Provide --grid-center X Y Z, or set --center-mode protein for apo fallback.")
        return 1

    # 4b. Verify the shared grid box is actually valid for every receptor.
    # One grid is passed to all structures, so they must share a coordinate frame.
    if len(ensemble_structures) > 1 and not args.skip_frame_check:
        comp_ids = None
        if 'Ligand' in df.columns:
            ens_rows = df[df['PDB_ID'].notna()]
            comp_ids = {}
            for _, row in ens_rows.iterrows():
                lig = row.get('Ligand')
                if pd.notna(lig) and str(lig).strip():
                    pdb_id = str(row['PDB_ID']).strip()
                    comp_ids[pdb_id] = str(lig).strip()
                    comp_ids[f"{pdb_id}_aligned"] = str(lig).strip()
                    comp_ids[f"{pdb_id}_prepared"] = str(lig).strip()
        if not check_ensemble_coframed(
            ensemble_structures, grid_center, args.grid_size, comp_ids=comp_ids
        ):
            return 1
    elif args.skip_frame_check:
        logger.warning("\n! Ensemble coordinate-frame check SKIPPED (--skip-frame-check).")
    
    # 5. Setup experiment
    logger.info("\n" + "="*70)
    logger.info("EXPERIMENT CONFIGURATION")
    logger.info("="*70)
    
    single_structure = ensemble_structures[0]  # Use first as single
    
    logger.info(f"Single structure: {Path(single_structure).name}")
    logger.info(f"Ensemble: {len(ensemble_structures)} structures")
    for i, struct in enumerate(ensemble_structures, 1):
        logger.info(f"  {i}. {Path(struct).name}")
    
    grid_params = {
        'center': grid_center,
        'size': args.grid_size
    }
    
    logger.info(f"\nGrid parameters:")
    logger.info(f"  Center: {grid_params['center']}")
    logger.info(f"  Size: {grid_params['size']}")
    logger.info(f"\nDocking settings:")
    logger.info(f"  Exhaustiveness: {args.exhaustiveness}")
    logger.info(f"  CPU cores: {args.cpu}")
    logger.info(f"  Ligand caching: {'enabled' if args.use_cache else 'disabled'}")
    logger.info(f"\nBootstrap iterations: {args.bootstrap}")
    
    # 6. Create docking manager
    logger.info("\nInitializing docking manager...")
    
    if args.use_cache:
        docking_manager = CachedDockingAdapter(
            cache_dir=str(output_dir / 'ligand_cache'),
            smina_binary=args.smina_binary,
            exhaustiveness=args.exhaustiveness,
            cpu=args.cpu
        )
    else:
        docking_manager = DockingManagerAdapter(
            smina_binary=args.smina_binary,
            exhaustiveness=args.exhaustiveness,
            cpu=args.cpu
        )
    
    # 7. Create experiment runner
    experiment = EnsembleValidationExperiment(
        docking_manager=docking_manager,
        output_dir=args.output
    )
    
    # 8. Load benchmark dataset
    logger.info("\nLoading benchmark dataset...")
    dataset = experiment.load_benchmark_dataset(
        actives_file=str(actives_file),
        decoys_file=str(decoys_file)
    )
    
    # 9. Run experiment
    logger.info("\n" + "="*70)
    logger.info("STARTING ENSEMBLE VALIDATION EXPERIMENT")
    logger.info("="*70)
    logger.info("⏰ This will take several hours...")
    logger.info(f"Total compounds: {len(dataset)}")
    logger.info(f"Ensemble structures: {len(ensemble_structures)}")
    logger.info(f"Total docking runs: {len(dataset) * (1 + len(ensemble_structures))}")
    logger.info(f"Progress will be logged to: {output_dir / 'experiment.log'}")
    
    try:
        results = experiment.run_full_experiment(
            benchmark_dataset=dataset,
            single_structure=single_structure,
            ensemble_structures=ensemble_structures,
            grid_params=grid_params,
            n_bootstrap=args.bootstrap
        )
        
        logger.info("\n" + "="*70)
        logger.info("✅ EXPERIMENT COMPLETED SUCCESSFULLY!")
        logger.info("="*70)
        logger.info(f"Results saved to: {args.output}")
        logger.info("\nKey files:")
        logger.info(f"  📄 {args.output}/SUMMARY_REPORT.txt  ← READ THIS FIRST")
        logger.info(f"  📊 {args.output}/ensemble_validation_plots.png")
        logger.info(f"  📈 {args.output}/metrics.json")
        logger.info(f"  📋 {args.output}/single_results.csv")
        logger.info(f"  📋 {args.output}/ensemble_results.csv")
        
        # Print recommendation
        rec = results['recommendation']
        logger.info("\n" + "="*70)
        logger.info("🎯 FINAL RECOMMENDATION")
        logger.info("="*70)
        logger.info(f"Decision: {rec['recommendation']}")
        logger.info(f"Reason: {rec['reason']}")
        logger.info(f"AUC improvement: {rec['auc_improvement_pct']:+.1f}%")
        logger.info(f"Statistically significant: {rec['auc_significant']}")

        if args.compare_all_singles:
            logger.info("\n" + "="*70)
            logger.info("COMPARING ENSEMBLE VS EACH SINGLE STRUCTURE")
            logger.info("="*70)
            summary_rows = []
            ensemble_results = experiment.ensemble_results
            ensemble_metrics = results['ensemble_metrics']

            for idx, struct in enumerate(ensemble_structures):
                pdb_id = Path(struct).stem.replace("_prepared", "")
                if struct == single_structure:
                    continue

                logger.info(f"\nRunning single-structure docking for: {Path(struct).name}")
                single_df = experiment._run_single_experiment(
                    dataset,
                    struct,
                    grid_params
                )

                # Save per-single results
                single_out = output_dir / f"single_{pdb_id}_results.csv"
                single_df.to_csv(single_out, index=False)
                logger.info(f"  Saved: {single_out}")

                single_analysis = EnrichmentAnalysis(single_df)
                single_metrics = single_analysis.calculate_all_metrics()

                boot = BootstrapValidation.compare_methods(
                    single_df,
                    ensemble_results,
                    n_iterations=args.bootstrap
                )

                summary_rows.append({
                    'single_structure': Path(struct).name,
                    'single_auc': single_metrics['AUC'],
                    'ensemble_auc': ensemble_metrics['AUC'],
                    'auc_diff': ensemble_metrics['AUC'] - single_metrics['AUC'],
                    'single_ef1': single_metrics['EF_1%'],
                    'ensemble_ef1': ensemble_metrics['EF_1%'],
                    'ef1_diff': ensemble_metrics['EF_1%'] - single_metrics['EF_1%'],
                    'auc_p_value': boot['auc']['p_value'],
                    'auc_significant': boot['auc']['significant'],
                    'ef1_p_value': boot['ef1']['p_value'],
                    'ef1_significant': boot['ef1']['significant']
                })

            if summary_rows:
                summary_df = pd.DataFrame(summary_rows)
                summary_file = output_dir / "single_vs_ensemble_summary.csv"
                summary_df.to_csv(summary_file, index=False)
                logger.info(f"\nSaved comparison summary: {summary_file}")

        return 0
        
    except Exception as e:
        logger.exception(f"❌ Experiment failed: {e}")
        return 1
    finally:
        docking_manager.cleanup()


if __name__ == "__main__":
    sys.exit(main())
