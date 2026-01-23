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
    for idx, row in df.iterrows():
        if pd.isna(row['decoy SMILES']):
            continue
        
        smiles = row['decoy SMILES']
        mol = Chem.MolFromSmiles(smiles)
        
        if mol is None:
            logger.warning(f"  Invalid decoy SMILES for row {idx}")
            continue
        
        mol = Chem.AddHs(mol)
        result = AllChem.EmbedMolecule(mol, randomSeed=42)
        if result != 0:
            continue
        
        AllChem.UFFOptimizeMolecule(mol)
        
        decoy_id = row.get('decoy compound', f'decoy_{idx}')
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


def calculate_consensus_center(structures):
    """Calculate consensus grid center from ligand positions in PDBQT files."""
    logger.info("\nCalculating consensus grid center...")

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
            name = pdbqt_file.name
            base = None
            for suffix in ("_prepared.pdbqt", "_prep.pdbqt", "_receptor.pdbqt"):
                if name.endswith(suffix):
                    base = name[:-len(suffix)]
                    break
            if base is None and name.endswith(".pdbqt"):
                base = name[:-6]

            pdb_candidates = []
            if base:
                pdb_candidates.append(pdbqt_file.with_name(f"{base}.pdb"))
            pdb_candidates.append(pdbqt_file.with_suffix(".pdb"))

            pdb_path = next((p for p in pdb_candidates if p.exists()), None)
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
        logger.error("Could not extract any ligand positions!")
        logger.error("Using default center [0, 0, 0] - VERIFY THIS IS CORRECT!")
        return [0.0, 0.0, 0.0]
    
    consensus = np.mean(centers, axis=0)
    std = np.std(centers, axis=0)
    
    logger.info(f"\n  Consensus: ({consensus[0]:.2f}, {consensus[1]:.2f}, {consensus[2]:.2f})")
    logger.info(f"  Std dev: ({std[0]:.2f}, {std[1]:.2f}, {std[2]:.2f})")
    
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
    grid_center = calculate_consensus_center(ensemble_structures)
    
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
        
        return 0
        
    except Exception as e:
        logger.exception(f"❌ Experiment failed: {e}")
        return 1
    finally:
        docking_manager.cleanup()


if __name__ == "__main__":
    sys.exit(main())
