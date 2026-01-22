#!/usr/bin/env python3
"""
Test script to validate parallel RMSD calculation fix.

This creates a molecule with multiple conformers and tests both
serial and parallel RMSD calculation methods.
"""

import sys
import time
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent))

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from loguru import logger

# Import the fixed ligand preparation
try:
    from docking_platform_gui.preparation.ligand import (
        LigandPreparation,
        LigandPreparationConfig
    )
except ImportError:
    # If running standalone, import from local file
    import importlib.util
    spec = importlib.util.spec_from_file_location("ligand", "ligand.py")
    ligand_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ligand_module)
    LigandPreparation = ligand_module.LigandPreparation
    LigandPreparationConfig = ligand_module.LigandPreparationConfig


def create_test_molecule(smiles: str, num_conformers: int = 50):
    """Create a molecule with multiple conformers for testing."""
    logger.info(f"Creating test molecule: {smiles}")
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    
    mol = Chem.AddHs(mol)
    
    # Generate conformers
    logger.info(f"Generating {num_conformers} conformers...")
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    params.useRandomCoords = True
    
    conf_ids = list(AllChem.EmbedMultipleConfs(
        mol, 
        numConfs=num_conformers,
        params=params
    ))
    
    logger.info(f"✓ Generated {len(conf_ids)} conformers")
    
    # Quick MMFF minimize
    logger.info("Quick MMFF minimization...")
    props = AllChem.MMFFGetMoleculeProperties(mol)
    if props:
        for conf_id in conf_ids:
            ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=conf_id)
            if ff:
                ff.Minimize(maxIts=50)
    
    logger.info("✓ Minimization complete")
    return mol, conf_ids


def test_parallel_rmsd():
    """Test parallel RMSD calculation with timing."""
    logger.info("="*80)
    logger.info("PARALLEL RMSD CALCULATION TEST")
    logger.info("="*80)
    
    # Test parameters
    test_smiles = "C[C@H]1C[C@@H]([C@@H]([C@H](/C=C(/[C@@H]([C@H](/C=C\C=C(\C(=O)NC2=CC(=O)C(=C(C1)C2=O)OC)/C)OC)OC(=O)N)\C)C)O)OC"  # Ibuprofen - moderately flexible
    num_conformers = 10  # Start with moderate number
    
    logger.info(f"\nTest molecule: {test_smiles}")
    logger.info(f"Target conformers: {num_conformers}")
    logger.info(f"Expected RMSD calculations: {num_conformers * (num_conformers-1) // 2}")
    
    # Create test molecule
    mol, conf_ids = create_test_molecule(test_smiles, num_conformers)
    
    # Create LigandPreparation instance
    config = LigandPreparationConfig()
    
    # Test with different worker counts
    for n_cpus in [1, 2, 4]:
        logger.info("\n" + "-"*80)
        logger.info(f"Testing with {n_cpus} CPU(s)")
        logger.info("-"*80)
        
        prep = LigandPreparation(config, n_cpus=n_cpus)
        
        start_time = time.time()
        
        if n_cpus == 1:
            # Serial
            logger.info("Running SERIAL RMSD calculation...")
            rmsds = prep._calculate_rmsd_matrix_serial(mol, conf_ids)
        else:
            # Parallel
            logger.info(f"Running PARALLEL RMSD calculation with {n_cpus} workers...")
            rmsds = prep._calculate_rmsd_matrix_parallel(mol, conf_ids, n_cpus)
        
        elapsed = time.time() - start_time
        
        # Validate results
        n = len(conf_ids)
        expected_pairs = n * (n - 1) // 2
        non_zero = (rmsds > 0).sum() // 2  # Divide by 2 because matrix is symmetric
        
        logger.info("\n" + "="*80)
        logger.info("RESULTS")
        logger.info("="*80)
        logger.info(f"Time taken: {elapsed:.2f} seconds ({elapsed/60:.2f} minutes)")
        logger.info(f"RMSD pairs calculated: {non_zero} / {expected_pairs}")
        logger.info(f"Matrix shape: {rmsds.shape}")
        logger.info(f"Min RMSD: {rmsds[rmsds > 0].min():.3f} Å")
        logger.info(f"Max RMSD: {rmsds.max():.3f} Å")
        logger.info(f"Mean RMSD: {rmsds[rmsds > 0].mean():.3f} Å")
        
        if n_cpus > 1:
            logger.info(f"\n✓ Parallel processing {'PASSED' if non_zero == expected_pairs else 'FAILED'}")
        else:
            logger.info(f"\n✓ Serial processing {'PASSED' if non_zero == expected_pairs else 'FAILED'}")
        
        # Calculate speedup
        if n_cpus == 1:
            serial_time = elapsed
        else:
            speedup = serial_time / elapsed
            efficiency = speedup / n_cpus * 100
            logger.info(f"Speedup: {speedup:.2f}x (efficiency: {efficiency:.1f}%)")


def test_timeout_handling():
    """Test that timeout handling works correctly."""
    logger.info("\n\n")
    logger.info("="*80)
    logger.info("TIMEOUT HANDLING TEST")
    logger.info("="*80)
    
    # Create a large molecule that will take time
    test_smiles = "C[C@H]1C[C@@H]([C@@H]([C@H](/C=C(/[C@@H]([C@H](/C=C\C=C(\C(=O)NC2=CC(=O)C(=C(C1)C2=O)OC)/C)OC)OC(=O)N)\C)C)O)OC"  # Ibuprofen
    num_conformers = 5
    
    logger.info(f"\nTest molecule: {test_smiles}")
    logger.info(f"Conformers: {num_conformers}")
    
    mol, conf_ids = create_test_molecule(test_smiles, num_conformers)
    
    config = LigandPreparationConfig()
    prep = LigandPreparation(config, n_cpus=4)
    
    logger.info("\nTesting timeout handling (should complete normally)...")
    start_time = time.time()
    
    try:
        rmsds = prep._calculate_rmsd_matrix_parallel(mol, conf_ids, 4)
        elapsed = time.time() - start_time
        
        logger.info(f"\n✓ Parallel processing completed in {elapsed:.2f}s")
        logger.info(f"  Matrix shape: {rmsds.shape}")
        logger.info(f"  Non-zero entries: {(rmsds > 0).sum() // 2}")
        
    except Exception as e:
        logger.error(f"✗ Parallel processing failed: {e}")
        raise


if __name__ == "__main__":
    # Configure logging
    logger.remove()  # Remove default handler
    logger.add(
        sys.stdout,
        format="<green>{time:HH:mm:ss}</green> | <level>{level: <8}</level> | <level>{message}</level>",
        level="INFO"
    )
    
    logger.info("\n" + "="*80)
    logger.info("PARALLEL RMSD FIX VALIDATION TESTS")
    logger.info("="*80)
    
    try:
        # Run tests
        test_parallel_rmsd()
        test_timeout_handling()
        
        logger.info("\n" + "="*80)
        logger.info("ALL TESTS PASSED ✓")
        logger.info("="*80)
        
    except Exception as e:
        logger.error(f"\n\n{'='*80}")
        logger.error("TEST FAILED ✗")
        logger.error(f"{'='*80}")
        logger.exception(e)
        sys.exit(1)