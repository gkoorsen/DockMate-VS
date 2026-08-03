"""
Enhanced RMSD calculation with active/decoy handling.

This module extends the base RMSD calculation to properly handle:
1. Self-docking validation (same molecule - calculate RMSD)
2. Enrichment studies (different molecules - skip RMSD)
"""

from pathlib import Path
from typing import Optional, Tuple
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from loguru import logger

from docking_platform_gui.utils.rmsd import (
    calculate_rmsd as _base_calculate_rmsd,
    _load_molecule,
    _load_molecules
)


def _flatten_for_comparison(mol: Chem.Mol) -> Optional[Chem.Mol]:
    """Reduce a molecule to its heavy-atom connectivity graph.

    Ligands extracted from a PDB file carry no CONECT bond orders, so RDKit
    reads every bond as single and every aromatic ring as saturated. Comparing
    such a molecule against a SMILES-derived reference by Morgan fingerprint
    gives a near-zero Tanimoto even when they are the same compound.

    Flattening both sides — single bonds, no aromaticity, no charges, no
    explicit hydrogens — makes the comparison depend on the atom/bond graph
    alone, which survives the round trip through PDB.
    """
    try:
        flat = Chem.RWMol(Chem.RemoveHs(Chem.Mol(mol), sanitize=False))
        for atom in flat.GetAtoms():
            atom.SetIsAromatic(False)
            atom.SetFormalCharge(0)
            atom.SetNoImplicit(True)
            atom.SetNumExplicitHs(0)
        for bond in flat.GetBonds():
            bond.SetIsAromatic(False)
            bond.SetBondType(Chem.BondType.SINGLE)
        out = flat.GetMol()
        Chem.SanitizeMol(out, sanitizeOps=Chem.SanitizeFlags.SANITIZE_SYMMRINGS)
        return out
    except Exception as exc:
        logger.debug(f"Could not flatten molecule for comparison: {exc}")
        return None


def are_same_molecule(
    mol1: Chem.Mol,
    mol2: Chem.Mol,
    similarity_threshold: float = 0.85
) -> bool:
    """
    Check if two molecules are the same (or very similar).
    
    Uses Tanimoto similarity of Morgan fingerprints.
    
    Args:
        mol1: First molecule
        mol2: Second molecule
        similarity_threshold: Tanimoto threshold (default 0.85)
        
    Returns:
        True if molecules are the same/similar enough for RMSD
    """
    try:
        # Generate Morgan fingerprints
        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=2048)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=2048)
        
        # Calculate Tanimoto similarity
        similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
        
        logger.debug(f"Molecular similarity: {similarity:.3f}")
        
        if similarity >= similarity_threshold:
            return True

        # A PDB-derived ligand has no bond orders, so the direct comparison
        # above under-reports badly for the very case we most need to detect:
        # a crystal ligand being redocked into its own structure. Retry on the
        # bond-order-agnostic connectivity graph before declaring a decoy.
        flat1, flat2 = _flatten_for_comparison(mol1), _flatten_for_comparison(mol2)
        if flat1 is None or flat2 is None:
            return False
        ffp1 = AllChem.GetMorganFingerprintAsBitVect(flat1, 2, nBits=2048)
        ffp2 = AllChem.GetMorganFingerprintAsBitVect(flat2, 2, nBits=2048)
        flat_similarity = DataStructs.TanimotoSimilarity(ffp1, ffp2)
        logger.debug(
            f"Connectivity-only similarity: {flat_similarity:.3f} "
            f"(bond-order-aware was {similarity:.3f})"
        )
        return flat_similarity >= similarity_threshold
        
    except Exception as e:
        logger.warning(f"Similarity check failed: {e}")
        # Fall back to atom count comparison
        return abs(mol1.GetNumAtoms() - mol2.GetNumAtoms()) <= 2


def calculate_rmsd_safe(
    reference_file: Path,
    docked_file: Path,
    check_similarity: bool = True,
    similarity_threshold: float = 0.85,
    use_symmetry: bool = True
) -> Optional[float]:
    """
    Safely calculate RMSD, returning None if molecules are different.
    
    This function first checks if the reference and docked molecules are
    similar enough to warrant RMSD calculation. This prevents errors when
    comparing actives (same molecule) vs decoys (different molecules).
    
    Args:
        reference_file: Crystal/reference structure (PDB, SDF, or PDBQT)
        docked_file: Docked structure (PDB, SDF, or PDBQT)
        check_similarity: Whether to check molecular similarity first
        similarity_threshold: Tanimoto threshold for considering molecules same
        use_symmetry: Whether to consider molecular symmetry in RMSD
        
    Returns:
        RMSD in Angstroms if molecules match, None if different or error
        
    Example:
        >>> rmsd = calculate_rmsd_safe(crystal_pdb, docked_pdbqt)
        >>> if rmsd is not None:
        ...     print(f"Self-docking RMSD: {rmsd:.2f}Å")
        ... else:
        ...     print("Decoy molecule - RMSD not applicable")
    """
    try:
        # Load molecules for similarity check
        if check_similarity:
            ref_mol = _load_molecule(reference_file)
            docked_mols = _load_molecules(docked_file)
            
            if ref_mol is None or not docked_mols:
                logger.warning(
                    f"Could not load molecules: ref={reference_file}, "
                    f"docked={docked_file}"
                )
                return None
            
            # Check if molecules are similar (use first pose)
            if not are_same_molecule(ref_mol, docked_mols[0], similarity_threshold):
                logger.info(
                    f"Molecules are different (decoy case) - skipping RMSD"
                )
                return None
        
        # Molecules are similar - calculate RMSD
        rmsd = _base_calculate_rmsd(reference_file, docked_file, use_symmetry)
        return rmsd
        
    except Exception as e:
        logger.warning(f"RMSD calculation failed: {e}")
        return None


def detect_molecule_type(
    ligand_smiles: str,
    crystal_ligand_file: Path
) -> str:
    """
    Detect if a ligand is an active (same as crystal) or decoy (different).
    
    Args:
        ligand_smiles: SMILES of the ligand to dock
        crystal_ligand_file: Path to crystal ligand structure
        
    Returns:
        'active' if same molecule, 'decoy' if different, 'unknown' if can't determine
    """
    try:
        # Parse input SMILES
        input_mol = Chem.MolFromSmiles(ligand_smiles)
        if input_mol is None:
            return 'unknown'
        
        # Load crystal ligand
        crystal_mol = _load_molecule(crystal_ligand_file)
        if crystal_mol is None:
            return 'unknown'
        
        # Check similarity
        if are_same_molecule(input_mol, crystal_mol, similarity_threshold=0.85):
            return 'active'
        else:
            return 'decoy'
            
    except Exception as e:
        logger.warning(f"Molecule type detection failed: {e}")
        return 'unknown'


class RedockResult:
    """Enhanced redock result with active/decoy awareness."""
    
    def __init__(
        self,
        pdb_id: str,
        ligand_name: str,
        molecule_type: str,  # 'active' or 'decoy'
        rmsd: Optional[float],
        score: float,
        success: Optional[bool],
        **kwargs
    ):
        self.pdb_id = pdb_id
        self.ligand_name = ligand_name
        self.molecule_type = molecule_type
        self.rmsd = rmsd
        self.score = score
        self.success = success
        self.metadata = kwargs
    
    def is_active(self) -> bool:
        """Check if this is an active (self-docking) case."""
        return self.molecule_type == 'active'
    
    def is_decoy(self) -> bool:
        """Check if this is a decoy (enrichment) case."""
        return self.molecule_type == 'decoy'
    
    def has_valid_rmsd(self) -> bool:
        """Check if RMSD is valid (not None)."""
        return self.rmsd is not None
    
    def to_dict(self) -> dict:
        """Convert to dictionary for export."""
        return {
            'pdb_id': self.pdb_id,
            'ligand_name': self.ligand_name,
            'molecule_type': self.molecule_type,
            'rmsd': self.rmsd,
            'score': self.score,
            'success': self.success,
            **self.metadata
        }


def calculate_enrichment_metrics(
    results: list,
    score_key: str = 'score',
    lower_is_better: bool = True
) -> dict:
    """
    Calculate enrichment metrics from redock results.
    
    Properly handles actives vs decoys:
    - Actives: Use RMSD for success definition
    - Decoys: Use score only for ranking
    - ROC AUC: Based on score ranking
    
    Args:
        results: List of RedockResult objects or dicts
        score_key: Key for docking score
        lower_is_better: True if lower scores are better
        
    Returns:
        Dictionary of enrichment metrics
    """
    from sklearn.metrics import roc_auc_score, roc_curve
    
    # Separate actives and decoys
    actives = []
    decoys = []
    
    for result in results:
        if isinstance(result, dict):
            mol_type = result.get('molecule_type', 'unknown')
            score = result.get(score_key, 0.0)
        else:
            mol_type = result.molecule_type
            score = getattr(result, score_key, 0.0)
        
        if mol_type == 'active':
            actives.append(score)
        elif mol_type == 'decoy':
            decoys.append(score)
    
    if not actives or not decoys:
        logger.warning("Need both actives and decoys for enrichment calculation")
        return {
            'roc_auc': None,
            'ef_1': None,
            'ef_5': None,
            'ef_10': None
        }
    
    # Create labels and scores for ROC
    labels = [1] * len(actives) + [0] * len(decoys)
    scores = actives + decoys
    
    # Invert scores if lower is better (for ROC curve)
    if lower_is_better:
        scores = [-s for s in scores]
    
    # Calculate ROC AUC
    try:
        roc_auc = roc_auc_score(labels, scores)
    except Exception as e:
        logger.warning(f"ROC AUC calculation failed: {e}")
        roc_auc = None
    
    # Calculate enrichment factors
    def calculate_ef(percent: float) -> Optional[float]:
        """Calculate enrichment factor at given percentage."""
        try:
            n_total = len(labels)
            n_actives = sum(labels)
            n_select = int(n_total * percent / 100)
            
            if n_select == 0:
                return None
            
            # Sort by score (descending for better scores)
            sorted_indices = sorted(
                range(len(scores)),
                key=lambda i: scores[i],
                reverse=True
            )
            
            # Count actives in top N
            actives_found = sum(labels[i] for i in sorted_indices[:n_select])
            
            # EF = (actives_found / n_select) / (n_actives / n_total)
            ef = (actives_found / n_select) / (n_actives / n_total)
            return ef
            
        except Exception as e:
            logger.warning(f"EF calculation failed: {e}")
            return None
    
    return {
        'roc_auc': roc_auc,
        'ef_1': calculate_ef(1.0),
        'ef_5': calculate_ef(5.0),
        'ef_10': calculate_ef(10.0),
        'n_actives': len(actives),
        'n_decoys': len(decoys)
    }


def calculate_redock_statistics(results: list) -> dict:
    """
    Calculate comprehensive redock statistics.
    
    Separates metrics for actives (RMSD-based) and enrichment (score-based).
    
    Args:
        results: List of RedockResult objects or dicts
        
    Returns:
        Dictionary with separate active and enrichment statistics
    """
    # Separate actives and decoys
    actives = []
    decoys = []
    
    for result in results:
        if isinstance(result, dict):
            mol_type = result.get('molecule_type', 'unknown')
        else:
            mol_type = result.molecule_type
        
        if mol_type == 'active':
            actives.append(result)
        elif mol_type == 'decoy':
            decoys.append(result)
    
    stats = {
        'total_cases': len(results),
        'n_actives': len(actives),
        'n_decoys': len(decoys)
    }
    
    # Active statistics (RMSD-based)
    if actives:
        valid_rmsds = []
        successes = []
        
        for result in actives:
            if isinstance(result, dict):
                rmsd = result.get('rmsd')
                success = result.get('success')
            else:
                rmsd = result.rmsd
                success = result.success
            
            if rmsd is not None:
                valid_rmsds.append(rmsd)
            if success is not None:
                successes.append(success)
        
        stats['active_stats'] = {
            'n_valid_rmsd': len(valid_rmsds),
            'mean_rmsd': np.mean(valid_rmsds) if valid_rmsds else None,
            'median_rmsd': np.median(valid_rmsds) if valid_rmsds else None,
            'success_rate': (
                sum(successes) / len(successes) * 100 
                if successes else None
            )
        }
    else:
        stats['active_stats'] = None
    
    # Enrichment statistics (score-based)
    if actives and decoys:
        stats['enrichment_stats'] = calculate_enrichment_metrics(results)
    else:
        stats['enrichment_stats'] = None
    
    return stats
