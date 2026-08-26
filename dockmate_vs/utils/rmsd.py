"""
RMSD calculation utilities.

Provides functions to calculate RMSD between crystal and docked ligand structures.
"""

from pathlib import Path
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFMCS
from typing import Optional


def calculate_rmsd(
    reference_file: Path,
    docked_file: Path,
    use_symmetry: bool = True
) -> float:
    """
    Calculate RMSD between reference and docked structures.

    Args:
        reference_file: Crystal/reference structure (PDB, SDF, or PDBQT)
        docked_file: Docked structure (PDB, SDF, or PDBQT)
        use_symmetry: Whether to consider molecular symmetry

    Returns:
        Best RMSD in Angstroms

    Raises:
        ValueError: If files cannot be read or molecules don't match
    """
    # Load reference molecule
    ref_mol = _load_molecule(reference_file)
    if ref_mol is None:
        raise ValueError(f"Could not load reference molecule: {reference_file}")

    # Load docked molecules (may have multiple poses)
    docked_mols = _load_molecules(docked_file)
    if not docked_mols:
        raise ValueError(f"Could not load docked molecules: {docked_file}")

    # Calculate RMSD for each pose
    rmsds = []
    for mol in docked_mols:
        try:
            rmsd = coordinate_rmsd(ref_mol, mol, use_symmetry=use_symmetry)
            rmsds.append(rmsd)
            continue
        except Exception:
            pass

        try:
            rmsd = _rmsd_via_mcs(ref_mol, mol)
            rmsds.append(rmsd)
        except Exception:
            # Skip poses that can't be aligned
            continue

    if not rmsds:
        raise ValueError("Could not calculate RMSD for any pose")

    return min(rmsds)


def _load_molecule(file_path: Path) -> Optional[Chem.Mol]:
    """Load a single molecule from file."""
    suffix = file_path.suffix.lower()

    if suffix == '.pdb':
        mol = Chem.MolFromPDBFile(str(file_path), removeHs=False)
        if mol is None:
            mol = Chem.MolFromPDBFile(str(file_path), removeHs=False, sanitize=False)
        return mol

    elif suffix == '.sdf' or suffix == '.sd':
        with _silence_rdkit():
            suppl = Chem.SDMolSupplier(str(file_path), removeHs=False)
            mol = suppl[0] if len(suppl) > 0 else None
            if mol is None:
                suppl = Chem.SDMolSupplier(str(file_path), removeHs=False, sanitize=False)
                mol = suppl[0] if len(suppl) > 0 else None
        return mol

    elif suffix == '.pdbqt':
        # PDBQT: Convert to PDB first, then load
        pdb_content = _pdbqt_to_pdb(file_path)
        if pdb_content:
            mol = Chem.MolFromPDBBlock(pdb_content, removeHs=False)
            if mol is None:
                mol = Chem.MolFromPDBBlock(pdb_content, removeHs=False, sanitize=False)
            return mol
        return None

    else:
        raise ValueError(f"Unsupported file format: {suffix}")


def _load_molecules(file_path: Path) -> list:
    """Load multiple molecules from file (for multi-pose docking results)."""
    suffix = file_path.suffix.lower()
    molecules = []

    if suffix == '.pdb':
        # PDB typically has one structure
        mol = Chem.MolFromPDBFile(str(file_path), removeHs=False)
        if mol is None:
            mol = Chem.MolFromPDBFile(str(file_path), removeHs=False, sanitize=False)
        if mol:
            molecules.append(mol)

    elif suffix == '.sdf' or suffix == '.sd':
        # SDF can have multiple structures
        with _silence_rdkit():
            suppl = Chem.SDMolSupplier(str(file_path), removeHs=False)
            for mol in suppl:
                if mol is not None:
                    molecules.append(mol)
            if not molecules:
                suppl = Chem.SDMolSupplier(str(file_path), removeHs=False, sanitize=False)
                for mol in suppl:
                    if mol is not None:
                        molecules.append(mol)

    elif suffix == '.pdbqt':
        # PDBQT: Parse models
        pdb_blocks = _parse_pdbqt_models(file_path)
        for pdb_block in pdb_blocks:
            mol = Chem.MolFromPDBBlock(pdb_block, removeHs=False)
            if mol is None:
                mol = Chem.MolFromPDBBlock(pdb_block, removeHs=False, sanitize=False)
            if mol:
                molecules.append(mol)

    else:
        raise ValueError(f"Unsupported file format: {suffix}")

    return molecules


def coordinate_rmsd(
    ref_mol: Chem.Mol,
    probe_mol: Chem.Mol,
    use_symmetry: bool = True,
    max_matches: int = 1000,
) -> float:
    """Calculate heavy-atom RMSD in the existing receptor coordinate frame.

    Docked and crystal ligands are already expressed in the same receptor
    frame.  Aligning one ligand onto the other would erase displacement from
    the binding site and report conformational RMSD instead of pose RMSD.
    Atom mappings are therefore evaluated without changing either conformer.
    """
    ref_noh = Chem.RemoveHs(ref_mol, sanitize=False)
    probe_noh = Chem.RemoveHs(probe_mol, sanitize=False)

    ref_atoms = ref_noh.GetNumAtoms()
    probe_atoms = probe_noh.GetNumAtoms()
    if ref_atoms != probe_atoms:
        raise ValueError(
            f"Heavy-atom counts differ: ref {ref_atoms}, probe {probe_atoms}"
        )
    if ref_atoms < 1:
        raise ValueError("RMSD requires at least one heavy atom")

    mcs = rdFMCS.FindMCS(
        [ref_noh, probe_noh],
        ringMatchesRingOnly=False,
        completeRingsOnly=False,
        matchValences=False,
        bondCompare=rdFMCS.BondCompare.CompareAny,
        atomCompare=rdFMCS.AtomCompare.CompareElements
    )

    if mcs.numAtoms != ref_atoms:
        raise ValueError(
            f"Molecules are not graph-equivalent: matched {mcs.numAtoms} of "
            f"{ref_atoms} heavy atoms"
        )

    pattern = Chem.MolFromSmarts(mcs.smartsString)
    if pattern is None:
        raise ValueError("MCS pattern generation failed")
    ref_matches = ref_noh.GetSubstructMatches(
        pattern, uniquify=not use_symmetry, maxMatches=max_matches
    )
    probe_matches = probe_noh.GetSubstructMatches(
        pattern, uniquify=not use_symmetry, maxMatches=max_matches
    )
    if not ref_matches or not probe_matches:
        raise ValueError("MCS substructure match failed")

    ref_conf = ref_noh.GetConformer()
    probe_conf = probe_noh.GetConformer()
    best = float("inf")
    for ref_match in ref_matches:
        ref_xyz = np.asarray([ref_conf.GetAtomPosition(i) for i in ref_match])
        for probe_match in probe_matches:
            probe_xyz = np.asarray([probe_conf.GetAtomPosition(i) for i in probe_match])
            value = float(np.sqrt(np.mean(np.sum((ref_xyz - probe_xyz) ** 2, axis=1))))
            best = min(best, value)
    return best


def _rmsd_via_mcs(ref_mol: Chem.Mol, probe_mol: Chem.Mol) -> float:
    """Compatibility wrapper for coordinate-frame, MCS-mapped RMSD."""
    return coordinate_rmsd(ref_mol, probe_mol, use_symmetry=True)


class _silence_rdkit:
    """Context manager to silence RDKit warning/error output."""

    def __enter__(self):
        self._disabled = ["rdApp.error", "rdApp.warning"]
        for name in self._disabled:
            RDLogger.DisableLog(name)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        for name in self._disabled:
            RDLogger.EnableLog(name)
        return False


def _pdbqt_to_pdb(pdbqt_file: Path) -> Optional[str]:
    """Convert PDBQT to PDB format (simple conversion, first model only)."""
    try:
        with open(pdbqt_file) as f:
            lines = []
            for line in f:
                if line.startswith(('ATOM', 'HETATM')):
                    # Remove PDBQT charge columns (last 2 fields)
                    lines.append(line[:66] + '\n')
                elif line.startswith(('END', 'ENDMDL')):
                    break
                elif line.startswith('MODEL'):
                    continue
            return ''.join(lines)
    except Exception:
        return None


def _parse_pdbqt_models(pdbqt_file: Path) -> list:
    """Parse all MODELs from PDBQT file."""
    models = []
    current_model = []

    with open(pdbqt_file) as f:
        in_model = False

        for line in f:
            if line.startswith('MODEL'):
                in_model = True
                current_model = []

            elif line.startswith('ENDMDL'):
                in_model = False
                if current_model:
                    models.append(''.join(current_model))
                current_model = []

            elif in_model and line.startswith(('ATOM', 'HETATM')):
                # Remove PDBQT charge columns
                current_model.append(line[:66] + '\n')

    # Handle case where there's no MODEL/ENDMDL tags (single structure)
    if not models:
        single_model = []
        with open(pdbqt_file) as f:
            for line in f:
                if line.startswith(('ATOM', 'HETATM')):
                    single_model.append(line[:66] + '\n')
                elif line.startswith('END'):
                    break
        if single_model:
            models.append(''.join(single_model))

    return models


def calculate_rmsd_matrix(
    reference_file: Path,
    docked_file: Path,
    max_poses: int = 20
) -> np.ndarray:
    """
    Calculate RMSD matrix for all poses.

    Args:
        reference_file: Crystal/reference structure
        docked_file: Docked structure with multiple poses
        max_poses: Maximum number of poses to calculate

    Returns:
        Array of RMSDs for each pose

    Example:
        rmsds = calculate_rmsd_matrix(crystal.pdb, docked.pdbqt)
        print(f"Best: {rmsds.min():.2f}Å")
        print(f"Worst: {rmsds.max():.2f}Å")
        print(f"Mean: {rmsds.mean():.2f}Å")
    """
    ref_mol = _load_molecule(reference_file)
    if ref_mol is None:
        raise ValueError(f"Could not load reference: {reference_file}")

    docked_mols = _load_molecules(docked_file)
    if not docked_mols:
        raise ValueError(f"Could not load docked poses: {docked_file}")

    rmsds = []
    for mol in docked_mols[:max_poses]:
        try:
            rmsd = coordinate_rmsd(ref_mol, mol, use_symmetry=True)
            rmsds.append(rmsd)
        except Exception:
            rmsds.append(999.9)  # Failed alignment

    return np.array(rmsds)


# Convenience function for quick RMSD calculation
def quick_rmsd(ref_pdb: str, docked_pdbqt: str) -> float:
    """
    Quick RMSD calculation with simple API.

    Example:
        rmsd = quick_rmsd("crystal.pdb", "docked.pdbqt")
        print(f"RMSD: {rmsd:.2f}Å")
    """
    return calculate_rmsd(Path(ref_pdb), Path(docked_pdbqt))
