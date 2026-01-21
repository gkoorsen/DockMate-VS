"""
Co-crystal ligand extraction and binding site definition.

This module implements binding site definition from co-crystal structures,
following best practices from the molecular docking literature.

Key capabilities:
- Extract ligand from PDB co-crystal structure
- Calculate binding box center and size
- RMSD calculation for self-docking validation
- Support for validation that RMSD < 2.0Å (success criterion)

References:
- Binding site definition: Eberhardt et al., J. Chem. Inf. Model., 2021
- RMSD validation: Réau et al., J. Chem. Inf. Model., 2018
"""

import numpy as np
from pathlib import Path
from typing import Optional, Tuple, List
from dataclasses import dataclass

try:
    import MDAnalysis as mda
    from MDAnalysis.analysis import align, rms
    HAS_MDANALYSIS = True
except ImportError:
    HAS_MDANALYSIS = False

from loguru import logger


class BindingSiteError(Exception):
    """Raised when binding site definition fails."""
    pass


@dataclass
class BindingSite:
    """
    Binding site definition with box parameters.

    Attributes:
        center: Box center coordinates (x, y, z) in Angstroms
        size: Box dimensions (x, y, z) in Angstroms
        ligand_resname: Residue name of co-crystal ligand
        ligand_chain: Chain ID of co-crystal ligand
        ligand_atoms: Number of atoms in co-crystal ligand
        source_pdb: Path to source PDB file
    """
    center: np.ndarray
    size: np.ndarray
    ligand_resname: str
    ligand_chain: str
    ligand_atoms: int
    source_pdb: str
    self_docking_result: Optional['SelfDockingResult'] = None

    def to_dict(self) -> dict:
        """Convert to dictionary for serialization."""
        result = {
            'center': self.center.tolist(),
            'size': self.size.tolist(),
            'ligand_resname': self.ligand_resname,
            'ligand_chain': self.ligand_chain,
            'ligand_atoms': self.ligand_atoms,
            'source_pdb': self.source_pdb
        }
        if self.self_docking_result:
            result['self_docking_rmsd'] = self.self_docking_result.rmsd
            result['self_docking_success'] = self.self_docking_result.success
        return result


@dataclass
class SelfDockingResult:
    """
    Results from self-docking validation.

    Attributes:
        rmsd: RMSD between docked and crystal pose (Angstroms)
        success: Whether RMSD < threshold (default 2.0Å)
        threshold: RMSD threshold for success
        docked_pose: Docked ligand coordinates
        crystal_pose: Crystal ligand coordinates
    """
    rmsd: float
    success: bool
    threshold: float
    docked_pose: Optional[np.ndarray] = None
    crystal_pose: Optional[np.ndarray] = None

    def to_dict(self) -> dict:
        """Convert to dictionary for serialization."""
        return {
            'rmsd': self.rmsd,
            'success': self.success,
            'threshold': self.threshold
        }


class BindingSiteDefinition:
    """
    Define binding site from co-crystal structure.

    This class extracts ligand information from a PDB file containing
    a protein-ligand complex and defines the docking box parameters.

    The box is centered on the ligand centroid with configurable margins.
    Standard practice is 12Å margin from ligand extremes.

    Examples:
        >>> bsd = BindingSiteDefinition(margin=12.0)
        >>> site = bsd.from_cocrystal(
        ...     pdb_path="1hsg.pdb",
        ...     ligand_resname="MK1",
        ...     ligand_chain="A"
        ... )
        >>> print(f"Center: {site.center}")
        >>> print(f"Size: {site.size}")
    """

    def __init__(self, margin: float = 12.0):
        """
        Initialize binding site definition.

        Args:
            margin: Distance (Angstroms) to extend box beyond ligand extremes.
                   Literature standard is 12Å to allow ligand rotation/translation.

        References:
            Margin selection: Eberhardt et al., J. Chem. Inf. Model., 2021
        """
        if not HAS_MDANALYSIS:
            raise ImportError(
                "MDAnalysis is required for binding site definition. "
                "Install with: pip install mdanalysis"
            )

        self.margin = margin
        logger.info(f"Initialized BindingSiteDefinition with {margin}Å margin")

    def from_cocrystal(
        self,
        pdb_path: str,
        ligand_resname: str,
        ligand_chain: Optional[str] = None
    ) -> BindingSite:
        """
        Extract binding site from co-crystal structure.

        Args:
            pdb_path: Path to PDB file containing protein-ligand complex
            ligand_resname: Residue name of ligand (e.g., "MK1", "ATP")
            ligand_chain: Chain ID of ligand (optional, auto-detect if None)

        Returns:
            BindingSite object with box parameters

        Raises:
            BindingSiteError: If ligand not found or extraction fails
        """
        pdb_path = Path(pdb_path)
        if not pdb_path.exists():
            raise BindingSiteError(f"PDB file not found: {pdb_path}")

        logger.info(f"Loading co-crystal structure: {pdb_path}")

        try:
            universe = mda.Universe(str(pdb_path))
        except Exception as e:
            raise BindingSiteError(f"Failed to load PDB: {e}")

        # Find ligand
        ligand = self._find_ligand(universe, ligand_resname, ligand_chain)

        # Extract actual ligand chain (may differ from specified if fallback was used)
        actual_chain = ligand.segments[0].segid
        if ligand_chain is None:
            ligand_chain = actual_chain
            logger.info(f"Auto-detected ligand chain: {ligand_chain}")
        elif ligand_chain != actual_chain:
            logger.warning(
                f"Specified chain '{ligand_chain}' does not contain ligand. "
                f"Using actual chain '{actual_chain}'"
            )
            ligand_chain = actual_chain

        # Calculate box parameters
        center = self._calculate_center(ligand)
        size = self._calculate_size(ligand)

        site = BindingSite(
            center=center,
            size=size,
            ligand_resname=ligand_resname,
            ligand_chain=ligand_chain,
            ligand_atoms=len(ligand.atoms),
            source_pdb=str(pdb_path)
        )

        logger.info(
            f"Binding site defined: center={center}, size={size}, "
            f"ligand_atoms={len(ligand.atoms)}"
        )

        return site

    def _find_ligand(
        self,
        universe: 'mda.Universe',
        ligand_resname: str,
        ligand_chain: Optional[str] = None
    ) -> 'mda.AtomGroup':
        """
        Find ligand in structure.

        Args:
            universe: MDAnalysis Universe
            ligand_resname: Ligand residue name
            ligand_chain: Optional chain ID

        Returns:
            AtomGroup containing ligand atoms

        Raises:
            BindingSiteError: If ligand not found
        """
        # Build selection string
        if ligand_chain:
            selection = f"resname {ligand_resname} and segid {ligand_chain}"
        else:
            selection = f"resname {ligand_resname}"

        logger.debug(f"Searching for ligand with selection: {selection}")

        try:
            ligand = universe.select_atoms(selection)
        except Exception as e:
            raise BindingSiteError(f"Failed to select ligand: {e}")

        if len(ligand.atoms) == 0:
            # Try alternative selection without segid
            selection = f"resname {ligand_resname}"
            ligand = universe.select_atoms(selection)

            if len(ligand.atoms) == 0:
                raise BindingSiteError(
                    f"Ligand '{ligand_resname}' not found in structure. "
                    f"Available residues: {set(universe.residues.resnames)}"
                )

        logger.info(f"Found ligand '{ligand_resname}' with {len(ligand.atoms)} atoms")
        return ligand

    def _calculate_center(self, ligand: 'mda.AtomGroup') -> np.ndarray:
        """
        Calculate box center from ligand centroid.

        Args:
            ligand: Ligand atom group

        Returns:
            Center coordinates (x, y, z)
        """
        center = ligand.center_of_geometry()
        logger.debug(f"Calculated box center: {center}")
        return center

    def _calculate_size(self, ligand: 'mda.AtomGroup') -> np.ndarray:
        """
        Calculate box size from ligand extent plus margin.

        The box extends `self.margin` Angstroms beyond ligand extremes
        in each dimension.

        Args:
            ligand: Ligand atom group

        Returns:
            Box dimensions (x, y, z)
        """
        positions = ligand.positions

        # Calculate extent in each dimension
        min_coords = positions.min(axis=0)
        max_coords = positions.max(axis=0)
        extent = max_coords - min_coords

        # Add margin on both sides (2 * margin total per dimension)
        size = extent + (2 * self.margin)

        logger.debug(f"Ligand extent: {extent}")
        logger.debug(f"Box size (with {self.margin}Å margin): {size}")

        return size

    def extract_ligand_coordinates(
        self,
        pdb_path: str,
        ligand_resname: str,
        ligand_chain: Optional[str] = None
    ) -> np.ndarray:
        """
        Extract ligand coordinates for RMSD calculation.

        Args:
            pdb_path: Path to PDB file
            ligand_resname: Ligand residue name
            ligand_chain: Optional chain ID

        Returns:
            Nx3 array of ligand atom coordinates
        """
        universe = mda.Universe(str(pdb_path))
        ligand = self._find_ligand(universe, ligand_resname, ligand_chain)
        return ligand.positions.copy()

    def calculate_rmsd(
        self,
        coords1: np.ndarray,
        coords2: np.ndarray,
        align: bool = True
    ) -> float:
        """
        Calculate RMSD between two sets of coordinates.

        Args:
            coords1: First coordinate set (Nx3)
            coords2: Second coordinate set (Nx3)
            align: Whether to align structures before RMSD calculation
                  (recommended for docking validation)

        Returns:
            RMSD in Angstroms

        Raises:
            BindingSiteError: If coordinate shapes don't match

        References:
            RMSD calculation: Coutsias et al., J. Comput. Chem., 2004
        """
        if coords1.shape != coords2.shape:
            raise BindingSiteError(
                f"Coordinate shape mismatch: {coords1.shape} vs {coords2.shape}"
            )

        if align:
            # Align coords2 to coords1 using Kabsch algorithm
            coords2_aligned = self._kabsch_align(coords1, coords2)
        else:
            coords2_aligned = coords2

        # Calculate RMSD
        diff = coords1 - coords2_aligned
        rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))

        return rmsd

    def _kabsch_align(
        self,
        coords_ref: np.ndarray,
        coords_mobile: np.ndarray
    ) -> np.ndarray:
        """
        Align mobile coordinates to reference using Kabsch algorithm.

        Args:
            coords_ref: Reference coordinates (Nx3)
            coords_mobile: Mobile coordinates to align (Nx3)

        Returns:
            Aligned mobile coordinates

        References:
            Kabsch algorithm: Kabsch, Acta Cryst., 1976
        """
        # Center coordinates
        center_ref = coords_ref.mean(axis=0)
        center_mobile = coords_mobile.mean(axis=0)

        coords_ref_centered = coords_ref - center_ref
        coords_mobile_centered = coords_mobile - center_mobile

        # Calculate rotation matrix using SVD
        correlation = np.dot(coords_mobile_centered.T, coords_ref_centered)
        U, S, Vt = np.linalg.svd(correlation)

        # Handle reflection case
        d = np.linalg.det(np.dot(Vt.T, U.T))
        if d < 0:
            Vt[-1, :] *= -1

        rotation = np.dot(Vt.T, U.T)

        # Apply rotation and translation
        coords_aligned = np.dot(coords_mobile_centered, rotation.T) + center_ref

        return coords_aligned

    def validate_self_docking(
        self,
        crystal_coords: np.ndarray,
        docked_coords: np.ndarray,
        threshold: float = 2.0
    ) -> SelfDockingResult:
        """
        Validate self-docking by comparing docked pose to crystal pose.

        Success criterion: RMSD < 2.0Å is standard in literature.

        Args:
            crystal_coords: Crystal ligand coordinates (Nx3)
            docked_coords: Docked ligand coordinates (Nx3)
            threshold: RMSD threshold for success (default 2.0Å)

        Returns:
            SelfDockingResult with RMSD and success status

        References:
            Success criterion: Réau et al., J. Chem. Inf. Model., 2018
            70% success rate is typical for modern docking engines
        """
        rmsd = self.calculate_rmsd(crystal_coords, docked_coords, align=True)
        success = rmsd < threshold

        result = SelfDockingResult(
            rmsd=rmsd,
            success=success,
            threshold=threshold,
            docked_pose=docked_coords.copy(),
            crystal_pose=crystal_coords.copy()
        )

        logger.info(
            f"Self-docking validation: RMSD={rmsd:.2f}Å, "
            f"Success={success} (threshold={threshold}Å)"
        )

        return result
