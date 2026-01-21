"""
Abstract base class for docking engines.

Defines the interface that all docking engines must implement,
supporting multiple flexibility modes (rigid, flexible side-chains,
ensemble, hybrid multi-phase).

References:
- AutoDock Vina: Eberhardt et al., J. Chem. Inf. Model., 2021
- Smina: Koes et al., J. Chem. Inf. Model., 2013
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import List, Optional, Dict, Any
from pathlib import Path
from enum import Enum
import numpy as np

from loguru import logger


class FlexibilityMode(str, Enum):
    """Docking flexibility modes."""
    RIGID = "rigid"
    FLEXIBLE_SIDECHAIN = "flexible_sidechain"
    ENSEMBLE = "ensemble"
    HYBRID = "hybrid"


@dataclass
class DockingPose:
    """
    Single docking pose result.

    Attributes:
        score: Binding affinity score (kcal/mol)
        rmsd_lb: Lower bound RMSD from best pose (Angstroms)
        rmsd_ub: Upper bound RMSD from best pose (Angstroms)
        coordinates: Ligand atom coordinates (Nx3 array)
        pdbqt: PDBQT format output
        mode_number: Pose number (1 = best)
    """
    score: float
    rmsd_lb: float = 0.0
    rmsd_ub: float = 0.0
    coordinates: Optional[np.ndarray] = None
    pdbqt: Optional[str] = None
    mode_number: int = 1


@dataclass
class DockingResult:
    """
    Complete docking result for a ligand.

    Attributes:
        ligand_id: Ligand identifier
        poses: List of docking poses (sorted by score)
        best_score: Best binding affinity (kcal/mol)
        success: Whether docking completed successfully
        runtime: Docking runtime (seconds)
        error_message: Error message if docking failed
        metadata: Additional metadata
    """
    ligand_id: str
    poses: List[DockingPose] = field(default_factory=list)
    best_score: Optional[float] = None
    success: bool = True
    runtime: Optional[float] = None
    error_message: Optional[str] = None
    metadata: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        """Calculate best_score from poses if not set."""
        if self.best_score is None and self.poses:
            self.best_score = min(pose.score for pose in self.poses)

    def to_dict(self) -> dict:
        """Convert to dictionary for serialization."""
        return {
            'ligand_id': self.ligand_id,
            'best_score': self.best_score,
            'num_poses': len(self.poses),
            'success': self.success,
            'runtime': self.runtime,
            'error_message': self.error_message,
            'metadata': self.metadata,
            'poses': [
                {
                    'mode': pose.mode_number,
                    'score': pose.score,
                    'rmsd_lb': pose.rmsd_lb,
                    'rmsd_ub': pose.rmsd_ub
                }
                for pose in self.poses
            ]
        }


class DockingEngine(ABC):
    """
    Abstract base class for docking engines.

    All docking engines must implement this interface to ensure
    compatibility with the workflow orchestrator.

    Subclasses should implement engine-specific command generation,
    execution, and output parsing.
    """

    def __init__(
        self,
        exhaustiveness: int = 8,
        num_modes: int = 9,
        energy_range: float = 3.0,
        cpu: int = 1,
        seed: Optional[int] = None
    ):
        """
        Initialize docking engine.

        Args:
            exhaustiveness: Thoroughness of search (higher = more thorough)
                           Default: 8 (standard)
                           Flexible/ensemble: recommend 16-32
            num_modes: Maximum number of binding modes to generate
            energy_range: Maximum energy difference from best mode (kcal/mol)
            cpu: Number of CPU cores to use
            seed: Random seed for reproducibility (None = random)
        """
        self.exhaustiveness = exhaustiveness
        self.num_modes = num_modes
        self.energy_range = energy_range
        self.cpu = cpu
        self.seed = seed

        logger.info(
            f"Initialized {self.__class__.__name__} with "
            f"exhaustiveness={exhaustiveness}, num_modes={num_modes}, cpu={cpu}"
        )

    @abstractmethod
    def dock(
        self,
        receptor_path: str,
        ligand_path: str,
        center: np.ndarray,
        size: np.ndarray,
        output_path: Optional[str] = None,
        **kwargs
    ) -> DockingResult:
        """
        Perform docking calculation.

        Args:
            receptor_path: Path to receptor PDBQT file
            ligand_path: Path to ligand PDBQT file
            center: Box center coordinates (x, y, z)
            size: Box dimensions (x, y, z)
            output_path: Optional path for output file
            **kwargs: Engine-specific parameters

        Returns:
            DockingResult with poses and scores

        Raises:
            DockingError: If docking fails
        """
        pass

    @abstractmethod
    def parse_output(self, output_file: str) -> DockingResult:
        """
        Parse docking output file.

        Args:
            output_file: Path to output file

        Returns:
            DockingResult parsed from output

        Raises:
            DockingError: If parsing fails
        """
        pass

    @abstractmethod
    def get_version(self) -> str:
        """
        Get engine version string.

        Returns:
            Version string
        """
        pass

    def validate_inputs(
        self,
        receptor_path: str,
        ligand_path: str
    ) -> bool:
        """
        Validate input files exist and are readable.

        Args:
            receptor_path: Path to receptor file
            ligand_path: Path to ligand file

        Returns:
            True if valid

        Raises:
            DockingError: If validation fails
        """
        receptor = Path(receptor_path)
        ligand = Path(ligand_path)

        if not receptor.exists():
            raise DockingError(f"Receptor file not found: {receptor_path}")

        if not ligand.exists():
            raise DockingError(f"Ligand file not found: {ligand_path}")

        if not receptor.is_file():
            raise DockingError(f"Receptor path is not a file: {receptor_path}")

        if not ligand.is_file():
            raise DockingError(f"Ligand path is not a file: {ligand_path}")

        return True

    def calculate_rmsd(
        self,
        coords1: np.ndarray,
        coords2: np.ndarray
    ) -> float:
        """
        Calculate RMSD between two coordinate sets.

        Args:
            coords1: First coordinate set (Nx3)
            coords2: Second coordinate set (Nx3)

        Returns:
            RMSD in Angstroms
        """
        if coords1.shape != coords2.shape:
            raise DockingError(
                f"Coordinate shape mismatch: {coords1.shape} vs {coords2.shape}"
            )

        diff = coords1 - coords2
        rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
        return rmsd


class DockingError(Exception):
    """Raised when docking operation fails."""
    pass
