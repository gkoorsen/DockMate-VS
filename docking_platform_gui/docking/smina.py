"""
Smina docking engine implementation.

Smina is a fork of AutoDock Vina with enhanced scoring functions
and additional features for molecular docking.

References:
- Smina: Koes et al., J. Chem. Inf. Model., 2013
- AutoDock Vina: Eberhardt et al., J. Chem. Inf. Model., 2021
- Scoring functions: Quiroga & Villarreal, PLoS ONE, 2016

Key advantages of Smina:
- Custom scoring functions (including Vinardo, AD4)
- Flexible side-chain docking support
- Atom typing compatible with Vina
- Open source and actively maintained
"""

import subprocess
import tempfile
import time
from pathlib import Path
from typing import Optional, List
import numpy as np
import re

from loguru import logger

from docking_platform_gui.docking.base import (
    DockingEngine,
    DockingResult,
    DockingPose,
    DockingError,
    FlexibilityMode
)


class SminaDockingEngine(DockingEngine):
    """
    Smina docking engine for rigid and flexible docking.

    Supports:
    - Rigid receptor docking (default)
    - Flexible side-chain docking (specify residues)
    - Custom scoring functions
    - All standard Vina parameters

    Examples:
        >>> engine = SminaDockingEngine(exhaustiveness=16, cpu=4)
        >>> result = engine.dock(
        ...     receptor_path="protein.pdbqt",
        ...     ligand_path="ligand.pdbqt",
        ...     center=np.array([10.0, 20.0, 30.0]),
        ...     size=np.array([20.0, 20.0, 20.0])
        ... )
        >>> print(f"Best score: {result.best_score} kcal/mol")
    """

    def __init__(
        self,
        smina_binary: str = "smina",
        scoring_function: str = "vina",
        exhaustiveness: int = 8,
        num_modes: int = 9,
        energy_range: float = 3.0,
        cpu: int = 1,
        seed: Optional[int] = None,
        timeout_sec: int = 1200
    ):
        """
        Initialize Smina docking engine.

        Args:
            smina_binary: Path to smina executable (default: "smina" in PATH)
            scoring_function: Scoring function to use
                            Options: vina (default), vinardo, ad4_scoring
            exhaustiveness: Search exhaustiveness (8=standard, 16-32=thorough)
            num_modes: Maximum binding modes to generate
            energy_range: Max energy difference from best mode (kcal/mol)
            cpu: Number of CPU cores
            seed: Random seed for reproducibility
            timeout_sec: Docking timeout in seconds

        References:
            Scoring functions: Quiroga & Villarreal, PLoS ONE, 2016
            - Vina: General-purpose, fast
            - Vinardo: Improved for hydrophobic interactions
            - AD4: AutoDock 4 scoring (more weight on H-bonds)
        """
        super().__init__(
            exhaustiveness=exhaustiveness,
            num_modes=num_modes,
            energy_range=energy_range,
            cpu=cpu,
            seed=seed
        )

        self.smina_binary = smina_binary
        self.scoring_function = scoring_function
        self.timeout_sec = timeout_sec

        # Verify binary exists
        self._verify_binary()

        logger.info(f"Using scoring function: {scoring_function}")

    def _verify_binary(self):
        """Verify smina binary is available."""
        try:
            result = subprocess.run(
                [self.smina_binary, "--version"],
                capture_output=True,
                text=True,
                timeout=5
            )
            version = result.stdout.strip() if result.stdout else "unknown"
            logger.info(f"Smina version: {version}")
        except FileNotFoundError:
            raise DockingError(
                f"Smina binary not found: {self.smina_binary}. "
                "Install smina or provide correct path."
            )
        except Exception as e:
            logger.warning(f"Could not verify smina version: {e}")

    def dock(
        self,
        receptor_path: str,
        ligand_path: str,
        center: np.ndarray,
        size: np.ndarray,
        output_path: Optional[str] = None,
        flexible_residues: Optional[str] = None,
        **kwargs
    ) -> DockingResult:
        """
        Perform docking with Smina.

        Args:
            receptor_path: Path to receptor PDBQT file
            ligand_path: Path to ligand PDBQT file
            center: Box center (x, y, z) in Angstroms
            size: Box size (x, y, z) in Angstroms
            output_path: Optional output file path (temp file if None)
            flexible_residues: Optional path to flexible residues PDBQT
            **kwargs: Additional smina parameters

        Returns:
            DockingResult with docking poses

        Raises:
            DockingError: If docking fails
        """
        # Validate inputs
        self.validate_inputs(receptor_path, ligand_path)

        # Extract ligand ID from path
        ligand_id = Path(ligand_path).stem

        # Create temp output if needed
        if output_path is None:
            temp_file = tempfile.NamedTemporaryFile(
                mode='w',
                suffix='.pdbqt',
                delete=False
            )
            output_path = temp_file.name
            temp_file.close()
        else:
            # Ensure output directory exists
            Path(output_path).parent.mkdir(parents=True, exist_ok=True)

        # Build command
        command = self._build_command(
            receptor_path=receptor_path,
            ligand_path=ligand_path,
            output_path=output_path,
            center=center,
            size=size,
            flexible_residues=flexible_residues,
            **kwargs
        )

        # Execute docking
        logger.info(f"Docking {ligand_id}...")
        start_time = time.time()

        try:
            result = subprocess.run(
                command,
                capture_output=True,
                text=True,
                timeout=self.timeout_sec
            )
            runtime = time.time() - start_time

            # Log stdout/stderr for debugging
            if result.stdout:
                logger.debug(f"Smina stdout: {result.stdout}")
            if result.stderr:
                logger.debug(f"Smina stderr: {result.stderr}")

            if result.returncode != 0:
                error_msg = result.stderr if result.stderr else "Unknown error"
                logger.error(f"Smina failed: {error_msg}")
                return DockingResult(
                    ligand_id=ligand_id,
                    success=False,
                    error_message=error_msg,
                    runtime=runtime
                )

            # Parse output
            docking_result = self.parse_output(output_path)
            docking_result.ligand_id = ligand_id
            docking_result.runtime = runtime

            logger.info(
                f"Docking completed: {ligand_id}, "
                f"best_score={docking_result.best_score:.2f} kcal/mol, "
                f"runtime={runtime:.1f}s"
            )

            return docking_result

        except subprocess.TimeoutExpired:
            runtime = time.time() - start_time
            logger.error(f"Docking timeout for {ligand_id}")
            return DockingResult(
                ligand_id=ligand_id,
                success=False,
                error_message="Docking timeout",
                runtime=runtime
            )

        except Exception as e:
            runtime = time.time() - start_time
            logger.error(f"Docking error for {ligand_id}: {e}")
            return DockingResult(
                ligand_id=ligand_id,
                success=False,
                error_message=str(e),
                runtime=runtime
            )

    def _build_command(
        self,
        receptor_path: str,
        ligand_path: str,
        output_path: str,
        center: np.ndarray,
        size: np.ndarray,
        flexible_residues: Optional[str] = None,
        **kwargs
    ) -> List[str]:
        """
        Build smina command line.

        Args:
            receptor_path: Receptor PDBQT path
            ligand_path: Ligand PDBQT path
            output_path: Output PDBQT path
            center: Box center (x, y, z)
            size: Box size (x, y, z)
            flexible_residues: Optional flexible residues PDBQT
            **kwargs: Additional parameters

        Returns:
            Command as list of strings
        """
        command = [
            self.smina_binary,
            "--receptor", receptor_path,
            "--ligand", ligand_path,
            "--out", output_path,
            "--center_x", str(center[0]),
            "--center_y", str(center[1]),
            "--center_z", str(center[2]),
            "--size_x", str(size[0]),
            "--size_y", str(size[1]),
            "--size_z", str(size[2]),
            "--exhaustiveness", str(self.exhaustiveness),
            "--num_modes", str(self.num_modes),
            "--energy_range", str(self.energy_range),
            "--cpu", str(self.cpu)
        ]

        # Add scoring function if not default
        if self.scoring_function != "vina":
            command.extend(["--scoring", self.scoring_function])

        # Add seed if specified
        if self.seed is not None:
            command.extend(["--seed", str(self.seed)])

        # Add flexible residues if specified
        if flexible_residues is not None:
            command.extend(["--flex", flexible_residues])

        # Add any additional parameters
        for key, value in kwargs.items():
            command.extend([f"--{key}", str(value)])

        logger.debug(f"Smina command: {' '.join(command)}")

        return command

    def parse_output(self, output_file: str) -> DockingResult:
        """
        Parse Smina output PDBQT file.

        Smina/Vina output format:
        - REMARK VINA RESULT: score rmsd_lb rmsd_ub
        - Multiple MODEL sections for each pose

        Args:
            output_file: Path to output PDBQT file

        Returns:
            DockingResult with parsed poses

        Raises:
            DockingError: If parsing fails
        """
        output_path = Path(output_file)

        if not output_path.exists():
            raise DockingError(f"Output file not found: {output_file}")

        try:
            with open(output_path, 'r') as f:
                content = f.read()

            poses = []
            mode_number = 0

            # Parse each MODEL section
            models = content.split('MODEL')[1:]  # Skip before first MODEL

            for model_text in models:
                mode_number += 1

                # Extract REMARK line with scores
                # Try Vina format: REMARK VINA RESULT:    -9.2      0.000      0.000
                match = re.search(
                    r'REMARK\s+VINA\s+RESULT:\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)',
                    model_text
                )

                if match:
                    score = float(match.group(1))
                    rmsd_lb = float(match.group(2))
                    rmsd_ub = float(match.group(3))
                else:
                    # Try Smina format: REMARK minimizedAffinity -3.83267069
                    match_smina = re.search(
                        r'REMARK\s+minimizedAffinity\s+([-\d.]+)',
                        model_text
                    )
                    if match_smina:
                        score = float(match_smina.group(1))
                        rmsd_lb = 0.0  # Smina doesn't provide RMSD in output
                        rmsd_ub = 0.0
                    else:
                        # Skip models without recognizable score format
                        continue

                # Extract coordinates (ATOM lines)
                coordinates = self._extract_coordinates(model_text)

                # Extract PDBQT for this model
                pdbqt = f"MODEL {mode_number}\n{model_text}"

                pose = DockingPose(
                    score=score,
                    rmsd_lb=rmsd_lb,
                    rmsd_ub=rmsd_ub,
                    coordinates=coordinates,
                    pdbqt=pdbqt,
                    mode_number=mode_number
                )
                poses.append(pose)

            if not poses:
                raise DockingError("No poses found in output file")

            # Create result
            result = DockingResult(
                ligand_id="parsed",  # Will be set by caller
                poses=poses,
                success=True
            )

            logger.debug(f"Parsed {len(poses)} poses from {output_file}")

            return result

        except Exception as e:
            logger.error(f"Failed to parse output file {output_file}: {e}")
            raise DockingError(f"Output parsing failed: {e}")

    def _extract_coordinates(self, model_text: str) -> np.ndarray:
        """
        Extract atom coordinates from MODEL text.

        Args:
            model_text: Text of MODEL section

        Returns:
            Nx3 array of coordinates
        """
        coords = []

        for line in model_text.split('\n'):
            if line.startswith('ATOM') or line.startswith('HETATM'):
                # PDBQT format: positions are columns 31-54
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    coords.append([x, y, z])
                except (ValueError, IndexError):
                    continue

        if coords:
            return np.array(coords)
        else:
            return np.array([]).reshape(0, 3)

    def dock_with_seed_variation(
        self,
        ligand_pdbqt: Path,
        receptor_pdbqt: Path,
        output_dir: Path,
        center: np.ndarray,
        box_size: np.ndarray,
        seeds: List[int],
        ligand_id: str = "ligand",
        flexible_residues: Optional[List[str]] = None,
        selection_mode: str = "best_score",
        **kwargs
    ) -> DockingResult:
        """
        Dock with multiple random seeds and select best result.

        This explores different search trajectories to find better poses.

        Args:
            ligand_pdbqt: Path to ligand PDBQT file
            receptor_pdbqt: Path to receptor PDBQT file
            output_dir: Directory for output files
            center: Box center coordinates
            box_size: Box dimensions
            seeds: List of random seeds to try
            ligand_id: Ligand identifier
            flexible_residues: Optional flexible residues
            selection_mode: How to select best result ("best_score" or "best_rmsd")
            **kwargs: Additional docking parameters

        Returns:
            Best DockingResult across all seeds
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        logger.info(f"Docking {ligand_id} with {len(seeds)} different seeds")

        all_results = []

        for i, seed in enumerate(seeds, 1):
            logger.debug(f"  Seed {i}/{len(seeds)}: {seed}")

            # Output file for this seed
            output_pdbqt = output_dir / f"{ligand_id}_seed{seed}.pdbqt"

            # Dock with this seed
            self.seed = seed  # Update seed

            try:
                result = self.dock(
                    ligand_path=str(ligand_pdbqt),
                    receptor_path=str(receptor_pdbqt),
                    output_path=str(output_pdbqt),
                    center=center,
                    size=box_size,
                    flexible_residues=flexible_residues,
                    **kwargs
                )

                if result.success:
                    result.seed_used = seed
                    all_results.append(result)
                    logger.debug(f"    Score: {result.best_score:.2f} kcal/mol")

            except Exception as e:
                logger.warning(f"    Failed with seed {seed}: {e}")
                continue

        if not all_results:
            raise DockingError("All seeds failed to produce results")

        # Select best result
        if selection_mode == "best_score":
            best_result = min(all_results, key=lambda r: r.best_score)
            logger.info(f"✅ Best result from seed {best_result.seed_used}: "
                       f"{best_result.best_score:.2f} kcal/mol")

        elif selection_mode == "best_rmsd":
            # Would need reference structure to calculate RMSD
            logger.warning("best_rmsd mode requires reference structure, using best_score")
            best_result = min(all_results, key=lambda r: r.best_score)

        else:
            raise ValueError(f"Unknown selection mode: {selection_mode}")

        # Add metadata
        best_result.ligand_id = ligand_id
        best_result.num_seeds_tried = len(seeds)
        best_result.num_seeds_successful = len(all_results)

        return best_result

    def get_version(self) -> str:
        """
        Get Smina version.

        Returns:
            Version string
        """
        try:
            result = subprocess.run(
                [self.smina_binary, "--version"],
                capture_output=True,
                text=True,
                timeout=5
            )
            return result.stdout.strip() if result.stdout else "unknown"
        except Exception:
            return "unknown"
