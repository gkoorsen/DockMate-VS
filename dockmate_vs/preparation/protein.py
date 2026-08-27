"""Configurable protein preparation for docking campaigns.

The pipeline uses PDBFixer/OpenMM when available, optional PROPKA and Reduce
steps, configurable water handling, and Open Babel PDBQT conversion. Flexible
receptor preparation additionally requires Meeko's ``mk_prepare_receptor.py``.
"""

from typing import List, Optional, Tuple, Set
from pathlib import Path
from dataclasses import dataclass
import tempfile
import subprocess
import copy
import shutil
import logging

import numpy as np
from loguru import logger

try:
    from Bio import PDB
    from Bio.PDB import PDBIO, Select
except ImportError:
    logger.error("Biopython not installed")
    raise

try:
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile
except ImportError:
    # Keep analysis and GUI modules importable when optional preparation tools
    # are unavailable; preparation methods degrade gracefully below.
    PDBFixer = None
    PDBFile = None
    logger.warning("PDBFixer or OpenMM not installed; protein fixing is unavailable")

from dockmate_vs.config.schema import (
    ProteinPreparationConfig,
    WaterHandling,
    WaterRetentionConfig
)


@dataclass
class PreparedProtein:
    """Container for prepared protein data."""

    protein_id: str
    structure: PDB.Structure.Structure
    pdbqt_path: Optional[str] = None
    pdb_path: Optional[str] = None
    structure_no_h: Optional[PDB.Structure.Structure] = None  # Clean structure without hydrogens
    waters_retained: int = 0
    waters_removed: int = 0
    missing_residues_added: List[str] = None
    pka_results: Optional[dict] = None

    def __post_init__(self):
        if self.missing_residues_added is None:
            self.missing_residues_added = []


class WaterSelector(Select):
    """BioPython selector for filtering water molecules."""

    def __init__(self, waters_to_keep: Set[Tuple[str, int]]):
        """
        Initialize water selector.

        Args:
            waters_to_keep: Set of (chain_id, residue_number) tuples
        """
        self.waters_to_keep = waters_to_keep

    def accept_residue(self, residue):
        """Accept residue if it's not water or if it's in keep list."""
        if residue.get_resname() in ['HOH', 'WAT']:
            chain_id = residue.get_parent().get_id()
            res_num = residue.get_id()[1]
            return (chain_id, res_num) in self.waters_to_keep
        return True  # Keep all non-water residues


class ProteinPreparationError(Exception):
    """Raised when protein preparation fails."""
    pass


class ProteinPreparation:
    """
    Protein preparation pipeline with customizable water handling.

    Provides three explicit water handling modes that must be validated for the
    target-specific protocol:
    - remove_all: Remove all crystallographic waters (default, fastest)
    - retain_all: Keep all crystallographic waters
    - selective: Retention based on configured structural criteria

    Example:
        >>> config = ProteinPreparationConfig(
        ...     water_handling=WaterHandling.SELECTIVE,
        ...     water_retention=WaterRetentionConfig(
        ...         distance_to_protein=3.0,
        ...         distance_to_ligand_site=3.0
        ...     )
        ... )
        >>> prep = ProteinPreparation(config)
        >>> protein = prep.prepare_receptor("protein.pdb")
    """

    def __init__(self, config: ProteinPreparationConfig):
        """
        Initialize protein preparation.

        Args:
            config: Configuration for protein preparation
        """
        self.config = config
        self.parser = PDB.PDBParser(QUIET=True)
        logger.info(f"ProteinPreparation initialized: pH {config.ph}, "
                   f"water_handling={config.water_handling.value}")

    def prepare_receptor(
        self,
        pdb_path: str,
        ligand_site_center: Optional[np.ndarray] = None,
        protein_id: Optional[str] = None
    ) -> PreparedProtein:
        """
        Full protein preparation pipeline.

        Args:
            pdb_path: Path to input PDB file
            ligand_site_center: Optional center of ligand binding site (for selective water retention)
            protein_id: Optional protein identifier

        Returns:
            PreparedProtein object

        Raises:
            ProteinPreparationError: If preparation fails
        """
        if protein_id is None:
            protein_id = Path(pdb_path).stem

        logger.info(f"Preparing protein {protein_id}: {pdb_path}")

        try:
            # 1. Load PDB structure
            structure = self._load_pdb(pdb_path)

            # 2. Fix structure with PDBFixer
            structure = self._fix_structure(pdb_path, structure)

            # 3. Handle water molecules
            structure, waters_retained, waters_removed = self._handle_waters(
                structure,
                ligand_site_center
            )

            # 3b. Strip requested heteroatoms (normally the co-crystal ligand
            # that defines the docking site). Must happen before hydrogens are
            # added and before PDBQT conversion, or the ligand ends up in the
            # receptor and blocks its own pocket.
            if getattr(self.config, "remove_hetero_residues", None):
                structure = self._remove_hetero_residues(
                    structure,
                    self.config.remove_hetero_residues
                )

            # Save clean structure without hydrogens (for PDBQT conversion)
            # Deep copy before hydrogen addition
            structure_no_h = copy.deepcopy(structure)

            # 4. Predict pKa with PROPKA (if available)
            pka_results = self._predict_pka(structure)

            # 5. Add hydrogens at target pH
            structure = self._add_hydrogens(structure, pka_results)

            # 6. Optimize H-bonds with Reduce (if available)
            if self.config.optimize_hbonds:
                structure = self._optimize_hbonds(structure)

            # 7. Validate structure (if enabled)
            if self.config.validate_structure:
                self._validate_structure(structure)

            # Create result
            result = PreparedProtein(
                protein_id=protein_id,
                structure=structure,
                structure_no_h=structure_no_h,  # Save clean structure for PDBQT conversion
                waters_retained=waters_retained,
                waters_removed=waters_removed,
                pka_results=pka_results
            )

            logger.info(f"Protein preparation complete: {protein_id}")
            logger.info(f"  Waters: {waters_retained} retained, {waters_removed} removed")

            return result

        except Exception as e:
            logger.error(f"Protein preparation failed for {protein_id}: {e}")
            raise ProteinPreparationError(f"Failed to prepare {protein_id}: {e}")

    def _load_pdb(self, pdb_path: str) -> PDB.Structure.Structure:
        """Load PDB structure."""
        if not Path(pdb_path).exists():
            raise ProteinPreparationError(f"PDB file not found: {pdb_path}")

        structure = self.parser.get_structure("protein", pdb_path)
        return structure

    def _fix_structure(
        self,
        pdb_path: str,
        structure: PDB.Structure.Structure
    ) -> PDB.Structure.Structure:
        """
        Fix structure using PDBFixer.

        Adds missing atoms and models missing loops.

        Args:
            pdb_path: Path to original PDB
            structure: BioPython structure

        Returns:
            Fixed structure
        """
        logger.info("Fixing structure with PDBFixer")

        try:
            # Initialize PDBFixer
            fixer = PDBFixer(filename=pdb_path)

            # Find missing residues
            fixer.findMissingResidues()

            # Add missing atoms
            if self.config.add_missing_atoms:
                fixer.findMissingAtoms()
                fixer.addMissingAtoms()
                logger.debug("Added missing atoms")

            # Model missing loops (can be slow)
            if self.config.model_missing_loops:
                missing_residues = fixer.missingResidues
                if missing_residues:
                    if not hasattr(fixer, "addMissingResidues"):
                        logger.info("PDBFixer does not support addMissingResidues; skipping loop modeling")
                    else:
                        logger.info(f"Modeling {len(missing_residues)} missing residues/loops")
                        # Note: This can be slow for large gaps
                        # In production, may want to add timeout or skip large gaps
                        try:
                            fixer.addMissingResidues()
                            logger.debug("Modeled missing loops")
                        except Exception as e:
                            logger.warning(f"Loop modeling failed: {e}")

            # Save fixed structure to temp file
            with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as tmp:
                PDBFile.writeFile(fixer.topology, fixer.positions, tmp)
                temp_path = tmp.name

            # Reload with BioPython
            structure = self.parser.get_structure("protein_fixed", temp_path)

            # Cleanup
            Path(temp_path).unlink()

            return structure

        except Exception as e:
            logger.error(f"PDBFixer failed: {e}")
            # Return original structure if fixing fails
            return structure

    def _handle_waters(
        self,
        structure: PDB.Structure.Structure,
        ligand_site_center: Optional[np.ndarray] = None
    ) -> Tuple[PDB.Structure.Structure, int, int]:
        """
        Handle water molecules according to configuration.

        Implements three modes:
        - remove_all: Remove all waters (default)
        - retain_all: Keep all waters
        - selective: Keep structural waters based on criteria

        Args:
            structure: Input structure
            ligand_site_center: Center of binding site (for selective mode)

        Returns:
            (modified_structure, waters_retained, waters_removed)
        """
        mode = self.config.water_handling

        logger.info(f"Water handling mode: {mode.value}")

        # Count initial waters
        initial_waters = sum(
            1 for residue in structure.get_residues()
            if residue.get_resname() in ['HOH', 'WAT']
        )

        if initial_waters == 0:
            logger.info("No water molecules found")
            return structure, 0, 0

        if mode == WaterHandling.REMOVE_ALL:
            # Remove all waters
            structure = self._remove_all_waters(structure)
            return structure, 0, initial_waters

        elif mode == WaterHandling.RETAIN_ALL:
            # Keep all waters
            logger.info(f"Retaining all {initial_waters} water molecules")
            return structure, initial_waters, 0

        elif mode == WaterHandling.SELECTIVE:
            # Selective retention
            if not self.config.water_retention:
                logger.warning("Selective mode but no retention config, using defaults")
                self.config.water_retention = WaterRetentionConfig()

            structure, retained = self._selective_water_retention(
                structure,
                ligand_site_center
            )
            removed = initial_waters - retained
            return structure, retained, removed

        else:
            raise ValueError(f"Unknown water handling mode: {mode}")

    def _remove_all_waters(
        self,
        structure: PDB.Structure.Structure
    ) -> PDB.Structure.Structure:
        """Remove all water molecules from structure."""
        logger.info("Removing all water molecules")

        for model in structure:
            for chain in model:
                residues_to_remove = []
                for residue in chain:
                    if residue.get_resname() in ['HOH', 'WAT']:
                        residues_to_remove.append(residue.get_id())

                for res_id in residues_to_remove:
                    chain.detach_child(res_id)

        return structure

    def _remove_hetero_residues(
        self,
        structure: PDB.Structure.Structure,
        resnames: List[str]
    ) -> PDB.Structure.Structure:
        """Remove named heteroatom residues (e.g. the co-crystal site ligand).

        Docking into a receptor that still contains its own co-crystal ligand
        silently produces nonsense: the pocket is occupied, so poses are pushed
        to the periphery and the crystal ligand cannot even reproduce its own
        geometry. Comparison is case-insensitive because PDB resnames vary in
        case and preparation may rewrite HETATM records as ATOM.
        """
        wanted = {str(r).strip().upper() for r in resnames if str(r).strip()}
        if not wanted:
            return structure

        removed = 0
        for model in structure:
            for chain in model:
                residues_to_remove = [
                    residue.get_id()
                    for residue in chain
                    if residue.get_resname().strip().upper() in wanted
                ]
                for res_id in residues_to_remove:
                    chain.detach_child(res_id)
                    removed += 1

        if removed:
            logger.info(
                f"Removed {removed} heteroatom residue(s) from receptor: "
                f"{sorted(wanted)}"
            )
        else:
            logger.warning(
                f"Requested removal of {sorted(wanted)} but no matching "
                f"residue was found in the receptor"
            )

        return structure

    def _selective_water_retention(
        self,
        structure: PDB.Structure.Structure,
        ligand_site_center: Optional[np.ndarray] = None
    ) -> Tuple[PDB.Structure.Structure, int]:
        """
        Selectively retain structural waters based on criteria.

        Keeps waters that:
        1. Are within distance threshold of protein atoms
        2. Are within distance threshold of ligand site (if provided)
        3. Have sufficient contacts
        4. Meet B-factor/occupancy criteria

        Args:
            structure: Input structure
            ligand_site_center: Optional center of binding site

        Returns:
            (modified_structure, number_of_retained_waters)
        """
        config = self.config.water_retention
        logger.info(f"Selective water retention: "
                   f"protein_dist={config.distance_to_protein}Å, "
                   f"site_dist={config.distance_to_ligand_site}Å")

        # Get all protein atoms (non-water heavy atoms)
        protein_atoms = []
        for atom in structure.get_atoms():
            residue = atom.get_parent()
            if residue.get_resname() not in ['HOH', 'WAT'] and atom.element != 'H':
                protein_atoms.append(atom)

        protein_coords = np.array([atom.get_coord() for atom in protein_atoms])

        # Identify waters to keep
        waters_to_keep = set()

        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_resname() not in ['HOH', 'WAT']:
                        continue

                    # Get water oxygen coordinate
                    try:
                        water_o = residue['O'].get_coord()
                    except KeyError:
                        continue  # Skip if no oxygen atom

                    # Check distance to protein
                    protein_distances = np.linalg.norm(
                        protein_coords - water_o,
                        axis=1
                    )
                    min_protein_dist = np.min(protein_distances)

                    # Count contacts within threshold
                    protein_contacts = np.sum(
                        protein_distances <= config.distance_to_protein
                    )

                    # Check criteria
                    keep_water = True

                    # Criterion 1: Distance to protein
                    if min_protein_dist > config.distance_to_protein:
                        keep_water = False

                    # Criterion 2: Sufficient contacts
                    if protein_contacts < config.min_contacts:
                        keep_water = False

                    # Criterion 3: Distance to ligand site (if provided)
                    if ligand_site_center is not None:
                        site_dist = np.linalg.norm(water_o - ligand_site_center)
                        if site_dist > config.distance_to_ligand_site:
                            keep_water = False

                    # Criterion 4: B-factor (if specified)
                    if config.max_bfactor is not None:
                        try:
                            bfactor = residue['O'].get_bfactor()
                            if bfactor > config.max_bfactor:
                                keep_water = False
                        except:
                            pass

                    # Criterion 5: Occupancy
                    try:
                        occupancy = residue['O'].get_occupancy()
                        if occupancy < config.min_occupancy:
                            keep_water = False
                    except:
                        pass

                    if keep_water:
                        chain_id = chain.get_id()
                        res_num = residue.get_id()[1]
                        waters_to_keep.add((chain_id, res_num))

        # Remove waters not in keep list
        for model in structure:
            for chain in model:
                residues_to_remove = []
                for residue in chain:
                    if residue.get_resname() in ['HOH', 'WAT']:
                        chain_id = chain.get_id()
                        res_num = residue.get_id()[1]
                        if (chain_id, res_num) not in waters_to_keep:
                            residues_to_remove.append(residue.get_id())

                for res_id in residues_to_remove:
                    chain.detach_child(res_id)

        retained_count = len(waters_to_keep)
        logger.info(f"Retained {retained_count} structural waters")

        return structure, retained_count

    def _predict_pka(
        self,
        structure: PDB.Structure.Structure
    ) -> Optional[dict]:
        """
        Predict pKa values using PROPKA.

        Args:
            structure: Input structure

        Returns:
            Dictionary of pKa predictions (or None if PROPKA unavailable)
        """
        try:
            import propka.run as propka_run
            from propka.molecular_container import MolecularContainer

            logger.info("Running PROPKA for pKa prediction")

            # Save structure to temp file (PROPKA handles single-model best)
            models = list(structure.get_models())
            with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as tmp:
                io = PDBIO()
                io.set_structure(structure)
                if len(models) > 1:
                    first_model_id = models[0].id

                    class _FirstModelSelect(PDB.Select):
                        def accept_model(self, model):  # noqa: D401
                            return model.id == first_model_id

                    logger.info(
                        "Multiple models detected (%d); using first model for PROPKA",
                        len(models)
                    )
                    io.save(tmp.name, select=_FirstModelSelect())
                else:
                    io.save(tmp.name)
                temp_path = tmp.name

            # Run PROPKA
            try:
                propka_logger = logging.getLogger("propka")
                previous_level = propka_logger.level
                propka_logger.setLevel(logging.ERROR)
                try:
                    import contextlib
                    import io
                    import signal

                    class _PropkaTimeout(Exception):
                        pass

                    def _alarm_handler(signum, frame):  # noqa: D401
                        raise _PropkaTimeout("PROPKA timed out")

                    alarm_set = False
                    if hasattr(signal, "SIGALRM"):
                        signal.signal(signal.SIGALRM, _alarm_handler)
                        signal.alarm(120)
                        alarm_set = True

                    stdout_buf = io.StringIO()
                    stderr_buf = io.StringIO()
                    with contextlib.redirect_stdout(stdout_buf), contextlib.redirect_stderr(stderr_buf):
                        molecule = propka_run.single(
                            temp_path,
                            optargs=['--quiet'],
                            write_pka=False
                        )
                    if stderr_buf.getvalue().strip():
                        logger.debug("PROPKA stderr: %s", stderr_buf.getvalue().strip())
                    if alarm_set:
                        signal.alarm(0)
                finally:
                    if "alarm_set" in locals() and alarm_set and hasattr(signal, "SIGALRM"):
                        signal.alarm(0)
                    propka_logger.setLevel(previous_level)

                # Extract pKa results
                pka_results = {}
                groups = []
                if hasattr(molecule, "get_titratable_groups"):
                    groups = molecule.get_titratable_groups()
                elif hasattr(molecule, "conformations"):
                    conf = None
                    if isinstance(molecule.conformations, dict):
                        conf = molecule.conformations.get("AVR")
                        if conf is None:
                            conf = next(iter(molecule.conformations.values()), None)
                    elif molecule.conformations:
                        conf = molecule.conformations[0]

                    if conf is not None and hasattr(conf, "get_titratable_groups"):
                        groups = conf.get_titratable_groups()

                def _group_key(group) -> str:
                    residue_type = getattr(group, "residue_type", None)
                    residue_number = getattr(group, "residue_number", None)
                    chain_id = getattr(group, "chain_id", None)
                    if residue_type and residue_number is not None and chain_id:
                        return f"{residue_type}{residue_number}{chain_id}"

                    label = getattr(group, "label", "")
                    if label:
                        parts = label.split()
                        if len(parts) >= 3:
                            return f"{parts[0]}{parts[1]}{parts[2]}"
                        return label.replace(" ", "")

                    return "UNKNOWN"

                for group in groups:
                    residue_key = _group_key(group)
                    pka_results[residue_key] = {
                        'pka': getattr(group, "pka_value", None),
                        'charge': getattr(group, "charge", None)
                    }

                logger.debug(f"PROPKA predicted {len(pka_results)} pKa values")

                # Cleanup
                Path(temp_path).unlink()

                return pka_results

            except Exception as e:
                error_msg = str(e)
                if "group_type" in error_msg and "NoneType" in error_msg:
                    logger.debug("PROPKA failed (known issue): %s", error_msg)
                else:
                    logger.warning(f"PROPKA execution failed: {e}")
                Path(temp_path).unlink()
                return None

        except ImportError:
            logger.warning("PROPKA not available, skipping pKa prediction")
            return None

    def _add_hydrogens(
        self,
        structure: PDB.Structure.Structure,
        pka_results: Optional[dict] = None
    ) -> PDB.Structure.Structure:
        """
        Add hydrogens at target pH.

        Uses PDBFixer to add hydrogens with appropriate protonation states.

        Args:
            structure: Input structure
            pka_results: Optional pKa predictions from PROPKA

        Returns:
            Structure with hydrogens added
        """
        logger.info(f"Adding hydrogens at pH {self.config.ph}")

        # Save structure to temp file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as tmp:
            io = PDBIO()
            io.set_structure(structure)
            io.save(tmp.name)
            temp_path = tmp.name

        try:
            # Use PDBFixer to add hydrogens
            fixer = PDBFixer(filename=temp_path)
            fixer.addMissingHydrogens(self.config.ph)

            # Save with hydrogens
            with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as tmp2:
                PDBFile.writeFile(fixer.topology, fixer.positions, tmp2)
                temp_path2 = tmp2.name

            # Reload structure
            structure = self.parser.get_structure("protein_h", temp_path2)

            # Cleanup
            Path(temp_path).unlink()
            Path(temp_path2).unlink()

            return structure

        except Exception as e:
            logger.error(f"Hydrogen addition failed: {e}")
            Path(temp_path).unlink()
            return structure

    def _optimize_hbonds(
        self,
        structure: PDB.Structure.Structure
    ) -> PDB.Structure.Structure:
        """
        Optimize hydrogen bonds using Reduce.

        Flips Asn/Gln/His sidechains for optimal H-bonding.

        Args:
            structure: Input structure

        Returns:
            Optimized structure
        """
        logger.info("Optimizing H-bonds with Reduce")

        # Check if Reduce is available
        reduce_path = shutil.which("reduce")
        if not reduce_path:
            logger.info("Reduce not available, skipping H-bond optimization")
            return structure

        # Save structure to temp file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as tmp:
            io = PDBIO()
            io.set_structure(structure)
            io.save(tmp.name)
            temp_path = tmp.name

        try:
            # Run Reduce
            with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as tmp_out:
                output_path = tmp_out.name

            result = subprocess.run(
                [reduce_path, '-FLIP', '-Quiet', temp_path],
                capture_output=True,
                text=True,
                timeout=60
            )

            if result.returncode == 0:
                # Write output
                with open(output_path, 'w') as f:
                    f.write(result.stdout)

                # Reload structure
                structure = self.parser.get_structure("protein_reduced", output_path)
                logger.debug("H-bond optimization complete")

                # Cleanup
                Path(output_path).unlink()
            else:
                logger.warning(f"Reduce failed: {result.stderr}")

            Path(temp_path).unlink()

        except Exception as e:
            logger.error(f"Reduce execution failed: {e}")
            Path(temp_path).unlink()

        return structure

    def _validate_structure(self, structure: PDB.Structure.Structure) -> None:
        """
        Validate structure quality.

        Logs warnings for potential issues.

        Args:
            structure: Structure to validate
        """
        logger.info("Validating structure")

        # Count atoms
        num_atoms = sum(1 for _ in structure.get_atoms())
        logger.debug(f"Total atoms: {num_atoms}")

        # Check for missing atoms (CA atoms as proxy)
        num_residues = sum(1 for _ in structure.get_residues())
        num_ca = sum(1 for atom in structure.get_atoms() if atom.get_name() == 'CA')

        if num_ca < num_residues * 0.9:
            logger.warning(f"Potential missing atoms: {num_ca} CA atoms for {num_residues} residues")

        # More validation could be added here (MolProbity, etc.)

    def to_pdbqt(
        self,
        prepared_protein: PreparedProtein,
        output_path: str
    ) -> str:
        """
        Convert a prepared rigid protein to PDBQT format using Open Babel.

        Args:
            prepared_protein: Prepared protein object
            output_path: Output path for PDBQT file

        Returns:
            Path to PDBQT file
        """
        logger.info(f"Converting to PDBQT: {output_path}")

        # Use the clean structure without hydrogens if available
        # This avoids valence issues from hydrogen manipulation
        structure_for_pdbqt = prepared_protein.structure_no_h if prepared_protein.structure_no_h else prepared_protein.structure

        # Create selector to filter alternate locations
        class CleanProteinSelect(Select):
            def accept_atom(self, atom):
                # Skip hydrogen atoms if they somehow got into structure_no_h
                if atom.element == 'H':
                    return False

                # Skip alternate locations (keep only first)
                if atom.is_disordered():
                    altloc = atom.get_altloc()
                    if altloc != ' ' and altloc != 'A':
                        return False

                return True

        # Save cleaned PDB to temp file
        # For NMR structures with multiple models, use only the first model
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as tmp:
            io = PDBIO()
            # If structure has multiple models, use only first model
            if len(list(structure_for_pdbqt.get_models())) > 1:
                logger.debug(f"Multiple models detected, using only first model")
                io.set_structure(structure_for_pdbqt[0])  # First model only
            else:
                io.set_structure(structure_for_pdbqt)
            io.save(tmp.name, CleanProteinSelect())
            temp_pdb = tmp.name

        logger.debug(f"Saved clean structure to: {temp_pdb}")

        try:
            # Use Open Babel for PDBQT conversion (more forgiving than Meeko)
            # The -xr flag preserves charges and makes it rigid
            result = subprocess.run(
                [
                    'obabel',
                    temp_pdb,
                    '-O', output_path,
                    '-xr'  # Make receptor rigid (required for docking)
                ],
                capture_output=True,
                timeout=30,
                text=True
            )

            # Check if conversion succeeded
            if result.returncode != 0:
                error_msg = f"Open Babel conversion failed: {result.stderr}"
                raise ProteinPreparationError(error_msg)

            if not Path(output_path).exists():
                raise ProteinPreparationError("PDBQT file was not created")

            # Verify the output has content
            if Path(output_path).stat().st_size == 0:
                raise ProteinPreparationError("PDBQT file is empty")

            logger.info(f"PDBQT saved: {output_path}")
            prepared_protein.pdbqt_path = output_path

            # Cleanup
            Path(temp_pdb).unlink()

            return output_path

        except subprocess.TimeoutExpired:
            logger.error("PDBQT conversion timed out")
            Path(temp_pdb).unlink()
            raise ProteinPreparationError("PDBQT conversion timed out after 30 seconds")
        except Exception as e:
            logger.error(f"PDBQT conversion failed: {e}")
            if Path(temp_pdb).exists():
                Path(temp_pdb).unlink()
            raise ProteinPreparationError(f"PDBQT conversion failed: {e}")

    def prepare_flexible_receptor(
        self,
        prepared_protein: PreparedProtein,
        flexible_residues: List[str],
        output_receptor_path: str,
        output_flex_path: str
    ) -> Tuple[str, str]:
        """
        Prepare flexible receptor with flexible side-chains.

        Uses Meeko to create two PDBQT files:
        1. Rigid receptor (protein minus flexible residues)
        2. Flexible residues (will move during docking)

        Args:
            prepared_protein: Prepared protein object
            flexible_residues: List of residue specifications (e.g., ["A:315", "A:318"])
                              Format: chain:residue_number
            output_receptor_path: Path for rigid receptor PDBQT
            output_flex_path: Path for flexible residues PDBQT

        Returns:
            Tuple of (receptor_pdbqt_path, flex_pdbqt_path)

        Raises:
            ProteinPreparationError: If flexible preparation fails

        Flexible-receptor performance is target dependent and must be validated
        against the campaign's own pose and enrichment controls.
        """
        logger.info(f"Preparing flexible receptor with {len(flexible_residues)} flexible residues")
        logger.info(f"Flexible residues: {', '.join(flexible_residues)}")

        # Validate residue format
        for res in flexible_residues:
            if ':' not in res:
                raise ProteinPreparationError(
                    f"Invalid residue format '{res}'. Expected format: 'chain:number' (e.g., 'A:315')"
                )

        # Save PDB to temp file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as tmp:
            io = PDBIO()
            io.set_structure(prepared_protein.structure)
            io.save(tmp.name)
            temp_pdb = tmp.name

        try:
            # Convert flexible residues list to Meeko format
            # Meeko expects: A:315 A:318 (space-separated)
            flex_string = ' '.join(flexible_residues)

            # Use mk_prepare_receptor.py with -f flag for flexible residues
            # This creates two files:
            # - output_receptor_path: rigid receptor
            # - output_flex_path: flexible residues
            command = [
                'mk_prepare_receptor.py',
                '-i', temp_pdb,
                '-o', output_receptor_path,
                '-f', flex_string,
                '-g', output_flex_path  # -g specifies the flexible residues output
            ]

            logger.debug(f"Running command: {' '.join(command)}")

            result = subprocess.run(
                command,
                capture_output=True,
                text=True,
                timeout=60
            )

            if result.returncode != 0:
                error_msg = result.stderr if result.stderr else "Unknown error"
                raise ProteinPreparationError(
                    f"Meeko flexible preparation failed: {error_msg}"
                )

            # Check both output files exist
            if not Path(output_receptor_path).exists():
                raise ProteinPreparationError(
                    f"Rigid receptor PDBQT not created: {output_receptor_path}"
                )

            if not Path(output_flex_path).exists():
                raise ProteinPreparationError(
                    f"Flexible residues PDBQT not created: {output_flex_path}"
                )

            logger.info(f"Flexible receptor prepared:")
            logger.info(f"  Rigid: {output_receptor_path}")
            logger.info(f"  Flex:  {output_flex_path}")

            # Cleanup
            Path(temp_pdb).unlink()

            return output_receptor_path, output_flex_path

        except subprocess.TimeoutExpired:
            logger.error("Meeko flexible preparation timed out")
            Path(temp_pdb).unlink()
            raise ProteinPreparationError("Flexible preparation timed out")

        except Exception as e:
            logger.error(f"Flexible preparation failed: {e}")
            if Path(temp_pdb).exists():
                Path(temp_pdb).unlink()
            raise ProteinPreparationError(f"Flexible preparation failed: {e}")
