"""Ligand preparation with RDKit and Open Babel.

The current pipeline performs largest-fragment selection and charge
normalization, RDKit tautomer and conformer generation, optional MMFF94s
minimization, and Open Babel PDBQT conversion. It does not claim pH-aware
ionization-state enumeration.
"""

from typing import List, Optional, Tuple
from pathlib import Path
from dataclasses import dataclass
from loguru import logger
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
from multiprocessing import Pool, cpu_count, TimeoutError as MPTimeoutError
from functools import partial
import os
import time

from dockmate_vs.config.schema import LigandPreparationConfig


# Helper functions for parallel RMSD calculation
def _calculate_rmsd_pair(args):
    """
    Calculate RMSD between two conformers.
    
    Args:
        args: Tuple of (mol_binary, conf_id_i, conf_id_j)
    
    Returns:
        Tuple of (i, j, rmsd)
    """
    mol_binary, i, j, conf_id_i, conf_id_j = args
    
    try:
        # Deserialize molecule from binary
        mol = Chem.Mol(mol_binary)
        
        # Calculate RMSD
        rmsd = AllChem.GetBestRMS(mol, mol, conf_id_i, conf_id_j)
        return (i, j, rmsd)
    except Exception as e:
        logger.warning(f"RMSD calculation failed for conformers {i},{j}: {e}")
        return (i, j, 999.9)


def _calculate_rmsd_pairs_batch(mol_binary, pairs_with_indices):
    """
    Calculate multiple RMSD pairs in a batch (for one worker).
    
    Args:
        mol_binary: Serialized RDKit molecule
        pairs_with_indices: List of (i, j, conf_id_i, conf_id_j) tuples
    
    Returns:
        List of (i, j, rmsd) tuples
    """
    results = []
    mol = Chem.Mol(mol_binary)
    
    for i, j, conf_id_i, conf_id_j in pairs_with_indices:
        try:
            rmsd = AllChem.GetBestRMS(mol, mol, conf_id_i, conf_id_j)
            results.append((i, j, rmsd))
        except Exception:
            results.append((i, j, 999.9))
    
    return results


@dataclass
class PreparedLigand:
    """Container for prepared ligand data."""

    mol_id: str
    smiles: str
    mol: Chem.Mol
    pdbqt: Optional[str] = None
    conformer_id: int = 0
    protonation_state: Optional[str] = None
    tautomer_id: int = 0
    energy: Optional[float] = None

    def to_pdbqt_file(self, output_path: str) -> None:
        """Write PDBQT to file."""
        if not self.pdbqt:
            raise ValueError("PDBQT string not generated yet")

        with open(output_path, 'w') as f:
            f.write(self.pdbqt)


class LigandPreparationError(Exception):
    """Raised when ligand preparation fails."""
    pass


class LigandPreparation:
    """Convert SMILES strings into docking-ready PDBQT variants.

    Users should inspect charge-sensitive compounds because the standardization
    step normalizes neutralizable formal charges and does not enumerate pH-aware
    ionization states.

    Example:
        >>> config = LigandPreparationConfig()
        >>> prep = LigandPreparation(config)
        >>> ligands = prep.prepare_from_smiles("CCO", "ethanol")
        >>> ligands[0].to_pdbqt_file("ethanol.pdbqt")
    """

    def __init__(self, config: LigandPreparationConfig, n_cpus: Optional[int] = None):
        """
        Initialize ligand preparation.

        Args:
            config: Configuration for ligand preparation
            n_cpus: Number of CPUs to use for parallel processing. 
                   If None, auto-detects and caps at 8.
                   Set to 1 to disable parallel processing.
        """
        self.config = config
        self.uncharger = rdMolStandardize.Uncharger()
        
        # Handle CPU count
        if n_cpus is None:
            self.n_cpus = min(cpu_count(), 8)
            logger.info(f"LigandPreparation initialized: pH {config.ph_range}, "
                       f"max_tautomers={config.max_tautomers}, "
                       f"CPUs={self.n_cpus} (auto-detected)")
        else:
            self.n_cpus = max(1, int(n_cpus))  # Ensure at least 1
            logger.info(f"LigandPreparation initialized: pH {config.ph_range}, "
                       f"max_tautomers={config.max_tautomers}, "
                       f"CPUs={self.n_cpus} (user-specified)")

    def prepare_from_smiles(
        self,
        smiles: str,
        mol_id: str,
        enumerate_states: bool = True
    ) -> List[PreparedLigand]:
        """
        Full preparation pipeline from SMILES.

        Args:
            smiles: Input SMILES string
            mol_id: Molecule identifier
            enumerate_states: Whether to enumerate normalized-state tautomers

        Returns:
            List of prepared ligands (multiple if enumeration enabled)

        Raises:
            LigandPreparationError: If preparation fails
        """
        logger.info(f"Preparing ligand {mol_id}: {smiles}")

        try:
            # 1. Parse SMILES
            mol = self._parse_smiles(smiles)
            if mol is None:
                raise LigandPreparationError(f"Invalid SMILES: {smiles}")

            logger.info(f"  Standardizing molecule for {mol_id}")
            # 2. Standardize molecule
            mol = self._standardize_molecule(mol)

            logger.info(f"  Selecting normalized charge state for {mol_id}")
            # 3. Retain one normalized charge state before tautomer enumeration.
            if enumerate_states:
                protonated_mols = self._enumerate_protonation(mol)
            else:
                protonated_mols = [(mol, "standardized")]

            results = []

            for prot_idx, (prot_mol, prot_state) in enumerate(protonated_mols):
                # 4. Enumerate tautomers
                if enumerate_states and self.config.max_tautomers > 1:
                    tautomers = self._enumerate_tautomers(prot_mol)
                else:
                    tautomers = [prot_mol]

                for taut_idx, taut_mol in enumerate(tautomers):
                    tautomer_id = f"{mol_id}_p{prot_idx}_t{taut_idx}"
                    logger.info(f"    Generating 3D conformers for {tautomer_id}")
                    # 5. Generate 3D conformers
                    conf_mols = self._generate_3d_conformers(taut_mol)
                    logger.info(f"    Generated {len(conf_mols)} conformers for {tautomer_id}")
                    logger.info(f"    Next: Will process each conformer through MMFF minimization → PDBQT conversion")

                    for conf_idx, conf_mol in enumerate(conf_mols):
                        # 6. Energy minimization
                        variant_id = f"{mol_id}_p{prot_idx}_t{taut_idx}_c{conf_idx}"
                        logger.info(f"")
                        logger.info(f"    {'='*70}")
                        logger.info(f"    CONFORMER {conf_idx+1}/{len(conf_mols)}: {variant_id}")
                        logger.info(f"    ={'='*70}")
                        
                        # Step 1: Energy minimization
                        logger.info(f"    [STEP 1/3] MMFF94s Energy Minimization")
                        if self.config.mmff_minimize:
                            logger.info(f"      → Starting MMFF minimization for {variant_id}...")
                            minimized_mol, energy = self._minimize_mmff94s(conf_mol)
                            energy_str = f"{energy:.2f} kcal/mol" if energy is not None else "N/A"
                            logger.info(f"      ✓ Energy minimization complete: energy={energy_str}")
                        else:
                            minimized_mol = conf_mol
                            energy = None
                            energy_str = "N/A"
                            logger.info(f"      → Skipping minimization (disabled in config)")

                        # Step 2: PDBQT conversion
                        logger.info(f"    [STEP 2/3] PDBQT Conversion")
                        logger.info(f"      → About to convert {variant_id} to PDBQT format using Open Babel...")
                        logger.info(f"      → Molecule info: {minimized_mol.GetNumHeavyAtoms()} heavy atoms, "
                                   f"{minimized_mol.GetNumAtoms()} total atoms")
                        pdbqt = self._to_pdbqt(minimized_mol, mol_id)

                        pdbqt_lines = len(pdbqt.splitlines())
                        logger.info(f"      ✓ PDBQT conversion complete ({pdbqt_lines} lines)")

                        # Step 3: Create result object
                        logger.info(f"    [STEP 3/3] Creating PreparedLigand Result")
                        logger.info(f"      → Building result object for {variant_id}...")
                        result = PreparedLigand(
                            mol_id=variant_id,
                            smiles=Chem.MolToSmiles(minimized_mol),
                            mol=minimized_mol,
                            pdbqt=pdbqt,
                            conformer_id=conf_idx,
                            protonation_state=prot_state,
                            tautomer_id=taut_idx,
                            energy=energy
                        )
                        results.append(result)
                        logger.info(f"      ✓ Result object created and added to results list")
                        logger.info(f"    {'='*70}")
                        logger.info(f"")

            logger.info(f"Prepared {len(results)} variants for {mol_id}")
            return results

        except Exception as e:
            logger.error(f"Ligand preparation failed for {mol_id}: {e}")
            raise LigandPreparationError(f"Failed to prepare {mol_id}: {e}")

    def _parse_smiles(self, smiles: str) -> Optional[Chem.Mol]:
        """Parse SMILES string to RDKit molecule."""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.error(f"Failed to parse SMILES: {smiles}")
        return mol

    def _standardize_molecule(self, mol: Chem.Mol) -> Chem.Mol:
        """
        Standardize molecule (remove salts, neutralize).

        Args:
            mol: Input molecule

        Returns:
            Standardized molecule
        """
        # Remove salts (keep largest fragment)
        remover = rdMolStandardize.LargestFragmentChooser()
        mol = remover.choose(mol)

        # Normalize removable formal charges. This is not pH-aware ionization.
        mol = self.uncharger.uncharge(mol)

        return mol

    def _enumerate_protonation(
        self,
        mol: Chem.Mol
    ) -> List[Tuple[Chem.Mol, str]]:
        """Return the single normalized charge state used by this release.

        Earlier versions labelled unchanged molecule copies as protonated and
        deprotonated variants. Those copies were chemically identical and only
        repeated docking work. The pH range remains in the configuration for
        backward compatibility but is not interpreted as ionization enumeration.

        Args:
            mol: Input molecule

        Returns:
            List of (molecule, state_description) tuples
        """
        return [(Chem.Mol(mol), "standardized")]

    def _enumerate_tautomers(self, mol: Chem.Mol) -> List[Chem.Mol]:
        """Enumerate and rank tautomers by chemical likelihood."""
        enumerator = rdMolStandardize.TautomerEnumerator()
        enumerator.SetMaxTautomers(self.config.max_tautomers * 2)  # Generate more
        
        tautomers = []
        for taut in enumerator.Enumerate(mol):
            tautomers.append(taut)
        
        if not tautomers:
            return [mol]
        
        # Score tautomers
        scored = []
        for taut in tautomers:
            score = self._score_tautomer(taut)
            scored.append((score, taut))
    
        # Sort by score (higher is better) and take top N
        scored.sort(reverse=True, key=lambda x: x[0])
        return [taut for score, taut in scored[:self.config.max_tautomers]]

    def _score_tautomer(self, mol: Chem.Mol) -> float:
        """Score tautomer by chemical likelihood."""
        score = 0.0
        
        # Prefer aromatic systems
        score += sum(1 for atom in mol.GetAromaticAtoms()) * 2.0
        
        # Penalize charge separation
        formal_charges = [atom.GetFormalCharge() for atom in mol.GetAtoms()]
        score -= sum(abs(c) for c in formal_charges) * 3.0
        
        # Prefer carbonyl over enol
        carbonyl_pattern = Chem.MolFromSmarts('[C]=[O]')
        enol_pattern = Chem.MolFromSmarts('[C]=[C]-[OH]')
        
        score += len(mol.GetSubstructMatches(carbonyl_pattern)) * 1.0
        score -= len(mol.GetSubstructMatches(enol_pattern)) * 0.5
        
        # Prefer amide over iminol
        amide_pattern = Chem.MolFromSmarts('[C](=[O])-[N]')
        iminol_pattern = Chem.MolFromSmarts('[C](-[OH])=[N]')
        
        score += len(mol.GetSubstructMatches(amide_pattern)) * 2.0
        score -= len(mol.GetSubstructMatches(iminol_pattern)) * 1.5
        
        return score


    def _generate_3d_conformers(self, mol: Chem.Mol) -> List[Chem.Mol]:
        """
        Generate 3D conformers using ETKDGv3 with RMSD clustering.
        
        Args:
            mol: Input 2D molecule
            
        Returns:
            List of molecules, each with a single diverse conformer
        """
        mol = Chem.AddHs(mol)
        
        # Check if conformer generation needed (skip for very rigid/small molecules)
        n_rotatable = Descriptors.NumRotatableBonds(mol)
        if n_rotatable == 0 or mol.GetNumHeavyAtoms() < 5:
            # Very rigid - just generate one conformer
            params = AllChem.ETKDGv3()
            params.randomSeed = 42
            conf_id = AllChem.EmbedMolecule(mol, params)
            if conf_id == -1:
                raise LigandPreparationError("Failed to generate 3D coordinates")
            return [mol]
        
        # Generate multiple conformers for flexible molecules
        num_confs_target = self.config.max_conformers
        num_confs_generate = max(num_confs_target * 3, 10)  # Generate extras for clustering
        
        # For very large molecules, reduce conformer generation slightly to help performance
        # Parallel RMSD processing makes clustering much faster, so we can be less aggressive
        num_atoms = mol.GetNumHeavyAtoms()
        if num_atoms > 60:
            # Only for VERY large molecules (>60 atoms), reduce conformer count
            logger.info(f"Very large molecule ({num_atoms} atoms): reducing conformers for performance")
            num_confs_target = min(num_confs_target, 8)
            num_confs_generate = num_confs_target * 2
        elif num_atoms > 50:
            # For large molecules (50-60 atoms), slight reduction
            logger.info(f"Large molecule ({num_atoms} atoms): slight conformer reduction")
            num_confs_target = min(num_confs_target, 10)
            num_confs_generate = num_confs_target * 2
        
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        params.useRandomCoords = True
        params.numThreads = 0  # Use all available cores
        
        conf_ids = list(AllChem.EmbedMultipleConfs(
            mol, 
            numConfs=num_confs_generate, 
            params=params
        ))
        
        if not conf_ids:
            logger.warning("ETKDGv3 multi-conformer embedding failed, using single conformer")
            result = AllChem.EmbedMolecule(mol, params)
            if result == -1:
                result = AllChem.EmbedMolecule(mol)
                if result == -1:
                    raise LigandPreparationError("Failed to generate 3D coordinates")
            conf_ids = [0]
        
        # Minimize all conformers
        if self.config.mmff_minimize:
            props = AllChem.MMFFGetMoleculeProperties(mol)
            if props:
                energies = []
                logger.info(f"Minimizing {len(conf_ids)} conformers...")  # ADD THIS
                for i, conf_id in enumerate(conf_ids, 1):
                    if i % 5 == 0:  # ADD THIS - log every 5 conformers
                        logger.info(f"  Progress: {i}/{len(conf_ids)} conformers minimized")
                    ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=conf_id)
                    if ff:
                        ff.Minimize(maxIts=200)
                        energies.append(ff.CalcEnergy())
                    else:
                        energies.append(999.9)
            else:
                energies = [0.0] * len(conf_ids)
        else:
            energies = [0.0] * len(conf_ids)
        
        # Cluster conformers by RMSD if we have multiple
        logger.info(f"Conformer minimization complete. Total conformers: {len(conf_ids)}, target: {num_confs_target}")
        if len(conf_ids) > num_confs_target:
            logger.info(f"Starting RMSD clustering to select {num_confs_target} diverse conformers from {len(conf_ids)} total...")
            logger.info(f"This will calculate {len(conf_ids)*(len(conf_ids)-1)//2} pairwise RMSD values...")
            selected_conf_ids = self._cluster_conformers(
                mol,
                conf_ids,
                energies,
                rmsd_threshold=0.5,  # Angstroms
                max_keep=num_confs_target
            )
            logger.info(f"RMSD clustering complete. Selected {len(selected_conf_ids)} conformers.")
        else:
            logger.info(f"Skipping clustering (have {len(conf_ids)} conformers, target is {num_confs_target})")
            selected_conf_ids = conf_ids
        
        # Create separate molecule for each conformer
        logger.info(f"Creating separate molecule objects for {len(selected_conf_ids)} conformers...")
        conf_mols = []
        for conf_id in selected_conf_ids:
            conf_mol = Chem.Mol(mol)
            # Remove all conformers except the one we want
            for other_id in reversed([c.GetId() for c in conf_mol.GetConformers()]):
                if other_id != conf_id:
                    conf_mol.RemoveConformer(other_id)
            conf_mols.append(conf_mol)
        
        logger.info(f"✓ Successfully generated {len(conf_mols)} diverse conformers from {len(conf_ids)} total")
        return conf_mols
    
    def _cluster_conformers(
        self,
        mol: Chem.Mol,
        conf_ids: List[int],
        energies: List[float],
        rmsd_threshold: float = 0.5,
        max_keep: int = 10
    ) -> List[int]:
        """
        Cluster conformers by RMSD and select diverse representatives.
        
        Strategy: Greedy selection starting with lowest energy, then pick
        conformers maximally different from already selected ones.
        
        Args:
            mol: Molecule with multiple conformers
            conf_ids: List of conformer IDs
            energies: Energy for each conformer
            rmsd_threshold: Minimum RMSD between selected conformers
            max_keep: Maximum number of conformers to keep
            
        Returns:
            List of selected conformer IDs
        """
        if len(conf_ids) <= max_keep:
            return conf_ids
        
        n = len(conf_ids)
        total_rmsds = n * (n - 1) // 2
        
        logger.info(f"  Calculating {total_rmsds} pairwise RMSDs for {n} conformers...")
        logger.info(f"  Molecule has {mol.GetNumHeavyAtoms()} heavy atoms - this may take a while...")
        
        # Check if parallel processing should be forced off (for troubleshooting)
        force_serial = os.environ.get('LIGAND_PREP_FORCE_SERIAL', '0') == '1'
        
        # Determine if we should use parallel processing
        use_parallel = total_rmsds > 50 and self.n_cpus > 1 and not force_serial
        
        if force_serial:
            logger.warning(f"  Parallel processing disabled by LIGAND_PREP_FORCE_SERIAL environment variable")
            logger.info(f"  Using serial processing (forced)")
            rmsds = self._calculate_rmsd_matrix_serial(mol, conf_ids)
        elif use_parallel:
            logger.info(f"  Using parallel processing with {self.n_cpus} workers...")
            rmsds = self._calculate_rmsd_matrix_parallel(mol, conf_ids, self.n_cpus)
        else:
            if self.n_cpus == 1:
                logger.info(f"  Using serial processing (single CPU mode)")
            else:
                logger.info(f"  Using serial processing (small molecule)")
            rmsds = self._calculate_rmsd_matrix_serial(mol, conf_ids)
        
        logger.info(f"  RMSD matrix calculation complete!")
        logger.info(f"  RMSD matrix calculation complete!")
        
        # Start with lowest energy conformer
        logger.info(f"  Selecting {max_keep} diverse conformers using greedy algorithm...")
        selected_indices = [int(np.argmin(energies))]
        logger.debug(f"    Starting with conformer {selected_indices[0]} (lowest energy: {min(energies):.2f} kcal/mol)")
        
        # Greedily add most different conformers
        while len(selected_indices) < max_keep:
            # For each unselected conformer, find minimum distance to any selected
            max_min_dist = -1
            best_idx = -1
            
            for i in range(n):
                if i in selected_indices:
                    continue
                
                # Distance to nearest selected conformer
                min_dist = min(rmsds[i, j] for j in selected_indices)
                
                # Pick the one with maximum minimum distance
                if min_dist > max_min_dist:
                    max_min_dist = min_dist
                    best_idx = i
            
            if best_idx == -1 or max_min_dist < rmsd_threshold:
                break  # No more diverse conformers available
            
            selected_indices.append(best_idx)
        
        logger.info(f"  Greedy selection complete. Selected {len(selected_indices)} conformers.")
        
        # Return conformer IDs in energy order
        selected_with_energy = [(energies[i], conf_ids[i]) for i in selected_indices]
        selected_with_energy.sort()  # Sort by energy
        
        logger.debug(f"  Final conformer selection (by energy): {[conf_id for _, conf_id in selected_with_energy]}")
        
        return [conf_id for energy, conf_id in selected_with_energy]
    
    def _calculate_rmsd_matrix_serial(
        self,
        mol: Chem.Mol,
        conf_ids: List[int]
    ) -> np.ndarray:
        """
        Calculate RMSD matrix using serial processing (with progress reporting).
        
        Args:
            mol: Molecule with multiple conformers
            conf_ids: List of conformer IDs
            
        Returns:
            Square RMSD matrix
        """
        n = len(conf_ids)
        rmsds = np.zeros((n, n))
        total_rmsds = n * (n - 1) // 2
        completed = 0
        report_interval = max(total_rmsds // 20, 1)  # Report every 5%
        
        for i in range(n):
            for j in range(i+1, n):
                try:
                    rmsd = AllChem.GetBestRMS(mol, mol, conf_ids[i], conf_ids[j])
                    rmsds[i, j] = rmsd
                    rmsds[j, i] = rmsd
                except:
                    rmsds[i, j] = 999.9
                    rmsds[j, i] = 999.9
                
                completed += 1
                if completed % report_interval == 0:
                    percent = (completed * 100) // total_rmsds
                    logger.info(f"    RMSD calculation progress: {completed}/{total_rmsds} ({percent}%)")
        
        return rmsds
    
    def _calculate_rmsd_matrix_parallel(
        self,
        mol: Chem.Mol,
        conf_ids: List[int],
        n_workers: int
    ) -> np.ndarray:
        """
        Calculate RMSD matrix using parallel processing with real-time progress.
        
        Uses streaming results (starmap) to:
        - Monitor progress in real-time
        - Detect and handle stalled workers
        - Provide better timeout handling
        
        Args:
            mol: Molecule with multiple conformers
            conf_ids: List of conformer IDs
            n_workers: Number of worker processes
            
        Returns:
            Square RMSD matrix
        """
        n = len(conf_ids)
        rmsds = np.zeros((n, n))
        
        # Serialize molecule once for all workers
        mol_binary = mol.ToBinary()
        
        # Generate all pairs
        pairs = []
        for i in range(n):
            for j in range(i+1, n):
                pairs.append((i, j, conf_ids[i], conf_ids[j]))
        
        total_pairs = len(pairs)
        logger.info(f"    Distributing {total_pairs} RMSD calculations across {n_workers} workers...")
        
        # Smart batching: scale with problem size
        if total_pairs < 50:
            # Very small - serial is likely faster due to overhead
            logger.info(f"    Small calculation count ({total_pairs}), considering serial processing...")
            logger.info(f"    Using serial processing for better efficiency")
            return self._calculate_rmsd_matrix_serial(mol, conf_ids)
        elif total_pairs < 200:
            # Small - use fewer batches (2 per worker) to minimize overhead
            num_batches = n_workers * 2
            batch_size = max(1, total_pairs // num_batches)
        else:
            # Large - use many batches for progress monitoring
            # Aim for ~100 batches or 10 per worker, whichever is larger
            num_batches = max(min(100, total_pairs // 10), n_workers * 10)
            batch_size = max(1, total_pairs // num_batches)
        
        batches = []
        for i in range(0, len(pairs), batch_size):
            batch = pairs[i:i+batch_size]
            batches.append((mol_binary, batch))
        
        logger.info(f"    Created {len(batches)} batches (~{batch_size} calculations per batch)")
        logger.info(f"    Expected total time: {total_pairs * 0.01:.1f}s (estimated)")
        
        # Fixed timeout: 30 minutes
        total_timeout = 1800  # 30 minutes in seconds
        
        logger.info(f"    Total timeout: {total_timeout/60:.1f} minutes ({total_timeout}s)")
        
        try:
            logger.info(f"    Creating multiprocessing Pool with {n_workers} workers...")
            
            # Use multiprocessing pool with starmap for streaming results
            with Pool(processes=n_workers) as pool:
                logger.info(f"    Pool created. Starting parallel RMSD calculations...")
                
                # Track progress
                completed_batches = 0
                last_progress_time = time.time()
                start_time = time.time()
                
                # Use starmap for streaming results
                try:
                    results_iter = pool.starmap(_calculate_rmsd_pairs_batch, batches, chunksize=1)
                    
                    for batch_results in results_iter:
                        # Process this batch's results
                        for i, j, rmsd in batch_results:
                            rmsds[i, j] = rmsd
                            rmsds[j, i] = rmsd
                        
                        completed_batches += 1
                        current_time = time.time()
                        
                        # Log progress every 10% or every 30 seconds
                        if (completed_batches % max(1, len(batches) // 10) == 0 or 
                            current_time - last_progress_time > 30):
                            percent = (completed_batches * 100) // len(batches)
                            elapsed = current_time - start_time
                            rate = completed_batches / elapsed if elapsed > 0 else 0
                            eta = (len(batches) - completed_batches) / rate if rate > 0 else 0
                            logger.info(
                                f"    Progress: {completed_batches}/{len(batches)} batches ({percent}%) - "
                                f"Elapsed: {elapsed/60:.1f} min, ETA: {eta/60:.1f} min"
                            )
                            last_progress_time = current_time
                        
                        # Check if we're over total timeout
                        if current_time - start_time > total_timeout:
                            logger.error(f"    Total timeout exceeded ({total_timeout/60:.1f} min)!")
                            logger.error(f"    Completed {completed_batches}/{len(batches)} batches ({completed_batches*100//len(batches)}%)")
                            pool.terminate()
                            pool.join()
                            logger.warning("    Falling back to serial processing...")
                            return self._calculate_rmsd_matrix_serial(mol, conf_ids)
                    
                    elapsed_total = time.time() - start_time
                    logger.info(f"    ✓ Parallel RMSD calculation complete in {elapsed_total/60:.2f} minutes!")
                    
                except StopIteration:
                    logger.info(f"    All batches processed successfully")
                    
            
        except MPTimeoutError as e:
            logger.error(f"    Multiprocessing timeout: {e}")
            logger.warning("    Falling back to serial processing...")
            return self._calculate_rmsd_matrix_serial(mol, conf_ids)
        except OSError as e:
            logger.error(f"    Multiprocessing OSError: {e}")
            logger.warning("    This can happen with GUI applications or certain systems")
            logger.warning("    Falling back to serial processing...")
            return self._calculate_rmsd_matrix_serial(mol, conf_ids)
        except Exception as e:
            logger.error(f"    Parallel processing failed: {type(e).__name__}: {e}")
            logger.warning("    Falling back to serial processing...")
            import traceback
            logger.debug(f"    Traceback: {traceback.format_exc()}")
            return self._calculate_rmsd_matrix_serial(mol, conf_ids)
        
        return rmsds
    
    def _minimize_mmff94s(
        self,
        mol: Chem.Mol
    ) -> Tuple[Chem.Mol, float]:
        """
        Energy minimization using MMFF94s force field.

        Args:
            mol: Input molecule with 3D coordinates

        Returns:
            (minimized_molecule, final_energy)
        """
        # Get MMFF properties
        props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94s')

        if props is None:
            logger.warning("MMFF94s properties could not be computed, skipping minimization")
            return mol, None

        # Set up force field
        ff = AllChem.MMFFGetMoleculeForceField(mol, props)

        if ff is None:
            logger.warning("MMFF94s force field could not be set up, skipping minimization")
            return mol, None

        # Minimize
        ff.Initialize()
        result = ff.Minimize(maxIts=200)
        if result != 0:
            result = ff.Minimize(maxIts=1000)
            if result != 0:
                logger.info(f"MMFF94s minimization did not converge (code: {result}); using best found")

        final_energy = ff.CalcEnergy()

        logger.debug(f"MMFF94s minimization: energy = {final_energy:.2f} kcal/mol")

        return mol, final_energy

    def _to_pdbqt(self, mol: Chem.Mol, mol_id: str) -> str:
        """
        Convert molecule to PDBQT format using Open Babel.

        Args:
            mol: Input molecule with 3D coordinates
            mol_id: Molecule identifier

        Returns:
            PDBQT string
        """
        logger.debug(f"      _to_pdbqt() called for {mol_id}")
        logger.debug(f"      Steps planned: SDF creation → Open Babel conversion → Cleanup")
        
        import tempfile
        import subprocess
        import os
        import signal
        from pathlib import Path

        def _run_obabel(command: list, output_path: str, timeout_secs: int = 60) -> Tuple[int, str, str, str]:
            logger.debug(f"        → Executing: {' '.join(command)}")
            logger.debug(f"        → Timeout: {timeout_secs}s, will kill process group if exceeded")
            logger.debug(f"        → Starting subprocess.run() NOW...")
            
            import time
            start_time = time.time()
            
            try:
                # Use preexec_fn to create a new process group for proper timeout handling
                result = subprocess.run(
                    command,
                    capture_output=True,
                    timeout=timeout_secs,
                    text=True,
                    preexec_fn=os.setsid if hasattr(os, 'setsid') else None
                )
                elapsed = time.time() - start_time
                logger.debug(f"        ✓ obabel completed in {elapsed:.1f}s with returncode: {result.returncode}")
            except subprocess.TimeoutExpired as e:
                elapsed = time.time() - start_time
                logger.error(f"        ✗ obabel TIMEOUT after {elapsed:.1f}s (limit: {timeout_secs}s)")
                # Try to kill the process group
                try:
                    if hasattr(os, 'killpg'):
                        os.killpg(os.getpgid(e.args[1].pid), signal.SIGKILL)
                        logger.debug(f"        → Process group killed")
                except Exception as kill_err:
                    logger.debug(f"        → Could not kill process group: {kill_err}")
                raise
            
            try:
                pdbqt_data = Path(output_path).read_text()
            except Exception as e:
                logger.warning(f"Could not read PDBQT output: {e}")
                pdbqt_data = ""
            return result.returncode, result.stdout, result.stderr, pdbqt_data

        try:
            # Write molecule to temporary SDF file
            logger.debug(f"      Creating temporary SDF file for {mol_id}...")
            with tempfile.NamedTemporaryFile(mode='w', suffix='.sdf', delete=False) as tmp_in:
                writer = Chem.SDWriter(tmp_in.name)
                writer.write(mol)
                writer.close()
                temp_sdf = tmp_in.name
            logger.debug(f"      ✓ SDF file created: {temp_sdf}")

            # Create temporary PDBQT output file
            logger.debug(f"      Creating temporary PDBQT output file...")
            with tempfile.NamedTemporaryFile(mode='w', suffix='.pdbqt', delete=False) as tmp_out:
                temp_pdbqt = tmp_out.name
            logger.debug(f"      ✓ PDBQT output file ready: {temp_pdbqt}")

            # Convert using Open Babel
            # -h flag adds hydrogens
            logger.debug(f"      Setting up Open Babel command...")
            base_cmd = [
                'obabel',
                temp_sdf,
                '-O', temp_pdbqt,
                '-h'
            ]
            logger.debug(f"      Base command: {' '.join(base_cmd)}")

            charge_models = self._charge_models_for(mol)
            
            # Adjust timeout based on molecule size
            num_atoms = mol.GetNumHeavyAtoms()
            timeout = 60 if num_atoms > 40 else 30  # Longer timeout for complex molecules
            
            logger.info(f"      Charge model strategy for {mol_id}:")
            logger.info(f"        - Models to try: {charge_models}")
            logger.info(f"        - Molecule size: {num_atoms} heavy atoms")
            logger.info(f"        - Timeout per attempt: {timeout}s")
            logger.info(f"        - Will try each model until one succeeds")
            logger.info(f"      Beginning Open Babel conversion attempts...")
            
            last_stderr = ""
            pdbqt_string = ""
            rc = 1

            for idx, charge_model in enumerate(charge_models, 1):
                charge_label = charge_model if charge_model else "none"
                logger.info(f"      → Attempt {idx}/{len(charge_models)}: charge model '{charge_label}'")
                cmd = list(base_cmd)
                if charge_model:
                    cmd += ['--partialcharge', charge_model]
                rc, stdout, stderr, pdbqt_string = _run_obabel(cmd, temp_pdbqt, timeout)
                last_stderr = stderr

                if rc == 0 and pdbqt_string.strip() and "0 molecules converted" not in stderr:
                    # Success!
                    if charge_model is None:
                        logger.warning(
                            f"        ⚠ Conversion succeeded without charge model (charges may be missing)"
                        )
                    elif charge_model != "gasteiger":
                        logger.info(
                            f"        ✓ SUCCESS with charge model '{charge_model}'"
                        )
                    else:
                        logger.info(f"        ✓ SUCCESS with standard Gasteiger charges")
                    logger.debug(f"        Output: {len(pdbqt_string)} chars, returncode={rc}")
                    break
                else:
                    # Failed
                    charge_label = charge_model or "none"
                    failure_reason = "empty output" if not pdbqt_string.strip() else \
                                   "zero molecules" if "0 molecules converted" in stderr else \
                                   f"returncode {rc}"
                    logger.warning(
                        f"        ✗ FAILED with charge model '{charge_label}': {failure_reason}"
                    )
                    if stderr.strip():
                        logger.debug(f"        stderr: {stderr.strip()[:200]}")  # First 200 chars

            # After the loop - check if any attempt succeeded
            logger.info(f"      Open Babel conversion result:")
            if rc != 0 or not pdbqt_string.strip():
                logger.error(f"        ✗ ALL {len(charge_models)} charge models FAILED")
                logger.error(f"        Last error: {last_stderr[:200]}")
                error_msg = f"Open Babel conversion failed: {last_stderr}"
                raise LigandPreparationError(error_msg)
            else:
                logger.info(f"        ✓ Conversion SUCCESSFUL ({len(pdbqt_string)} chars)")

            # Cleanup temporary files
            logger.debug(f"      Cleaning up temporary files...")
            Path(temp_sdf).unlink()
            Path(temp_pdbqt).unlink()
            logger.debug(f"      ✓ Temporary files deleted")

            if not pdbqt_string.strip():
                raise LigandPreparationError("Empty PDBQT output")

            return pdbqt_string

        except subprocess.TimeoutExpired:
            logger.error("PDBQT conversion timed out")
            if Path(temp_sdf).exists():
                Path(temp_sdf).unlink()
            if Path(temp_pdbqt).exists():
                Path(temp_pdbqt).unlink()
            raise LigandPreparationError("PDBQT conversion timed out")
        except Exception as e:
            logger.error(f"PDBQT conversion failed: {e}")
            if 'temp_sdf' in locals() and Path(temp_sdf).exists():
                Path(temp_sdf).unlink()
            if 'temp_pdbqt' in locals() and Path(temp_pdbqt).exists():
                Path(temp_pdbqt).unlink()
            raise LigandPreparationError(f"PDBQT conversion failed: {e}")

    def prepare_batch(
        self,
        smiles_list: List[Tuple[str, str]],
        max_variants_per_ligand: int = 8
    ) -> List[PreparedLigand]:
        """
        Prepare multiple ligands in batch.

        Args:
            smiles_list: List of (smiles, mol_id) tuples
            max_variants_per_ligand: Maximum variants to keep per ligand

        Returns:
            List of all prepared ligands
        """
        all_results = []

        for smiles, mol_id in smiles_list:
            try:
                results = self.prepare_from_smiles(smiles, mol_id)
                # Limit variants
                results = results[:max_variants_per_ligand]
                all_results.extend(results)

            except LigandPreparationError as e:
                logger.error(f"Skipping {mol_id}: {e}")
                continue

        logger.info(f"Batch preparation: {len(all_results)} ligands from {len(smiles_list)} inputs")
        return all_results

    def _charge_models_for(self, mol: Chem.Mol) -> List[Optional[str]]:
        """Select Open Babel charge models, with metal-aware fallbacks."""
        if self._contains_metal(mol):
            return ["qeq", "qtpie", "eem", None]
        return ["gasteiger", "mmff94", None]

    def _contains_metal(self, mol: Chem.Mol) -> bool:
        """Return True if molecule contains common metal ions."""
        metal_symbols = {
            "LI", "NA", "K", "RB", "CS",
            "MG", "CA", "SR", "BA",
            "MN", "FE", "CO", "NI", "CU", "ZN",
            "MO", "W", "AG", "AU", "CD", "HG",
            "PT", "PD", "AL", "GA", "IN", "SN", "PB", "BI"
        }
        for atom in mol.GetAtoms():
            if atom.GetSymbol().upper() in metal_symbols:
                return True
        return False

    def validate_ligand(self, mol: Chem.Mol) -> Tuple[bool, List[str]]:
        """
        Validate ligand meets basic criteria.

        Args:
            mol: Molecule to validate

        Returns:
            (is_valid, list_of_issues)
        """
        issues = []

        # Check has 3D coordinates
        if mol.GetNumConformers() == 0:
            issues.append("No 3D conformer")

        # Check reasonable size
        num_atoms = mol.GetNumHeavyAtoms()
        if num_atoms < 3:
            issues.append(f"Too few atoms: {num_atoms}")
        elif num_atoms > 100:
            issues.append(f"Too many atoms: {num_atoms}")

        # Check molecular weight
        mw = Descriptors.MolWt(mol)
        if mw > 900:
            issues.append(f"Molecular weight too high: {mw:.1f}")

        # Check has reasonable geometry
        try:
            conf = mol.GetConformer()
            positions = conf.GetPositions()
            if positions.shape[0] == 0:
                issues.append("Empty conformer")
        except:
            issues.append("Invalid conformer")

        return len(issues) == 0, issues
