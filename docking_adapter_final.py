"""
Adapter to integrate ensemble validation with existing docking platform infrastructure.

Uses your existing modules:
- docking_platform_gui.preparation.ligand (LigandPreparation, PreparedLigand)
- docking_platform_gui.preparation.protein (ProteinPreparation)
- docking_platform_gui.docking.smina (SminaDockingEngine)
- docking_platform_gui.docking.base (DockingResult, DockingPose, DockingEngine)
- docking_platform_gui.config.schema (LigandPreparationConfig)
"""

import tempfile
from pathlib import Path
from loguru import logger
import numpy as np
from rdkit import Chem
from typing import Optional, List, Dict

from docking_platform_gui.config.schema import LigandPreparationConfig
from docking_platform_gui.preparation.ligand import LigandPreparation, PreparedLigand
from docking_platform_gui.docking.smina import SminaDockingEngine
from docking_platform_gui.docking.base import DockingResult


class DockingManagerAdapter:
    """
    Adapter between ensemble validation experiment and your existing docking platform.
    
    This adapter:
    1. Uses your LigandPreparation for ligand prep (SMILES → PDBQT)
    2. Uses your SminaDockingEngine for docking
    3. Handles ensemble docking (dock to multiple structures, take best score)
    4. Returns results in format expected by ensemble validation experiment
    
    Features:
    - Automatic ligand preparation from SMILES
    - Proper pH 7.4 protonation states
    - 3D conformer generation with ETKDGv3
    - MMFF94s energy minimization
    - PDBQT conversion via Meeko/OpenBabel
    """
    
    def __init__(
        self,
        use_vina: bool = False,
        vina_binary: Optional[str] = None,
        smina_binary: Optional[str] = None,
        exhaustiveness: int = 16,
        num_modes: int = 20,
        cpu: int = 4,
        ligand_config: Optional[LigandPreparationConfig] = None
    ):
        """
        Initialize docking manager.
        
        Parameters
        ----------
        use_vina : bool
            Use Vina (True) or Smina (False) - currently only Smina supported
        vina_binary : str, optional
            Path to vina executable (not used currently)
        smina_binary : str, optional
            Path to smina executable (uses system PATH if None)
        exhaustiveness : int
            Docking exhaustiveness (16 recommended for ensemble)
        num_modes : int
            Number of binding modes to generate
        cpu : int
            Number of CPU cores for parallel processing
        ligand_config : LigandPreparationConfig, optional
            Custom ligand prep config (uses default if None)
        """
        if use_vina:
            logger.warning("Vina not currently supported - using Smina instead")
        
        # Initialize ligand preparation
        if ligand_config is None:
            ligand_config = LigandPreparationConfig(
                ph_range=(7.0, 7.4),  # Physiological pH
                max_tautomers=1,      # Keep simple for benchmarking speed
                max_conformers=1,     # Single conformer for benchmarking
                mmff_minimize=True    # Energy minimization
            )
        
        self.ligand_config = ligand_config
        self.ligand_prep = LigandPreparation(ligand_config, n_cpus=cpu)
        
        logger.info(f"Ligand preparation initialized with pH {ligand_config.ph_range}")
        
        # Initialize docking engine
        self.engine = SminaDockingEngine(
            binary_path=smina_binary,
            exhaustiveness=exhaustiveness,
            num_modes=num_modes,
            cpu=cpu
        )
        
        logger.info(f"Smina docking engine initialized (exhaustiveness={exhaustiveness}, "
                   f"modes={num_modes}, cpu={cpu})")
        
        # Temp directory for intermediate files
        self.temp_dir = Path(tempfile.mkdtemp(prefix='ensemble_docking_'))
        logger.info(f"Temp directory: {self.temp_dir}")
    
    def dock_compound(
        self,
        compound,
        protein: str,
        grid_params: Dict
    ) -> Dict:
        """
        Dock a single compound to a single protein.
        
        Parameters
        ----------
        compound : rdkit.Chem.Mol or str
            Molecule to dock (RDKit mol object or SMILES string)
        protein : str or Path
            Path to prepared protein PDBQT file
        grid_params : dict
            Grid box parameters with keys:
            - center: [x, y, z]
            - size: [x, y, z]
        
        Returns
        -------
        dict : Docking results with keys:
            - best_score: float (kcal/mol)
            - best_pose_path: str (path to output PDBQT)
            - all_scores: list of float
            - docking_result: DockingResult object
        """
        # 1. Prepare ligand (SMILES → PDBQT)
        ligand_pdbqt_path = self._prepare_ligand(compound)
        
        # 2. Validate protein is PDBQT
        protein_path = Path(protein)
        if not protein_path.suffix == '.pdbqt':
            logger.warning(f"Protein file is not .pdbqt: {protein_path}")
            logger.warning("Assuming it's already in PDBQT format")
        
        # 3. Extract grid parameters
        center = np.array(grid_params['center'])
        size = np.array(grid_params.get('size', [20, 20, 20]))
        
        # 4. Run docking
        output_pdbqt = self.temp_dir / f'docked_{ligand_pdbqt_path.stem}.pdbqt'
        
        try:
            # Call your SminaDockingEngine.dock()
            result: DockingResult = self.engine.dock(
                receptor_path=str(protein_path),
                ligand_path=str(ligand_pdbqt_path),
                center=center,
                size=size,
                output_path=str(output_pdbqt)
            )
            
            if not result.success:
                logger.error(f"Docking failed: {result.error_message}")
                raise Exception(result.error_message)
            
            # Extract all scores from poses
            all_scores = [pose.score for pose in result.poses]
            best_score = result.best_score
            
            logger.debug(f"Docking complete: {len(result.poses)} poses, "
                        f"best = {best_score:.2f} kcal/mol")
            
            return {
                'best_score': best_score,
                'best_pose_path': str(output_pdbqt),
                'all_scores': all_scores,
                'docking_result': result
            }
            
        except Exception as e:
            logger.error(f"Docking failed: {e}")
            # Return worst possible score
            return {
                'best_score': 0.0,
                'best_pose_path': None,
                'all_scores': [0.0],
                'docking_result': None,
                'error': str(e)
            }
    
    def dock_ensemble(
        self,
        compound,
        proteins: List[str],
        grid_params: Dict
    ) -> Dict:
        """
        Dock a single compound to ensemble of proteins.
        
        Uses "best score" strategy: dock to all structures, return minimum score.
        
        Parameters
        ----------
        compound : rdkit.Chem.Mol or str
            Molecule to dock
        proteins : list of str/Path
            Paths to prepared protein PDBQT files
        grid_params : dict
            Grid box parameters (same for all structures)
        
        Returns
        -------
        dict : Ensemble docking results with keys:
            - best_score: float (best across ensemble)
            - best_structure_id: int (which structure, 1-indexed)
            - individual_scores: list of float
            - individual_results: list of dict
        """
        individual_results = []
        individual_scores = []
        
        logger.debug(f"Ensemble docking to {len(proteins)} structures...")
        
        for i, protein in enumerate(proteins, 1):
            logger.debug(f"  Structure {i}/{len(proteins)}: {Path(protein).name}")
            
            # Dock to this ensemble member
            result = self.dock_compound(
                compound=compound,
                protein=protein,
                grid_params=grid_params
            )
            
            individual_results.append(result)
            individual_scores.append(result['best_score'])
        
        # Find best score (most negative)
        best_idx = np.argmin(individual_scores)
        best_score = individual_scores[best_idx]
        
        logger.debug(f"Ensemble complete: best is structure {best_idx+1} "
                    f"with {best_score:.2f} kcal/mol")
        
        return {
            'best_score': best_score,
            'best_structure_id': best_idx + 1,  # 1-indexed
            'individual_scores': individual_scores,
            'individual_results': individual_results
        }
    
    def _prepare_ligand(self, compound) -> Path:
        """
        Prepare ligand for docking using your LigandPreparation pipeline.
        
        Pipeline: SMILES → Dimorphite-DL → Tautomers → ETKDGv3 → MMFF94s → PDBQT
        
        Parameters
        ----------
        compound : rdkit.Chem.Mol or str
            Molecule or SMILES
        
        Returns
        -------
        Path : Path to prepared ligand PDBQT file
        """
        # Convert to SMILES if needed
        if isinstance(compound, str):
            smiles = compound
            mol = Chem.MolFromSmiles(smiles)
        elif isinstance(compound, Chem.Mol):
            mol = compound
            smiles = Chem.MolToSmiles(mol)
        else:
            raise ValueError(f"Invalid compound type: {type(compound)}")
        
        if mol is None:
            raise ValueError(f"Invalid compound: {compound}")
        
        # Get compound ID
        if mol.HasProp('_Name'):
            mol_id = mol.GetProp('_Name')
        else:
            mol_id = 'ligand'
        
        # Prepare ligand using your LigandPreparation
        logger.debug(f"Preparing ligand {mol_id} from SMILES: {smiles[:50]}...")
        
        prepared_ligands: List[PreparedLigand] = self.ligand_prep.prepare_from_smiles(
            smiles=smiles,
            mol_id=mol_id,
            enumerate_states=False  # Single state for benchmarking speed
        )
        
        if not prepared_ligands:
            raise ValueError(f"Ligand preparation failed for {mol_id}")
        
        # Use first prepared variant
        prepared_ligand: PreparedLigand = prepared_ligands[0]
        
        if not prepared_ligand.pdbqt:
            raise ValueError(f"PDBQT generation failed for {mol_id}")
        
        # Save PDBQT to temp file
        ligand_pdbqt = self.temp_dir / f'{mol_id}_prepared.pdbqt'
        prepared_ligand.to_pdbqt_file(str(ligand_pdbqt))
        
        logger.debug(f"  ✓ Ligand prepared: {ligand_pdbqt}")
        
        return ligand_pdbqt
    
    def cleanup(self):
        """Clean up temporary files."""
        import shutil
        if self.temp_dir.exists():
            try:
                shutil.rmtree(self.temp_dir)
                logger.info(f"Cleaned up temp directory: {self.temp_dir}")
            except Exception as e:
                logger.warning(f"Failed to clean up temp directory: {e}")


class CachedDockingAdapter(DockingManagerAdapter):
    """
    Enhanced adapter with ligand caching.
    
    Caches prepared ligands to avoid expensive re-computation when
    screening the same compounds multiple times.
    """
    
    def __init__(
        self,
        cache_dir: str = "ligand_cache",
        **kwargs
    ):
        """
        Initialize with ligand caching.
        
        Parameters
        ----------
        cache_dir : str
            Directory for ligand cache
        **kwargs : dict
            All other arguments passed to DockingManagerAdapter
        """
        super().__init__(**kwargs)
        
        from docking_platform_gui.preparation.ligand_cache import LigandCache
        
        self.cache = LigandCache(cache_dir=cache_dir)
        logger.info(f"Ligand caching enabled: {cache_dir}")
    
    def _prepare_ligand(self, compound) -> Path:
        """
        Prepare ligand with caching.
        
        First checks cache. If found, loads instantly.
        If not found, prepares and saves to cache.
        """
        # Convert to SMILES
        if isinstance(compound, str):
            smiles = compound
            mol = Chem.MolFromSmiles(smiles)
        else:
            mol = compound
            smiles = Chem.MolToSmiles(mol)
        
        if mol is None:
            raise ValueError(f"Invalid compound")
        
        # Get mol ID
        mol_id = mol.GetProp('_Name') if mol.HasProp('_Name') else 'ligand'
        
        # Try to get from cache (or prepare if not cached)
        prepared_ligands = self.cache.get(
            smiles=smiles,
            mol_id=mol_id,
            config=self.ligand_config,
            n_cpus=self.ligand_prep.n_cpus,
            enumerate_states=False
        )
        
        if not prepared_ligands:
            raise ValueError(f"Ligand preparation failed for {mol_id}")
        
        prepared_ligand = prepared_ligands[0]
        
        # Save to temp file
        ligand_pdbqt = self.temp_dir / f'{mol_id}_prepared.pdbqt'
        prepared_ligand.to_pdbqt_file(str(ligand_pdbqt))
        
        return ligand_pdbqt
