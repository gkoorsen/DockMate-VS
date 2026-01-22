"""Caching system for prepared ligands to avoid expensive re-computation."""

import hashlib
import pickle
from pathlib import Path
from typing import List, Optional
from loguru import logger

from docking_platform_gui.preparation.ligand import (
    PreparedLigand,
    LigandPreparation,
    LigandPreparationConfig
)


class LigandCache:
    """
    Cache for prepared ligands to avoid re-computation.
    
    Caches prepared ligands to disk using SMILES + configuration hash as key.
    Automatically invalidates cache when preparation parameters change.
    
    Example:
        >>> cache = LigandCache(cache_dir="ligand_cache")
        >>> config = LigandPreparationConfig(ph_range=(7.0, 7.4))
        >>> 
        >>> # First call: prepares and caches (slow)
        >>> ligands = cache.get("CCO", "ethanol", config, n_cpus=4)
        >>> 
        >>> # Second call: loads from cache (instant)
        >>> ligands = cache.get("CCO", "ethanol", config, n_cpus=4)
    """
    
    def __init__(self, cache_dir: str = "ligand_cache"):
        """
        Initialize ligand cache.
        
        Args:
            cache_dir: Directory to store cached ligands
        """
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(exist_ok=True)
        logger.info(f"Ligand cache initialized: {self.cache_dir.absolute()}")
    
    def _get_cache_key(self, smiles: str, config: LigandPreparationConfig) -> str:
        """
        Generate unique cache key from SMILES and configuration.
        
        Args:
            smiles: SMILES string
            config: Preparation configuration
            
        Returns:
            12-character hash string
        """
        # Include all relevant config parameters that affect output
        cache_input = (
            f"{smiles}|"
            f"{config.ph_range}|"
            f"{config.max_tautomers}|"
            f"{config.max_conformers}|"
            f"{config.mmff_minimize}"
        )
        return hashlib.md5(cache_input.encode()).hexdigest()[:12]
    
    def get(
        self,
        smiles: str,
        mol_id: str,
        config: LigandPreparationConfig,
        n_cpus: int = 4,
        enumerate_states: bool = True
    ) -> List[PreparedLigand]:
        """
        Get prepared ligands, using cache if available.
        
        Checks cache first. If found, loads from disk (fast).
        If not found, prepares ligands and saves to cache for next time.
        
        Args:
            smiles: SMILES string
            mol_id: Molecule identifier
            config: Preparation configuration
            n_cpus: Number of CPUs for parallel processing
            enumerate_states: Whether to enumerate protonation states
        
        Returns:
            List of prepared ligands
        """
        # Generate cache key
        cache_key = self._get_cache_key(smiles, config)
        cache_file = self.cache_dir / f"{mol_id}_{cache_key}.pkl"
        
        # Try to load from cache
        if cache_file.exists():
            logger.info(f"✓ Loading cached ligands for {mol_id}")
            try:
                with open(cache_file, 'rb') as f:
                    ligands = pickle.load(f)
                logger.info(f"  Loaded {len(ligands)} variants from cache (instant)")
                return ligands
            except Exception as e:
                logger.warning(f"Cache load failed for {mol_id}: {e}, will re-prepare")
        
        # Not cached or load failed - prepare it
        logger.info(f"⚙ Preparing ligands for {mol_id} (not in cache)...")
        prep = LigandPreparation(config, n_cpus=n_cpus)
        ligands = prep.prepare_from_smiles(smiles, mol_id, enumerate_states)
        
        # Save to cache
        try:
            with open(cache_file, 'wb') as f:
                pickle.dump(ligands, f)
            logger.info(f"✓ Cached {len(ligands)} variants for {mol_id}")
        except Exception as e:
            logger.warning(f"Failed to cache {mol_id}: {e}")
        
        return ligands
    
    def clear(self, mol_id: Optional[str] = None):
        """
        Clear cache for specific molecule or all molecules.
        
        Args:
            mol_id: Molecule ID to clear, or None to clear all
        """
        if mol_id:
            # Clear specific molecule (all cached variants)
            cleared = 0
            for cache_file in self.cache_dir.glob(f"{mol_id}_*.pkl"):
                cache_file.unlink()
                cleared += 1
            if cleared:
                logger.info(f"Cleared {cleared} cached variant(s) for {mol_id}")
            else:
                logger.info(f"No cached variants found for {mol_id}")
        else:
            # Clear all cached ligands
            cache_files = list(self.cache_dir.glob("*.pkl"))
            for cache_file in cache_files:
                cache_file.unlink()
            logger.info(f"Cleared all {len(cache_files)} cached ligands")
    
    def stats(self):
        """Print cache statistics."""
        cache_files = list(self.cache_dir.glob("*.pkl"))
        
        if not cache_files:
            logger.info("Cache is empty")
            return
        
        total_size = sum(f.stat().st_size for f in cache_files)
        
        # Count unique molecules
        unique_mols = set()
        for f in cache_files:
            mol_id = f.stem.rsplit('_', 1)[0]  # Remove hash suffix
            unique_mols.add(mol_id)
        
        logger.info(f"Cache statistics:")
        logger.info(f"  Unique molecules: {len(unique_mols)}")
        logger.info(f"  Total variants: {len(cache_files)}")
        logger.info(f"  Total size: {total_size / 1024 / 1024:.2f} MB")
        logger.info(f"  Average size: {total_size / len(cache_files) / 1024:.1f} KB per variant")
    
    def list_cached(self) -> List[str]:
        """
        List all cached molecule IDs.
        
        Returns:
            List of molecule IDs that have cached variants
        """
        cache_files = list(self.cache_dir.glob("*.pkl"))
        
        # Extract unique molecule IDs
        unique_mols = set()
        for f in cache_files:
            mol_id = f.stem.rsplit('_', 1)[0]  # Remove hash suffix
            unique_mols.add(mol_id)
        
        return sorted(unique_mols)
