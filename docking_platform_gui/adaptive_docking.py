"""
Adaptive Docking Pipeline

Implements intelligent cascading docking protocol based on investigation findings.
Starts with "best guess" configuration, then tries alternatives if initial result
doesn't meet RMSD threshold.

Key Learnings Applied:
1. Box size is THE critical parameter (4Å margin best)
2. rDock outperforms Smina (86% improvement)
3. Waters hurt pose prediction (remove all)
4. Exhaustiveness plateaus at 16
5. Scoring function has minimal impact

Strategy:
- Protocol 1: Smina + 4Å margin (fast, usually works)
- Protocol 2: rDock + 4Å margin (if Smina fails, try better algorithm)
- Protocol 3: Explore box sizes (if both fail, maybe wrong margin)
- Protocol 4: Last resort options (flexible, higher exhaustiveness)
"""

import logging
import time
import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple

from docking_platform_gui.preparation.protein import (
    ProteinPreparation,
    ProteinPreparationConfig,
    WaterHandling
)
from docking_platform_gui.preparation.ligand import (
    LigandPreparation,
    LigandPreparationConfig
)
from docking_platform_gui.binding_site.cocrystal import BindingSiteDefinition
from docking_platform_gui.docking.smina import SminaDockingEngine
from docking_platform_gui.utils.rmsd import calculate_rmsd


logger = logging.getLogger(__name__)


@dataclass
class DockingProtocol:
    """Configuration for a specific docking protocol."""
    name: str
    description: str
    engine: str  # 'smina' or 'rdock'
    margin: float
    scoring: str
    exhaustiveness: int
    flexible_residues: Optional[List[str]] = None
    expected_runtime_min: float = 2.0  # minutes


@dataclass
class DockingResult:
    """Result from a docking attempt."""
    protocol_name: str
    rmsd: float
    score: float
    num_poses: int
    runtime_sec: float
    output_file: Path
    success: bool  # < 2.0Å threshold
    error_message: Optional[str] = None
    variant_label: Optional[str] = None


class AdaptiveDockingPipeline:
    """
    Adaptive docking pipeline that tries multiple protocols in order.

    Cascading strategy:
    1. Best guess (Smina + 4Å margin) - 90% success expected
    2. rDock algorithm (if Smina fails) - 86% better
    3. Explore margins (if both fail) - maybe wrong box size
    4. Last resort (flexible, high exhaustiveness) - edge cases
    """

    def __init__(
        self,
        output_dir: Path,
        rmsd_threshold: float = 2.0,
        max_protocols: int = 5,
        rdock_available: bool = False,
        rdock_root: Optional[Path] = None,
        smina_binary: str = "smina",
        vina_binary: str = "vina",
        smina_timeout_sec: int = 1200,
        smina_exhaustiveness: int = 16,
        use_vina: bool = False,
        ligand_variant_mode: str = "first",
        variant_select_by: str = "rmsd",
        max_tautomers: int = 8,
        max_conformers: int = 10
    ):
        """
        Initialize adaptive pipeline.

        Args:
            output_dir: Directory for outputs
            rmsd_threshold: Success threshold (default 2.0Å)
            max_protocols: Maximum protocols to try before giving up
            rdock_available: Whether rDock is installed
            rdock_root: Path to rDock installation
            smina_binary: Path to smina executable (default: "smina" in PATH)
            vina_binary: Path to vina executable (default: "vina" in PATH)
            smina_timeout_sec: Smina timeout in seconds
            smina_exhaustiveness: Base Smina exhaustiveness for protocols
            use_vina: Use Vina binary instead of Smina (vina scoring only)
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.rmsd_threshold = rmsd_threshold
        self.max_protocols = max_protocols
        self.rdock_available = rdock_available
        self.rdock_root = rdock_root
        self.smina_binary = smina_binary
        self.vina_binary = vina_binary
        self.smina_timeout_sec = smina_timeout_sec
        self.smina_exhaustiveness = smina_exhaustiveness
        self.use_vina = use_vina
        self.ligand_variant_mode = ligand_variant_mode
        self.variant_select_by = variant_select_by
        self.max_tautomers = max_tautomers
        self.max_conformers = max_conformers
        self.docking_binary = self.vina_binary if self.use_vina else self.smina_binary
        self.allow_custom_scoring = not self.use_vina

        # Initialize protocol cascade
        self.protocols = self._build_protocol_cascade()

        logger.info(f"Adaptive pipeline initialized with {len(self.protocols)} protocols")
        logger.info(f"RMSD threshold: {rmsd_threshold}Å")
        logger.info(f"rDock available: {rdock_available}")
        logger.info(f"Docking binary: {self.docking_binary}")

    def _build_protocol_cascade(self) -> List[DockingProtocol]:
        """
        Build ordered list of protocols to try.

        Based on investigation findings:
        - Protocol 1: Smina + 4Å (best guess, fastest)
        - Protocol 2: rDock + 4Å (better algorithm)
        - Protocol 3: Smina + 6Å (maybe need more space)
        - Protocol 4: Smina + 8Å (traditional default)
        - Protocol 5: Smina + 4Å + higher exhaustiveness (last resort)
        """
        protocols = []
        base_exhaustiveness = self.smina_exhaustiveness
        high_exhaustiveness = max(32, base_exhaustiveness)
        base_scoring = "vina"
        alt_scoring = "vinardo" if self.allow_custom_scoring else "vina"

        # PROTOCOL 1: Best Guess (Smina + optimal margin)
        protocols.append(DockingProtocol(
            name="Protocol1_BestGuess",
            description="Smina + 4Å margin + vina scoring (optimal configuration)",
            engine="smina",
            margin=4.0,
            scoring=base_scoring,
            exhaustiveness=base_exhaustiveness,
            expected_runtime_min=1.5
        ))

        # PROTOCOL 2: rDock (if available)
        if self.rdock_available:
            protocols.append(DockingProtocol(
                name="Protocol2_rDock",
                description="rDock + 4Å margin (superior algorithm, 86% better)",
                engine="rdock",
                margin=4.0,
                scoring=base_scoring,  # Not used by rDock, but keep for consistency
                exhaustiveness=base_exhaustiveness,
                expected_runtime_min=2.5
            ))

        # PROTOCOL 3: Larger margin (maybe constrained too much)
        protocols.append(DockingProtocol(
            name="Protocol3_LargerMargin",
            description="Smina + 6Å margin (less constrained search)",
            engine="smina",
            margin=6.0,
            scoring=base_scoring,
            exhaustiveness=base_exhaustiveness,
            expected_runtime_min=1.8
        ))

        # PROTOCOL 4: Traditional default (8Å)
        protocols.append(DockingProtocol(
            name="Protocol4_TraditionalDefault",
            description="Smina + 8Å margin + vinardo scoring",
            engine="smina",
            margin=8.0,
            scoring=alt_scoring,  # Vinardo better with larger boxes
            exhaustiveness=base_exhaustiveness,
            expected_runtime_min=2.0
        ))

        # PROTOCOL 5: High exhaustiveness (last resort)
        protocols.append(DockingProtocol(
            name="Protocol5_HighExhaustiveness",
            description="Smina + 4Å margin + exhaustiveness=32 (thorough search)",
            engine="smina",
            margin=4.0,
            scoring=base_scoring,
            exhaustiveness=high_exhaustiveness,
            expected_runtime_min=3.0
        ))

        return protocols[:self.max_protocols]

    def run_adaptive_docking(
        self,
        pdb_file: Path,
        ligand_smiles: str,
        ligand_name: str,
        ligand_resname: str,
        ligand_chain: str,
        crystal_ligand_pdb: Optional[Path] = None
    ) -> Tuple[DockingResult, List[DockingResult]]:
        """
        Run adaptive docking cascade until success or exhaustion.

        Args:
            pdb_file: Receptor PDB file
            ligand_smiles: Ligand SMILES string
            ligand_name: Ligand identifier
            ligand_resname: Residue name in PDB (for binding site)
            ligand_chain: Chain ID in PDB (for binding site)
            crystal_ligand_pdb: Optional crystal structure for RMSD (if None, extract from pdb_file)

        Returns:
            (best_result, all_results) - Best result and list of all attempts
        """
        logger.info("=" * 60)
        logger.info("ADAPTIVE DOCKING PIPELINE")
        logger.info("=" * 60)
        logger.info(f"Receptor: {pdb_file}")
        logger.info(f"Ligand: {ligand_name} ({ligand_smiles})")
        logger.info(f"Threshold: {self.rmsd_threshold}Å")
        logger.info("")

        logger.info("Preparing receptor once for all protocols...")
        prepared_receptor_pdbqt, prepared_receptor_pdb = self._prepare_receptor(
            pdb_file=pdb_file
        )

        # Extract crystal ligand if not provided
        if crystal_ligand_pdb is None:
            crystal_ligand_pdb = self._extract_crystal_ligand(
                pdb_file,
                ligand_resname,
                ligand_chain
            )

        enumerate_states = not self._contains_metal(ligand_smiles)
        if not enumerate_states:
            logger.info("Skipping protonation/tautomer enumeration for metal-containing ligand")

        variants = self._prepare_ligand_variants(
            ligand_smiles=ligand_smiles,
            ligand_name=ligand_name,
            enumerate_states=enumerate_states
        )

        mode = self.ligand_variant_mode
        if mode == "best":
            variants = [self._select_best_variant(variants)]
        elif mode == "first":
            variants = variants[:1]

        all_results: List[DockingResult] = []
        best_result: Optional[DockingResult] = None

        for variant in variants:
            variant_label = variant["label"]
            ligand_pdbqt = variant["pdbqt"]
            variant_smiles = variant.get("smiles") or ligand_smiles
            variant_root = self.output_dir / "variants" / variant_label
            variant_root.mkdir(parents=True, exist_ok=True)
            if crystal_ligand_pdb and crystal_ligand_pdb.exists():
                dest = variant_root / "crystal_ligand.pdb"
                if not dest.exists():
                    try:
                        shutil.copy2(crystal_ligand_pdb, dest)
                    except Exception:
                        pass

            variant_results: List[DockingResult] = []
            variant_best: Optional[DockingResult] = None

            for i, protocol in enumerate(self.protocols, 1):
                logger.info(f"Attempting Protocol {i}/{len(self.protocols)}: {protocol.name}")
                logger.info(f"  Description: {protocol.description}")
                logger.info(f"  Expected runtime: ~{protocol.expected_runtime_min:.1f} min")
                logger.info("")

                try:
                    result = self._run_single_protocol(
                        protocol=protocol,
                        pdb_file=pdb_file,
                        ligand_smiles=variant_smiles,
                        ligand_name=ligand_name,
                        ligand_resname=ligand_resname,
                        ligand_chain=ligand_chain,
                        crystal_ligand_pdb=crystal_ligand_pdb,
                        prepared_ligand_pdbqt=ligand_pdbqt,
                        prepared_receptor_pdbqt=prepared_receptor_pdbqt,
                        prepared_receptor_pdb=prepared_receptor_pdb,
                        protocol_root=variant_root
                    )
                    result.variant_label = variant_label
                    variant_results.append(result)
                    all_results.append(result)

                    if variant_best is None or result.rmsd < variant_best.rmsd:
                        variant_best = result

                    if result.success:
                        logger.info("")
                        logger.info("=" * 60)
                        logger.info(f"✅ SUCCESS! Protocol {i} achieved RMSD < {self.rmsd_threshold}Å")
                        logger.info(f"   RMSD: {result.rmsd:.2f}Å")
                        logger.info(f"   Score: {result.score:.2f}")
                        logger.info(f"   Runtime: {result.runtime_sec:.1f}s")
                        logger.info("=" * 60)
                        break
                    else:
                        if result.error_message:
                            logger.info(f"⚠️  Protocol {i} failed: {result.error_message}")
                        else:
                            logger.info(f"⚠️  Protocol {i} failed (RMSD: {result.rmsd:.2f}Å)")
                        logger.info(f"   Continuing to next protocol...")
                        logger.info("")

                except Exception as e:
                    logger.error(f"❌ Protocol {i} crashed: {e}")
                    logger.info(f"   Continuing to next protocol...")
                    logger.info("")
                    continue

            if variant_best is None and variant_results:
                variant_best = self._select_best_result(variant_results, by="rmsd")

            if variant_best:
                if best_result is None:
                    best_result = variant_best
                else:
                    best_result = self._select_best_result(
                        [best_result, variant_best],
                        by=self.variant_select_by
                    )

        # Final summary
        if best_result and best_result.success:
            logger.info(f"🎉 Pipeline succeeded! Best RMSD: {best_result.rmsd:.2f}Å")
        elif best_result:
            logger.warning(f"⚠️  Pipeline exhausted. Best RMSD: {best_result.rmsd:.2f}Å (threshold: {self.rmsd_threshold}Å)")
            logger.warning(f"   Tried {len(all_results)} protocols.")
            logger.warning(f"   Consider:")
            logger.warning(f"   - Verify binding site definition")
            logger.warning(f"   - Check ligand SMILES correctness")
            logger.warning(f"   - Try flexible docking")
            logger.warning(f"   - Manual inspection of results")
        else:
            logger.error("❌ All protocols failed or crashed. No results available.")
            logger.warning(f"   Tried {len(all_results)} protocols.")
            logger.warning("   Check errors above for details.")

        return best_result, all_results

    def _run_single_protocol(
        self,
        protocol: DockingProtocol,
        pdb_file: Path,
        ligand_smiles: str,
        ligand_name: str,
        ligand_resname: str,
        ligand_chain: str,
        crystal_ligand_pdb: Path,
        prepared_ligand_pdbqt: Optional[Path] = None,
        prepared_receptor_pdbqt: Optional[Path] = None,
        prepared_receptor_pdb: Optional[Path] = None,
        protocol_root: Optional[Path] = None
    ) -> DockingResult:
        """Run a single docking protocol."""
        start_time = time.time()

        root_dir = protocol_root or self.output_dir
        protocol_dir = root_dir / protocol.name
        protocol_dir.mkdir(parents=True, exist_ok=True)

        # Step 1: Prepare receptor
        logger.info("  [1/5] Preparing receptor...")
        if prepared_receptor_pdbqt and prepared_receptor_pdbqt.exists():
            receptor_pdbqt = prepared_receptor_pdbqt
            logger.info("  Using cached receptor PDBQT: %s", receptor_pdbqt)
        else:
            prep_protein = ProteinPreparation(ProteinPreparationConfig(
                ph=7.4,
                water_handling=WaterHandling.REMOVE_ALL  # Critical: no waters!
            ))
            prepared_protein = prep_protein.prepare_receptor(str(pdb_file))
            receptor_pdbqt = protocol_dir / "receptor.pdbqt"
            prep_protein.to_pdbqt(prepared_protein, str(receptor_pdbqt))

        # Step 2: Prepare ligand
        logger.info("  [2/5] Preparing ligand from SMILES...")
        if prepared_ligand_pdbqt and prepared_ligand_pdbqt.exists():
            ligand_pdbqt = prepared_ligand_pdbqt
            logger.info("  Using cached ligand PDBQT: %s", ligand_pdbqt)
        else:
            prep_ligand = LigandPreparation(LigandPreparationConfig(
                max_tautomers=self.max_tautomers,
                max_conformers=self.max_conformers,
                use_etkdg_v3=True,
                mmff_minimize=True
            ))
            ligand_pdbqt = protocol_dir / f"{ligand_name}.pdbqt"
            prepared_ligands = prep_ligand.prepare_from_smiles(ligand_smiles, ligand_name)
            if not prepared_ligands:
                raise ValueError(f"Failed to prepare ligand from SMILES: {ligand_smiles}")
            # Use the first (lowest energy) variant
            prepared_ligands[0].to_pdbqt_file(str(ligand_pdbqt))

        # Step 3: Define binding site
        logger.info(f"  [3/5] Defining binding site (margin={protocol.margin}Å)...")
        binding_site_def = BindingSiteDefinition(margin=protocol.margin)
        binding_site = binding_site_def.from_cocrystal(
            pdb_file,
            ligand_resname=ligand_resname,
            ligand_chain=ligand_chain
        )

        # Step 4: Run docking
        logger.info(f"  [4/5] Docking with {protocol.engine}...")
        if protocol.engine == "smina":
            result = self._dock_with_smina(
                receptor_pdbqt=receptor_pdbqt,
                ligand_pdbqt=ligand_pdbqt,
                binding_site=binding_site,
                output_dir=protocol_dir,
                scoring=protocol.scoring,
                exhaustiveness=protocol.exhaustiveness
            )
        elif protocol.engine == "rdock":
            result = self._dock_with_rdock(
                receptor_pdb=pdb_file,
                ligand_smiles=ligand_smiles,
                ligand_name=ligand_name,
                binding_site=binding_site,
                output_dir=protocol_dir,
                prepared_receptor_pdb=prepared_receptor_pdb
            )
        else:
            raise ValueError(f"Unknown engine: {protocol.engine}")

        # Step 5: Calculate RMSD
        logger.info("  [5/5] Calculating RMSD...")
        runtime_sec = time.time() - start_time

        if not result.get('success'):
            error_message = result.get('error_message') or "unknown error"
            logger.warning(
                "  Docking failed for %s: %s",
                ligand_name,
                error_message
            )
            return DockingResult(
                protocol_name=protocol.name,
                rmsd=999.9,
                score=result.get('score', 0.0),
                num_poses=result.get('num_poses', 0),
                runtime_sec=runtime_sec,
                output_file=result['output_file'],
                success=False,
                error_message=f"Docking failed: {error_message}"
            )

        try:
            rmsd = calculate_rmsd(crystal_ligand_pdb, result['output_file'])
        except Exception as e:
            logger.warning("  RMSD calculation failed: %s", e)
            return DockingResult(
                protocol_name=protocol.name,
                rmsd=999.9,
                score=result.get('score', 0.0),
                num_poses=result.get('num_poses', 0),
                runtime_sec=runtime_sec,
                output_file=result['output_file'],
                success=False,
                error_message=f"RMSD calculation failed: {e}"
            )

        return DockingResult(
            protocol_name=protocol.name,
            rmsd=rmsd,
            score=result['score'],
            num_poses=result['num_poses'],
            runtime_sec=runtime_sec,
            output_file=result['output_file'],
            success=(rmsd < self.rmsd_threshold),
            error_message=None
        )

    def _prepare_ligand(self, ligand_smiles: str, ligand_name: str) -> Path:
        """Prepare ligand once per pipeline run and return PDBQT path."""
        prep_ligand = LigandPreparation(LigandPreparationConfig(
            max_tautomers=self.max_tautomers,
            max_conformers=self.max_conformers,
            use_etkdg_v3=True,
            mmff_minimize=True
        ))
        enumerate_states = not self._contains_metal(ligand_smiles)
        if not enumerate_states:
            logger.info(
                "Skipping protonation/tautomer enumeration for metal-containing ligand"
            )

        prepared_ligands = prep_ligand.prepare_from_smiles(
            ligand_smiles,
            ligand_name,
            enumerate_states=enumerate_states
        )
        if not prepared_ligands:
            raise ValueError(f"Failed to prepare ligand from SMILES: {ligand_smiles}")

        ligand_pdbqt = self.output_dir / f"{ligand_name}_prepared.pdbqt"
        prepared_ligands[0].to_pdbqt_file(str(ligand_pdbqt))
        return ligand_pdbqt

    def _prepare_ligand_variants(
        self,
        ligand_smiles: str,
        ligand_name: str,
        enumerate_states: bool = True
    ) -> List[dict]:
        prep_ligand = LigandPreparation(LigandPreparationConfig(
            max_tautomers=self.max_tautomers,
            max_conformers=self.max_conformers,
            use_etkdg_v3=True,
            mmff_minimize=True
        ))
        prepared_ligands = prep_ligand.prepare_from_smiles(
            ligand_smiles,
            ligand_name,
            enumerate_states=enumerate_states
        )
        if not prepared_ligands:
            raise ValueError(f"Failed to prepare ligand from SMILES: {ligand_smiles}")

        variants_dir = self.output_dir / "ligand_variants"
        variants_dir.mkdir(parents=True, exist_ok=True)
        variants = []
        for idx, lig in enumerate(prepared_ligands, 1):
            variant_label = f"{ligand_name}_v{idx}"
            ligand_pdbqt = variants_dir / f"{variant_label}.pdbqt"
            lig.to_pdbqt_file(str(ligand_pdbqt))
            variants.append({
                "label": variant_label,
                "pdbqt": ligand_pdbqt,
                "energy": lig.energy,
                "smiles": lig.smiles
            })
        return variants

    def _select_best_variant(self, variants: List[dict]) -> dict:
        if not variants:
            raise ValueError("No ligand variants available")
        def _key(item: dict) -> tuple:
            energy = item.get("energy")
            return (energy is None, energy if energy is not None else 0.0)
        return sorted(variants, key=_key)[0]

    def _select_best_result(self, results: List[DockingResult], by: str) -> Optional[DockingResult]:
        if not results:
            return None
        metric = (lambda r: r.score) if by == "score" else (lambda r: r.rmsd)
        def _key(r: DockingResult) -> float:
            value = metric(r)
            return value if value is not None else float("inf")
        successes = [r for r in results if r.success]
        if successes:
            return min(successes, key=_key)
        return min(results, key=_key)

    def _prepare_receptor(
        self,
        pdb_file: Path,
        water_handling: Optional[WaterHandling] = None
    ) -> Tuple[Path, Path]:
        """Prepare receptor once per pipeline run and return PDBQT and PDB paths."""
        if water_handling is None:
            water_handling = WaterHandling.REMOVE_ALL
        elif isinstance(water_handling, str):
            water_handling = WaterHandling(water_handling)

        prep_protein = ProteinPreparation(ProteinPreparationConfig(
            ph=7.4,
            water_handling=water_handling
        ))
        prepared_protein = prep_protein.prepare_receptor(str(pdb_file))
        receptor_pdbqt = self.output_dir / "receptor_prepared.pdbqt"
        prep_protein.to_pdbqt(prepared_protein, str(receptor_pdbqt))

        receptor_pdb = self.output_dir / "receptor_prepared.pdb"
        with open(receptor_pdbqt) as f_in, open(receptor_pdb, "w") as f_out:
            for line in f_in:
                if line.startswith(("ATOM", "HETATM", "TER", "END")):
                    if line.startswith(("ATOM", "HETATM")):
                        f_out.write(line[:66] + "\n")
                    else:
                        f_out.write(line)

        return receptor_pdbqt, receptor_pdb

    def _contains_metal(self, ligand_smiles: str) -> bool:
        """Return True if SMILES contains common metal ions."""
        try:
            from rdkit import Chem
        except Exception:
            return False

        mol = Chem.MolFromSmiles(ligand_smiles)
        if mol is None:
            return False

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

    def _dock_with_smina(
        self,
        receptor_pdbqt: Path,
        ligand_pdbqt: Path,
        binding_site: 'BindingSite',
        output_dir: Path,
        scoring: str,
        exhaustiveness: int
    ) -> dict:
        """Run Smina docking."""
        from docking_platform_gui.docking.smina import SminaDockingEngine

        engine = SminaDockingEngine(
            smina_binary=self.docking_binary,
            scoring_function=scoring,
            exhaustiveness=exhaustiveness,
            num_modes=20,
            cpu=4,
            seed=42,
            timeout_sec=self.smina_timeout_sec
        )

        output_pdbqt = output_dir / "docked.pdbqt"

        docking_result = engine.dock(
            receptor_path=str(receptor_pdbqt),
            ligand_path=str(ligand_pdbqt),
            center=binding_site.center,
            size=binding_site.size,
            output_path=str(output_pdbqt)
        )

        return {
            'output_file': output_pdbqt,
            'score': docking_result.poses[0].score if docking_result.poses else 0.0,
            'num_poses': len(docking_result.poses),
            'success': docking_result.success,
            'error_message': docking_result.error_message,
            'runtime_sec': docking_result.runtime
        }

    def _dock_with_rdock(
        self,
        receptor_pdb: Path,
        ligand_smiles: str,
        ligand_name: str,
        binding_site: 'BindingSite',
        output_dir: Path,
        prepared_receptor_pdb: Optional[Path] = None,
        radius_override: Optional[float] = None,
        runs: int = 20,
        seed: int = 42
    ) -> dict:
        """
        Run rDock docking.
        """
        start_time = time.time()

        if not self.rdock_root:
            return {
                "success": False,
                "error_message": "rDock root not configured",
                "output_file": output_dir / "docked.sd",
                "score": 0.0,
                "num_poses": 0
            }

        rbdock = self.rdock_root / "bin" / "rbdock"
        rbcavity = self.rdock_root / "bin" / "rbcavity"
        if not rbdock.exists() or not rbcavity.exists():
            return {
                "success": False,
                "error_message": f"rDock binaries not found under {self.rdock_root}",
                "output_file": output_dir / "docked.sd",
                "score": 0.0,
                "num_poses": 0
            }

        obabel = shutil.which("obabel")
        if not obabel:
            return {
                "success": False,
                "error_message": "Open Babel (obabel) not found in PATH",
                "output_file": output_dir / "docked.sd",
                "score": 0.0,
                "num_poses": 0
            }

        env = os.environ.copy()
        env["RBT_ROOT"] = str(self.rdock_root)
        env["RBT_HOME"] = str(Path.home())
        lib_path = str(self.rdock_root / "lib")
        if "DYLD_LIBRARY_PATH" in env:
            env["DYLD_LIBRARY_PATH"] = f"{lib_path}:{env['DYLD_LIBRARY_PATH']}"
        else:
            env["DYLD_LIBRARY_PATH"] = lib_path
        if "LD_LIBRARY_PATH" in env:
            env["LD_LIBRARY_PATH"] = f"{lib_path}:{env['LD_LIBRARY_PATH']}"
        else:
            env["LD_LIBRARY_PATH"] = lib_path

        try:
            # Prepare receptor (PDB -> MOL2)
            if prepared_receptor_pdb and prepared_receptor_pdb.exists():
                receptor_pdb_path = prepared_receptor_pdb
                logger.info("Using cached receptor PDB for rDock: %s", receptor_pdb_path)
            else:
                prep = ProteinPreparation(ProteinPreparationConfig(
                    ph=7.4,
                    water_handling=WaterHandling.REMOVE_ALL
                ))
                prepared = prep.prepare_receptor(str(receptor_pdb))
                receptor_pdbqt = output_dir / "receptor_prepared.pdbqt"
                prep.to_pdbqt(prepared, str(receptor_pdbqt))

                receptor_pdb_path = output_dir / "receptor_prepared.pdb"
                with open(receptor_pdbqt) as f_in, open(receptor_pdb_path, "w") as f_out:
                    for line in f_in:
                        if line.startswith(("ATOM", "HETATM", "TER", "END")):
                            if line.startswith(("ATOM", "HETATM")):
                                f_out.write(line[:66] + "\n")
                            else:
                                f_out.write(line)

            receptor_mol2 = output_dir / "receptor.mol2"
            result = subprocess.run(
                [
                    obabel,
                    str(receptor_pdb_path),
                    "-O", str(receptor_mol2),
                    "-p", "7.4",
                    "-h"
                ],
                capture_output=True,
                text=True,
                timeout=60,
                env=env
            )
            if result.returncode != 0 or not receptor_mol2.exists():
                raise RuntimeError(f"Open Babel receptor conversion failed: {result.stderr}")

            # Prepare ligand SDF from SMILES
            ligand_sdf = output_dir / f"{ligand_name}.sdf"
            from rdkit import Chem, RDLogger
            from rdkit.Chem import AllChem
            from contextlib import contextmanager, nullcontext

            def _write_sdf_with_obabel(smiles: str, output_path: Path) -> None:
                result = subprocess.run(
                    [
                        obabel,
                        f"-:{smiles}",
                        "--gen3d",
                        "-O", str(output_path),
                        "-p", "7.4",
                        "-h"
                    ],
                    capture_output=True,
                    text=True,
                    env=env
                )
                if result.returncode != 0 or not output_path.exists():
                    raise RuntimeError(f"Open Babel ligand conversion failed: {result.stderr}")

            @contextmanager
            def _silence_rdkit_warnings():
                RDLogger.DisableLog("rdApp.warning")
                try:
                    yield
                finally:
                    RDLogger.EnableLog("rdApp.warning")

            metal_ligand = self._contains_metal(ligand_smiles)
            if metal_ligand:
                logger.info("Metal detected in ligand; skipping MMFF minimization for rDock.")

            with _silence_rdkit_warnings() if metal_ligand else nullcontext():
                mol = Chem.MolFromSmiles(ligand_smiles)
                if mol is None:
                    raise RuntimeError("Invalid ligand SMILES")
                mol = Chem.AddHs(mol)

                params = AllChem.ETKDGv3()
                params.randomSeed = int(seed)
                params.useRandomCoords = True
                conf_ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=10, params=params))
                if not conf_ids:
                    conf_id = AllChem.EmbedMolecule(mol, params)
                    conf_ids = [conf_id] if conf_id != -1 else []

                if not conf_ids:
                    _write_sdf_with_obabel(ligand_smiles, ligand_sdf)
                else:
                    best_conf_id = conf_ids[0]
                    if not metal_ligand:
                        try:
                            props = AllChem.MMFFGetMoleculeProperties(mol)
                            energies = []
                            if props:
                                for conf_id in conf_ids:
                                    ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=conf_id)
                                    if ff:
                                        ff.Minimize(maxIts=200)
                                        energies.append((conf_id, ff.CalcEnergy()))
                            if energies:
                                best_conf_id = min(energies, key=lambda x: x[1])[0]
                        except Exception:
                            pass

                    writer = Chem.SDWriter(str(ligand_sdf))
                    writer.write(mol, confId=best_conf_id)
                    writer.close()

            # Create receptor parameter file
            if radius_override is not None:
                radius = float(radius_override)
            else:
                radius = max(binding_site.size) / 2.0
                radius = max(6.0, min(20.0, float(radius)))

            receptor_prm = output_dir / "receptor.prm"
            receptor_prm.write_text(
                f"""RBT_PARAMETER_FILE_V1.00
TITLE {receptor_pdb.stem} {ligand_name}

RECEPTOR_FILE {receptor_mol2.name}
RECEPTOR_FLEX 3.0

##################################################################
### CAVITY DEFINITION: REFERENCE LIGAND METHOD
##################################################################
SECTION MAPPER
    SITE_MAPPER RbtLigandSiteMapper
    REF_MOL {ligand_sdf.name}
    RADIUS {radius:.1f}
    SMALL_SPHERE 1.0
    MIN_VOLUME 100
    MAX_CAVITIES 1
    VOL_INCR 0.0
    GRIDSTEP 0.5
END_SECTION

##################################################################
### CAVITY RESTRAINT PENALTY
##################################################################
SECTION CAVITY
    SCORING_FUNCTION RbtCavityGridSF
    WEIGHT 1.0
END_SECTION
"""
            )

            # Generate cavity
            cavity_file = output_dir / "receptor.as"
            result = subprocess.run(
                [str(rbcavity), "-r", receptor_prm.name, "-was"],
                cwd=str(output_dir),
                capture_output=True,
                text=True,
                env=env
            )
            if result.returncode != 0 or not cavity_file.exists():
                raise RuntimeError(f"rbcavity failed: {result.stderr}")

            # Run docking
            dock_prm = self.rdock_root / "data" / "scripts" / "dock.prm"
            if not dock_prm.exists():
                raise RuntimeError(f"Docking protocol not found: {dock_prm}")

            output_root = output_dir / "docked"
            result = subprocess.run(
                [
                    str(rbdock),
                    "-i", ligand_sdf.name,
                    "-o", output_root.stem,
                    "-r", receptor_prm.name,
                    "-p", str(dock_prm),
                    "-n", str(runs),
                    "-s", str(seed)
                ],
                cwd=str(output_dir),
                capture_output=True,
                text=True,
                env=env
            )
            output_sd = output_dir / "docked.sd"
            if result.returncode != 0 or not output_sd.exists():
                raise RuntimeError(f"rbdock failed: {result.stderr}")

            # Parse scores
            from rdkit import Chem, RDLogger
            RDLogger.DisableLog("rdApp.error")
            RDLogger.DisableLog("rdApp.warning")
            try:
                poses = [
                    mol for mol in Chem.SDMolSupplier(str(output_sd), sanitize=False)
                    if mol is not None
                ]
            finally:
                RDLogger.EnableLog("rdApp.error")
                RDLogger.EnableLog("rdApp.warning")
            scores = []
            for mol in poses:
                props = mol.GetPropsAsDict()
                score = props.get("SCORE", props.get("SCORE.INTER"))
                if score is not None:
                    scores.append(float(score))

            best_score = min(scores) if scores else 0.0

            return {
                "output_file": output_sd,
                "score": best_score,
                "num_poses": len(poses),
                "success": True,
                "error_message": None,
                "runtime_sec": time.time() - start_time
            }
        except Exception as exc:
            return {
                "success": False,
                "error_message": str(exc),
                "output_file": output_dir / "docked.sd",
                "score": 0.0,
                "num_poses": 0,
                "runtime_sec": time.time() - start_time
            }

    def _extract_crystal_ligand(
        self,
        pdb_file: Path,
        ligand_resname: str,
        ligand_chain: str
    ) -> Path:
        """Extract crystal ligand from PDB for RMSD reference."""
        output_pdb = self.output_dir / "crystal_ligand.pdb"
        target_res_id = None
        saw_model = False
        in_first_model = True

        with open(pdb_file) as f_in, open(output_pdb, 'w') as f_out:
            for line in f_in:
                if line.startswith('MODEL'):
                    if saw_model:
                        in_first_model = False
                    else:
                        saw_model = True
                        in_first_model = True
                    continue
                if line.startswith('ENDMDL') and saw_model:
                    break
                if not in_first_model:
                    continue
                if not line.startswith(('HETATM', 'ATOM')):
                    continue

                res = line[17:20].strip()
                if res != ligand_resname:
                    continue
                chain = line[21:22].strip()
                if ligand_chain is not None and chain != ligand_chain:
                    continue
                res_seq = line[22:26].strip()
                ins_code = line[26:27].strip()
                res_id = (res_seq, ins_code)
                if target_res_id is None:
                    target_res_id = res_id
                if res_id == target_res_id:
                    f_out.write(line)

        if output_pdb.stat().st_size == 0:
            raise ValueError(
                f"No ligand found: resname={ligand_resname}, chain={ligand_chain}"
            )

        return output_pdb

    def generate_report(self, all_results: List[DockingResult]) -> str:
        """Generate markdown report of all attempts."""
        report = []
        report.append("# Adaptive Docking Pipeline Report\n")
        report.append(f"**RMSD Threshold:** {self.rmsd_threshold}Å\n")
        report.append(f"**Protocols Attempted:** {len(all_results)}\n")
        report.append("\n---\n")

        for i, result in enumerate(all_results, 1):
            status = "✅ SUCCESS" if result.success else "❌ FAILED"
            report.append(f"\n## Protocol {i}: {result.protocol_name}\n")
            report.append(f"**Status:** {status}\n")
            report.append(f"**RMSD:** {result.rmsd:.2f}Å\n")
            report.append(f"**Score:** {result.score:.2f}\n")
            report.append(f"**Poses:** {result.num_poses}\n")
            report.append(f"**Runtime:** {result.runtime_sec:.1f}s\n")
            report.append(f"**Output:** {result.output_file}\n")

        # Best result
        best = min(all_results, key=lambda r: r.rmsd)
        report.append("\n---\n")
        report.append(f"\n## Best Result: {best.protocol_name}\n")
        report.append(f"**RMSD:** {best.rmsd:.2f}Å\n")
        report.append(f"**Success:** {'Yes' if best.success else 'No'}\n")

        return "".join(report)


# Convenience function
def run_adaptive_docking(
    pdb_file: Path,
    ligand_smiles: str,
    ligand_name: str,
    ligand_resname: str,
    ligand_chain: str,
    output_dir: Path,
    rmsd_threshold: float = 2.0,
    rdock_available: bool = False,
    rdock_root: Optional[Path] = None
) -> Tuple[DockingResult, List[DockingResult]]:
    """
    Convenience function to run adaptive docking pipeline.

    Example:
        best_result, all_results = run_adaptive_docking(
            pdb_file=Path("1gkd.pdb"),
            ligand_smiles="CNC(=O)[C@@H](N)C(C)(C)C",
            ligand_name="BUM",
            ligand_resname="BUM",
            ligand_chain="A",
            output_dir=Path("output/adaptive_docking")
        )

        print(f"Best RMSD: {best_result.rmsd:.2f}Å")
        print(f"Success: {best_result.success}")
    """
    pipeline = AdaptiveDockingPipeline(
        output_dir=output_dir,
        rmsd_threshold=rmsd_threshold,
        rdock_available=rdock_available,
        rdock_root=rdock_root
    )

    return pipeline.run_adaptive_docking(
        pdb_file=pdb_file,
        ligand_smiles=ligand_smiles,
        ligand_name=ligand_name,
        ligand_resname=ligand_resname,
        ligand_chain=ligand_chain
    )
