"""rDock adapter, including its reference-centred cavity preparation."""

import os
import shutil
import subprocess
import time
from pathlib import Path
from typing import Optional

from loguru import logger
from rdkit import Chem

from .base import DockingEngine, DockingResult, DockingPose
from dockmate_vs.preparation.protein import ProteinPreparation
from dockmate_vs.config.schema import ProteinPreparationConfig, WaterHandling
from dockmate_vs.binding_site.cocrystal import BindingSite


class RDockEngine(DockingEngine):
    def __init__(self, root=None, runs=20, seed=42, radius=None, timeout_sec=1200,
                 obabel_timeout_sec=30.0, binaries=None):
        super().__init__(seed=seed)
        self.rdock_root = Path(root).expanduser() if root else None
        self.runs = runs
        self.radius = radius
        self.timeout_sec = timeout_sec
        self.obabel_timeout_sec = obabel_timeout_sec
        self.binaries = binaries or {}

    @classmethod
    def from_config(cls, config):
        options = config.get("engine_options") or {}
        unexpected = set(options) - {"binaries", "runs", "radius"}
        if unexpected:
            raise ValueError(f"Unsupported rDock engine options: {', '.join(sorted(unexpected))}")
        return cls(
            root=config.get("rdock_root"), runs=options.get("runs", config.get("rdock_runs", 20)),
            seed=config.get("rdock_seed", config.get("seed", 42)),
            radius=options.get("radius", config.get("rdock_radius")),
            timeout_sec=config.get("timeout", 1200),
            obabel_timeout_sec=config.get("obabel_timeout_sec", 30.0),
            binaries=options.get("binaries"),
        )

    def dock(self, receptor_path, ligand_path, center, size, output_path=None, **kwargs):
        output = Path(output_path or Path(ligand_path).with_name("docked.sd"))
        output.parent.mkdir(parents=True, exist_ok=True)
        result = self._run(
            receptor_pdb=Path(receptor_path),
            ligand_smiles=kwargs["ligand_smiles"],
            ligand_name=kwargs.get("ligand_name", Path(ligand_path).stem),
            binding_site=BindingSite(center=center, size=size, ligand_resname="",
                                     ligand_chain="", ligand_atoms=0, source_pdb=str(receptor_path)),
            output_dir=output.parent, prepared_receptor_pdb=kwargs.get("prepared_receptor_pdb"),
            reference_ligand_pdb=kwargs.get("reference_ligand_pdb"),
            radius_override=self.radius, runs=self.runs, seed=self.seed,
        )
        if not result["success"]:
            return DockingResult(ligand_id=Path(ligand_path).stem, success=False,
                                 error_message=result["error_message"], runtime=result.get("runtime_sec"))
        actual_output = Path(result["output_file"])
        if actual_output != output:
            actual_output.replace(output)
        parsed = self.parse_output(str(output))
        parsed.runtime = result.get("runtime_sec")
        return parsed

    def parse_output(self, output_file):
        poses = []
        for index, mol in enumerate(Chem.SDMolSupplier(str(output_file), removeHs=False, sanitize=False), 1):
            if mol is None:
                raise ValueError(f"Unreadable rDock pose {index}")
            score = mol.GetPropsAsDict().get("SCORE", mol.GetPropsAsDict().get("SCORE.INTER"))
            if score is None:
                raise ValueError(f"rDock pose {index} has no score")
            poses.append(DockingPose(score=float(score), mode_number=index,
                                     coordinates=mol.GetConformer().GetPositions()))
        return DockingResult(ligand_id=Path(output_file).stem, poses=poses, success=bool(poses))

    def get_version(self):
        binary = self.binaries.get("rbdock") or (str(self.rdock_root / "bin/rbdock") if self.rdock_root else "rbdock")
        result = subprocess.run([binary, "-h"], capture_output=True, text=True, timeout=5)
        return (result.stdout or result.stderr).strip()

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

    @staticmethod
    def _translate_sdf_to_center(source: Path, destination: Path, center) -> None:
        """Copy an SDF while translating its molecular centroid to ``center``."""
        from rdkit import Chem

        molecule = Chem.SDMolSupplier(
            str(source), removeHs=False, sanitize=False
        )[0]
        if molecule is None or molecule.GetNumConformers() == 0:
            raise RuntimeError("Could not build rDock cavity reference molecule")
        conformer = molecule.GetConformer()
        coordinates = [
            conformer.GetAtomPosition(index)
            for index in range(molecule.GetNumAtoms())
        ]
        centroid = tuple(
            sum(getattr(point, axis) for point in coordinates) / len(coordinates)
            for axis in ("x", "y", "z")
        )
        shift = tuple(float(center[index]) - centroid[index] for index in range(3))
        for index, point in enumerate(coordinates):
            conformer.SetAtomPosition(
                index,
                (point.x + shift[0], point.y + shift[1], point.z + shift[2]),
            )
        writer = Chem.SDWriter(str(destination))
        writer.write(molecule)
        writer.close()

    def _run(
        self,
        receptor_pdb: Path,
        ligand_smiles: str,
        ligand_name: str,
        binding_site: 'BindingSite',
        output_dir: Path,
        prepared_receptor_pdb: Optional[Path] = None,
        reference_ligand_pdb: Optional[Path] = None,
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

        rbdock = Path(self.binaries.get("rbdock") or self.rdock_root / "bin" / "rbdock")
        rbcavity = Path(self.binaries.get("rbcavity") or self.rdock_root / "bin" / "rbcavity")
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
                    obabel_timeout_sec=getattr(self, "obabel_timeout_sec", 30.0),
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
                timeout=self.obabel_timeout_sec,
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
                    env=env, timeout=self.obabel_timeout_sec,
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

            # RbtLigandSiteMapper uses the coordinates of REF_MOL. The docking
            # input conformer is generated around the origin, so using it
            # directly creates a cavity far from most crystallographic sites.
            cavity_reference_sdf = output_dir / "cavity_reference.sdf"
            if reference_ligand_pdb and reference_ligand_pdb.exists():
                result = subprocess.run(
                    [
                        obabel, str(reference_ligand_pdb),
                        "-O", str(cavity_reference_sdf), "-h",
                    ],
                    capture_output=True,
                    text=True,
                    timeout=self.obabel_timeout_sec,
                    env=env,
                )
                if result.returncode != 0 or not cavity_reference_sdf.exists():
                    raise RuntimeError(
                        "Open Babel crystal-ligand conversion failed: " + result.stderr
                    )
            else:
                self._translate_sdf_to_center(
                    ligand_sdf, cavity_reference_sdf, binding_site.center
                )

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
    REF_MOL {cavity_reference_sdf.name}
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
                env=env, timeout=self.timeout_sec,
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
                env=env, timeout=self.timeout_sec,
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
