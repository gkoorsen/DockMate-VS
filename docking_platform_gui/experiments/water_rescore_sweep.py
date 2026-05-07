#!/usr/bin/env python3
"""
Sweep water handling modes and rescoring methods for a single active+decoy set.

Uses the same core components as the GUI (ProteinPreparation, BindingSiteDefinition,
LigandPreparation, Smina docking/rescore) to isolate effects of:
  - Water handling: remove_all, retain_all, selective
  - Rescoring: vina, vinardo, ad4_scoring, dkoes_fast, dkoes_scoring
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Optional
import subprocess
import tempfile
from datetime import datetime

import numpy as np
import pandas as pd
from loguru import logger
from rdkit import Chem

from docking_platform_gui.adaptive_docking import AdaptiveDockingPipeline
from docking_platform_gui.binding_site.cocrystal import BindingSiteDefinition
from docking_platform_gui.config.schema import (
    ProteinPreparationConfig,
    WaterHandling,
    WaterRetentionConfig,
    LigandPreparationConfig
)
from docking_platform_gui.docking.smina import SminaDockingEngine
from docking_platform_gui.preparation.ligand import LigandPreparation


RESCORE_METHODS = ["vina", "vinardo", "ad4_scoring", "dkoes_fast", "dkoes_scoring"]
WATER_MODES = [WaterHandling.REMOVE_ALL, WaterHandling.RETAIN_ALL, WaterHandling.SELECTIVE]


def _load_active_smiles(crystal_ligand_pdb: Path) -> str:
    mol = Chem.MolFromPDBFile(str(crystal_ligand_pdb), removeHs=False)
    if mol is None:
        mol = Chem.MolFromPDBFile(str(crystal_ligand_pdb), removeHs=False, sanitize=False)
    if mol is None:
        raise RuntimeError(f"Failed to read crystal ligand: {crystal_ligand_pdb}")
    mol = Chem.RemoveHs(mol)
    return Chem.MolToSmiles(mol)


def _read_template_pairs(template_path: Path) -> List[Dict[str, str]]:
    df = pd.read_excel(template_path)
    col_map = {c.strip().lower(): c for c in df.columns}
    pdb_col = col_map.get("pdb_id")
    ligand_col = col_map.get("ligand")
    decoy_name_col = col_map.get("decoy compound")
    decoy_smiles_col = col_map.get("decoy smiles")
    if not pdb_col or not ligand_col:
        raise ValueError("Template missing PDB_ID or Ligand columns.")
    if not decoy_name_col or not decoy_smiles_col:
        raise ValueError("Template missing decoy compound or decoy SMILES columns.")

    pairs = []
    for _, row in df.iterrows():
        pdb_id = str(row[pdb_col]).strip().upper()
        ligand = str(row[ligand_col]).strip().upper()
        decoy_name = str(row[decoy_name_col]).strip()
        decoy_smiles = str(row[decoy_smiles_col]).strip()
        if not pdb_id or len(pdb_id) != 4 or not ligand:
            continue
        if not decoy_name or decoy_name.lower() == "nan":
            continue
        if not decoy_smiles or decoy_smiles.lower() == "nan":
            continue
        pairs.append(
            {
                "pdb_id": pdb_id,
                "ligand": ligand,
                "decoy_name": decoy_name,
                "decoy_smiles": decoy_smiles
            }
        )
    if not pairs:
        raise ValueError("No decoy rows found in template.")
    # assume single PDB/ligand in this template
    return pairs


def _rescore_with_smina(
    smina_binary: str,
    receptor_pdbqt: Path,
    ligand_pdbqt: Path,
    scoring: str
) -> Optional[float]:
    # smina --score_only expects a single model without MODEL/ENDMDL wrappers
    ligand_to_score = ligand_pdbqt
    try:
        first_model_lines: List[str] = []
        in_model = False
        with open(ligand_pdbqt, "r") as handle:
            for line in handle:
                if line.startswith("MODEL"):
                    in_model = True
                    first_model_lines = []
                    continue
                if in_model:
                    if line.startswith("ENDMDL"):
                        break
                    first_model_lines.append(line)
                elif line.startswith("END"):
                    break
        if not first_model_lines:
            with open(ligand_pdbqt, "r") as handle:
                for line in handle:
                    if line.startswith("END"):
                        break
                    first_model_lines.append(line)
        if first_model_lines:
            with tempfile.NamedTemporaryFile(mode="w", suffix=".pdbqt", delete=False) as tmp:
                tmp.write("".join(first_model_lines))
                ligand_to_score = Path(tmp.name)
    except Exception:
        ligand_to_score = ligand_pdbqt

    cmd = [
        smina_binary,
        "--score_only",
        "--receptor", str(receptor_pdbqt),
        "--ligand", str(ligand_to_score),
        "--scoring", scoring
    ]
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=False,
            timeout=300
        )
    except Exception as exc:
        logger.warning("Rescore failed: {}", exc)
        return None
    if result.returncode != 0:
        logger.warning("Rescore failed ({}): {}", scoring, result.stderr.strip())
        return None
    for line in result.stdout.splitlines():
        if line.startswith("REMARK VINA RESULT"):
            parts = line.split()
            if len(parts) >= 4:
                try:
                    return float(parts[3])
                except ValueError:
                    return None
    return None


def _detect_supported_scoring(smina_binary: str) -> Optional[List[str]]:
    cmd = [smina_binary, "--help"]
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=False,
            timeout=10
        )
    except Exception as exc:
        logger.warning("Unable to query smina scoring methods: {}", exc)
        return None
    help_text = (result.stdout or "") + "\n" + (result.stderr or "")
    for line in help_text.splitlines():
        if "--scoring" in line and ":" in line:
            # Example: --scoring (default = vina): vina, vinardo, ad4_scoring, ...
            parts = line.split(":")
            if len(parts) >= 2:
                candidates = [p.strip() for p in parts[-1].split(",")]
                values = [c for c in candidates if c]
                if values:
                    return values
    return None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb", required=True, help="Input PDB (e.g., 2WEQ.pdb)")
    parser.add_argument("--crystal-ligand", required=True, help="Crystal ligand PDB")
    parser.add_argument("--template", required=True, help="Template Excel file")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--smina", default="smina", help="Path to smina binary")
    parser.add_argument("--vina", default="vina", help="Path to vina binary")
    parser.add_argument("--box-margin", type=float, default=4.0, help="Binding site margin (Å)")
    parser.add_argument("--exhaustiveness", type=int, default=8)
    parser.add_argument("--max-tautomers", type=int, default=2)
    parser.add_argument("--max-conformers", type=int, default=2)
    parser.add_argument("--num-modes", type=int, default=9)
    parser.add_argument("--energy-range", type=float, default=3.0)
    parser.add_argument("--cpu", type=int, default=1)
    args = parser.parse_args()

    pdb_path = Path(args.pdb)
    crystal_ligand = Path(args.crystal_ligand)
    template_path = Path(args.template)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out_csv = outdir / "water_rescore_sweep.csv"

    supported = _detect_supported_scoring(args.smina)
    if supported:
        logger.info("Smina supported scoring: {}", ", ".join(supported))
        rescore_methods = [m for m in RESCORE_METHODS if m in supported]
    else:
        logger.warning("Unable to detect smina scoring methods; using default list")
        rescore_methods = RESCORE_METHODS

    if not rescore_methods:
        raise RuntimeError("No supported rescoring methods available in this smina build.")

    pairs = _read_template_pairs(template_path)
    pdb_id = pairs[0]["pdb_id"]
    ligand_name = pairs[0]["ligand"]
    active_smiles = _load_active_smiles(crystal_ligand)

    ligand_prep = LigandPreparation(
        LigandPreparationConfig(
            max_tautomers=args.max_tautomers,
            max_conformers=args.max_conformers
        ),
        n_cpus=args.cpu
    )

    results = []
    existing_rows: Dict[tuple, Dict[str, object]] = {}
    if out_csv.exists():
        existing_df = pd.read_csv(out_csv)
        for _, row in existing_df.iterrows():
            key = (row.get("water_handling"), row.get("compound"))
            existing_rows[key] = row.to_dict()
        logger.info("Loaded existing results: {} rows", len(existing_rows))

    for water_mode in WATER_MODES:
        logger.info("Water handling: {}", water_mode.value)
        prep_cfg = ProteinPreparationConfig(
            water_handling=water_mode,
            water_retention=WaterRetentionConfig() if water_mode == WaterHandling.SELECTIVE else None
        )

        pipeline = AdaptiveDockingPipeline(
            output_dir=outdir,
            smina_binary=args.smina,
            vina_binary=args.vina
        )

        receptor_pdbqt, _ = pipeline._prepare_receptor(
            str(pdb_path),
            water_handling=water_mode
        )

        binding_site = BindingSiteDefinition(margin=args.box_margin).from_cocrystal(
            str(pdb_path),
            ligand_resname=ligand_name
        )

        # Build compound list: active + decoys
        compounds = [
            {"name": ligand_name, "smiles": active_smiles, "is_active": True}
        ]
        for row in pairs:
            compounds.append(
                {
                    "name": row["decoy_name"],
                    "smiles": row["decoy_smiles"],
                    "is_active": False
                }
            )

        for comp in compounds:
            label = f"{pdb_id}_{ligand_name}_{comp['name']}_{water_mode.value}"
            case_dir = outdir / label
            case_dir.mkdir(parents=True, exist_ok=True)
            key = (water_mode.value, comp["name"])
            existing = existing_rows.get(key, {})

            output_pdbqt = case_dir / "docked.pdbqt"
            if not output_pdbqt.exists():
                prepared = ligand_prep.prepare_from_smiles(comp["smiles"], comp["name"], enumerate_states=True)
                if not prepared:
                    logger.warning("Ligand prep failed: {}", comp["name"])
                    continue
                ligand_pdbqt = case_dir / f"{comp['name']}.pdbqt"
                prepared[0].to_pdbqt_file(str(ligand_pdbqt))

                engine = SminaDockingEngine(
                    smina_binary=args.smina,
                    scoring_function="vina",
                    exhaustiveness=args.exhaustiveness,
                    num_modes=args.num_modes,
                    energy_range=args.energy_range,
                    cpu=args.cpu
                )

                docking_result = engine.dock(
                    receptor_path=str(receptor_pdbqt),
                    ligand_path=str(ligand_pdbqt),
                    center=binding_site.center,
                    size=binding_site.size,
                    output_path=str(output_pdbqt)
                )
            else:
                docking_result = None

            row = {
                "water_handling": water_mode.value,
                "compound": comp["name"],
                "is_active": comp["is_active"],
                "dock_score": docking_result.best_score if docking_result else existing.get("dock_score"),
                "output_file": str(output_pdbqt) if output_pdbqt.exists() else existing.get("output_file"),
                "updated_at": datetime.utcnow().isoformat()
            }

            if output_pdbqt.exists():
                for method in rescore_methods:
                    existing_val = existing.get(f"rescore_{method}")
                    if existing_val is not None and not pd.isna(existing_val):
                        row[f"rescore_{method}"] = existing_val
                        continue
                    score = _rescore_with_smina(
                        smina_binary=args.smina,
                        receptor_pdbqt=receptor_pdbqt,
                        ligand_pdbqt=output_pdbqt,
                        scoring=method
                    )
                    row[f"rescore_{method}"] = score

            results.append(row)

    pd.DataFrame(results).to_csv(out_csv, index=False)
    logger.info("Saved results: {}", out_csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
