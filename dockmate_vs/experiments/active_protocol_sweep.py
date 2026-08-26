#!/usr/bin/env python3
"""Resumable active-only sweep of water handling and search exhaustiveness."""

from __future__ import annotations

import argparse
import json
import queue
from dataclasses import asdict
from pathlib import Path

import pandas as pd
from loguru import logger

from dockmate_vs.gui.app import DockMateVSApp


WATER_MODES = ("remove_all", "retain_all", "selective")
EXHAUSTIVENESS = (8, 16, 32)
PROTOCOL_VERSION = 2


def _active_rows(template: Path) -> list[dict]:
    frame = pd.read_excel(template)
    controls = frame[frame["decoy compound"].notna()].copy()
    rows = []
    seen = set()
    for _, row in controls.iterrows():
        key = (str(row["PDB_ID"]).upper(), str(row["Ligand"]).upper())
        if key in seen:
            continue
        seen.add(key)
        rows.append({
            "pdb_id": key[0],
            "ligand": key[1],
            "name": str(row["Target_Ligand"]),
            "smiles": str(row["SMILES"]),
        })
    return rows


def _write(rows: list[dict], path: Path) -> None:
    temporary = path.with_suffix(".tmp")
    pd.DataFrame(rows).to_csv(temporary, index=False)
    temporary.replace(path)


def _ligand_chain(pdb_file: Path, resname: str) -> str:
    """Choose one co-crystal copy when a biological assembly has duplicates."""
    for line in pdb_file.read_text().splitlines():
        if line.startswith("HETATM") and line[17:20].strip().upper() == resname:
            return line[21].strip() or "A"
    raise ValueError(f"Ligand {resname} not found in {pdb_file}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--template", type=Path, required=True)
    parser.add_argument("--pdb-dir", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--smina", required=True)
    parser.add_argument("--vina", default="vina")
    parser.add_argument("--cpu", type=int, default=4)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    csv_path = args.outdir / "active_protocol_sweep.csv"
    rows = pd.read_csv(csv_path).to_dict("records") if csv_path.exists() else []
    completed = {
        (str(row["pdb_id"]), str(row["water_handling"]), int(row["exhaustiveness"]))
        for row in rows
        if row.get("status") in ("complete", "unsupported")
        and int(row.get("protocol_version", -1)) == PROTOCOL_VERSION
    }
    app = object.__new__(DockMateVSApp)
    app._queue = queue.Queue()

    manifest = {
        "template": str(args.template),
        "pdb_dir": str(args.pdb_dir),
        "water_modes": WATER_MODES,
        "exhaustiveness": EXHAUSTIVENESS,
        "num_modes": 20,
        "max_tautomers": 2,
        "max_conformers": 3,
        "seed": 42,
        "variant_mode": "adaptive",
        "variant_selection": "rmsd",
        "protocol_version": PROTOCOL_VERSION,
    }
    (args.outdir / "manifest.json").write_text(json.dumps(manifest, indent=2))

    for active in _active_rows(args.template):
        pdb_file = args.pdb_dir / f"{active['pdb_id']}.pdb"
        if not pdb_file.exists():
            raise FileNotFoundError(pdb_file)
        ligand_chain = _ligand_chain(pdb_file, active["ligand"])
        for water in WATER_MODES:
            for exhaustiveness in EXHAUSTIVENESS:
                key = (active["pdb_id"], water, exhaustiveness)
                if key in completed:
                    logger.info("Skipping completed condition: {}", key)
                    continue
                case_dir = args.outdir / f"{active['pdb_id']}_{active['ligand']}_{water}_e{exhaustiveness}"
                config = {
                    "engine": "smina", "rdock_root": ".",
                    "smina_binary": args.smina, "vina_binary": args.vina,
                    "ligand_variant_mode": "adaptive", "variant_select_by": "rmsd",
                    "max_tautomers": 2, "max_conformers": 3, "n_cpus": args.cpu,
                    "water_handling": water, "size_override": None, "box_margin": 4.0,
                    "exhaustiveness": exhaustiveness, "num_modes": 20,
                    "energy_range": 3.0, "cpu": args.cpu, "seed": 42,
                    "timeout": 3600, "scoring": "vina", "apo_site_mode": "auto",
                    "rdock_radius": None, "rdock_runs": 20, "rdock_seed": 42,
                }
                logger.info("Running {}", key)
                row = {
                    **active, "water_handling": water,
                    "exhaustiveness": exhaustiveness, "status": "failed",
                    "protocol_version": PROTOCOL_VERSION,
                }
                try:
                    if app._has_covalent_ligand_link(
                        pdb_file, active["ligand"], ligand_chain
                    ):
                        row["status"] = "unsupported"
                        row["error"] = (
                            "Covalent receptor-ligand LINK: standard Vina/Smina docking is invalid"
                        )
                        raise StopIteration
                    result = app._run_single_case(
                        pdb_file=pdb_file, ligand_name=active["name"], ligand_chain=ligand_chain,
                        smiles=active["smiles"], case_dir=case_dir, threshold=2.0,
                        single_cfg=config, ligand_resname=active["ligand"],
                        site_mode="cocrystal", run_mode="single",
                    )
                    row.update(asdict(result))
                    row["status"] = "complete"
                except StopIteration:
                    logger.warning("Unsupported covalent complex: {}", key)
                except Exception as exc:
                    logger.exception("Condition failed: {}", key)
                    row["error"] = str(exc)
                rows = [
                    old for old in rows
                    if (str(old.get("pdb_id")), str(old.get("water_handling")), int(old.get("exhaustiveness", -1))) != key
                ]
                rows.append(row)
                _write(rows, csv_path)

    logger.info("Sweep complete: {}", csv_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
