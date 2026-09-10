#!/usr/bin/env python3
"""Build the structure-free ESR1 supplement from verified frozen local evidence."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import zipfile
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
NAME = "DockMate-VS_ESR1_Supplement_S1"
IDENTIFIERS = {
    "pdb_id", "ligand_resname", "ligand_chain", "mode", "engine", "protocol",
    "dock_name", "control_label", "case_id", "target_name", "status",
    "box_definition", "rescore_method", "water_handling", "site_method",
}
NUMBERS = {
    "best_rmsd", "success", "runtime_sec", "pose_count", "best_score",
    "top1_rmsd", "top5_rmsd", "top10_rmsd", "best_rmsd_rank", "rmsd_best_score",
    "rmsd_mean", "rmsd_median", "rmsd_std", "near_native_fraction",
    "score_rmsd_pearson", "score_rmsd_spearman", "ligand_charge", "rescore_score",
    "rescore_pose_count", "rescore_top1_rmsd", "rescore_top5_rmsd", "rescore_top10_rmsd",
    "rescore_best_rmsd_rank", "rescore_rmsd_best_score", "rescore_score_rmsd_pearson",
    "rescore_score_rmsd_spearman", "docking_completed", "variants_prepared",
    "variants_docked", "molecular_weight", "logp", "tpsa", "rotatable_bonds",
    "exhaustiveness", "seed",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2) + "\n")


def write_csv(path: Path, columns: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def result_table(source: Path, destination: Path) -> dict:
    with source.open(newline="") as handle:
        reader = csv.DictReader(handle)
        columns = [name for name in reader.fieldnames if name in IDENTIFIERS | NUMBERS]
        rows = [{key: row[key] for key in columns} for row in reader]
    for row in rows:
        for key in set(columns) & NUMBERS:
            if row[key] not in {"", "True", "False", "N/A"}:
                float(row[key])
        if any("/Users/" in value or "\\Users\\" in value for value in row.values()):
            raise ValueError(f"Host path survived projection: {source.name}")
    write_csv(destination, columns, rows)
    return {"source_sha256": sha256(source), "rows": len(rows), "columns": columns}


def portable(value):
    if isinstance(value, dict):
        result = {}
        for key, item in value.items():
            if key in {"input_file", "output_dir", "planned_cases", "cases", "resume_signature"}:
                continue
            if key.endswith("_binary"):
                result[key] = key.removesuffix("_binary")
            elif key == "rdock_root":
                result[key] = "<set locally>"
            else:
                result[key] = portable(item)
        return result
    if isinstance(value, list):
        return [portable(item) for item in value]
    return value


def build(frozen: Path, output: Path) -> Path:
    template = REPO / "supplementary/esr1"
    root = output / NAME
    if root.exists() or (output / f"{NAME}.zip").exists():
        raise FileExistsError("Choose a new output directory; an existing supplement is not overwritten.")
    original_hashes = {}
    for line in (frozen / "SHA256SUMS.txt").read_text().splitlines():
        checksum, name = line.split("  ", 1)
        path = frozen / name
        if not path.resolve().is_relative_to(frozen.resolve()):
            raise ValueError(f"Unexpected checksum path: {name}")
        if sha256(path) != checksum:
            raise ValueError(f"Frozen evidence checksum mismatch: {name}")
        original_hashes[name] = checksum
    root.mkdir(parents=True)
    for source in template.rglob("*"):
        if source.is_file() and source.suffix in {".py", ".md", ".txt"}:
            target = root / source.relative_to(template)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
    for name in ("analyze_esr1_seed42.py", "prepare_esr1_ant.py"):
        if sha256(root / "provenance" / name) != original_hashes[f"provenance/{name}"]:
            raise ValueError(f"The archived analysis/preparation code changed: {name}")
    shutil.copy2(REPO / "LICENSE", root / "LICENSE")

    with (frozen / "input/screening_subset80_800_seed42_compounds.csv").open(newline="") as handle:
        compounds = list(csv.DictReader(handle))
    if len(compounds) != 880 or sum(row["label"] == "1" for row in compounds) != 80:
        raise ValueError("Unexpected source subset counts")
    selection = [{"screening_order": i + 1, "source_id": row["Source_ID"]}
                 for i, row in enumerate(compounds)]
    write_csv(root / "selection.csv", ["screening_order", "source_id"], selection)
    cleaning = json.loads((frozen / "input/cleaning_manifest.json").read_text())
    cleaning.update({
        "reconstructed_compounds_csv_sha256": original_hashes[
            "input/screening_subset80_800_seed42_compounds.csv"],
        "structural_inputs": {
            name: {"sha256": original_hashes[f"input/{name}"], "redistributed": False}
            for name in ("1XP1.pdb", "AIH_ideal.sdf")
        },
    })
    write_json(root / "source.json", cleaning)

    projections = {}
    for source_name, target_name in (
        ("results/redock_results.csv", "results/redock_results.csv"),
        ("protocol_development/sweep_results.csv", "protocol_development/protocol_development_results.csv"),
        ("protocol_development/seed_results.csv", "protocol_development/seed_results.csv"),
    ):
        projections[target_name] = result_table(frozen / source_name, root / target_name)

    for filename in ("esr1_seed42_screening_metrics.csv", "esr1_seed42_property_diagnostics.csv",
                     "esr1_seed42_property_model.json"):
        target = root / "publication/tables" / filename
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(frozen / "publication/tables" / filename, target)

    manifest = json.loads((frozen / "results/run_manifest.json").read_text())
    config = portable(manifest["config"])
    config["input_file"] = "reconstructed/screening_subset80_800_seed42.xlsx"
    write_json(root / "campaigns/screening.json", {
        "recorded_software": manifest["software"], "config": config,
        "note": "Historical settings; adaptive mode is unsupported in DockMate-VS 0.1.1.",
    })
    for source_name, name in (("sweep_manifest.json", "protocol-sweep.json"),
                              ("seed_manifest.json", "protocol-seeds.json")):
        config = portable(json.loads((frozen / "protocol_development" / source_name).read_text()))
        config["input_file"] = "reconstructed/protocol_development_1XP1.xlsx"
        write_json(root / "campaigns" / name, config)
    provenance = json.loads((frozen / "provenance/provenance.json").read_text())
    provenance["dataset_source"]["redistribution_note"] = (
        "Compound structures and original input workbooks are not redistributed; "
        "obtain the source archive separately and reconstruct locally."
    )
    provenance["archive_scope"] = {
        "included": "Identifiers, selected result columns, portable settings, statistics and scripts",
        "excluded": "SMILES, coordinate files, workbooks, poses, source archive and external binaries",
    }
    write_json(root / "provenance/provenance.json", provenance)
    write_json(root / "provenance/frozen-evidence-checksums.json", original_hashes)
    write_json(root / "provenance/packaging.json", {
        "result_projections": projections,
        "transformations": [
            "Only allowlisted identifier and numerical result columns retained; values not recomputed.",
            "Host paths and unrelated case lists removed from campaign manifests.",
            "Reference statistical outputs and two archived scripts retained byte-for-byte.",
            "No original input workbooks, structures, poses or analytical figure bitmaps included.",
        ],
    })
    members = sorted(path for path in root.rglob("*") if path.is_file())
    (root / "SHA256SUMS.txt").write_text("".join(
        f"{sha256(path)}  {path.relative_to(root).as_posix()}\n" for path in members
    ))
    archive = output / f"{NAME}.zip"
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as handle:
        for path in sorted(root.rglob("*")):
            if path.is_file():
                handle.write(path, arcname=f"{NAME}/{path.relative_to(root).as_posix()}")
    (output / f"{NAME}.zip.sha256").write_text(f"{sha256(archive)}  {archive.name}\n")
    return archive


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--frozen", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    print(build(args.frozen.expanduser(), args.output.expanduser()))
