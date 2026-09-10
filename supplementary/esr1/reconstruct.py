#!/usr/bin/env python3
"""Reconstruct the exact ESR1 subset from a separately obtained source archive."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import tarfile
import tempfile
from pathlib import Path

import pandas as pd

from provenance import prepare_esr1_ant as preparation


ROOT = Path(__file__).resolve().parent
MEMBERS = {
    "actives.smi": "LIT-PCBA_full/ESR1_ant/actives.smi",
    "inactives.smi": "LIT-PCBA_full/ESR1_ant/inactives.smi",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_sources(archive: Path, destination: Path) -> None:
    # Read only exact, regular-file members; never extract archive paths or links.
    with tarfile.open(archive, "r:gz") as handle:
        members = handle.getmembers()
        for name, expected in MEMBERS.items():
            matches = [member for member in members if member.name == expected]
            if len(matches) != 1:
                raise ValueError(f"Expected one archive member: {expected}")
            member = matches[0]
            if not member.isfile() or member.size > 10 * 1024 * 1024:
                raise ValueError(f"Invalid source member: {expected}")
            with handle.extractfile(member) as stream:
                (destination / name).write_bytes(stream.read())


def reconstruct(archive: Path, output: Path, root: Path = ROOT) -> dict:
    metadata = json.loads((root / "source.json").read_text())
    if sha256(archive) != metadata["archive_sha256"]:
        raise ValueError("Source archive SHA-256 differs from the historical archive.")
    if output.exists():
        raise FileExistsError(f"Choose a new output directory: {output}")

    with tempfile.TemporaryDirectory(prefix="esr1-reconstruct-") as directory:
        work = Path(directory)
        read_sources(archive, work)
        active_raw = preparation.load_smiles(work / "actives.smi", "active", 1)
        inactive_raw = preparation.load_smiles(work / "inactives.smi", "inactive", 0)
        actives, active_duplicates = preparation.deduplicate(active_raw)
        inactives, inactive_duplicates = preparation.deduplicate(inactive_raw)
        conflicts = {r.canonical_isomeric_smiles for r in actives} & {
            r.canonical_isomeric_smiles for r in inactives
        }
        active_clean = [r for r in actives if r.canonical_isomeric_smiles not in conflicts]
        inactive_clean = [r for r in inactives if r.canonical_isomeric_smiles not in conflicts]
        counts = {
            "raw_actives": len(active_raw), "raw_inactives": len(inactive_raw),
            "duplicate_active_records": len(active_duplicates),
            "duplicate_inactive_records": len(inactive_duplicates),
            "unique_active_structures_before_conflict_removal": len(actives),
            "unique_inactive_structures_before_conflict_removal": len(inactives),
            "label_conflicting_structures": len(conflicts),
            "clean_actives": len(active_clean), "clean_inactives": len(inactive_clean),
            "screening_jobs": len(active_clean) + len(inactive_clean),
        }
        if counts != metadata["counts"]:
            raise ValueError(f"Cleaning counts differ from the historical values: {counts}")
        active_frame = pd.DataFrame(preparation.screening_rows(active_clean))
        inactive_frame = pd.DataFrame(preparation.screening_rows(inactive_clean))
        subset = pd.concat([
            active_frame,
            inactive_frame.sample(n=800, replace=False, random_state=42),
        ], ignore_index=True).sample(frac=1, random_state=42).reset_index(drop=True)
        with (root / "selection.csv").open(newline="") as handle:
            expected = list(csv.DictReader(handle))
        selection = [
            {"screening_order": str(index + 1), "source_id": str(source_id)}
            for index, source_id in enumerate(subset["Source_ID"])
        ]
        if selection != expected:
            raise ValueError("Reconstructed subset/order differs from selection.csv.")
        compound_csv = work / "screening_subset80_800_seed42_compounds.csv"
        subset.to_csv(compound_csv, index=False)
        checksum = sha256(compound_csv)
        if checksum != metadata["reconstructed_compounds_csv_sha256"]:
            raise ValueError("Reconstructed CSV differs; check the recorded RDKit/pandas versions.")

        # Publish outputs only after the structures, order and checksum are verified.
        output.mkdir(parents=True)
        (output / compound_csv.name).write_bytes(compound_csv.read_bytes())
        workbook_metadata = [
            ("Dataset", "LIT-PCBA ESR1 antagonist; reconstructed manuscript subset"),
            ("Source archive SHA-256", metadata["archive_sha256"]),
            ("Compound CSV SHA-256", checksum),
            ("Sampling seed", 42),
            ("Campaign settings", "See campaigns/; historical adaptive mode is not supported in 0.1.1"),
        ]
        preparation.write_workbook(
            output / "screening_subset80_800_seed42.xlsx", subset, workbook_metadata
        )
        protocol = pd.DataFrame([{
            "Protein": "ESR1 antagonist", "Target_Ligand": "AIH_native_redock",
            "PDB_ID": "1XP1", "Ligand": "AIH", "Chain": "A", "SMILES": "",
            "label": 1,
            "Notes": "Native 1XP1/AIH pose-recovery control; SMILES resolved from RCSB",
        }])
        preparation.write_workbook(output / "protocol_development_1XP1.xlsx", protocol, workbook_metadata)
        report = {"counts": counts, "subset_actives": 80, "subset_inactives": 800,
                  "compound_csv_sha256": checksum, "selection_order_verified": True}
        (output / "reconstruction_report.json").write_text(json.dumps(report, indent=2) + "\n")
        return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    print(json.dumps(reconstruct(args.archive.expanduser(), args.output.expanduser()), indent=2))


if __name__ == "__main__":
    main()
