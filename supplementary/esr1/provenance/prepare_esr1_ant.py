#!/usr/bin/env python3
"""Prepare the local LIT-PCBA ESR1 antagonist SoftwareX benchmark."""

from __future__ import annotations

import csv
import hashlib
import json
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import pandas as pd
from rdkit import Chem


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "source" / "LIT-PCBA_full" / "ESR1_ant"
ARCHIVE = ROOT / "LIT-PCBA_full.tar.gz"
OUTPUT = ROOT / "prepared"

SOURCE_URL = "https://drugdesign.unistra.fr/downloads/datasets/LIT-PCBA_full.tar.gz"
DATASET_PAGE = "https://lab.drugdesign.unistra.fr/datasets/lit-pcba/"
PAPER_DOI = "10.1021/acs.jcim.0c00155"

RECEPTOR = {
    "pdb_id": "1XP1",
    "ligand": "AIH",
    "chain": "A",
    "resolution_angstrom": 1.8,
    "r_free": 0.269,
    "selection": (
        "Prespecified structural-QC choice: wild-type, noncovalent ESR1 ligand-binding "
        "domain; 1.8 A resolution; only terminal residues 552-554 reported missing."
    ),
}


@dataclass(frozen=True)
class Record:
    activity_class: str
    label: int
    source_id: str
    source_line: int
    smiles: str
    canonical_isomeric_smiles: str


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def source_id_key(value: str) -> Tuple[int, str]:
    try:
        return int(value), value
    except ValueError:
        return 10**30, value


def load_smiles(path: Path, activity_class: str, label: int) -> List[Record]:
    records = []
    for line_number, line in enumerate(path.read_text().splitlines(), 1):
        fields = line.split()
        if len(fields) < 2:
            raise ValueError(f"Malformed record at {path}:{line_number}")
        smiles, source_id = fields[0], fields[1]
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            raise ValueError(f"Invalid SMILES at {path}:{line_number}: {smiles}")
        records.append(
            Record(
                activity_class=activity_class,
                label=label,
                source_id=source_id,
                source_line=line_number,
                smiles=smiles,
                canonical_isomeric_smiles=Chem.MolToSmiles(
                    molecule, canonical=True, isomericSmiles=True
                ),
            )
        )
    return records


def deduplicate(records: Iterable[Record]) -> Tuple[List[Record], List[dict]]:
    groups: Dict[str, List[Record]] = defaultdict(list)
    for record in records:
        groups[record.canonical_isomeric_smiles].append(record)

    retained = []
    excluded = []
    for canonical, members in groups.items():
        ordered = sorted(members, key=lambda item: source_id_key(item.source_id))
        keep = ordered[0]
        retained.append(keep)
        for duplicate in ordered[1:]:
            excluded.append(
                {
                    **asdict(duplicate),
                    "retained_source_id": keep.source_id,
                    "exclusion_reason": "duplicate canonical isomeric structure within class",
                }
            )
    retained.sort(key=lambda item: source_id_key(item.source_id))
    return retained, excluded


def write_csv(path: Path, rows: List[dict]) -> None:
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def screening_rows(records: Iterable[Record]) -> List[dict]:
    rows = []
    for record in records:
        rows.append(
            {
                "Protein": "ESR1 antagonist",
                "Target_Ligand": f"LITPCBA_{record.activity_class}_{record.source_id}",
                "PDB_ID": RECEPTOR["pdb_id"],
                "Ligand": RECEPTOR["ligand"],
                "Chain": RECEPTOR["chain"],
                "SMILES": record.smiles,
                "label": record.label,
                "Source_ID": record.source_id,
                "Activity_Class": record.activity_class,
                "Canonical_Isomeric_SMILES": record.canonical_isomeric_smiles,
                "Notes": "LIT-PCBA ESR1_ant; cleaned by prepare_esr1_ant.py",
            }
        )
    return rows


def write_workbook(path: Path, jobs: pd.DataFrame, metadata: List[Tuple[str, object]]) -> None:
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        jobs.to_excel(writer, sheet_name="Docking_Jobs", index=False)
        pd.DataFrame(metadata, columns=["Field", "Value"]).to_excel(
            writer, sheet_name="Metadata", index=False
        )
        sheet = writer.book["Docking_Jobs"]
        sheet.freeze_panes = "A2"
        sheet.auto_filter.ref = sheet.dimensions


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    active_raw = load_smiles(SOURCE / "actives.smi", "active", 1)
    inactive_raw = load_smiles(SOURCE / "inactives.smi", "inactive", 0)
    active_unique, active_duplicates = deduplicate(active_raw)
    inactive_unique, inactive_duplicates = deduplicate(inactive_raw)

    active_by_structure = {item.canonical_isomeric_smiles: item for item in active_unique}
    inactive_by_structure = {item.canonical_isomeric_smiles: item for item in inactive_unique}
    conflicts = sorted(set(active_by_structure) & set(inactive_by_structure))

    conflict_rows = []
    for canonical in conflicts:
        for record in (active_by_structure[canonical], inactive_by_structure[canonical]):
            conflict_rows.append(
                {
                    **asdict(record),
                    "exclusion_reason": "same canonical isomeric structure has both labels",
                }
            )

    active_clean = [
        item for item in active_unique if item.canonical_isomeric_smiles not in conflicts
    ]
    inactive_clean = [
        item for item in inactive_unique if item.canonical_isomeric_smiles not in conflicts
    ]
    cleaned = active_clean + inactive_clean

    duplicate_rows = active_duplicates + inactive_duplicates
    write_csv(OUTPUT / "excluded_duplicates.csv", duplicate_rows)
    write_csv(OUTPUT / "excluded_label_conflicts.csv", conflict_rows)
    write_csv(OUTPUT / "cleaned_compounds.csv", screening_rows(cleaned))

    archive_hash = sha256(ARCHIVE)
    audit = {
        "dataset": "LIT-PCBA ESR1_ant",
        "source_url": SOURCE_URL,
        "dataset_page": DATASET_PAGE,
        "paper_doi": PAPER_DOI,
        "archive_sha256": archive_hash,
        "canonicalization": "RDKit canonical isomeric SMILES",
        "deduplication": "within each class; lowest numeric source ID retained",
        "label_conflicts": "excluded from both classes",
        "counts": {
            "raw_actives": len(active_raw),
            "raw_inactives": len(inactive_raw),
            "duplicate_active_records": len(active_duplicates),
            "duplicate_inactive_records": len(inactive_duplicates),
            "unique_active_structures_before_conflict_removal": len(active_unique),
            "unique_inactive_structures_before_conflict_removal": len(inactive_unique),
            "label_conflicting_structures": len(conflicts),
            "clean_actives": len(active_clean),
            "clean_inactives": len(inactive_clean),
            "screening_jobs": len(cleaned),
        },
        "receptor": RECEPTOR,
    }
    (OUTPUT / "benchmark_manifest.json").write_text(json.dumps(audit, indent=2) + "\n")

    common_metadata = [
        ("Dataset", "LIT-PCBA ESR1_ant"),
        ("Dataset page", DATASET_PAGE),
        ("Source archive", SOURCE_URL),
        ("Archive SHA-256", archive_hash),
        ("Paper DOI", PAPER_DOI),
        ("Receptor", f"{RECEPTOR['pdb_id']} chain {RECEPTOR['chain']}"),
        ("Co-crystal ligand", RECEPTOR["ligand"]),
        ("Receptor selection", RECEPTOR["selection"]),
        ("Cleaning", audit["canonicalization"] + "; " + audit["label_conflicts"]),
    ]

    screening = pd.DataFrame(screening_rows(cleaned))
    write_workbook(
        OUTPUT / "LIT-PCBA_ESR1_ant_screening_1XP1.xlsx",
        screening,
        common_metadata
        + [
            ("Clean actives", len(active_clean)),
            ("Clean inactives", len(inactive_clean)),
            ("Total docking jobs", len(cleaned)),
            ("Analysis", "Use label 1/0 for enrichment; assay-active RMSD is undefined"),
        ],
    )

    protocol = pd.DataFrame(
        [
            {
                "Protein": "ESR1 antagonist",
                "Target_Ligand": "AIH_native_redock",
                "PDB_ID": RECEPTOR["pdb_id"],
                "Ligand": RECEPTOR["ligand"],
                "Chain": RECEPTOR["chain"],
                "SMILES": "",
                "label": 1,
                "Notes": "Native 1XP1/AIH pose-recovery control; SMILES resolved from RCSB",
            }
        ]
    )
    write_workbook(
        OUTPUT / "LIT-PCBA_ESR1_ant_protocol_development_1XP1.xlsx",
        protocol,
        common_metadata
        + [
            ("Purpose", "Native-ligand pose recovery before enrichment screening"),
            ("Expected jobs", 1),
        ],
    )

    readme = f"""# LIT-PCBA ESR1 antagonist benchmark

This local benchmark uses the downloaded official LIT-PCBA archive and receptor
`1XP1`, chain `A`, co-crystal ligand `AIH`. The prepared screening workbook is
a deterministic cleaned derivative of the archive, not an unchanged copy of the
published benchmark table.

## Prepared inputs

- `LIT-PCBA_ESR1_ant_protocol_development_1XP1.xlsx`: native AIH redocking.
- `LIT-PCBA_ESR1_ant_screening_1XP1.xlsx`: {len(active_clean)} cleaned actives and
  {len(inactive_clean)} cleaned inactives ({len(cleaned)} jobs).
- `benchmark_manifest.json`: provenance, checksum, receptor choice, and counts.
- `excluded_duplicates.csv`: duplicate source records removed within a class.
- `excluded_label_conflicts.csv`: exact structures excluded from both classes.

## Audit counts

- Raw archive records: {len(active_raw)} actives and {len(inactive_raw)} inactives.
- Unique structures before conflict removal: {len(active_unique)} actives and
  {len(inactive_unique)} inactives.
- Exact active/inactive label conflicts removed from both classes: {len(conflicts)}.
- Cleaned screening set: {len(active_clean)} actives and {len(inactive_clean)}
  inactives.

The LIT-PCBA dataset page reports ESR_antago summary counts separately from this
downloaded archive. Use `benchmark_manifest.json` and the archive SHA-256 when
describing this local reproducibility example.

## GUI workflow

1. Run the protocol-development workbook in **Protocol Development** mode.
2. Lock the protocol before examining enrichment.
3. Run the screening workbook in **Screen compounds** mode with **Use SMILES**.
4. Report ROC AUC and early enrichment. RMSD applies only to the native AIH
   redock, not to LIT-PCBA assay actives.

The activity `label` column controls enrichment membership. It does not imply
that an assay active has a crystallographic reference pose.
"""
    (OUTPUT / "README.md").write_text(readme)

    print(json.dumps(audit["counts"], indent=2))
    print(f"Wrote benchmark files to {OUTPUT}")


if __name__ == "__main__":
    main()
