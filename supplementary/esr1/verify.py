#!/usr/bin/env python3
"""Verify supplement integrity and compare locally reconstructed analysis outputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_integrity(root: Path = ROOT) -> None:
    for line in (root / "SHA256SUMS.txt").read_text().splitlines():
        expected, name = line.split("  ", 1)
        path = root / name
        if not path.resolve().is_relative_to(root.resolve()):
            raise ValueError(f"Unexpected checksum path: {name}")
        if sha256(path) != expected:
            raise ValueError(f"Checksum mismatch: {name}")
    with (root / "selection.csv").open(newline="") as handle:
        selection = list(csv.DictReader(handle))
    with (root / "results/redock_results.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(selection) != 880 or len(rows) != 880:
        raise ValueError("Expected 880 selected compounds and screening result rows")
    if len({row["source_id"] for row in selection}) != 880:
        raise ValueError("Duplicate selected source identifiers")
    labels = [int(row["control_label"]) for row in rows]
    if labels.count(1) != 80 or labels.count(0) != 800:
        raise ValueError("Unexpected active/inactive counts")
    if {row["source_id"] for row in selection} != {row["dock_name"].rsplit("_", 1)[-1] for row in rows}:
        raise ValueError("Selection and result identifiers differ")
    projections = json.loads((root / "provenance/packaging.json").read_text())["result_projections"]
    for name, metadata in projections.items():
        with (root / name).open(newline="") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames != metadata["columns"] or len(list(reader)) != metadata["rows"]:
                raise ValueError(f"Result projection changed: {name}")
    print("Package checksums, 880 identifiers and 80/800 result labels verified.")


def verify_reconstruction(output: Path, root: Path = ROOT) -> None:
    source = json.loads((root / "source.json").read_text())
    if sha256(output / "screening_subset80_800_seed42_compounds.csv") != source["reconstructed_compounds_csv_sha256"]:
        raise ValueError("Reconstructed compound CSV checksum differs")
    import pandas as pd
    expected = pd.read_csv(output / "screening_subset80_800_seed42_compounds.csv")
    actual = pd.read_excel(output / "screening_subset80_800_seed42.xlsx", sheet_name="Docking_Jobs")
    pd.testing.assert_frame_equal(actual, expected)
    print("Reconstructed CSV is byte-identical; workbook job rows match.")


def verify_analysis(output: Path, root: Path = ROOT) -> None:
    import numpy as np
    import pandas as pd

    for filename in ("esr1_seed42_screening_metrics.csv", "esr1_seed42_property_diagnostics.csv"):
        expected = pd.read_csv(root / "publication/tables" / filename)
        actual = pd.read_csv(output / "tables" / filename)
        pd.testing.assert_frame_equal(actual, expected, check_exact=False, rtol=1e-12, atol=1e-12)
    filename = "esr1_seed42_property_model.json"
    expected = json.loads((root / "publication/tables" / filename).read_text())
    actual = json.loads((output / "tables" / filename).read_text())
    if set(actual) != set(expected):
        raise ValueError("Property-model fields differ")
    for key in expected:
        np.testing.assert_allclose(actual[key], expected[key], rtol=1e-12, atol=1e-12)
    print("Screening metrics, bootstrap intervals, property diagnostics and all ten model folds match.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reconstructed", type=Path)
    parser.add_argument("--reanalysis", type=Path)
    args = parser.parse_args()
    verify_integrity()
    if args.reconstructed:
        verify_reconstruction(args.reconstructed)
    if args.reanalysis:
        verify_analysis(args.reanalysis)
