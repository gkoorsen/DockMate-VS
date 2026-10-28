"""Verify the molecular-weight-controlled diagnostics (Figure 5 and Table 3).

Regenerates the diagnostics from this package's own results/redock_results.csv and
compares every value against the figures printed in the manuscript.

Usage (from the package root):
    python verification/verify_mw_controlled.py

Exit status 0 means every value matched at the printed precision.
"""
import json
import math
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
EXPECTED = json.load(open(ROOT / "verification/mw_controlled_expected.json"))
TOL = 5e-5  # values are printed to four decimals


def main():
    out = ROOT / "verification/_regenerated_mw_controlled.csv"
    subprocess.run([sys.executable, str(ROOT / "analysis/esr1_mw_controlled.py"),
                    str(ROOT / "results/redock_results.csv")],
                   cwd=out.parent, check=True, capture_output=True)
    (out.parent / "esr1_mw_controlled.csv").replace(out)
    got = pd.read_csv(out).set_index("analysis")
    packaged = pd.read_csv(ROOT / "publication/tables/esr1_mw_controlled.csv").set_index("analysis")

    failures = []
    if list(got.index) != list(EXPECTED):
        failures.append("regenerated row identities/order differ from manuscript expectations")
    if not got.index.is_unique:
        print("FAIL duplicate analysis rows in regenerated output")
        return 1
    for name, (auc, lo, hi) in EXPECTED.items():
        if name not in got.index:
            failures.append(f"{name}: absent from regenerated output")
            continue
        row = got.loc[name]
        for label, expected in (("auc", auc), ("ci_low", lo), ("ci_high", hi)):
            if expected is None:
                if not pd.isna(row[label]):
                    failures.append(f"{name}.{label}: expected a missing value, got {row[label]}")
                continue
            try:
                value = float(row[label])
            except (TypeError, ValueError):
                failures.append(f"{name}.{label}: expected a finite number, got {row[label]}")
                continue
            if not math.isfinite(value):
                failures.append(f"{name}.{label}: expected a finite number, got {value}")
            elif abs(value - expected) > TOL:
                failures.append(f"{name}.{label}: regenerated {value} vs manuscript {expected}")
    if not got.index.equals(packaged.index):
        failures.append("row order differs between regenerated and packaged tables")

    packaged_matches = got.round(4).equals(packaged.round(4))
    if not packaged_matches:
        failures.append("packaged table differs from regenerated table")
    print(f"rows checked: {len(EXPECTED)}")
    print("packaged table matches regenerated table:", packaged_matches)
    for f in failures:
        print("FAIL", f)
    print("RESULT:", "all values match the manuscript" if not failures else f"{len(failures)} mismatches")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
