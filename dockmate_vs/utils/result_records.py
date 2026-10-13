"""Portable RMSD records, including compatibility with historical sentinels."""

import math


RMSD_FIELDS = (
    "best_rmsd", "top1_rmsd", "top5_rmsd", "top10_rmsd", "rmsd_best_score",
    "rmsd_mean", "rmsd_median", "rmsd_std", "rescore_top1_rmsd",
    "rescore_top5_rmsd", "rescore_top10_rmsd", "rescore_rmsd_best_score", "rmsd",
)


def rmsd_value(value):
    """Return an angstrom distance, or None for absent/legacy failed values."""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) and 0 <= number < 900 else None


def normalize_rmsd_record(record):
    """Copy a result for export without representing a failed RMSD as a distance.

    New screening records declare applicability explicitly. Historical screening
    controls lack this information; retain their native-control interpretation.
    """
    row = dict(record)
    for key in RMSD_FIELDS:
        if key in row:
            row[key] = rmsd_value(row[key])
    status = row.get("rmsd_status")
    label = str(row.get("control_label", "")).lower()
    best = row.get("best_rmsd", row.get("rmsd"))
    if status == "not_applicable" or label in {"0", "0.0", "decoy"}:
        status = "not_applicable"
        for key in RMSD_FIELDS:
            if key in row:
                row[key] = None
    elif best is not None:
        status = "ok"
    elif status != "not_computable":
        is_screen = str(row.get("mode", "")).lower() == "screening"
        is_active = label in {"1", "1.0", "active"}
        status = "not_applicable" if is_screen and not is_active else "not_computable"
    row["rmsd_status"] = status
    if status != "ok":
        row["success"] = False
    if status == "not_computable":
        reason = "Native-pose RMSD unavailable: no valid reference-to-pose atom mapping was obtained."
        error = row.get("error_message")
        if not isinstance(error, str) or not error.strip():
            row["error_message"] = reason
    if isinstance(row.get("protocol_attempts"), list):
        row["protocol_attempts"] = [
            normalize_rmsd_record(attempt) for attempt in row["protocol_attempts"]
        ]
    return row
