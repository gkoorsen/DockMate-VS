import json

import pandas as pd
import pytest

from dockmate_vs.gui.app import RedockResult
from dockmate_vs.headless import HeadlessDockMateRunner, _results_from_rows, regenerate_report
from dockmate_vs.utils.result_records import normalize_rmsd_record


def _result(**kwargs):
    values = dict(
        pdb_id="1ABC", ligand_resname="LIG", ligand_chain="A", mode="single",
        engine="smina", protocol="single", best_rmsd=1.2, top1_rmsd=1.2,
        success=True, runtime_sec=1.0, control_label=1, docking_completed=True,
    )
    values.update(kwargs)
    return RedockResult(**values)


def test_failed_rmsd_export_reload_report_and_resume(tmp_path):
    runner = HeadlessDockMateRunner()
    pose = tmp_path / "docked.pdbqt"
    pose.write_text("MODEL 1\nENDMDL\n")
    good = _result(case_id="good", output_file=str(pose))
    bad = _result(
        case_id="bad", output_file=str(pose), best_rmsd=999.9, top1_rmsd=999.9,
        success=False, protocol_attempts=[{"rmsd": 999.9, "success": False}],
    )
    results_path = tmp_path / "redock_results.json"
    csv_path = tmp_path / "redock_results.csv"
    progress_path = tmp_path / "redock_progress.json"
    runner._write_results(results_path, csv_path, [good, bad], 2.0)
    runner._write_progress(progress_path, [good, bad])

    for path in (results_path, csv_path, progress_path):
        assert "999.9" not in path.read_text()
    records = json.loads(results_path.read_text())["results"]
    assert [row["rmsd_status"] for row in records] == ["ok", "not_computable"]
    assert records[1]["best_rmsd"] is None
    assert "RMSD" in records[1]["error_message"]
    assert records[1]["protocol_attempts"][0]["rmsd"] is None
    reloaded = _results_from_rows(pd.read_csv(csv_path).to_dict("records"))
    assert reloaded[1].best_rmsd is None
    assert reloaded[1].docking_completed
    regenerate_report(results_path)
    summary = json.loads((tmp_path / "redock_summary.json").read_text())
    assert summary["mean_best_rmsd"] == pytest.approx(1.2)
    assert summary["success_rate_best"] == pytest.approx(50.0)
    assert summary["success_rate_top1"] == pytest.approx(50.0)
    assert summary["per_target_pose_recovery"][0]["cases"] == 2
    assert summary["by_engine"]["smina"]["count"] == 2

    manifest = {"resume_version": "test", "config": {"single": {"charge_handling": "preserve"}}, "cases": []}
    manifest_path = tmp_path / "run_manifest.json"
    manifest_path.write_text(json.dumps(manifest))
    resumed = runner._load_resumable_results(manifest_path, progress_path, manifest)
    assert {row.case_id for row in resumed} == {"good", "bad"}


def test_screening_without_reference_is_not_applicable():
    row = _result(
        mode="screening", control_label=1, best_rmsd=999.9, top1_rmsd=None,
        rmsd_status="not_applicable", success=False,
    ).to_record()
    assert row["rmsd_status"] == "not_applicable"
    assert row["best_rmsd"] is None
    assert row["error_message"] is None


@pytest.mark.parametrize("mode", ["single", "adaptive", "screening"])
def test_decoy_has_no_native_pose_in_any_workflow(mode):
    row = _result(mode=mode, control_label=0, best_rmsd=999.9, top1_rmsd=None).to_record()
    assert row["rmsd_status"] == "not_applicable"
    assert row["error_message"] is None
    assert row["best_rmsd"] is None


@pytest.mark.parametrize("value", [None, 999.9, float("nan"), float("inf")])
def test_legacy_failed_rmsd_and_null_have_same_summary(value):
    runner = HeadlessDockMateRunner()
    summary = runner._build_summary([_result(best_rmsd=value, top1_rmsd=None, success=False)], 2.0)
    assert summary["mean_best_rmsd"] is None
    assert summary["success_rate_best"] == 0.0
    assert summary["per_target_pose_recovery"][0]["success_rate_best"] == 0.0


def test_export_preserves_score_and_existing_error():
    row = normalize_rmsd_record({
        "best_rmsd": 999.9, "best_score": -8.0, "error_message": "Ligand atom mapping failed",
    })
    assert row["best_score"] == -8.0
    assert "Ligand atom mapping failed" in row["error_message"]
    assert row == normalize_rmsd_record(row)
