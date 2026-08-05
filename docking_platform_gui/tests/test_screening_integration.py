import queue
import json
from dataclasses import asdict
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

import docking_platform_gui.gui.redock_analysis as redock_module
from docking_platform_gui.binding_site.cocrystal import BindingSite
from docking_platform_gui.docking.base import DockingPose, DockingResult
from docking_platform_gui.gui.redock_analysis import RedockAnalysisApp
from docking_platform_gui.gui.redock_analysis import RedockResult


class FakePipeline:
    variants = []
    adaptive_limit = 2
    receptor_site_ligand = None

    def __init__(self, *args, **kwargs):
        pass

    def _prepare_receptor(self, pdb_file, water_handling, site_ligand_resname=None):
        type(self).receptor_site_ligand = site_ligand_resname
        root = Path(pdb_file).parent
        receptor_pdbqt = root / "receptor.pdbqt"
        receptor_pdb = root / "receptor.pdb"
        receptor_pdbqt.write_text("RECEPTOR\n")
        receptor_pdb.write_text("RECEPTOR\n")
        return receptor_pdbqt, receptor_pdb

    def _contains_metal(self, smiles):
        return False

    def _prepare_ligand_variants(self, **kwargs):
        return self.variants

    def _adaptive_variant_selection(self, variants, **kwargs):
        return variants[:self.adaptive_limit]

    def _select_diverse_variants(self, variants, count):
        return variants[:count]


class FakeSminaEngine:
    scores = {}
    failed = set()

    def __init__(self, *args, **kwargs):
        pass

    def dock(self, ligand_path, output_path, **kwargs):
        label = Path(ligand_path).stem
        if label in self.failed:
            return DockingResult(ligand_id=label, success=False, error_message="failed")
        Path(output_path).write_text("ATOM\n")
        score = self.scores[label]
        return DockingResult(ligand_id=label, poses=[DockingPose(score=score)])


def _single_config(selection="score"):
    return {
        "engine": "smina",
        "rdock_root": ".",
        "smina_binary": "smina",
        "vina_binary": "vina",
        "ligand_variant_mode": "all",
        "variant_select_by": selection,
        "max_tautomers": 2,
        "max_conformers": 2,
        "n_cpus": 1,
        "water_handling": "remove_all",
        "size_override": None,
        "box_margin": 4.0,
        "exhaustiveness": 8,
        "num_modes": 9,
        "energy_range": 3.0,
        "cpu": 1,
        "seed": 42,
        "timeout": 30,
        "scoring": "vina",
        "apo_site_mode": "auto",
        "rdock_radius": None,
        "rdock_runs": 10,
        "rdock_seed": 42,
    }


def _mock_case(monkeypatch, tmp_path, scores, failed=()):
    app = object.__new__(RedockAnalysisApp)
    app._queue = queue.Queue()
    pdb_file = tmp_path / "1abc.pdb"
    pdb_file.write_text("HEADER\n")
    variants = []
    for label in scores:
        ligand = tmp_path / f"{label}.pdbqt"
        ligand.write_text("LIGAND\n")
        variants.append({"label": label, "pdbqt": ligand, "smiles": "CC"})
    FakePipeline.variants = variants
    FakeSminaEngine.scores = scores
    FakeSminaEngine.failed = set(failed)
    monkeypatch.setattr(redock_module, "AdaptiveDockingPipeline", FakePipeline)
    monkeypatch.setattr(redock_module, "SminaDockingEngine", FakeSminaEngine)
    site = BindingSite(
        center=np.zeros(3),
        size=np.ones(3) * 20,
        ligand_resname="LIG",
        ligand_chain="A",
        ligand_atoms=10,
        source_pdb=str(pdb_file),
    )
    monkeypatch.setattr(app, "_predict_binding_site", lambda **kwargs: (site, "mock"))

    # Deliberately make the weakest-scoring variant have the best RMSD. A
    # screening run must still select by score.
    rmsds = {label: float(i + 1) for i, label in enumerate(reversed(list(scores)))}
    monkeypatch.setattr(
        app,
        "_compute_pose_metrics",
        lambda docked_file, **kwargs: {
            "best_score": scores[docked_file.parent.name],
            "best_rmsd": rmsds[docked_file.parent.name],
            "pose_count": 3,
        },
    )
    return app, pdb_file


def test_screening_case_selects_best_score_not_best_rmsd(monkeypatch, tmp_path):
    scores = {"v1": -6.0, "v2": -9.0, "v3": -7.0}
    app, pdb_file = _mock_case(monkeypatch, tmp_path, scores)

    result = app._run_single_case(
        pdb_file=pdb_file,
        ligand_name="sample",
        ligand_chain="A",
        smiles="CC",
        case_dir=tmp_path / "case",
        threshold=2.0,
        single_cfg=_single_config("score"),
        ligand_resname="LIG",
        site_mode="prediction",
        run_mode="screening",
    )

    assert result.best_score == -9.0
    assert Path(result.output_file).parent.name == "v2"
    assert result.docking_completed is True
    assert result.variants_prepared == 3
    assert result.variants_docked == 3


def test_adaptive_variant_mode_docks_only_selected_subset(monkeypatch, tmp_path):
    scores = {f"v{i}": -float(i) for i in range(1, 7)}
    app, pdb_file = _mock_case(monkeypatch, tmp_path, scores)
    config = _single_config("score")
    config["ligand_variant_mode"] = "adaptive"

    result = app._run_single_case(
        pdb_file, "sample", "A", "CC", tmp_path / "case", 2.0,
        config, ligand_resname="LIG", site_mode="prediction", run_mode="screening",
    )

    assert result.variants_prepared == 6
    assert result.variants_docked == 2
    assert result.best_score == -2.0


def test_failed_variant_is_skipped_when_another_variant_succeeds(monkeypatch, tmp_path):
    scores = {"v1": -10.0, "v2": -8.0}
    app, pdb_file = _mock_case(monkeypatch, tmp_path, scores, failed={"v1"})

    result = app._run_single_case(
        pdb_file, "sample", "A", "CC", tmp_path / "case", 2.0,
        _single_config("score"), ligand_resname="LIG", site_mode="prediction",
        run_mode="screening",
    )

    assert result.best_score == -8.0
    assert Path(result.output_file).parent.name == "v2"


def test_all_failed_variants_raise_clear_error(monkeypatch, tmp_path):
    scores = {"v1": -10.0, "v2": -8.0}
    app, pdb_file = _mock_case(monkeypatch, tmp_path, scores, failed=scores)

    with pytest.raises(ValueError, match="failed for all ligand variants"):
        app._run_single_case(
            pdb_file, "sample", "A", "CC", tmp_path / "case", 2.0,
            _single_config("score"), ligand_resname="LIG", site_mode="prediction",
            run_mode="screening",
        )


def test_rescoring_extracts_first_pose_from_multimodel_pdbqt(monkeypatch, tmp_path):
    app = object.__new__(RedockAnalysisApp)
    case_dir = tmp_path / "case"
    case_dir.mkdir()
    (case_dir / "receptor_prepared.pdbqt").write_text("RECEPTOR\n")
    docked = case_dir / "docked.pdbqt"
    docked.write_text(
        "MODEL 1\nREMARK pose-one\nATOM one\nENDMDL\n"
        "MODEL 2\nREMARK pose-two\nATOM two\nENDMDL\n"
    )
    seen = {}

    def fake_run(command, **kwargs):
        ligand = Path(command[command.index("--ligand") + 1])
        seen["ligand"] = ligand
        seen["content"] = ligand.read_text()
        return SimpleNamespace(returncode=0, stdout="Affinity: -8.25 (kcal/mol)\n", stderr="")

    monkeypatch.setattr(redock_module.subprocess, "run", fake_run)
    result = app._rescore_with_smina(docked, case_dir, "smina", "vinardo")

    assert result["score"] == -8.25
    assert result["method"] == "smina_score_only:vinardo"
    assert seen["ligand"].name == "rescore_pose_1.pdbqt"
    assert "pose-one" in seen["content"]
    assert "pose-two" not in seen["content"]
    assert "MODEL" not in seen["content"]


def test_rescoring_returns_subprocess_error_instead_of_raising(monkeypatch, tmp_path):
    app = object.__new__(RedockAnalysisApp)
    case_dir = tmp_path / "case"
    case_dir.mkdir()
    (case_dir / "receptor_prepared.pdbqt").write_text("RECEPTOR\n")
    docked = case_dir / "docked.pdbqt"
    docked.write_text("ATOM\n")
    monkeypatch.setattr(
        redock_module.subprocess,
        "run",
        lambda *args, **kwargs: SimpleNamespace(
            returncode=1, stdout="", stderr="unsupported scoring function"
        ),
    )

    result = app._rescore_with_smina(docked, case_dir, "smina", "not_supported")

    assert result == {"error": "unsupported scoring function"}


def _resume_manifest(tmp_path, scoring="vina"):
    return {
        "created_at": "ignored-for-compatibility",
        "config": {
            "mode": "screening",
            "output_dir": tmp_path,
            "single": {"variant_select_by": "score", "scoring": scoring},
        },
        "cases": [
            {
                "pdb_id": "1ABC",
                "site_ligand": "LIG",
                "dock_name": "sample",
                "control_label": None,
                "case_id": "1ABC_LIG_sample",
            }
        ],
    }


def test_compatible_progress_resumes_completed_case(tmp_path):
    app = object.__new__(RedockAnalysisApp)
    manifest_path = tmp_path / "run_manifest.json"
    progress_path = tmp_path / "redock_progress.json"
    output = tmp_path / "docked.pdbqt"
    output.write_text("ATOM\n")
    manifest = _resume_manifest(tmp_path)
    app._write_json_atomic(manifest_path, manifest)
    result = RedockResult(
        pdb_id="1ABC",
        ligand_resname="LIG",
        ligand_chain="A",
        mode="screening",
        engine="smina",
        protocol="single",
        best_rmsd=999.9,
        success=False,
        runtime_sec=1.0,
        output_file=str(output),
        dock_name="sample",
        docking_completed=True,
        case_id="1ABC_LIG_sample",
    )
    progress_path.write_text(json.dumps({"results": [asdict(result)]}))

    resumed = app._load_resumable_results(manifest_path, progress_path, manifest)

    assert [item.case_id for item in resumed] == ["1ABC_LIG_sample"]


def test_changed_scoring_configuration_disables_resume(tmp_path):
    app = object.__new__(RedockAnalysisApp)
    manifest_path = tmp_path / "run_manifest.json"
    progress_path = tmp_path / "redock_progress.json"
    app._write_json_atomic(manifest_path, _resume_manifest(tmp_path, scoring="vina"))
    progress_path.write_text('{"results": []}')

    resumed = app._load_resumable_results(
        manifest_path, progress_path, _resume_manifest(tmp_path, scoring="vinardo")
    )

    assert resumed == []


@pytest.mark.parametrize("completed,has_output", [(False, True), (True, False)])
def test_failed_or_missing_output_is_retried(tmp_path, completed, has_output):
    app = object.__new__(RedockAnalysisApp)
    manifest_path = tmp_path / "run_manifest.json"
    progress_path = tmp_path / "redock_progress.json"
    output = tmp_path / "docked.pdbqt"
    if has_output:
        output.write_text("ATOM\n")
    manifest = _resume_manifest(tmp_path)
    app._write_json_atomic(manifest_path, manifest)
    result = RedockResult(
        pdb_id="1ABC",
        ligand_resname="LIG",
        ligand_chain="A",
        mode="screening",
        engine="smina",
        protocol="single",
        best_rmsd=999.9,
        success=False,
        runtime_sec=0.0,
        output_file=str(output),
        dock_name="sample",
        docking_completed=completed,
        case_id="1ABC_LIG_sample",
        error_message=None if completed else "failed",
    )
    progress_path.write_text(json.dumps({"results": [asdict(result)]}))

    assert app._load_resumable_results(manifest_path, progress_path, manifest) == []


def test_legacy_progress_infers_case_id(tmp_path):
    app = object.__new__(RedockAnalysisApp)
    manifest_path = tmp_path / "run_manifest.json"
    progress_path = tmp_path / "redock_progress.json"
    output = tmp_path / "docked.pdbqt"
    output.write_text("ATOM\n")
    manifest = _resume_manifest(tmp_path)
    app._write_json_atomic(manifest_path, manifest)
    result = RedockResult(
        pdb_id="1ABC",
        ligand_resname="LIG",
        ligand_chain="A",
        mode="screening",
        engine="smina",
        protocol="single",
        best_rmsd=999.9,
        success=False,
        runtime_sec=1.0,
        output_file=str(output),
        dock_name="sample",
        control_label=None,
        docking_completed=None,
    )
    progress_path.write_text(json.dumps({"results": [asdict(result)]}))

    resumed = app._load_resumable_results(manifest_path, progress_path, manifest)

    assert resumed[0].case_id == "1ABC_LIG_sample"


def test_atomic_progress_write_leaves_no_temporary_file(tmp_path):
    app = object.__new__(RedockAnalysisApp)
    path = tmp_path / "redock_progress.json"

    app._write_progress(path, [])

    assert json.loads(path.read_text()) == {"results": []}
    assert not (tmp_path / ".redock_progress.json.tmp").exists()


def test_worker_skips_compatible_completed_case(monkeypatch, tmp_path):
    app = object.__new__(RedockAnalysisApp)
    app._queue = queue.Queue()
    app.progress_dialog = None
    output = tmp_path / "case" / "docked.pdbqt"
    output.parent.mkdir()
    output.write_text("ATOM\n")
    pair = {
        "pdb_id": "1ABC",
        "ligand": "LIG",
        "site_ligand": "LIG",
        "dock_name": "sample",
        "control_label": None,
        "case_id": "1ABC_LIG_sample",
    }
    config = {
        "output_dir": tmp_path,
        "mode": "screening",
        "threshold": 2.0,
        "single": {"engine": "smina"},
        "rescore": {"enable": False},
    }
    manifest = {
        "config": config,
        "cases": [{
            "pdb_id": "1ABC",
            "site_ligand": "LIG",
            "dock_name": "sample",
            "control_label": None,
            "case_id": "1ABC_LIG_sample",
        }],
    }
    app._write_json_atomic(tmp_path / "run_manifest.json", manifest)
    completed = RedockResult(
        pdb_id="1ABC",
        ligand_resname="LIG",
        ligand_chain="A",
        mode="screening",
        engine="smina",
        protocol="single",
        best_rmsd=999.9,
        success=False,
        runtime_sec=1.0,
        output_file=str(output),
        dock_name="sample",
        docking_completed=True,
        case_id="1ABC_LIG_sample",
    )
    app._write_progress(tmp_path / "redock_progress.json", [completed])
    monkeypatch.setattr(
        app, "_download_pdb", lambda *args, **kwargs: pytest.fail("case was recalculated")
    )
    written = {}
    monkeypatch.setattr(
        app,
        "_write_results",
        lambda json_path, csv_path, results, threshold: written.update(results=results),
    )

    app._run_worker([pair], config)

    assert [result.case_id for result in written["results"]] == ["1ABC_LIG_sample"]
    messages = list(app._queue.queue)
    assert any(message[0] == "log" and "Skipping completed" in message[1] for message in messages)


def test_cancelled_worker_saves_partial_results_without_reporting_done(monkeypatch, tmp_path):
    app = object.__new__(RedockAnalysisApp)
    app._queue = queue.Queue()
    app.progress_dialog = SimpleNamespace(cancelled=True)
    pair = {
        "pdb_id": "1ABC",
        "ligand": "LIG",
        "site_ligand": "LIG",
        "dock_name": "sample",
        "control_label": None,
        "case_id": "1ABC_LIG_sample",
    }
    config = {
        "output_dir": tmp_path,
        "mode": "screening",
        "threshold": 2.0,
        "single": {"engine": "smina"},
        "rescore": {"enable": False},
    }
    monkeypatch.setattr(
        app, "_download_pdb", lambda *args, **kwargs: pytest.fail("cancelled case ran")
    )
    written = {}
    monkeypatch.setattr(
        app,
        "_write_results",
        lambda json_path, csv_path, results, threshold: written.update(results=results),
    )

    app._run_worker([pair], config)

    assert written["results"] == []
    message_types = [message[0] for message in app._queue.queue]
    assert "cancelled" in message_types
    assert "done" not in message_types
