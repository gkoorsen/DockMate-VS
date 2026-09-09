"""Explicit variant policies and removal of property-based sampling budgets."""

import json
from types import SimpleNamespace

import pytest

from dockmate_vs.adaptive_docking import AdaptiveDockingPipeline, DockingResult
from dockmate_vs.gui.app import DockMateVSApp
from dockmate_vs.headless import load_campaign_config


@pytest.mark.parametrize("mode, count", [("all", 20), ("first", 1), ("best", 1), ("thorough", 15)])
def test_explicit_variant_selection(mode, count):
    pipeline = object.__new__(AdaptiveDockingPipeline)
    pipeline.ligand_variant_mode = mode
    variants = [{"label": f"ligand_p0_t{i}_c0", "energy": -i} for i in range(20)]
    selected = pipeline._select_ligand_variants(variants)
    assert len(selected) == count
    assert pipeline._select_ligand_variants([]) == []
    if mode == "all":
        assert selected == variants
    elif mode == "first":
        assert selected == variants[:1]
    else:
        assert selected[0] == variants[-1]


@pytest.mark.parametrize("mode", ["adaptive", "unknown"])
def test_removed_or_unknown_modes_rejected_before_pipeline_writes(tmp_path, mode):
    output = tmp_path / "run"
    with pytest.raises(ValueError, match="Unsupported ligand_variant_mode"):
        AdaptiveDockingPipeline(output, ligand_variant_mode=mode)
    assert not output.exists()


@pytest.mark.parametrize("mode", [None, "all", "best", "first", "thorough", "adaptive"])
def test_headless_variant_selection(tmp_path, mode):
    (tmp_path / "input.xlsx").touch()
    payload = {"input_file": "input.xlsx", "output_dir": "run", "single": {}}
    if mode is not None:
        payload["single"]["ligand_variant_mode"] = mode
    config_file = tmp_path / "config.json"
    config_file.write_text(json.dumps(payload))
    if mode == "adaptive":
        with pytest.raises(ValueError, match="selection has been removed"):
            load_campaign_config(config_file, "screen")
    else:
        assert load_campaign_config(config_file, "screen")["single"]["ligand_variant_mode"] == (mode or "all")


@pytest.mark.parametrize("mode, expected", [
    ("all_score", ("all", "score")), ("all_rmsd", ("all", "rmsd")),
    ("best", ("best", "rmsd")), ("thorough", ("thorough", "rmsd")),
])
def test_gui_explicit_variant_mapping(mode, expected):
    app = object.__new__(DockMateVSApp)
    app.variant_mode_var = SimpleNamespace(get=lambda: mode)
    assert app._variant_config() == expected


def test_gui_does_not_fall_back_to_smart_selection():
    app = object.__new__(DockMateVSApp)
    app.variant_mode_var = SimpleNamespace(get=lambda: "adaptive")
    with pytest.raises(ValueError, match="Unsupported ligand variant mode"):
        app._variant_config()


@pytest.mark.parametrize("worker, mode, key", [
    ("_run_worker", "screening", "single"),
    ("_run_worker", "adaptive", "adaptive"),
    ("_run_protocol_worker", "protocol_development", "single"),
])
def test_removed_mode_cannot_overwrite_campaign(tmp_path, worker, mode, key):
    app = object.__new__(DockMateVSApp)
    config = {"mode": mode, "output_dir": tmp_path / "run",
              key: {"ligand_variant_mode": "adaptive"}}
    with pytest.raises(ValueError, match="selection has been removed"):
        getattr(app, worker)([], config)
    assert not config["output_dir"].exists()


@pytest.mark.parametrize("mode, count", [("all", 3), ("first", 1), ("best", 1)])
def test_protocol_cascade_honours_explicit_variant_policy(monkeypatch, tmp_path, mode, count):
    monkeypatch.chdir(tmp_path)
    pipeline = AdaptiveDockingPipeline(tmp_path / "run", ligand_variant_mode=mode,
                                      max_tautomers=7, max_conformers=9, max_protocols=1)
    monkeypatch.setattr(pipeline, "_prepare_receptor", lambda **kw: (tmp_path / "rec.pdbqt", tmp_path / "rec.pdb"))
    preparations = []
    variants = [{"label": f"ligand_p0_t{i}_c0", "energy": -i,
                 "pdbqt": tmp_path / f"v{i}.pdbqt"} for i in range(3)]
    def prepare(**kwargs):
        preparations.append(kwargs)
        return variants
    monkeypatch.setattr(pipeline, "_prepare_ligand_variants", prepare)
    docked = []
    def dock(**kwargs):
        docked.append(kwargs["prepared_ligand_pdbqt"])
        return DockingResult(protocol_name="test", success=True, score=-7., rmsd=1.,
                             runtime_sec=0., num_poses=1, output_file=tmp_path / "pose.pdbqt")
    monkeypatch.setattr(pipeline, "_run_single_protocol", dock)
    pipeline.run_adaptive_docking(tmp_path / "rec.pdb", "CC", "ligand", "LIG", "A",
                                  crystal_ligand_pdb=tmp_path / "crystal.pdb")
    assert len(docked) == count
    assert preparations == [{"ligand_smiles": "CC", "ligand_name": "ligand", "enumerate_states": True}]
    assert (pipeline.max_tautomers, pipeline.max_conformers) == (7, 9)
