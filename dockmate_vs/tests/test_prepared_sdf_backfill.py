import json
from dataclasses import asdict
from types import SimpleNamespace

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

import dockmate_vs.preparation.backfill as backfill_module
from dockmate_vs.config.schema import LigandPreparationConfig
from dockmate_vs.gui.app import DockMateVSApp, RedockResult
from dockmate_vs.preparation.backfill import (
    PreparedSdfBackfillTarget,
    backfill_prepared_ligand_sdfs,
    validate_prepared_variant,
)


def _pdbqt(x: float = 0.0) -> str:
    return (
        "ROOT\n"
        f"ATOM      1  C   UNL     1    {x:8.3f}{1.0:8.3f}{2.0:8.3f}  0.00  0.00    +0.000 C \n"
        "ENDROOT\n"
        "TORSDOF 0\n"
    )


def _prepared_ligand(pdbqt: str):
    molecule = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    assert AllChem.EmbedMolecule(molecule, randomSeed=7) == 0

    def _write_sdf(path):
        writer = Chem.SDWriter(str(path))
        try:
            writer.write(molecule)
        finally:
            writer.close()

    return SimpleNamespace(pdbqt=pdbqt, to_sdf_file=_write_sdf)


def test_backfill_writes_only_exact_saved_preparation_variants(monkeypatch, tmp_path):
    saved_pdbqt = tmp_path / "hit_v1.pdbqt"
    saved_pdbqt.write_text(_pdbqt())
    first_sdf = tmp_path / "case_a" / "hit_v1.sdf"
    second_sdf = tmp_path / "case_b" / "hit_v1.sdf"
    cache_calls = []

    class FakeCache:
        def __init__(self, _cache_dir):
            pass

        def get(self, **kwargs):
            cache_calls.append(kwargs)
            return [_prepared_ligand(_pdbqt())]

    monkeypatch.setattr(backfill_module, "LigandCache", FakeCache)
    targets = [
        PreparedSdfBackfillTarget(
            compound="hit",
            smiles="CCO",
            variant_index=1,
            pdbqt_path=saved_pdbqt,
            sdf_path=path,
            case_id=path.parent.name,
        )
        for path in (first_sdf, second_sdf)
    ]

    result = backfill_prepared_ligand_sdfs(
        targets,
        LigandPreparationConfig(max_tautomers=1, max_conformers=1),
        cache_dir=tmp_path / "cache",
    )

    assert len(cache_calls) == 1
    assert result.errors == []
    assert result.written == [first_sdf, second_sdf]
    assert all(Chem.SDMolSupplier(str(path), removeHs=False)[0] is not None for path in result.written)


def test_backfill_rejects_changed_preparation_coordinates(monkeypatch, tmp_path):
    saved_pdbqt = tmp_path / "hit_v1.pdbqt"
    saved_pdbqt.write_text(_pdbqt(0.0))
    output_sdf = tmp_path / "hit_v1.sdf"

    class FakeCache:
        def __init__(self, _cache_dir):
            pass

        def get(self, **_kwargs):
            return [_prepared_ligand(_pdbqt(0.5))]

    monkeypatch.setattr(backfill_module, "LigandCache", FakeCache)
    result = backfill_prepared_ligand_sdfs(
        [PreparedSdfBackfillTarget(
            compound="hit",
            smiles="CCO",
            variant_index=1,
            pdbqt_path=saved_pdbqt,
            sdf_path=output_sdf,
            case_id="1ABC_LIG_hit",
        )],
        LigandPreparationConfig(max_tautomers=1, max_conformers=1),
        cache_dir=tmp_path / "cache",
    )

    assert not output_sdf.exists()
    assert len(result.errors) == 1
    assert "coordinates differ" in result.errors[0]


def test_validate_prepared_variant_compares_all_three_coordinates(tmp_path):
    saved_pdbqt = tmp_path / "saved.pdbqt"
    saved_pdbqt.write_text(_pdbqt(0.0))

    valid, reason = validate_prepared_variant(_pdbqt(0.5), saved_pdbqt)

    assert valid is False
    assert "0.500 A" in reason


def test_results_backfill_plan_finds_transferred_workbook_and_legacy_policy(tmp_path):
    transfer = tmp_path / "transfer"
    run_root = transfer / "results" / "screening_run"
    templates = transfer / "templates"
    templates.mkdir(parents=True)
    run_root.mkdir(parents=True)
    workbook = templates / "screen.xlsx"
    pd.DataFrame([{
        "PDB_ID": "1ABC",
        "Ligand": "LIG",
        "Target_Ligand": "hit",
        "SMILES": "CCO",
    }]).to_excel(workbook, index=False)

    case_dir = run_root / "1ABC_LIG_hit"
    output_file = case_dir / "variants" / "hit_v1" / "docked.pdbqt"
    output_file.parent.mkdir(parents=True)
    output_file.write_text("MODEL 1\nENDMDL\n")
    prepared_pdbqt = case_dir / "ligand_variants" / "hit_v1.pdbqt"
    prepared_pdbqt.parent.mkdir(parents=True)
    prepared_pdbqt.write_text(_pdbqt())

    result = RedockResult(
        pdb_id="1ABC",
        ligand_resname="LIG",
        ligand_chain="",
        mode="screening",
        engine="smina",
        protocol="single",
        best_rmsd=999.9,
        success=False,
        runtime_sec=1.0,
        output_file=str(output_file),
        best_score=-7.5,
        dock_name="hit",
        control_label=None,
        case_id="1ABC_LIG_hit",
        docking_completed=True,
    )
    results_path = run_root / "redock_results.json"
    results_path.write_text(json.dumps({"results": [asdict(result)]}))
    (run_root / "run_manifest.json").write_text(json.dumps({
        "config": {
            "mode": "screening",
            "input_file": "/another/computer/screen.xlsx",
            "filters": {},
            "single": {"max_tautomers": 5, "max_conformers": 5, "n_cpus": 3},
        },
        "cases": [{
            "pdb_id": "1ABC",
            "site_ligand": "LIG",
            "dock_name": "hit",
            "control_label": None,
            "case_id": "1ABC_LIG_hit",
        }],
    }))

    app = object.__new__(DockMateVSApp)
    plan = app._prepared_sdf_backfill_plan(results_path)

    assert plan is not None
    assert plan["workbook"] == workbook.resolve()
    assert plan["errors"] == []
    assert plan["n_cpus"] == 3
    assert plan["config"].charge_handling == "neutralize"
    assert len(plan["targets"]) == 1
    assert plan["targets"][0].sdf_path == prepared_pdbqt.with_suffix(".sdf")


def test_results_selection_waits_for_automatic_backfill(monkeypatch, tmp_path):
    results_path = tmp_path / "redock_results.json"
    results_path.write_text('{"results": []}')
    app = object.__new__(DockMateVSApp)
    app.mode_var = SimpleNamespace(get=lambda: "screening")
    started = []
    displayed = []
    monkeypatch.setattr(
        app,
        "_start_automatic_prepared_sdf_backfill",
        lambda path: started.append(path) or True,
    )
    monkeypatch.setattr(app, "_display_results_selection", displayed.append)

    app._load_results_selection(tmp_path)

    assert started == [results_path]
    assert displayed == []


def test_results_selection_displays_immediately_when_no_backfill_is_needed(
    monkeypatch, tmp_path
):
    results_path = tmp_path / "redock_results.json"
    results_path.write_text('{"results": []}')
    app = object.__new__(DockMateVSApp)
    app.mode_var = SimpleNamespace(get=lambda: "screening")
    displayed = []
    monkeypatch.setattr(app, "_start_automatic_prepared_sdf_backfill", lambda _path: False)
    monkeypatch.setattr(app, "_display_results_selection", displayed.append)

    app._load_results_selection(tmp_path)

    assert displayed == [results_path]
