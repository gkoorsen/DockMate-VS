import json
import queue
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from dockmate_vs.adaptive_docking import AdaptiveDockingPipeline
from dockmate_vs.binding_site.cocrystal import BindingSiteDefinition
import dockmate_vs.gui.app as redock_analysis_module
from dockmate_vs.gui.app import DockMateVSApp, RedockResult
from dockmate_vs.utils.rmsd import coordinate_rmsd


def _app_without_tk() -> DockMateVSApp:
    return object.__new__(DockMateVSApp)


def test_executable_defaults_are_portable(monkeypatch, tmp_path):
    discovered = tmp_path / "smina"

    monkeypatch.setattr(
        redock_analysis_module.shutil,
        "which",
        lambda command: str(discovered) if command == "smina" else None,
    )

    assert redock_analysis_module._executable_default(None, "smina") == str(discovered)
    assert redock_analysis_module._executable_default(None, "vina") == "vina"
    assert redock_analysis_module._executable_default("~/bin/vina", "vina") == str(
        Path("~/bin/vina").expanduser()
    )


def test_rdock_root_default_prefers_environment_then_path(monkeypatch, tmp_path):
    configured = tmp_path / "configured-rdock"
    executable = tmp_path / "path-rdock" / "bin" / "rbdock"

    monkeypatch.setenv("RBT_ROOT", str(configured))
    monkeypatch.setattr(redock_analysis_module.shutil, "which", lambda _: str(executable))
    assert redock_analysis_module._rdock_root_default() == str(configured)

    monkeypatch.delenv("RBT_ROOT")
    assert redock_analysis_module._rdock_root_default() == str(executable.parent.parent)

    monkeypatch.setattr(redock_analysis_module.shutil, "which", lambda _: None)
    assert redock_analysis_module._rdock_root_default() == ""


def test_decoy_rows_expand_and_blank_decoy_rows_remain_samples(tmp_path: Path):
    template = tmp_path / "screen.xlsx"
    pd.DataFrame(
        [
            {
                "PDB_ID": "1ABC",
                "Ligand": "LIG",
                "Target_Ligand": "active",
                "SMILES": "CCO",
                "decoy compound": "d1,d2",
                "decoy SMILES": "CCC,CCN",
            },
            {
                "PDB_ID": "1ABC",
                "Ligand": "LIG",
                "Target_Ligand": "sample",
                "SMILES": "COC",
            },
        ]
    ).to_excel(template, index=False)

    pairs, _ = _app_without_tk()._load_pairs_from_excel(template)

    assert [(p["dock_name"], p["control_label"]) for p in pairs] == [
        ("active", 1),
        ("d1", 0),
        ("d2", 0),
        ("sample", None),
    ]
    assert pairs[0]["smiles"] == "CCO"
    assert pairs[-1]["smiles"] == "COC"


def test_screening_always_selects_ligand_variants_by_score():
    app = _app_without_tk()

    assert app._variant_selection_for_mode("screening", "rmsd") == "score"
    assert app._variant_selection_for_mode("screening", "score") == "score"
    assert app._variant_selection_for_mode("single", "rmsd") == "rmsd"


def test_filter_config_migrates_legacy_docking_platform_settings(
    monkeypatch, tmp_path: Path
):
    app = _app_without_tk()
    current = tmp_path / ".dockmate-vs" / "filters.json"
    legacy = tmp_path / ".docking_platform_gui" / "redock_filters.json"
    legacy.parent.mkdir(parents=True)
    legacy.write_text(json.dumps({"additives": ["GOL"], "cofactors": ["FAD"]}))

    monkeypatch.setattr(redock_analysis_module, "FILTERS_PATH", current)
    monkeypatch.setattr(redock_analysis_module, "LEGACY_FILTERS_PATH", legacy)
    monkeypatch.setattr(redock_analysis_module, "KNOWN_ADDITIVES", set())
    monkeypatch.setattr(redock_analysis_module, "COFACTORS", set())
    monkeypatch.setattr(redock_analysis_module, "ADDITIVES_ONLY", set())

    app._load_filter_config()

    assert json.loads(current.read_text()) == {"additives": ["GOL"], "cofactors": ["FAD"]}
    assert redock_analysis_module.KNOWN_ADDITIVES == {"GOL"}
    assert redock_analysis_module.COFACTORS == {"FAD"}


def test_pymol_overlay_preserves_receptor_frame_pose_coordinates(
    monkeypatch, tmp_path: Path
):
    app = _app_without_tk()
    receptor = tmp_path / "receptor_prepared.pdb"
    native = tmp_path / "crystal_ligand.pdb"
    receptor.write_text("END\n")
    native.write_text("END\n")
    poses = [object(), object()]
    reference = object()
    written = {}

    monkeypatch.setattr(app, "_ensure_receptor_pdb", lambda _case_dir: receptor)
    monkeypatch.setattr(
        app, "_load_poses_and_scores", lambda _path: (poses, [-9.0, -8.0])
    )
    monkeypatch.setattr(app, "_load_reference_mol", lambda _path: reference)
    monkeypatch.setattr(
        app, "_pose_rmsd", lambda _reference, pose: 4.0 if pose is poses[0] else 1.0
    )

    def _record_pose(molecule, _name, output_path, *args, **kwargs):
        output_path = Path(output_path)
        written[output_path.name] = molecule
        output_path.write_text("END\n")

    monkeypatch.setattr(app, "_write_ligand_pdb", _record_pose)
    case = {
        "output_file": tmp_path / "docked.pdbqt",
        "case_dir": tmp_path,
        "viewer_dir": tmp_path / "viewer",
        "crystal_ligand_pdb": native,
        "ligand": "AIH",
        "display_name": "AIH",
    }

    assert app._prepare_pymol_overlay(case) is not None
    assert written["best_score.pdb"] is poses[0]
    assert written["best_rmsd.pdb"] is poses[1]


def test_protocol_integer_lists_are_validated_and_deduplicated():
    assert DockMateVSApp._parse_positive_int_list("8, 16, 8", "Values") == [8, 16]

    with pytest.raises(ValueError, match="comma-separated"):
        DockMateVSApp._parse_positive_int_list("8, high", "Values")
    with pytest.raises(ValueError, match="greater than zero"):
        DockMateVSApp._parse_positive_int_list("8, 0", "Values")


def test_protocol_box_definitions_accept_margins_and_fixed_sizes():
    boxes = DockMateVSApp._parse_protocol_box_definitions(
        "margin:4; margin:6; 20x22x24; margin:4"
    )

    assert boxes == [
        {"label": "margin:4", "box_margin": 4.0, "size_override": None},
        {"label": "margin:6", "box_margin": 6.0, "size_override": None},
        {"label": "20x22x24", "box_margin": None, "size_override": (20.0, 22.0, 24.0)},
    ]

    with pytest.raises(ValueError, match="Use margin:4 or 20x20x20"):
        DockMateVSApp._parse_protocol_box_definitions("large")
    with pytest.raises(ValueError, match="greater than zero"):
        DockMateVSApp._parse_protocol_box_definitions("20x0x20")


def test_protocol_rescore_methods_are_validated_and_deduplicated():
    assert DockMateVSApp._parse_protocol_rescore_methods(
        "none, vina, vinardo, vina"
    ) == ["none", "vina", "vinardo"]

    with pytest.raises(ValueError, match="Unsupported rescoring"):
        DockMateVSApp._parse_protocol_rescore_methods("vina, imaginary")


def test_rescored_pose_metrics_follow_rescored_ranking(monkeypatch, tmp_path: Path):
    app = _app_without_tk()
    poses = [Chem.MolFromSmiles("C"), Chem.MolFromSmiles("CC"), Chem.MolFromSmiles("CCC")]
    rmsd_by_pose = {id(poses[0]): 0.8, id(poses[1]): 3.0, id(poses[2]): 1.5}
    monkeypatch.setattr(app, "_load_poses_and_scores", lambda _path: (poses, [None] * 3))
    monkeypatch.setattr(app, "_load_reference_mol", lambda _path: Chem.MolFromSmiles("C"))
    monkeypatch.setattr(app, "_pose_rmsd", lambda _reference, pose: rmsd_by_pose[id(pose)])

    metrics = app._compute_rescored_pose_metrics(
        tmp_path / "reference.pdb", tmp_path / "poses.pdbqt",
        scores=[-7.0, -9.0, -8.0], threshold=2.0,
        has_reference_pose=True,
    )

    assert metrics["rescore_score"] == -9.0
    assert metrics["rescore_top1_rmsd"] == 3.0
    assert metrics["rescore_top5_rmsd"] == 0.8
    assert metrics["rescore_best_rmsd_rank"] == 3
    assert metrics["rescore_rmsd_best_score"] == 3.0


def test_rdock_fallback_cavity_reference_is_translated_to_site_center(tmp_path: Path):
    molecule = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    AllChem.EmbedMolecule(molecule, randomSeed=7)
    source = tmp_path / "input.sdf"
    destination = tmp_path / "cavity_reference.sdf"
    writer = Chem.SDWriter(str(source))
    writer.write(molecule)
    writer.close()

    target = (12.5, -3.25, 8.75)
    AdaptiveDockingPipeline._translate_sdf_to_center(source, destination, target)

    translated = Chem.SDMolSupplier(
        str(destination), removeHs=False, sanitize=False
    )[0]
    conformer = translated.GetConformer()
    centroid = np.mean([
        list(conformer.GetAtomPosition(index))
        for index in range(translated.GetNumAtoms())
    ], axis=0)
    assert centroid == pytest.approx(target, abs=1e-3)


def test_protocol_report_groups_engine_and_box_definition(tmp_path: Path):
    rows = [
        {
            "status": "complete", "engine": "smina", "box_definition": "margin:4",
            "water_handling": "remove_all", "exhaustiveness": 8,
            "best_rmsd": 1.2, "success": True, "runtime_sec": 10.0,
        },
        {
            "status": "complete", "engine": "vina", "box_definition": "20x20x20",
            "water_handling": "retain_all", "exhaustiveness": 16,
            "best_rmsd": 2.4, "success": False, "runtime_sec": 20.0,
        },
    ]

    report = DockMateVSApp._write_protocol_report(rows, tmp_path).read_text()

    assert "| Engine | Box | Rescoring | Water |" in report
    assert "| smina | margin:4 | none | remove_all | 8 |" in report
    assert "| vina | 20x20x20 | none | retain_all | 16 |" in report
    assert "## Pose recovery across conditions" in report
    assert "## Recommended protocol per target" in report


def test_protocol_report_marks_disabled_rescoring_as_not_applicable(tmp_path: Path):
    rows = [{
        "status": "complete", "pdb_id": "1ABC", "ligand_resname": "LIG",
        "engine": "smina", "box_definition": "margin:4",
        "rescore_method": "none", "water_handling": "remove_all",
        "exhaustiveness": 8, "seed": 42, "best_rmsd": 1.0,
        "top1_rmsd": 1.2, "top5_rmsd": 1.0, "top10_rmsd": 1.0,
        "best_rmsd_rank": 2, "best_score": -8.0,
        "score_rmsd_spearman": 0.4, "success": True, "runtime_sec": 10.0,
        # Legacy files copied these baseline values into rescoring columns.
        "rescore_top1_rmsd": 1.2, "rescore_top5_rmsd": 1.0,
        "rescore_top10_rmsd": 1.0, "rescore_best_rmsd_rank": 2,
        "rescore_score": -8.0, "rescore_score_rmsd_spearman": 0.4,
    }]

    report = DockMateVSApp._write_protocol_report(rows, tmp_path).read_text()
    comparison_row = next(
        line for line in report.splitlines()
        if line.startswith("| smina | margin:4 | none |")
    )

    assert "| 1.20/1.00/1.00 | N/A/N/A/N/A |" in comparison_row
    assert "| 2.00 / N/A | -8.00 / N/A | 0.40 / N/A |" in comparison_row
    assert "| remove_all | none | 1 | 100.0% | 100.0% | N/A | 100.0% | N/A |" in report


def test_protocol_conditions_do_not_repeat_rdock_for_exhaustiveness():
    actives = [{"pdb_id": "1ABC", "site_ligand": "LIG"}]
    sweep = {
        "engines": ["smina", "rdock"],
        "box_definitions": [{"label": "margin:4"}],
        "water_modes": ["remove_all"],
        "exhaustiveness": [8, 16, 32],
        "seeds": [11, 42],
        "rescore_methods": ["none", "vinardo"],
    }

    conditions = DockMateVSApp._expand_protocol_conditions(actives, sweep)
    smina = [condition for condition in conditions if condition[1] == "smina"]
    rdock = [condition for condition in conditions if condition[1] == "rdock"]

    assert len(smina) == 12
    assert len(rdock) == 4
    assert {condition[5] for condition in rdock} == {None}


def test_protocol_report_collapses_legacy_rdock_exhaustiveness_rows(tmp_path: Path):
    rows = [
        {
            "status": "complete", "pdb_id": "1ABC", "ligand_resname": "LIG",
            "engine": "rdock", "box_definition": "margin:4",
            "rescore_method": "none", "water_handling": "remove_all",
            "exhaustiveness": exhaustiveness, "seed": 42,
            "best_rmsd": 1.8, "success": True, "runtime_sec": 10.0,
        }
        for exhaustiveness in (8, 16, 32)
    ]

    report = DockMateVSApp._write_protocol_report(rows, tmp_path).read_text()

    assert "- Unique crystal complexes evaluated: 1" in report
    assert "- Completed protocol-ranking conditions: 1" in report
    assert report.count("| rdock | margin:4 | none | remove_all | N/A |") == 1
    assert "| Conditions (with RMSD) |" in report


def test_protocol_recommendation_prefers_native_pose_rank_over_failed_top1_delta(
    tmp_path: Path,
):
    common = {
        "status": "complete", "pdb_id": "1XP1", "ligand_resname": "AIH",
        "target_name": "ESR1", "engine": "vina", "box_definition": "margin:6",
        "rescore_method": "none", "exhaustiveness": 8, "seed": 42,
        "success": True, "runtime_sec": 10.0, "top5_rmsd": 1.9,
    }
    rows = [
        {
            **common, "water_handling": "remove_all", "best_rmsd": 1.87,
            "top1_rmsd": 4.851, "best_rmsd_rank": 5,
        },
        {
            **common, "water_handling": "retain_all", "best_rmsd": 1.88,
            "top1_rmsd": 4.859, "best_rmsd_rank": 2,
        },
    ]

    report = DockMateVSApp._write_protocol_report(rows, tmp_path).read_text()
    recommendation = report.split("## Recommended protocol per target", 1)[1]

    assert "| ESR1 | vina | margin:6 | retain_all | 8 | baseline |" in recommendation


def test_protocol_markdown_parser_extracts_renderable_table():
    report = """# Protocol Development Summary

- Unsupported covalent conditions skipped: 3

Comparison explanation.

| Engine | Water | Baseline / rescored rank |
| --- | --- | ---: |
| smina | remove_all | 4.67 / 2.00 |
| rdock | selective | 10.50 / 12.50 |
"""

    prose, headers, rows = DockMateVSApp._parse_protocol_markdown(report)

    assert prose == [
        "Protocol Development Summary",
        "- Unsupported covalent conditions skipped: 3",
        "Comparison explanation.",
    ]
    assert headers == ["Engine", "Water", "Baseline / rescored rank"]
    assert rows == [
        ["smina", "remove_all", "4.67 / 2.00"],
        ["rdock", "selective", "10.50 / 12.50"],
    ]


def test_protocol_section_parser_keeps_each_summary_table_separate():
    report = """# Protocol Development Summary

| Engine | Water |
| --- | --- |
| smina | remove_all |

## Overall pose recovery

| Water | Top-1 |
| --- | ---: |
| remove_all | 75.0% |

## Recommended protocol per target

| Target | Ranking |
| --- | --- |
| NQO2 | rescored |
"""

    prose, tables = DockMateVSApp._parse_protocol_markdown_sections(report)

    assert prose[0] == "Protocol Development Summary"
    assert [table[0] for table in tables] == [
        "Protocol comparison", "Overall pose recovery", "Recommended protocol per target"
    ]
    assert tables[2][2] == [["NQO2", "rescored"]]


def test_protocol_chart_data_uses_completed_csv_rows_and_excludes_sentinels():
    frame = pd.DataFrame([
        {
            "status": "complete", "target_name": "Mpro", "engine": "smina",
            "box_definition": "margin:4", "rescore_method": "vinardo",
            "water_handling": "remove_all", "exhaustiveness": 8, "seed": 42,
            "best_rmsd": 1.0, "top1_rmsd": 3.0, "top5_rmsd": 1.0,
            "rescore_top1_rmsd": 1.5, "rescore_top5_rmsd": 1.0,
            "runtime_sec": 20.0,
        },
        {
            "status": "complete", "target_name": "Mpro", "engine": "smina",
            "box_definition": "margin:4", "rescore_method": "vinardo",
            "water_handling": "remove_all", "exhaustiveness": 8, "seed": 42,
            "best_rmsd": 999.9, "top1_rmsd": 999.9, "top5_rmsd": 999.9,
            "rescore_top1_rmsd": 999.9, "rescore_top5_rmsd": 999.9,
            "runtime_sec": 30.0,
        },
        {
            "status": "failed", "target_name": "Mpro", "engine": "smina",
            "box_definition": "margin:4", "rescore_method": "vinardo",
            "water_handling": "remove_all", "exhaustiveness": 8, "seed": 42,
            "best_rmsd": 0.1, "top1_rmsd": 0.1, "top5_rmsd": 0.1,
            "rescore_top1_rmsd": 0.1, "rescore_top5_rmsd": 0.1,
            "runtime_sec": 1.0,
        },
    ])

    data = DockMateVSApp._protocol_chart_data(frame)

    label = "smina | margin:4 | remove all | e8 | s42 | rescore:vinardo"
    assert data["pose_ranking_points"] == [{
        "engine": "smina", "rescore": True, "label": label,
        "x": 1.0, "y": 1.5,
        "tooltip": (
            f"{label}\nBest pose: 1.00 A; selected Top-1: 1.50 A"
        ),
    }]
    assert data["rescore_points"][0]["x"] == 3.0
    assert data["rescore_points"][0]["y"] == 1.5
    assert data["rescore_points"][0]["delta"] == -1.5
    assert data["runtime_points"][0]["x"] == 25.0
    assert data["top1_success"] == 1
    assert data["top5_success"] == 1
    assert data["best_pose_success"] == 1
    assert data["rescore_improved"] == 1
    assert data["rescore_unchanged"] == 0
    assert data["rescore_worse"] == 0


def test_protocol_chart_data_treats_seeds_as_replicates():
    frame = pd.DataFrame([
        {
            "status": "complete", "engine": "smina", "box_definition": "margin:4",
            "rescore_method": "none", "water_handling": "retain_all",
            "exhaustiveness": 8, "seed": seed, "best_rmsd": best,
            "top1_rmsd": top1, "top5_rmsd": best, "runtime_sec": 10.0,
        }
        for seed, best, top1 in ((1, 1.0, 1.5), (2, 1.2, 2.5), (3, 1.4, 3.5))
    ])

    data = DockMateVSApp._protocol_chart_data(frame)

    assert data["condition_count"] == 1
    assert "3 seeds" in data["pose_ranking_points"][0]["label"]
    assert data["pose_ranking_points"][0]["x"] == 1.2
    assert data["pose_ranking_points"][0]["y"] == 2.5
    assert data["top1_success"] == 0
    assert data["best_pose_success"] == 1


def test_protocol_chart_data_uses_every_rescored_condition():
    frame = pd.DataFrame([
        {
            "status": "complete", "pdb_id": "1ABC", "ligand_resname": "LIG",
            "engine": "smina", "box_definition": f"box:{index}",
            "rescore_method": "vinardo", "water_handling": "retain_all",
            "exhaustiveness": 8, "seed": 42, "best_rmsd": 1.0,
            "top1_rmsd": 3.0 + index / 10, "top5_rmsd": 1.0,
            "rescore_top1_rmsd": 2.5 + index / 10,
            "rescore_top5_rmsd": 1.0, "runtime_sec": 10.0 + index,
        }
        for index in range(20)
    ])

    data = DockMateVSApp._protocol_chart_data(frame)

    assert data["condition_count"] == 20
    assert len(data["pose_ranking_points"]) == 20
    assert len(data["rescore_points"]) == 20


def test_protocol_chart_data_switches_all_metrics_to_selected_top_n():
    frame = pd.DataFrame([{
        "status": "complete", "pdb_id": "1ABC", "ligand_resname": "LIG",
        "engine": "smina", "box_definition": "margin:4",
        "rescore_method": "vinardo", "water_handling": "remove_all",
        "exhaustiveness": 8, "seed": 42, "best_rmsd": 0.7,
        "top1_rmsd": 4.0, "top5_rmsd": 2.5, "top10_rmsd": 1.1,
        "rescore_top1_rmsd": 5.0, "rescore_top5_rmsd": 1.8,
        "rescore_top10_rmsd": 0.8, "runtime_sec": 20.0,
    }])

    top1 = DockMateVSApp._protocol_chart_data(frame, top_n=1)
    top5 = DockMateVSApp._protocol_chart_data(frame, top_n=5)
    top10 = DockMateVSApp._protocol_chart_data(frame, top_n=10)

    assert top1["pose_ranking_points"][0]["y"] == 5.0
    assert top5["pose_ranking_points"][0]["y"] == 1.8
    assert top10["pose_ranking_points"][0]["y"] == 0.8
    assert top10["rescore_points"][0]["x"] == 1.1
    assert top10["rescore_points"][0]["y"] == 0.8
    assert top10["runtime_points"][0]["y"] == 0.8
    assert {row["median"] for row in top10["factor_effects"]} == {0.8}
    assert top10["selected_success"] == 1
    assert top10["top1_success"] == 0
    assert top10["top5_success"] == 1
    assert top10["top10_success"] == 1

    with pytest.raises(ValueError, match="Top-1, Top-5, or Top-10"):
        DockMateVSApp._protocol_chart_data(frame, top_n=20)


def test_protocol_chart_data_collapses_legacy_rdock_exhaustiveness_duplicates():
    frame = pd.DataFrame([
        {
            "status": "complete", "pdb_id": "1ABC", "ligand_resname": "LIG",
            "engine": "rdock", "box_definition": "radius:6",
            "rescore_method": "none", "water_handling": "retain_all",
            "exhaustiveness": exhaustiveness, "seed": 42, "best_rmsd": 1.2,
            "top1_rmsd": 2.5, "top5_rmsd": 1.2, "runtime_sec": 50.0,
        }
        for exhaustiveness in (8, 16, 32)
    ])

    data = DockMateVSApp._protocol_chart_data(frame)

    assert data["source_rows"] == 3
    assert data["condition_count"] == 1
    assert len(data["pose_ranking_points"]) == 1
    engine_effect = next(
        row for row in data["factor_effects"] if row["label"] == "Engine: rdock"
    )
    assert engine_effect["n"] == 1


def test_screening_chart_data_prioritizes_enrichment_and_cross_structure_hits():
    summary = {
        "screening_score_direction": "lower",
        "per_structure_enrichment": [
            {
                "target_name": "Mpro", "pdb_id": "7GCO", "ligand": "LO0",
                "roc_auc": 0.93, "active_rank": 3, "decoys": 30,
                "score_margin": -0.22,
            },
            {
                "target_name": "Mpro", "pdb_id": "5REE", "ligand": "T1M",
                "roc_auc": 1.0, "active_rank": 1, "decoys": 30,
                "score_margin": 0.40,
            },
        ],
        "per_target_enrichment": [
            {
                "target_name": "Mpro", "structures": 2,
                "macro_roc_auc": 0.965, "pooled_roc_auc": 0.80,
            }
        ],
        "ef_1_percent": 4.0,
        "per_structure_screening": [
            {
                "target_name": "Mpro", "pdb_id": "7GCO", "ligand": "LO0",
                "cases": 20, "completed": 20, "completion_rate": 100.0,
                "scored": 20, "best_score": -10.0, "median_score": -7.0,
            }
        ],
        "screening_top_hits": [
            {"compound": "Compound A", "rank": 1},
            {"compound": "Compound A", "rank": 1},
            {"compound": "Compound B", "rank": 1},
            {"compound": "Compound A", "rank": 3},
        ],
    }

    data = DockMateVSApp._screening_chart_data(summary)

    assert data["structure_auc"][0][1] == 1.0
    assert data["score_margin"][0][1] == 0.40
    assert data["target_auc"] == [("Mpro (2 structures)", [0.965, 0.80])]
    assert data["top_hit_recurrence"][0] == ("Compound A (top-5 hits 3)", 2.0)
    assert data["score_advantage"] == [("Mpro | 7GCO/LO0", 3.0)]
    assert data["early_enrichment"] == [("EF 1%", 4.0)]


def test_mixed_screening_with_control_rmsd_uses_screening_charts():
    summary = {
        "n_samples": 40,
        "mean_best_rmsd": 1.2,
        "per_structure_enrichment": [{"pdb_id": "1ABC", "ligand": "LIG"}],
    }

    assert DockMateVSApp._is_screening_summary(summary) is True


def test_screening_chart_data_flags_incomplete_campaign_coverage():
    summary = {
        "total_cases": 100,
        "docking_completed": 80,
        "n_actives": 10,
        "n_decoys": 70,
        "control_actives": 8,
        "control_decoys": 52,
        "n_samples": 20,
        "screening_score_count": 15,
    }

    data = DockMateVSApp._screening_chart_data(summary)

    assert data["campaign_coverage"] == [
        ("Docking completed", 80.0),
        ("Controls with ranking scores", 75.0),
        ("Unknowns with ranking scores", 75.0),
    ]
    assert data["has_coverage_gap"] is True


def test_protocol_development_uses_unique_control_actives_only():
    pairs = [
        {"pdb_id": "1ABC", "site_ligand": "LIG", "chain": "A", "control_label": 1},
        {"pdb_id": "1ABC", "site_ligand": "LIG", "chain": "A", "control_label": 1},
        {"pdb_id": "1ABC", "site_ligand": "LIG", "chain": "A", "control_label": 0},
        {"pdb_id": "2DEF", "site_ligand": "DRG", "chain": None, "control_label": None},
        {"pdb_id": "3GHI", "site_ligand": "ACT", "chain": "B", "control_label": 1},
    ]

    actives = DockMateVSApp._protocol_active_pairs(pairs)

    assert [(item["pdb_id"], item["site_ligand"]) for item in actives] == [
        ("1ABC", "LIG"),
        ("3GHI", "ACT"),
    ]


def test_campaign_inputs_are_prefetched_before_offline_phase(monkeypatch, tmp_path):
    app = _app_without_tk()
    app._queue = queue.Queue()
    downloads = []
    monkeypatch.setattr(
        app, "_download_pdb",
        lambda pdb_id, _directory: downloads.append(pdb_id) or tmp_path / f"{pdb_id}.pdb"
    )
    monkeypatch.setattr(app, "_detect_ligand_chain", lambda _path, _ligand: "A")
    monkeypatch.setattr(app, "_get_ligand_smiles", lambda *_args: "CCO")
    pairs = [
        {"pdb_id": "1abc", "site_ligand": "LIG", "smiles": None, "chain": None},
        {"pdb_id": "1abc", "site_ligand": "LIG", "smiles": "CCC", "chain": None},
        {"pdb_id": "2def", "site_ligand": None, "smiles": "CCN", "chain": None},
    ]

    app._prefetch_remote_inputs(pairs, tmp_path)

    assert downloads == ["1ABC", "2DEF"]
    assert pairs[0]["smiles"] == "CCO"
    assert pairs[0]["chain"] == "A"
    assert pairs[0]["is_reference_ligand"] is True
    assert pairs[1]["is_reference_ligand"] is False
    assert pairs[2]["is_reference_ligand"] is False
    assert app._network_phase_complete is True


def test_activity_label_does_not_enable_reference_pose_metrics(monkeypatch, tmp_path):
    app = _app_without_tk()
    docked = tmp_path / "poses.pdbqt"
    docked.write_text("MODEL 1\nENDMDL\n")
    pose = Chem.MolFromSmiles("CC")
    monkeypatch.setattr(app, "_load_poses_and_scores", lambda _path: ([pose], [-8.0]))
    monkeypatch.setattr(
        app, "_load_reference_mol",
        lambda _path: pytest.fail("A score-only assay active must not load an RMSD reference"),
    )

    metrics = app._compute_pose_metrics(
        crystal_ligand_pdb=tmp_path / "native.pdb",
        docked_file=docked,
        threshold=2.0,
        has_reference_pose=False,
    )

    assert metrics == {"pose_count": 1, "best_score": -8.0}


def test_reference_identity_is_exact_not_analogue_similarity():
    app = _app_without_tk()

    assert app._same_ligand_smiles("c1ccccc1", "c1ccccc1") is True
    assert app._same_ligand_smiles("Cc1ccccc1", "c1ccccc1") is False
    assert app._same_ligand_smiles(None, "c1ccccc1") is None


def test_covalent_site_ligand_link_is_detected(tmp_path):
    pdb = tmp_path / "linked.pdb"
    pdb.write_text(
        "LINK         SG  CYS A  61                 C10 JMR A 301     1555   1555  1.68\n"
    )

    assert DockMateVSApp._has_covalent_ligand_link(pdb, "JMR", "A") is True
    assert DockMateVSApp._has_covalent_ligand_link(pdb, "JMR", "B") is False
    assert DockMateVSApp._has_covalent_ligand_link(pdb, "STI", "A") is False


def test_control_preserving_sample_is_balanced_and_keeps_sheet_order():
    pairs = []
    for pdb_id in ("1AAA", "2BBB", "3CCC", "4DDD"):
        pairs.append({"pdb_id": pdb_id, "ligand": "LIG", "control_label": 1})
        pairs.append({"pdb_id": pdb_id, "ligand": "LIG", "control_label": 0})
        pairs.append(
            {"pdb_id": pdb_id, "ligand": "LIG", "control_label": None, "dock_name": "s1"}
        )
        pairs.append(
            {"pdb_id": pdb_id, "ligand": "LIG", "control_label": None, "dock_name": "s2"}
        )

    sampled = _app_without_tk()._apply_random_sample(
        pairs, size=4, seed=42
    )
    non_controls = [p for p in sampled if p["control_label"] is None]

    assert len(sampled) == 12
    assert {p["pdb_id"] for p in non_controls} == {"1AAA", "2BBB", "3CCC", "4DDD"}
    assert [pairs.index(p) for p in sampled] == sorted(pairs.index(p) for p in sampled)


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        (1, 1),
        (0, 0),
        ("active", 1),
        ("positive", 1),
        ("decoy", 0),
        ("negative", 0),
        ("", None),
        (None, None),
        ("unknown", None),
        (2, None),
    ],
)
def test_control_label_parser(value, expected):
    assert _app_without_tk()._parse_control_label(value) == expected


def test_sampling_always_preserves_controls():
    pairs = [
        {"pdb_id": "1AAA", "ligand": "LIG", "control_label": 1},
        {"pdb_id": "1AAA", "ligand": "LIG", "control_label": 0},
        {"pdb_id": "1AAA", "ligand": "LIG", "control_label": None},
        {"pdb_id": "1AAA", "ligand": "LIG", "control_label": None},
    ]

    sampled = _app_without_tk()._apply_random_sample(
        pairs, size=1, seed=1
    )

    assert len(sampled) == 3
    assert pairs[0] in sampled
    assert pairs[1] in sampled
    assert sum(pair["control_label"] is None for pair in sampled) == 1


@pytest.mark.parametrize(
    ("scores", "expected_auc"),
    [
        ([(3.0, 1), (2.0, 0), (1.0, 0)], 1.0),
        ([(1.0, 1), (2.0, 0), (3.0, 0)], 0.0),
        ([(1.0, 1), (1.0, 0)], 0.5),
    ],
)
def test_enrichment_auc_handles_perfect_reversed_and_tied_scores(scores, expected_auc):
    metrics = _app_without_tk()._compute_enrichment_metrics(scores)

    assert metrics is not None
    assert metrics["roc_auc"] == pytest.approx(expected_auc)


def test_enrichment_requires_both_classes():
    app = _app_without_tk()

    assert app._compute_enrichment_metrics([]) is None
    assert app._compute_enrichment_metrics([(1.0, 1)]) is None
    assert app._compute_enrichment_metrics([(1.0, 0)]) is None


def test_enrichment_reports_tie_independent_average_precision():
    metrics = _app_without_tk()._compute_enrichment_metrics(
        [(2.0, 1), (2.0, 0), (1.0, 1), (0.0, 0)]
    )

    assert metrics is not None
    assert metrics["average_precision"] == pytest.approx(7.0 / 12.0)
    assert metrics["active_prevalence"] == pytest.approx(0.5)


def test_enrichment_factor_averages_a_tie_split_by_the_cutoff():
    metrics = _app_without_tk()._tie_aware_enrichment_factor(
        [(3.0, 1), (3.0, 0)] + [(2.0, 0)] * 8,
        percent=1.0,
    )

    assert metrics is not None
    assert metrics["selected"] == 1
    assert metrics["cutoff_tie_size"] == 2
    assert metrics["enrichment_factor"] == pytest.approx(5.0)
    assert metrics["enrichment_factor_min"] == pytest.approx(0.0)
    assert metrics["enrichment_factor_max"] == pytest.approx(10.0)


def _result(**overrides) -> RedockResult:
    values = {
        "pdb_id": "1ABC",
        "ligand_resname": "LIG",
        "ligand_chain": "A",
        "mode": "screening",
        "engine": "smina",
        "protocol": "single",
        "best_rmsd": 999.9,
        "success": False,
        "runtime_sec": 1.0,
    }
    values.update(overrides)
    return RedockResult(**values)


def test_rank_score_uses_rescore_precedence_and_correct_direction():
    app = _app_without_tk()

    assert app._rank_score_value(_result(best_score=-7.0)) == 7.0
    assert app._rank_score_value(_result(best_score=-7.0, rescore_score=-8.0)) == 8.0
    assert app._rank_score_value(
        _result(rescore_score=-8.0, rescore_cnn_score=0.7)
    ) == 0.7
    assert app._rank_score_value(
        _result(rescore_cnn_score=0.7, rescore_cnn_affinity=9.0)
    ) == 9.0
    assert app._rank_score_value(_result()) is None


def test_progress_file_distinguishes_docking_completion_from_rmsd_success(tmp_path: Path):
    path = tmp_path / "progress.json"
    results = [
        _result(best_score=-7.0, docking_completed=True),
        _result(error_message="failed", docking_completed=False),
    ]

    _app_without_tk()._write_progress(path, results)
    payload = json.loads(path.read_text())

    assert payload["results"][0]["docking_completed"] is True
    assert payload["results"][0]["success"] is False
    assert payload["results"][1]["docking_completed"] is False


def test_manifest_json_converter_handles_paths_and_numpy_values():
    app = _app_without_tk()
    payload = json.dumps(
        {"path": Path("run"), "array": np.array([1, 2]), "number": np.int64(3)},
        default=app._json_default,
    )

    assert json.loads(payload) == {"path": "run", "array": [1, 2], "number": 3}


def test_software_provenance_records_commit_dependencies_and_binary_versions(monkeypatch):
    def fake_command(command):
        if "rev-parse" in command:
            return "abc123\n"
        if "status" in command:
            return " M file.py\n"
        return "AutoDock Vina test build\n"

    monkeypatch.setattr(DockMateVSApp, "_command_version", staticmethod(fake_command))
    provenance = DockMateVSApp._software_provenance({
        "single": {
            "vina_binary": "/opt/vina",
            "smina_binary": "/opt/smina",
        }
    })

    assert provenance["git_commit"] == "abc123"
    assert provenance["git_dirty"] is True
    assert "numpy" in provenance["dependencies"]
    assert provenance["external_binaries"]["vina"]["path"] == "/opt/vina"
    assert provenance["external_binaries"]["vina"]["version"] == "AutoDock Vina test build"


def test_ligplot_and_dictionary_resolve_from_portable_environment(monkeypatch, tmp_path):
    root = tmp_path / "LigPlus"
    ligplot = root / "lib" / "exe_mac64" / "ligplot"
    linux_ligplot = root / "lib" / "exe_linux64" / "ligplot"
    components = root / "lib" / "data" / "components.cif"
    ligplus_jar = root / "LigPlus.jar"
    ligplot.parent.mkdir(parents=True)
    linux_ligplot.parent.mkdir(parents=True)
    components.parent.mkdir(parents=True)
    ligplot.write_text("#!/bin/sh\n")
    linux_ligplot.write_text("#!/bin/sh\n")
    components.write_text("data_components\n")
    ligplus_jar.write_bytes(b"jar")

    monkeypatch.setenv("LIGPLUS_ROOT", str(root))
    monkeypatch.delenv("LIGPLOT_BIN", raising=False)
    monkeypatch.delenv("LIGPLOT_HOME", raising=False)
    monkeypatch.delenv("HET_GROUP_DICTIONARY", raising=False)
    monkeypatch.setattr(redock_analysis_module.shutil, "which", lambda _: None)
    monkeypatch.setattr(redock_analysis_module.sys, "platform", "darwin")

    app = _app_without_tk()
    assert app._resolve_ligplot_bin() == str(ligplot.resolve())
    assert app._resolve_ligplus_jar(str(ligplot)) == ligplus_jar.resolve()
    assert app._resolve_components_cif(str(ligplot)) == components.resolve()


def test_java_resolver_and_ligplus_launcher(monkeypatch, tmp_path):
    java = tmp_path / "java"
    jar = tmp_path / "LigPlus.jar"
    drawing = tmp_path / "ligplot.drw"
    missing = tmp_path / "missing.drw"
    java.write_text("#!/bin/sh\n")
    jar.write_bytes(b"jar")
    drawing.write_text("drawing\n")
    monkeypatch.setenv("JAVA_BIN", str(java))

    calls = []

    def fake_popen(command, **kwargs):
        calls.append((command, kwargs))
        return object()

    monkeypatch.setattr(redock_analysis_module.subprocess, "Popen", fake_popen)

    assert DockMateVSApp._resolve_java_bin() == str(java.resolve())
    launched = DockMateVSApp._launch_ligplus_drawings(
        str(java.resolve()), jar, [("Native", drawing), ("Missing", missing)]
    )

    assert launched == 1
    assert calls[0][0] == [str(java.resolve()), "-jar", str(jar), str(drawing)]
    assert calls[0][1]["cwd"] == str(tmp_path)
    assert calls[0][1]["start_new_session"] is True


def test_viewer_receptor_removes_residual_site_ligand_and_connectivity(tmp_path):
    receptor = tmp_path / "receptor.pdb"
    output = tmp_path / "viewer" / "receptor_viewer.pdb"
    receptor.write_text(
        "\n".join(
            [
                "ATOM      1  CA  ALA A  10       0.000   0.000   0.000  1.00  0.00           C  ",
                "ATOM      2  C1  AIH L   1       1.000   0.000   0.000  1.00  0.00           C  ",
                "HETATM    3  O   HOH W   1       2.000   0.000   0.000  1.00  0.00           O  ",
                "CONECT    2    3",
                "CONECT    1    3",
                "END",
            ]
        )
        + "\n"
    )

    DockMateVSApp._prepare_viewer_receptor(receptor, "AIH", output)
    rendered = output.read_text()

    assert " AIH " not in rendered
    assert "HOH" in rendered
    assert "CONECT    2    3" not in rendered
    assert "CONECT    1    3" in rendered
    assert rendered.endswith("END\n")


def test_best_score_selection_reports_its_coordinate_rmsd(monkeypatch, tmp_path):
    app = _app_without_tk()
    poses = [object(), object()]
    reference = object()
    monkeypatch.setattr(app, "_load_poses_and_scores", lambda _path: (poses, [-9.0, -8.0]))
    monkeypatch.setattr(app, "_load_reference_mol", lambda _path: reference)
    monkeypatch.setattr(
        app, "_pose_rmsd", lambda _reference, pose: 3.25 if pose is poses[0] else 0.5
    )

    selected = app._select_best_pose(
        tmp_path / "native.pdb",
        tmp_path / "docked.pdbqt",
        selection_mode="best_score",
    )

    assert selected is not None
    assert selected[1] is poses[0]
    assert selected[2] == pytest.approx(3.25)
    assert selected[3] == pytest.approx(-9.0)


def test_ligplus_preparation_generates_native_and_distinct_selected_drawings(
    monkeypatch, tmp_path
):
    app = _app_without_tk()
    receptor = tmp_path / "receptor_prepared.pdb"
    receptor.write_text("END\n")
    reference = object()
    best_score_pose = object()
    best_rmsd_pose = object()

    monkeypatch.setattr(app, "_ensure_receptor_pdb", lambda _case_dir: receptor)

    def fake_select(_native, _docked, selection_mode="best_rmsd"):
        if selection_mode == "best_score":
            return reference, best_score_pose, 3.0, -9.0, 0, 2
        return reference, best_rmsd_pose, 0.5, -8.0, 1, 2

    monkeypatch.setattr(app, "_select_best_pose", fake_select)
    monkeypatch.setattr(
        app,
        "_prepare_ligplot_pair",
        lambda _case, ref, docked: (ref, docked, None, None, "AIH", "AIH"),
    )
    monkeypatch.setattr(
        app,
        "_write_ligand_pdb",
        lambda _mol, _name, path, **_kwargs: Path(path).write_text("END\n"),
    )
    monkeypatch.setattr(
        app,
        "_combine_complex",
        lambda _receptor, _ligand, path: Path(path).write_text("END\n"),
    )

    def fake_run(_bin, _complex, _resnum, _chain, out_dir):
        drawing = Path(out_dir) / "ligplot.drw"
        drawing.write_text("drawing\n")
        return None

    monkeypatch.setattr(app, "_run_ligplot", fake_run)
    case = {
        "case_dir": tmp_path,
        "viewer_dir": tmp_path / "viewer",
        "crystal_ligand_pdb": tmp_path / "native.pdb",
        "output_file": tmp_path / "docked.pdbqt",
        "ligand": "AIH",
    }

    drawings = app._prepare_ligplus_drawings(case, "/opt/ligplot")

    assert [label for label, _ in drawings] == ["Native", "Best score", "Best RMSD"]
    assert all(path.is_file() for _, path in drawings)


def test_ligplot_resolves_legacy_desktop_install(monkeypatch, tmp_path):
    home = tmp_path / "home"
    root = home / "Desktop" / "docking_platform" / "tools" / "LigPlus"
    ligplot = root / "lib" / "exe_mac64" / "ligplot"
    ligplot.parent.mkdir(parents=True)
    ligplot.write_text("#!/bin/sh\n")

    monkeypatch.delenv("LIGPLOT_BIN", raising=False)
    monkeypatch.delenv("LIGPLUS_ROOT", raising=False)
    monkeypatch.delenv("LIGPLOT_HOME", raising=False)
    monkeypatch.setattr(redock_analysis_module.shutil, "which", lambda _: None)
    monkeypatch.setattr(redock_analysis_module.sys, "platform", "darwin")
    monkeypatch.setattr(redock_analysis_module.Path, "home", staticmethod(lambda: home))

    app = _app_without_tk()
    assert app._resolve_ligplot_bin() == str(ligplot.resolve())


def test_ligplot_ligand_writer_replaces_duplicate_pdb_atom_names(tmp_path):
    mol = Chem.MolFromPDBBlock(
        "\n".join(
            [
                "HETATM    1  C   LIG L   1       0.000   0.000   0.000  1.00  0.00           C  ",
                "HETATM    2  C   LIG L   1       1.400   0.000   0.000  1.00  0.00           C  ",
                "HETATM    3  O   LIG L   1       2.600   0.000   0.000  1.00  0.00           O  ",
                "CONECT    1    2",
                "CONECT    2    3",
                "END",
            ]
        ),
        removeHs=False,
    )
    assert mol is not None

    output_pdb = tmp_path / "ligand.pdb"
    _app_without_tk()._write_ligand_pdb(mol, "LIG", output_pdb)

    atom_names = [
        line[12:16].strip()
        for line in output_pdb.read_text().splitlines()
        if line.startswith("HETATM")
    ]
    assert atom_names == ["C1", "C2", "O1"]
    conect = [
        line
        for line in output_pdb.read_text().splitlines()
        if line.startswith("CONECT")
    ]
    assert conect == ["CONECT    1    2", "CONECT    2    3"]


def test_ligplot_restores_dictionary_double_bonds_without_changing_labels(tmp_path):
    ps_file = tmp_path / "ligplot.ps"
    ps_file.write_text(
        "\n".join(
            [
                "/Ligbond_width   {    3.000 } def",
                "% Bond lines",
                "Ligbond_colour",
                "Ligbond_width W",
                " 10.00 20.00 30.00 20.00 L",
                " 30.00 20.00 40.00 30.00 L",
                "Nligbond_colour",
                " 50.00 50.00 60.00 60.00 L",
                "% Atoms",
                " 10.00 20.00 Ligatom_radius Sphere",
                "( C1 ) Ligatmname_size Center",
                " 30.00 20.00 Ligatom_radius Sphere",
                "( C2 ) Ligatmname_size Center",
                " 40.00 30.00 Ligatom_radius Sphere",
                "( O1 ) Ligatmname_size Center",
                "% Hydrophobic interactions",
                "(2.81) HBtext_size Print",
            ]
        )
        + "\n"
    )
    bond_orders = tmp_path / "hbadd.bonds"
    bond_orders.write_text(
        "    1  C1  LIG    1  L ->      2  C2  LIG    1  L   DOUB\n"
        "    2  C2  LIG    1  L ->      3  O1  LIG    1  L   SING\n"
    )

    restored = DockMateVSApp._restore_ligplot_double_bonds(
        ps_file,
        bond_orders,
        1,
        "L",
    )
    rendered = ps_file.read_text()

    assert restored == 1
    assert " 10.00 20.00 30.00 20.00 L" not in rendered
    assert " 10.00 22.25 30.00 22.25 L" in rendered
    assert " 10.00 17.75 30.00 17.75 L" in rendered
    assert " 30.00 20.00 40.00 30.00 L" in rendered
    assert "(2.81) HBtext_size Print" in rendered
    assert DockMateVSApp._restore_ligplot_double_bonds(
        ps_file,
        bond_orders,
        1,
        "L",
    ) == 0


def test_run_ligplot_uses_hbplus_and_ligplot_sequence(monkeypatch, tmp_path):
    root = tmp_path / "LigPlus"
    ligplot = root / "lib" / "exe_mac64" / "ligplot"
    hbadd = root / "lib" / "exe_mac64" / "hbadd"
    hbplus = root / "lib" / "exe_mac64" / "hbplus"
    prm = root / "lib" / "params" / "ligplot.prm"
    components = root / "lib" / "data" / "components.cif"
    for tool in (ligplot, hbadd, hbplus):
        tool.parent.mkdir(parents=True, exist_ok=True)
        tool.write_text("#!/bin/sh\n")
    prm.parent.mkdir(parents=True, exist_ok=True)
    prm.write_text("params\n")
    components.parent.mkdir(parents=True, exist_ok=True)
    components.write_text("data_components\n")

    complex_pdb = tmp_path / "complexes" / "complex.pdb"
    complex_pdb.parent.mkdir(parents=True, exist_ok=True)
    complex_pdb.write_text("END\n")
    out_dir = tmp_path / "ligplot_out"
    calls = []

    def fake_run(cmd, **kwargs):
        calls.append((cmd, kwargs))
        cwd = Path(kwargs["cwd"])
        if Path(cmd[0]).name == "ligplot":
            (cwd / "ligplot.ps").write_text("%!PS\n")
        return type("Result", (), {"returncode": 0, "stdout": "", "stderr": ""})()

    def fake_ps_to_png(self, ps_file, png_file):
        png_file.write_bytes(b"png")
        return True

    monkeypatch.setattr(redock_analysis_module.subprocess, "run", fake_run)
    monkeypatch.setattr(DockMateVSApp, "_ligplot_ps_to_png", fake_ps_to_png)

    app = _app_without_tk()
    png_file = app._run_ligplot(str(ligplot), complex_pdb, 1, "L", out_dir)

    assert png_file == out_dir / "ligplot.png"
    assert (out_dir / "ligplot.prm").exists()
    assert [Path(cmd[0]).name for cmd, _ in calls] == ["hbadd", "hbplus", "hbplus", "ligplot"]
    assert calls[0][0] == [str(hbadd), "complex.pdb", str(components.resolve())]
    assert calls[1][0] == [str(hbplus), "-L", "-h", "2.7", "-d", "3.35", "complex.pdb"]
    assert calls[2][0] == [str(hbplus), "-L", "-N", "complex.pdb"]
    assert calls[3][0] == [
        str(ligplot),
        "complex.pdb",
        "1",
        "1",
        "L",
        "-prm",
        "ligplot.prm",
        "-ctype",
        "1",
    ]


def test_pose_viewer_uses_protocol_development_csv_directly(tmp_path):
    protocol_csv = tmp_path / "protocol_development_results.csv"
    redock_json = tmp_path / "redock_results.json"

    assert DockMateVSApp._pose_results_csv(protocol_csv) == protocol_csv
    assert DockMateVSApp._pose_results_csv(redock_json) == tmp_path / "redock_results.csv"


def test_results_loader_resolves_screening_run_folder(tmp_path):
    run_dir = tmp_path / "screening_run"
    run_dir.mkdir()
    results = run_dir / "redock_results.json"
    results.write_text('{"results": []}')

    assert DockMateVSApp._result_file_for_selection(run_dir) == results


def test_results_loader_resolves_protocol_run_and_results_folders(tmp_path):
    run_dir = tmp_path / "protocol_run"
    protocol_dir = run_dir / "protocol_development"
    protocol_dir.mkdir(parents=True)
    results = protocol_dir / "protocol_development_results.csv"
    results.write_text("status\n")

    assert DockMateVSApp._result_file_for_selection(run_dir) == results
    assert DockMateVSApp._result_file_for_selection(protocol_dir) == results


def test_results_loader_accepts_exact_supported_file(tmp_path):
    results = tmp_path / "redock_results.csv"
    results.write_text("best_score\n")

    assert DockMateVSApp._result_file_for_selection(results) == results


def test_results_loader_rejects_parent_output_folder_and_unrelated_file(tmp_path):
    nested_run = tmp_path / "run_1"
    nested_run.mkdir()
    (nested_run / "redock_results.json").write_text('{"results": []}')
    unrelated = tmp_path / "notes.csv"
    unrelated.write_text("notes\n")

    assert DockMateVSApp._result_file_for_selection(tmp_path) is None
    assert DockMateVSApp._result_file_for_selection(unrelated) is None


def test_results_loader_prefers_current_workflow_when_both_exist(tmp_path):
    screening = tmp_path / "redock_results.json"
    screening.write_text('{"results": []}')
    protocol_dir = tmp_path / "protocol_development"
    protocol_dir.mkdir()
    protocol = protocol_dir / "protocol_development_results.csv"
    protocol.write_text("status\n")

    assert DockMateVSApp._result_file_for_selection(
        tmp_path, "screening"
    ) == screening
    assert DockMateVSApp._result_file_for_selection(
        tmp_path, "protocol_development"
    ) == protocol


def test_summary_counts_samples_separately_and_uses_explicit_controls():
    app = _app_without_tk()
    results = [
        _result(
            dock_name="active", control_label=1, best_score=-8.0,
            best_rmsd=1.0, success=True, docking_completed=True,
        ),
        _result(
            dock_name="decoy", control_label=0, best_score=-6.0,
            docking_completed=True,
        ),
        _result(
            dock_name="sample", control_label=None, best_score=-10.0,
            docking_completed=True,
        ),
    ]

    summary = app._build_summary(results, threshold=2.0)

    assert summary["n_actives"] == 1
    assert summary["n_decoys"] == 1
    assert summary["n_samples"] == 1
    assert summary["roc_auc"] == 1.0
    assert "1 actives, 1 decoys" in summary["interpretation"]["enrichment_message"]


def test_summary_reports_macro_and_pooled_auc_per_target():
    results = [
        _result(pdb_id="1AAA", ligand_resname="A01", target_name="NQO2",
                dock_name="active_a", control_label=1, best_score=-10.0,
                best_rmsd=1.0, top1_rmsd=3.0, top5_rmsd=1.0, top10_rmsd=1.0),
        _result(pdb_id="1AAA", ligand_resname="A01", target_name="NQO2",
                dock_name="decoy_a", control_label=0, best_score=-5.0),
        _result(pdb_id="2BBB", ligand_resname="B01", target_name="NQO2",
                dock_name="active_b", control_label=1, best_score=-8.0,
                best_rmsd=3.0, top1_rmsd=4.0, top5_rmsd=3.0, top10_rmsd=3.0),
        _result(pdb_id="2BBB", ligand_resname="B01", target_name="NQO2",
                dock_name="decoy_b", control_label=0, best_score=-9.0),
    ]

    summary = _app_without_tk()._build_summary(results, threshold=2.0)
    target = summary["per_target_enrichment"][0]

    assert target["target_name"] == "NQO2"
    assert target["structures"] == 2
    assert target["macro_roc_auc"] == pytest.approx(0.5)
    assert target["pooled_roc_auc"] == pytest.approx(0.75)
    pose = summary["per_target_pose_recovery"][0]
    assert pose["success_rate_best"] == pytest.approx(50.0)
    assert pose["success_rate_top1"] == pytest.approx(0.0)
    assert summary["screening_validation"] == "needs_review"
    markdown = _app_without_tk()._summary_to_markdown(summary)
    assert "Control Enrichment (per target)" in markdown
    assert "Pose Recovery (per target)" in markdown


def test_multi_active_assay_benchmark_is_not_treated_as_matched_decoys():
    results = [
        _result(
            pdb_id="1XP1", ligand_resname="AIH", target_name="ESR1",
            dock_name="active_1", control_label=1, best_score=-10.0,
            molecular_weight=500.0, logp=5.0, tpsa=20.0,
            rotatable_bonds=1, ligand_charge=0,
        ),
        _result(
            pdb_id="1XP1", ligand_resname="AIH", target_name="ESR1",
            dock_name="active_2", control_label=1, best_score=-8.0,
            molecular_weight=450.0, logp=4.0, tpsa=30.0,
            rotatable_bonds=2, ligand_charge=0,
        ),
        _result(
            pdb_id="1XP1", ligand_resname="AIH", target_name="ESR1",
            dock_name="inactive_1", control_label=0, best_score=-9.0,
            molecular_weight=150.0, logp=0.0, tpsa=100.0,
            rotatable_bonds=8, ligand_charge=1,
        ),
        _result(
            pdb_id="1XP1", ligand_resname="AIH", target_name="ESR1",
            dock_name="inactive_2", control_label=0, best_score=-7.0,
            molecular_weight=160.0, logp=0.0, tpsa=100.0,
            rotatable_bonds=8, ligand_charge=1,
        ),
    ]

    app = _app_without_tk()
    summary = app._build_summary(results, threshold=2.0)
    markdown = app._summary_to_markdown(summary)

    assert summary["enrichment_dataset_type"] == "assay_benchmark"
    assert summary["negative_class_label"] == "inactives"
    assert summary["screening_validation"] == "benchmark_result"
    assert summary["control_property_match_passed"] is None
    assert summary["control_property_diagnostics"] == []
    assert summary["per_structure_enrichment"][0]["best_active"] == "active_1"
    assert "Assay actives: 2" in markdown
    assert "Assay inactives: 2" in markdown
    assert "BENCHMARK RESULT" in markdown
    assert "Best assay active ranks first" in markdown
    assert "Best inactive" in markdown
    assert "INVALID DECOY MATCHING" not in markdown
    assert "Crystal ligand ranks first" not in markdown
    chart_data = summary["assay_benchmark_charts"]
    assert chart_data["actives"] == 2
    assert chart_data["inactives"] == 2
    assert chart_data["roc_curve"][0] == [0.0, 0.0]
    assert chart_data["roc_curve"][-1] == [1.0, 1.0]
    assert sum(chart_data["score_histogram"]["active_percent"]) == pytest.approx(100.0)
    assert sum(chart_data["score_histogram"]["inactive_percent"]) == pytest.approx(100.0)
    assert app._screening_chart_data(summary)["assay_benchmark"] == chart_data


def test_assay_benchmark_chart_curves_advance_score_ties_together():
    data = DockMateVSApp._assay_benchmark_chart_data(
        [(3.0, 1), (2.0, 1), (2.0, 0), (1.0, 0)], bins=2
    )

    assert data is not None
    assert data["roc_curve"] == [
        [0.0, 0.0], [0.0, 0.5], [0.5, 1.0], [1.0, 1.0]
    ]
    assert data["precision_recall_curve"] == [
        [0.0, 1.0], [0.5, 1.0], [1.0, pytest.approx(2.0 / 3.0)], [1.0, 0.5]
    ]
    assert data["cumulative_recovery_curve"] == [
        [0.0, 0.0], [0.25, 0.5], [0.75, 1.0], [1.0, 1.0]
    ]
    assert data["prevalence"] == pytest.approx(0.5)


def test_score_only_screening_summary_rebuilds_from_results_json(tmp_path: Path):
    app = _app_without_tk()
    results_path = tmp_path / "redock_results.json"
    summary_path = tmp_path / "redock_summary.json"
    results = [
        asdict(
            _result(
                pdb_id="1AAA",
                ligand_resname="L1",
                target_name="Mpro",
                dock_name="hit_1",
                best_score=-7.0,
                rescore_method="smina_score_only:vinardo",
                rescore_score=-8.0,
                pose_count=10,
                runtime_sec=10.0,
                docking_completed=True,
            )
        ),
        asdict(
            _result(
                pdb_id="1AAA",
                ligand_resname="L1",
                target_name="Mpro",
                dock_name="hit_2",
                best_score=-6.5,
                rescore_method="smina_score_only:vinardo",
                rescore_score=-6.0,
                pose_count=12,
                runtime_sec=20.0,
                docking_completed=True,
            )
        ),
        asdict(
            _result(
                pdb_id="1AAA",
                ligand_resname="L1",
                target_name="Mpro",
                dock_name="failed_hit",
                runtime_sec=0.0,
                docking_completed=False,
                error_message="Docking failed",
            )
        ),
    ]
    results_path.write_text(json.dumps({"results": results}))
    summary_path.write_text(
        json.dumps(
            {
                "threshold": 2.0,
                "mean_runtime_sec": 0.0,
                "success_rate_best": 0.0,
            }
        )
    )

    summary = app._summary_for_display(results_path, summary_path)

    assert summary["docking_completed"] == 2
    assert summary["docking_failed"] == 1
    assert summary["mean_runtime_sec"] == pytest.approx(15.0)
    assert summary["success_rate_best"] is None
    assert summary["screening_score_count"] == 2
    assert summary["screening_unscored_count"] == 1
    assert summary["screening_score_methods"] == ["vinardo (Smina score-only)"]
    saved_summary = json.loads(summary_path.read_text())
    assert saved_summary["docking_completed"] == 2
    assert saved_summary["docking_failed"] == 1
    saved_markdown = (tmp_path / "redock_summary.md").read_text()
    assert "## Screening score summary" in saved_markdown
    assert "Scored samples: 2/3" in saved_markdown
    assert summary["per_structure_screening"] == [
        {
            "target_name": "Mpro",
            "pdb_id": "1AAA",
            "ligand": "L1",
            "cases": 3,
            "completed": 2,
            "completion_rate": pytest.approx(66.6666666667),
            "scored": 2,
            "score_source": "vinardo (Smina score-only)",
            "best_compound": "hit_1",
            "best_score": -8.0,
            "median_score": -7.0,
        }
    ]


def test_score_only_screening_recreates_missing_summary_files(tmp_path: Path):
    app = _app_without_tk()
    results_path = tmp_path / "redock_results.json"
    summary_path = tmp_path / "redock_summary.json"
    results_path.write_text(
        json.dumps(
            {
                "results": [
                    asdict(
                        _result(
                            dock_name="screening_hit",
                            best_score=-7.5,
                            docking_completed=True,
                        )
                    )
                ]
            }
        )
    )

    summary = app._summary_for_display(results_path, summary_path)

    assert summary["docking_completed"] == 1
    assert summary_path.exists()
    assert (tmp_path / "redock_summary.md").exists()


def test_score_only_screening_markdown_omits_pose_metrics_and_reports_hits():
    app = _app_without_tk()
    results = [
        _result(
            pdb_id="1AAA",
            ligand_resname="L1",
            target_name="Mpro",
            dock_name="hit_1",
            best_score=-7.0,
            rescore_method="smina_score_only:vinardo",
            rescore_score=-8.0,
            pose_count=10,
            runtime_sec=10.0,
            docking_completed=True,
        ),
        _result(
            pdb_id="1AAA",
            ligand_resname="L1",
            target_name="Mpro",
            dock_name="failed_hit",
            runtime_sec=0.0,
            docking_completed=False,
            error_message="Docking failed",
        ),
    ]

    summary = app._build_summary(results, threshold=2.0)
    markdown = app._summary_to_markdown(summary)

    assert "RMSD and pose-recovery metrics are omitted" in markdown
    assert "## Screening score summary" in markdown
    assert "## Top-ranked compounds per structure" in markdown
    assert "## Failed screening cases" in markdown
    assert "Success rate (best pose)" not in markdown
    assert "Mean best RMSD" not in markdown
    assert "N/A" not in markdown


def test_screening_chart_data_uses_per_structure_completion_and_score_coverage():
    summary = {
        "per_structure_screening": [
            {
                "pdb_id": "1AAA",
                "ligand": "L1",
                "cases": 4,
                "completed": 3,
                "completion_rate": 75.0,
                "scored": 2,
            },
            {
                "pdb_id": "2BBB",
                "ligand": "L2",
                "cases": 2,
                "completed": 2,
                "completion_rate": 100.0,
                "scored": 2,
            },
        ],
        "mean_pose_count": 17.5,
    }

    data = DockMateVSApp._screening_chart_data(summary)

    assert data["completion"] == [("1AAA/L1", 75.0), ("2BBB/L2", 100.0)]
    assert data["score_coverage"] == [("1AAA/L1", 50.0), ("2BBB/L2", 100.0)]
    assert data["failures"] == [("1AAA/L1", 1.0), ("2BBB/L2", 0.0)]


def test_strong_auc_passes_enrichment_without_requiring_active_rank_one():
    results = [
        _result(dock_name="active", control_label=1, best_score=-9.5),
        _result(dock_name="better_decoy", control_label=0, best_score=-10.0),
    ]
    results.extend(
        _result(dock_name=f"decoy_{index}", control_label=0, best_score=-float(index))
        for index in range(1, 10)
    )

    summary = _app_without_tk()._build_summary(results, threshold=2.0)
    markdown = _app_without_tk()._summary_to_markdown(summary)

    assert summary["per_structure_enrichment"][0]["roc_auc"] == pytest.approx(0.9)
    assert summary["screening_validation"] == "passed_enrichment"
    assert "PASSED (ENRICHMENT)" in markdown
    assert "does not by itself invalidate useful enrichment" in markdown


def test_unlabelled_redock_is_not_counted_as_screening_sample(tmp_path):
    output = tmp_path / "docked.pdbqt"
    output.write_text("ATOM\n")
    result = _result(
        mode="adaptive", control_label=None, output_file=str(output),
        best_rmsd=1.0, success=True, docking_completed=None,
    )

    summary = _app_without_tk()._build_summary([result], threshold=2.0)

    assert summary["n_samples"] == 0
    assert summary["docking_completed"] == 1


def test_property_matching_rejects_charge_and_size_mismatches():
    app = _app_without_tk()
    active = _result(
        molecular_weight=400.0, logp=2.0, tpsa=80.0,
        rotatable_bonds=5, ligand_charge=1,
    )
    matched = _result(
        molecular_weight=380.0, logp=2.5, tpsa=70.0,
        rotatable_bonds=4, ligand_charge=1,
    )
    too_small = _result(
        molecular_weight=190.0, logp=2.0, tpsa=80.0,
        rotatable_bonds=5, ligand_charge=1,
    )
    wrong_charge = _result(
        molecular_weight=380.0, logp=2.5, tpsa=70.0,
        rotatable_bonds=4, ligand_charge=0,
    )

    assert app._property_matched(active, matched) is True
    assert app._property_matched(active, too_small) is False
    assert app._property_matched(active, wrong_charge) is False


def test_coordinate_rmsd_does_not_superimpose_displaced_pose():
    reference = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    AllChem.EmbedMolecule(reference, randomSeed=1)
    pose = Chem.Mol(reference)
    conformer = pose.GetConformer()
    for atom_index in range(pose.GetNumAtoms()):
        point = conformer.GetAtomPosition(atom_index)
        conformer.SetAtomPosition(atom_index, (point.x + 6.0, point.y, point.z))

    assert coordinate_rmsd(reference, pose) == pytest.approx(6.0)


def test_binding_site_self_docking_keeps_receptor_coordinate_frame():
    crystal = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    docked = crystal + np.array([5.0, 0.0, 0.0])

    result = BindingSiteDefinition().validate_self_docking(crystal, docked)

    assert result.rmsd == pytest.approx(5.0)
    assert result.success is False


def test_coordinate_rmsd_handles_symmetric_atom_permutations():
    reference = Chem.MolFromSmiles("ClCCl")
    pose = Chem.RenumberAtoms(reference, [2, 1, 0])
    reference_conf = Chem.Conformer(reference.GetNumAtoms())
    pose_conf = Chem.Conformer(pose.GetNumAtoms())
    coordinates = [(-1.0, 0.0, 0.0), (0.0, 0.0, 0.0), (1.0, 0.0, 0.0)]
    for index, xyz in enumerate(coordinates):
        reference_conf.SetAtomPosition(index, xyz)
        pose_conf.SetAtomPosition(index, xyz)
    reference.AddConformer(reference_conf)
    pose.AddConformer(pose_conf)

    assert coordinate_rmsd(reference, pose, use_symmetry=True) == pytest.approx(0.0)
