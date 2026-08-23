import json
import queue
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from docking_platform_gui.adaptive_docking import AdaptiveDockingPipeline
from docking_platform_gui.gui.redock_analysis import RedockAnalysisApp, RedockResult
from docking_platform_gui.utils.rmsd import coordinate_rmsd


def _app_without_tk() -> RedockAnalysisApp:
    return object.__new__(RedockAnalysisApp)


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

    pairs, _ = _app_without_tk()._load_pairs_from_excel(
        template, use_smiles=True, include_controls=True
    )

    assert [(p["dock_name"], p["control_label"]) for p in pairs] == [
        ("active", 1),
        ("d1", 0),
        ("d2", 0),
        ("sample", None),
    ]


def test_screening_always_selects_ligand_variants_by_score():
    app = _app_without_tk()

    assert app._variant_selection_for_mode("screening", "rmsd") == "score"
    assert app._variant_selection_for_mode("screening", "score") == "score"
    assert app._variant_selection_for_mode("single", "rmsd") == "rmsd"


def test_protocol_integer_lists_are_validated_and_deduplicated():
    assert RedockAnalysisApp._parse_positive_int_list("8, 16, 8", "Values") == [8, 16]

    with pytest.raises(ValueError, match="comma-separated"):
        RedockAnalysisApp._parse_positive_int_list("8, high", "Values")
    with pytest.raises(ValueError, match="greater than zero"):
        RedockAnalysisApp._parse_positive_int_list("8, 0", "Values")


def test_protocol_box_definitions_accept_margins_and_fixed_sizes():
    boxes = RedockAnalysisApp._parse_protocol_box_definitions(
        "margin:4; margin:6; 20x22x24; margin:4"
    )

    assert boxes == [
        {"label": "margin:4", "box_margin": 4.0, "size_override": None},
        {"label": "margin:6", "box_margin": 6.0, "size_override": None},
        {"label": "20x22x24", "box_margin": None, "size_override": (20.0, 22.0, 24.0)},
    ]

    with pytest.raises(ValueError, match="Use margin:4 or 20x20x20"):
        RedockAnalysisApp._parse_protocol_box_definitions("large")
    with pytest.raises(ValueError, match="greater than zero"):
        RedockAnalysisApp._parse_protocol_box_definitions("20x0x20")


def test_protocol_rescore_methods_are_validated_and_deduplicated():
    assert RedockAnalysisApp._parse_protocol_rescore_methods(
        "none, vina, vinardo, vina"
    ) == ["none", "vina", "vinardo"]

    with pytest.raises(ValueError, match="Unsupported rescoring"):
        RedockAnalysisApp._parse_protocol_rescore_methods("vina, imaginary")


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

    report = RedockAnalysisApp._write_protocol_report(rows, tmp_path).read_text()

    assert "| Engine | Box | Rescoring | Water |" in report
    assert "| smina | margin:4 | none | remove_all | 8 |" in report
    assert "| vina | 20x20x20 | none | retain_all | 16 |" in report
    assert "## Pose recovery across conditions" in report
    assert "## Recommended protocol per target" in report


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

    conditions = RedockAnalysisApp._expand_protocol_conditions(actives, sweep)
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

    report = RedockAnalysisApp._write_protocol_report(rows, tmp_path).read_text()

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

    report = RedockAnalysisApp._write_protocol_report(rows, tmp_path).read_text()
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

    prose, headers, rows = RedockAnalysisApp._parse_protocol_markdown(report)

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

    prose, tables = RedockAnalysisApp._parse_protocol_markdown_sections(report)

    assert prose[0] == "Protocol Development Summary"
    assert [table[0] for table in tables] == [
        "Protocol comparison", "Overall pose recovery", "Recommended protocol per target"
    ]
    assert tables[2][2] == [["NQO2", "rescored"]]


def test_protocol_chart_data_uses_completed_csv_rows_and_excludes_sentinels():
    frame = pd.DataFrame([
        {
            "status": "complete", "water_handling": "remove_all",
            "best_rmsd": 1.0, "top1_rmsd": 3.0, "top5_rmsd": 1.0,
            "rescore_top1_rmsd": 1.5, "rescore_top5_rmsd": 1.0,
        },
        {
            "status": "complete", "water_handling": "remove_all",
            "best_rmsd": 999.9, "top1_rmsd": 999.9, "top5_rmsd": 999.9,
            "rescore_top1_rmsd": 999.9, "rescore_top5_rmsd": 999.9,
        },
        {
            "status": "failed", "water_handling": "remove_all",
            "best_rmsd": 0.1, "top1_rmsd": 0.1, "top5_rmsd": 0.1,
            "rescore_top1_rmsd": 0.1, "rescore_top5_rmsd": 0.1,
        },
    ])

    data = RedockAnalysisApp._protocol_chart_data(frame)

    assert data["best"] == [("remove all", 100.0)]
    assert data["top1"] == [
        ("remove all\nbaseline", 0.0), ("remove all\nrescored", 100.0)
    ]
    assert data["top5"] == [
        ("remove all\nbaseline", 100.0), ("remove all\nrescored", 100.0)
    ]
    assert data["ranking_change"] == [
        ("Improved", 1.0), ("Unchanged", 0.0), ("Worse", 0.0)
    ]


def test_protocol_development_uses_unique_control_actives_only():
    pairs = [
        {"pdb_id": "1ABC", "site_ligand": "LIG", "chain": "A", "control_label": 1},
        {"pdb_id": "1ABC", "site_ligand": "LIG", "chain": "A", "control_label": 1},
        {"pdb_id": "1ABC", "site_ligand": "LIG", "chain": "A", "control_label": 0},
        {"pdb_id": "2DEF", "site_ligand": "DRG", "chain": None, "control_label": None},
        {"pdb_id": "3GHI", "site_ligand": "ACT", "chain": "B", "control_label": 1},
    ]

    actives = RedockAnalysisApp._protocol_active_pairs(pairs)

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

    assert RedockAnalysisApp._has_covalent_ligand_link(pdb, "JMR", "A") is True
    assert RedockAnalysisApp._has_covalent_ligand_link(pdb, "JMR", "B") is False
    assert RedockAnalysisApp._has_covalent_ligand_link(pdb, "STI", "A") is False


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
        pairs, size=4, seed=42, include_controls=True
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


def test_sampling_without_control_preservation_can_sample_any_case():
    pairs = [
        {"pdb_id": "1AAA", "ligand": "LIG", "control_label": 1},
        {"pdb_id": "1AAA", "ligand": "LIG", "control_label": 0},
        {"pdb_id": "1AAA", "ligand": "LIG", "control_label": None},
    ]

    sampled = _app_without_tk()._apply_random_sample(
        pairs, size=1, seed=1, include_controls=False
    )

    assert len(sampled) == 1
    assert sampled[0] in pairs


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


def test_pose_viewer_uses_protocol_development_csv_directly(tmp_path):
    protocol_csv = tmp_path / "protocol_development_results.csv"
    redock_json = tmp_path / "redock_results.json"

    assert RedockAnalysisApp._pose_results_csv(protocol_csv) == protocol_csv
    assert RedockAnalysisApp._pose_results_csv(redock_json) == tmp_path / "redock_results.csv"


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

    data = RedockAnalysisApp._screening_chart_data(summary)

    assert data["completion"] == [("1AAA\nL1", 75.0), ("2BBB\nL2", 100.0)]
    assert data["score_coverage"] == [("1AAA\nL1", 50.0), ("2BBB\nL2", 100.0)]
    assert data["failures"] == [("1AAA\nL1", 1.0), ("2BBB\nL2", 0.0)]
    assert data["pose_counts"] == [("Mean poses", 17.5)]


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
