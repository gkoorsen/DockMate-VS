import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

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
