from pathlib import Path

import pandas as pd

from docking_platform_gui.gui.redock_analysis import RedockAnalysisApp


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

