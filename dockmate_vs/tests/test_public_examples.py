from pathlib import Path

import pandas as pd
import yaml
from rdkit import Chem

from dockmate_vs.gui.app import DockMateVSApp
from dockmate_vs.headless import load_campaign_config


ROOT = Path(__file__).resolve().parents[2]
EXAMPLES = ROOT / "examples"


def test_campaign_configs_reference_public_aces_workbooks():
    protocol = yaml.safe_load((EXAMPLES / "campaign.protocol.yml").read_text())
    screen = yaml.safe_load((EXAMPLES / "campaign.screen.yml").read_text())

    assert protocol["input_file"] == "dude_aces_protocol_development_1E66.xlsx"
    assert screen["input_file"] == "dude_aces_screening_subset_1E66_seed42.xlsx"
    assert (EXAMPLES / protocol["input_file"]).is_file()
    assert (EXAMPLES / screen["input_file"]).is_file()
    condition_count = (
        len(protocol["protocol_sweep"]["engines"])
        * len(protocol["protocol_sweep"]["box_definitions"])
        * len(protocol["protocol_sweep"]["water_modes"])
        * len(protocol["protocol_sweep"]["exhaustiveness"])
        * len(protocol["protocol_sweep"]["seeds"])
        * len(protocol["protocol_sweep"]["rescore_methods"])
    )
    assert condition_count == 24
    assert screen["sampling"]["enabled"] is False
    assert screen["single"]["variant_select_by"] == "score"

    loaded_protocol = load_campaign_config(
        EXAMPLES / "campaign.protocol.yml", "protocol"
    )
    loaded_screen = load_campaign_config(EXAMPLES / "campaign.screen.yml", "screen")
    assert Path(loaded_protocol["input_file"]).name == protocol["input_file"]
    assert Path(loaded_screen["input_file"]).name == screen["input_file"]
    assert loaded_protocol["mode"] == "protocol_development"
    assert loaded_screen["mode"] == "screening"


def test_protocol_example_contains_native_aces_control():
    workbook = EXAMPLES / "dude_aces_protocol_development_1E66.xlsx"
    frame = pd.read_excel(workbook, sheet_name="Docking_Jobs", dtype=object)

    assert len(frame) == 1
    assert frame.iloc[0]["Target_Ligand"] == "HUX_native_redock"
    assert frame.iloc[0]["PDB_ID"] == "1E66"
    assert frame.iloc[0]["Ligand"] == "HUX"
    assert frame.iloc[0]["Chain"] == "A"
    assert frame.iloc[0]["label"] == 1
    assert pd.isna(frame.iloc[0]["SMILES"])

    app = object.__new__(DockMateVSApp)
    pairs, _ = app._load_pairs_from_excel(workbook)
    assert len(pairs) == 1
    assert pairs[0]["pdb_id"] == "1E66"
    assert pairs[0]["control_label"] == 1


def test_screening_example_contains_valid_dude_aces_benchmark():
    workbook = EXAMPLES / "dude_aces_screening_subset_1E66_seed42.xlsx"
    frame = pd.read_excel(workbook, sheet_name="Docking_Jobs", dtype=object)

    assert len(frame) == 220
    assert frame["label"].value_counts().to_dict() == {0: 200, 1: 20}
    assert set(frame["PDB_ID"]) == {"1E66"}
    assert set(frame["Ligand"]) == {"HUX"}
    assert set(frame["Chain"]) == {"A"}
    assert set(frame["Protein"]) == {"DUD-E ACES"}
    assert frame["Target_Ligand"].is_unique
    assert frame["Source_ID"].is_unique
    assert frame["SMILES"].notna().all()
    assert all(Chem.MolFromSmiles(value) is not None for value in frame["SMILES"])

    app = object.__new__(DockMateVSApp)
    pairs, _ = app._load_pairs_from_excel(workbook)
    assert len(pairs) == 220
    assert {pair["pdb_id"] for pair in pairs} == {"1E66"}
    assert sum(pair["control_label"] == 1 for pair in pairs) == 20
    assert sum(pair["control_label"] == 0 for pair in pairs) == 200

    metadata = pd.read_excel(workbook, sheet_name="Metadata")
    values = dict(zip(metadata["Field"], metadata["Value"]))
    assert values["Subset"] == "20 actives and 200 decoys"
    assert "random.Random(42)" in values["Selection"]
    assert values["Active source SHA-256"] == (
        "d3c04077ab78028fa104a585b91f8b1e80361d71923f66800a9386866a3d1168"
    )
    assert values["Decoy source SHA-256"] == (
        "da06b6b9d0a4ab51c5dfad2938242dc877ef453c67bc433d10568cd88fd1280d"
    )
