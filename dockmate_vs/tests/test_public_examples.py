from pathlib import Path

import pandas as pd
import yaml
from rdkit import Chem


ROOT = Path(__file__).resolve().parents[2]
EXAMPLES = ROOT / "examples"


def test_screening_campaign_uses_screening_workbook():
    config = yaml.safe_load((EXAMPLES / "campaign.screen.yml").read_text())

    assert config["input_file"] == "esr1_enrichment_smoke_test_1XP1.xlsx"
    assert (EXAMPLES / config["input_file"]).is_file()


def test_screening_example_contains_valid_enrichment_controls():
    workbook = EXAMPLES / "esr1_enrichment_smoke_test_1XP1.xlsx"
    frame = pd.read_excel(workbook, sheet_name="Docking_Jobs")

    assert len(frame) == 7
    assert frame["label"].value_counts().to_dict() == {0: 6, 1: 1}
    assert set(frame["PDB_ID"]) == {"1XP1"}
    assert set(frame["Ligand"]) == {"AIH"}
    assert set(frame["Chain"]) == {"A"}
    assert frame["Target_Ligand"].is_unique
    smiles = frame.loc[frame["SMILES"].notna(), "SMILES"]
    assert len(smiles) == 6
    assert all(Chem.MolFromSmiles(value) is not None for value in smiles)
