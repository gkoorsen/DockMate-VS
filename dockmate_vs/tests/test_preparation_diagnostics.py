import io
import subprocess

import pytest
from Bio.PDB import PDBParser
from rdkit import Chem

from dockmate_vs.adaptive_docking import AdaptiveDockingPipeline
from dockmate_vs.config.schema import LigandPreparationConfig, ProteinPreparationConfig
from dockmate_vs.preparation.ligand import LigandPreparation, LigandPreparationError
from dockmate_vs.preparation.protein import PreparedProtein, ProteinPreparation, ProteinPreparationError


def test_receptor_conversion_uses_configured_timeout_and_reports_tool(monkeypatch, tmp_path):
    atom = "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C\nEND\n"
    structure = PDBParser(QUIET=True).get_structure("test", io.StringIO(atom))
    prepared = PreparedProtein("test", structure)
    seen = []

    def timeout(command, **kwargs):
        seen.append(kwargs["timeout"])
        raise subprocess.TimeoutExpired(command, kwargs["timeout"])

    monkeypatch.setattr(subprocess, "run", timeout)
    preparation = ProteinPreparation(ProteinPreparationConfig(obabel_timeout_sec=95))
    with pytest.raises(ProteinPreparationError, match=r"obabel.*after [\d.]+s.*limit 95s"):
        preparation.to_pdbqt(prepared, str(tmp_path / "receptor.pdbqt"))
    assert seen == [95]


def test_ligand_conversion_uses_configured_timeout_and_reports_tool(monkeypatch):
    seen = []

    def timeout(command, **kwargs):
        seen.append(kwargs["timeout"])
        raise subprocess.TimeoutExpired(command, kwargs["timeout"])

    monkeypatch.setattr(subprocess, "run", timeout)
    preparation = LigandPreparation(LigandPreparationConfig(obabel_timeout_sec=95))
    with pytest.raises(LigandPreparationError, match=r"obabel.*after [\d.]+s.*limit 95s"):
        preparation._to_pdbqt(Chem.MolFromSmiles("CCO"), "ethanol")
    assert seen == [95]


@pytest.mark.parametrize("timeout", [0, -1, float("inf"), float("nan")])
def test_preparation_rejects_invalid_timeouts(timeout):
    for config in (LigandPreparationConfig, ProteinPreparationConfig):
        with pytest.raises(ValueError):
            config(obabel_timeout_sec=timeout)


def test_ligand_selection_error_lists_chains_and_requested_chain_residues(tmp_path):
    source = tmp_path / "input.pdb"
    source.write_text(
        "HETATM    1  C1  LIG B 101       0.000   0.000   0.000  1.00 20.00           C\n"
        "HETATM    2  C1  LIG C 102       1.000   0.000   0.000  1.00 20.00           C\n"
        "HETATM    3  O   HOH A 201       2.000   0.000   0.000  1.00 20.00           O\n"
        "END\n"
    )
    pipeline = object.__new__(AdaptiveDockingPipeline)
    pipeline.output_dir = tmp_path
    with pytest.raises(ValueError) as exc:
        pipeline._extract_crystal_ligand(source, "LIG", "A")
    assert "Available chains for LIG: B, C" in str(exc.value)
    assert "HOH 201" in str(exc.value)
    assert "Chain column" in str(exc.value)
    selected = pipeline._extract_crystal_ligand(source, "LIG", "B")
    assert "LIG B 101" in selected.read_text()
