import subprocess
from types import SimpleNamespace

import numpy as np

from dockmate_vs.docking.smina import SminaDockingEngine


def test_unsupported_vanadium_is_structured_and_never_substituted(monkeypatch, tmp_path):
    receptor = tmp_path / "receptor.pdbqt"
    receptor.write_text('ATOM      1  V   VO4 A 101       0.000   0.000   0.000  1.00  0.00    +0.000 V\n')
    ligand = tmp_path / "ligand.pdbqt"
    ligand.write_text("ROOT\nENDROOT\n")
    original = receptor.read_bytes()
    engine = object.__new__(SminaDockingEngine)
    engine.timeout_sec = 10
    monkeypatch.setattr(engine, "_build_command", lambda **kw: ["smina"])
    monkeypatch.setattr(subprocess, "run", lambda *a, **kw: SimpleNamespace(
        returncode=1, stdout="", stderr='Parse error: "V" is not a valid AutoDock type.'))
    result = engine.dock(str(receptor), str(ligand), np.zeros(3), np.ones(3),
                         output_path=str(tmp_path / "docked.pdbqt"))
    assert result.success is False
    assert result.poses == []
    assert result.metadata["failure_category"] == "unsupported_atom_type"
    assert result.metadata["unsupported_atom_type"] == "V"
    assert result.metadata["atom_substitution_policy"] == "none"
    assert "retained unchanged" in result.error_message
    assert receptor.read_bytes() == original


def test_other_engine_errors_are_not_classified_as_unsupported_atoms():
    message, metadata = SminaDockingEngine.failure_details("Could not open input")
    assert message == "Could not open input"
    assert metadata == {}
