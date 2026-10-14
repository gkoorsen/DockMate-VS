import subprocess
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
import yaml
from rdkit import Chem
from rdkit.Chem import AllChem

from dockmate_vs.adaptive_docking import AdaptiveDockingPipeline, DockingProtocol
from dockmate_vs.config.schema import DockingConfig
from dockmate_vs.docking import registry
from dockmate_vs.docking.base import DockingEngine, DockingPose, DockingResult
from dockmate_vs.docking.registry import EngineCapabilities, create_engine, register_engine
from dockmate_vs.gui.app import DockMateVSApp
from dockmate_vs.headless import HeadlessDockMateRunner, doctor, run_campaign


class ExampleEngine(DockingEngine):
    @classmethod
    def from_config(cls, config):
        instance = cls(seed=config["seed"])
        instance.binary = config["engine_options"]["binaries"]["example-dock"]
        instance.timeout = config["timeout"]
        return instance

    def dock(self, receptor_path, ligand_path, center, size, output_path=None, **kwargs):
        subprocess.run([self.binary, "--receptor", receptor_path, "--ligand", ligand_path,
                        "--out", output_path, "--seed", str(self.seed)],
                       check=True, capture_output=True, timeout=self.timeout)
        return self.parse_output(output_path)

    def parse_output(self, output_file):
        mols = Chem.SDMolSupplier(str(output_file), removeHs=False)
        return DockingResult("example", [DockingPose(float(m.GetProp("SCORE")), mode_number=i)
                                          for i, m in enumerate(mols, 1)])

    def get_version(self):
        return "example 1"


@pytest.fixture
def example(monkeypatch):
    monkeypatch.setattr(registry, "_ENGINES", dict(registry._ENGINES))
    return register_engine("example", ExampleEngine, EngineCapabilities(
        supports_exhaustiveness=False, supports_score_only_rescoring=False,
        binary_names=("example-dock",), receptor_format="pdb", ligand_format="sdf", output_suffix=".sdf",
    ))


def test_registration_validation_config_and_capabilities(example):
    assert DockingConfig(engine="example").engine == "example"
    engine = create_engine("example", {"seed": 7, "timeout": 12, "engine_options": {
        "example": {"binaries": {"example-dock": "/example/program"}},
    }})
    assert engine.binary == "/example/program"
    assert engine.seed == 7
    with pytest.raises(ValueError, match="already registered"):
        register_engine("example", ExampleEngine)
    with pytest.raises(ValueError, match="Unknown docking engine"):
        create_engine("absent", {})
    with pytest.raises(ValueError, match="must be a mapping"):
        registry.validate_engine_options({"example": []})
    sweep = {"engines": ["example"], "water_modes": ["remove_all"],
             "box_definitions": [{"label": "margin:4", "box_margin": 4}],
             "rescore_methods": ["none"], "exhaustiveness": [8, 16], "seeds": [42]}
    assert len(DockMateVSApp._expand_protocol_conditions([{}], sweep)) == 1
    sweep["rescore_methods"] = ["vinardo"]
    with pytest.raises(ValueError, match="rescor"):
        DockMateVSApp._expand_protocol_conditions([{}], sweep)


def test_doctor_lists_registered_binary_and_status(example, monkeypatch, capsys):
    monkeypatch.setattr(registry.shutil, "which", lambda name: "/bin/example-dock" if "example-dock" in name else None)
    assert doctor() == 0
    output = capsys.readouterr().out
    assert "example:" in output and "example-dock: /bin/example-dock" in output
    assert "MISSING" in output
    assert doctor(strict=True) == 1


def test_registered_engine_runs_and_resumes_real_protocol_campaign(example, monkeypatch, tmp_path):
    molecule = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    assert AllChem.EmbedMolecule(molecule, randomSeed=42) == 0
    reference = tmp_path / "1abc.pdb"
    for index, atom in enumerate(molecule.GetAtoms(), 1):
        info = Chem.AtomPDBResidueInfo(f"{atom.GetSymbol()}{index}".rjust(4), index)
        info.SetResidueName("UNL")
        info.SetResidueNumber(1)
        info.SetChainId("A")
        info.SetIsHeteroAtom(True)
        atom.SetMonomerInfo(info)
    Chem.MolToPDBFile(molecule, str(reference))
    ligand = tmp_path / "sample_v1.sdf"
    with Chem.SDWriter(str(ligand)) as writer:
        writer.write(molecule)
    input_pdbqt = ligand.with_suffix(".pdbqt")
    input_pdbqt.write_text("placeholder\n")
    monkeypatch.setattr(HeadlessDockMateRunner, "_prefetch_remote_inputs", lambda *a: None)
    monkeypatch.setattr(HeadlessDockMateRunner, "_download_pdb", lambda *a: reference)
    monkeypatch.setattr(HeadlessDockMateRunner, "_software_provenance", lambda *a: {})
    monkeypatch.setattr(HeadlessDockMateRunner, "_prepare_receptor_with_cache", lambda *a, **k: (input_pdbqt, reference))
    monkeypatch.setattr(AdaptiveDockingPipeline, "_prepare_ligand_variants", lambda *a, **k: [
        {"label": "sample_v1", "pdbqt": input_pdbqt, "smiles": "CCO"},
    ])
    calls = []

    def execute(command, **kwargs):
        if command[-1] in {"--version", "-V"}:
            return SimpleNamespace(returncode=0, stdout="test version", stderr="")
        assert "--out" in command, command
        calls.append((command, kwargs))
        output = Path(command[command.index("--out") + 1])
        posed = Chem.Mol(molecule)
        posed.SetProp("SCORE", "-7.5")
        with Chem.SDWriter(str(output)) as writer:
            writer.write(posed)
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", execute)
    workbook = tmp_path / "cases.xlsx"
    pd.DataFrame([{"PDB_ID": "1ABC", "Ligand": "UNL", "Target_Ligand": "control",
                   "SMILES": "CCO", "decoy compound": "decoy", "decoy SMILES": "CCC"}]).to_excel(workbook, index=False)
    config = tmp_path / "campaign.yml"
    payload = {"input_file": str(workbook), "output_dir": str(tmp_path / "run"),
               "single": {"engine": "example", "engine_options": {"example": {"binaries": {"example-dock": "/fake/dock"}}}},
               "protocol_sweep": {"engines": ["example"], "water_modes": ["remove_all"],
                                  "rescore_methods": ["none"], "exhaustiveness": [8, 16], "seeds": [23]}}
    config.write_text(yaml.safe_dump(payload))
    result_path = run_campaign(config, "protocol")
    rows = pd.read_csv(result_path)
    assert len(calls) == 1
    assert calls[0][0][0] == "/fake/dock"
    assert calls[0][0][-1] == "23"
    assert calls[0][1]["timeout"] == 1200
    assert rows.loc[0, "engine"] == "example"
    assert rows.loc[0, "best_score"] == -7.5
    assert rows.loc[0, "status"] == "complete"
    assert rows.loc[0, "best_rmsd"] < 0.01
    assert str(rows.loc[0, "case_id"]).startswith("1ABC_UNL_example_")
    assert "tk" not in HeadlessDockMateRunner().__dict__
    assert run_campaign(config, "protocol") == result_path
    assert len(calls) == 1
    payload["single"]["engine_options"]["example"]["binaries"]["example-dock"] = "/different/dock"
    config.write_text(yaml.safe_dump(payload))
    with pytest.raises(RuntimeError, match="settings have changed"):
        run_campaign(config, "protocol")


def test_adaptive_vina_uses_vina_registration(tmp_path):
    pipeline = AdaptiveDockingPipeline(tmp_path, use_vina=True, rdock_available=False)
    assert {protocol.engine for protocol in pipeline.protocols} == {"vina"}


def test_custom_engine_reached_by_adaptive_dispatch(example, monkeypatch, tmp_path):
    from dockmate_vs.binding_site.cocrystal import BindingSiteDefinition

    mol = Chem.MolFromSmiles("CCO")
    AllChem.Compute2DCoords(mol)
    ligand = tmp_path / "prepared.pdbqt"
    ligand.write_text("prepared placeholder")
    with Chem.SDWriter(str(ligand.with_suffix(".sdf"))) as writer:
        writer.write(mol)
    monkeypatch.setattr(BindingSiteDefinition, "from_cocrystal", lambda *a, **k: SimpleNamespace(center=np.zeros(3), size=np.ones(3) * 20))
    calls = []

    def execute(command, **kwargs):
        if "--out" not in command:
            return SimpleNamespace(returncode=0, stdout="version", stderr="")
        calls.append(command)
        mol.SetProp("SCORE", "-8.0")
        with Chem.SDWriter(command[command.index("--out") + 1]) as writer:
            writer.write(mol)
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", execute)
    pipeline = AdaptiveDockingPipeline(tmp_path / "adaptive", molecule_type="decoy",
                                      engine_options={"example": {"binaries": {"example-dock": "/fake/example"}}})
    protocol = DockingProtocol(name="custom", description="test", engine="example", margin=4,
                               scoring="unused", exhaustiveness=8, expected_runtime_min=1)
    result = pipeline._run_single_protocol(protocol, ligand, "CCO", "sample", "UNL", "A", None,
                                           prepared_ligand_pdbqt=ligand, prepared_receptor_pdbqt=ligand,
                                           prepared_receptor_pdb=ligand)
    assert calls[0][0] == "/fake/example"
    assert result.engine == "example"
    assert result.score == -8.0
    assert result.output_file.suffix == ".sdf"


def test_gui_controls_use_registered_capabilities(example):
    app = object.__new__(DockMateVSApp)
    changes = {}

    def widget(name):
        return SimpleNamespace(configure=lambda **kwargs: changes.setdefault(name, {}).update(kwargs))

    app.engine_var = SimpleNamespace(get=lambda: "example")
    app.rdock_single_frame = SimpleNamespace(grid=lambda: None, grid_remove=lambda: None)
    app._scoring_menu = widget("scoring")
    app._exhaustiveness_entry = widget("exhaustiveness")
    app._rescore_checkbox = widget("rescore")
    app._rescore_menu = widget("rescore_menu")
    app.enable_rescore_var = SimpleNamespace(set=lambda value: changes.update(enabled=value))
    app._update_engine()
    assert changes["exhaustiveness"]["state"] == "disabled"
    assert changes["rescore"]["state"] == "disabled"
    assert changes["enabled"] is False


@pytest.mark.parametrize("name", ["smina", "vina"])
@pytest.mark.parametrize("failure", [None, "exit", "timeout"])
def test_builtin_vina_family_dispatch_and_failure_handling(name, failure, monkeypatch, tmp_path):
    commands = []

    def execute(command, **kwargs):
        if "--version" in command:
            return SimpleNamespace(returncode=0, stdout="test version", stderr="")
        commands.append(command)
        assert kwargs["timeout"] == 19
        assert "--ligand_smiles" not in command and "--prepared_receptor_pdb" not in command
        if failure == "timeout":
            raise subprocess.TimeoutExpired(command, 19)
        if failure == "exit":
            return SimpleNamespace(returncode=1, stdout="", stderr="bad input")
        Path(command[command.index("--out") + 1]).write_text(
            "MODEL 1\nREMARK VINA RESULT: -8.1 0.0 0.0\nENDMDL\n"
            "MODEL 2\nREMARK VINA RESULT: -7.0 1.0 2.0\nENDMDL\n"
        )
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", execute)
    inputs = tmp_path / "prepared.pdbqt"
    inputs.write_text("ATOM\n")
    adapter = create_engine(name, {f"{name}_binary": f"/tools/{name}", "scoring": "vinardo", "timeout": 19})
    result = adapter.dock(str(inputs), str(inputs), np.zeros(3), np.ones(3) * 20,
                          str(tmp_path / "docked.pdbqt"), ligand_smiles="CCO", ligand_name="sample",
                          prepared_receptor_pdb=tmp_path / "receptor.pdb", reference_ligand_pdb=None)
    assert commands[0][0] == f"/tools/{name}"
    assert ("--scoring" in commands[0]) == (name == "smina")
    if failure:
        assert not result.success and result.error_message
    else:
        assert result.success
        assert [p.score for p in result.poses] == [-8.1, -7.0]


def test_rdock_adapter_keeps_cavity_commands_and_pose_order(monkeypatch, tmp_path):
    root = tmp_path / "rdock"
    (root / "bin").mkdir(parents=True)
    (root / "data/scripts").mkdir(parents=True)
    (root / "data/scripts/dock.prm").write_text("protocol\n")
    for binary in ("rbdock", "rbcavity"):
        (root / "bin" / binary).write_text("placeholder\n")
    receptor = tmp_path / "receptor.pdb"
    receptor.write_text("ATOM\n")
    monkeypatch.setattr(registry.shutil, "which", lambda name: str(name))
    commands = []

    def execute(command, **kwargs):
        commands.append((command, kwargs))
        binary = Path(command[0]).name
        if binary == "obabel":
            Path(command[command.index("-O") + 1]).write_text("MOL2\n")
            assert kwargs["timeout"] == 61
        elif binary == "rbcavity":
            (Path(kwargs["cwd"]) / "receptor.as").write_text("cavity\n")
        else:
            mol = Chem.SDMolSupplier(str(Path(kwargs["cwd"]) / command[command.index("-i") + 1]), removeHs=False)[0]
            with Chem.SDWriter(str(Path(kwargs["cwd"]) / "docked.sd")) as writer:
                for score in (-6.0, -9.0):
                    mol.SetProp("SCORE", str(score))
                    writer.write(mol)
            assert command[command.index("-n") + 1] == "3"
            assert command[command.index("-s") + 1] == "13"
            assert kwargs["timeout"] == 27
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", execute)
    adapter = create_engine("rdock", {"rdock_root": str(root), "rdock_runs": 3,
                                      "rdock_seed": 13, "timeout": 27, "obabel_timeout_sec": 61})
    output = tmp_path / "case/docked.sd"
    result = adapter.dock(str(receptor), "sample.pdbqt", np.ones(3) * 10, np.ones(3) * 20,
                          str(output), ligand_smiles="CCO", prepared_receptor_pdb=receptor)
    assert result.success, result.error_message
    assert [p.score for p in result.poses] == [-6.0, -9.0]
    assert result.best_score == -9.0
    assert [Path(command[0]).name for command, _ in commands] == ["obabel", "rbcavity", "rbdock"]
    cavity = Chem.SDMolSupplier(str(output.parent / "cavity_reference.sdf"), removeHs=False)[0]
    np.testing.assert_allclose(cavity.GetConformer().GetPositions().mean(axis=0), [10, 10, 10], atol=0.001)
