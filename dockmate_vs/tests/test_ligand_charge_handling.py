"""Charge preservation, legacy reproduction, and campaign isolation."""

import copy
from dataclasses import asdict
import json
import queue
from types import SimpleNamespace

import pandas as pd
import pytest
from rdkit import Chem

from dockmate_vs.adaptive_docking import AdaptiveDockingPipeline
from dockmate_vs.config.schema import LigandPreparationConfig
from dockmate_vs.gui.app import DockMateVSApp
from dockmate_vs.headless import load_campaign_config
from dockmate_vs.preparation.ligand import LigandPreparation
from dockmate_vs.preparation.ligand_cache import LigandCache
from dockmate_vs.tests.test_screening_integration import _protocol_resume_config, _resume_manifest
from dockmate_vs.tests.test_screening_pipeline import _result


def _charges(mol):
    return [atom.GetFormalCharge() for atom in mol.GetAtoms()]


@pytest.mark.parametrize("smiles, neutral", [
    ("CC[NH3+]", "CCN"),
    ("CC(=O)[O-]", "CC(=O)O"),
    ("[NH3+]CC(=O)[O-]", "NCC(=O)O"),
    ("C[N+](C)(C)C", "C[N+](C)(C)C"),
    ("CCO", "CCO"),
])
def test_input_charges_preserved_and_legacy_neutralization_reproduced(smiles, neutral):
    molecule = Chem.MolFromSmiles(smiles)
    prep = LigandPreparation(LigandPreparationConfig(), n_cpus=1)
    preserved = prep._standardize_molecule(molecule)
    assert Chem.MolToSmiles(preserved) == Chem.MolToSmiles(molecule)
    legacy = LigandPreparation(LigandPreparationConfig(charge_handling="neutralize"), n_cpus=1)
    assert Chem.MolToSmiles(legacy._standardize_molecule(molecule)) == (
        Chem.MolToSmiles(Chem.MolFromSmiles(neutral))
    )


def test_salt_removal_keeps_parent_charge():
    prep = LigandPreparation(LigandPreparationConfig(), n_cpus=1)
    result = prep._standardize_molecule(Chem.MolFromSmiles("CC[NH3+].[Cl-]"))
    assert Chem.MolToSmiles(result) == "CC[NH3+]"
    assert Chem.GetFormalCharge(result) == 1


def test_tautomers_cannot_neutralize_a_zwitterion(monkeypatch):
    molecule = Chem.MolFromSmiles("[NH3+]CC(=O)[O-]")
    neutral = Chem.MolFromSmiles("NCC(=O)O")
    assert Chem.GetFormalCharge(molecule) == Chem.GetFormalCharge(neutral) == 0
    enumerator = SimpleNamespace(SetMaxTautomers=lambda n: None,
                                 Enumerate=lambda mol: [neutral, molecule])
    monkeypatch.setattr("dockmate_vs.preparation.ligand.rdMolStandardize.TautomerEnumerator",
                        lambda: enumerator)
    prep = LigandPreparation(LigandPreparationConfig(), n_cpus=1)
    variants = prep._enumerate_tautomers(molecule)
    assert len(variants) == 1
    assert _charges(variants[0]) == _charges(molecule)


@pytest.mark.parametrize("policy, expected", [("preserve", [-1, 1]), ("neutralize", [])])
@pytest.mark.parametrize("enumerate_states", [False, True])
def test_charges_survive_through_3d_preparation_to_pdbqt_boundary(
    monkeypatch, policy, expected, enumerate_states
):
    prep = LigandPreparation(LigandPreparationConfig(
        charge_handling=policy, max_tautomers=2, max_conformers=1
    ), n_cpus=1)
    observed = []
    def convert(molecule, name):
        observed.append(sorted(charge for charge in _charges(molecule) if charge))
        assert molecule.GetNumConformers() == 1
        return "REMARK test PDBQT boundary\n"
    monkeypatch.setattr(prep, "_to_pdbqt", convert)
    results = prep.prepare_from_smiles("[NH3+]CC(=O)[O-]", "glycine", enumerate_states)
    assert results
    assert all(charges == expected for charges in observed)
    assert all(sorted(q for q in _charges(Chem.MolFromSmiles(r.smiles)) if q) == expected
               for r in results)
    assert {r.protonation_state for r in results} == {
        "input_charge" if policy == "preserve" else "standardized"
    }


def test_cache_does_not_reuse_a_different_charge_policy(tmp_path):
    cache = LigandCache(str(tmp_path))
    preserve = cache._get_cache_key("CC[NH3+]", LigandPreparationConfig())
    neutralize = cache._get_cache_key("CC[NH3+]", LigandPreparationConfig(charge_handling="neutralize"))
    assert preserve != neutralize


@pytest.mark.parametrize("policy", ["preserve", "neutralize"])
def test_pipeline_forwards_charge_policy_to_preparation(monkeypatch, tmp_path, policy):
    monkeypatch.chdir(tmp_path)
    pipeline = AdaptiveDockingPipeline(tmp_path / "run", charge_handling=policy)
    configs = []
    def get(**kwargs):
        configs.append(kwargs["config"])
        return [SimpleNamespace(energy=0., smiles="CC[NH3+]",
                                to_pdbqt_file=lambda path: None)]
    monkeypatch.setattr(pipeline.ligand_cache, "get", get)
    pipeline._prepare_ligand("CC[NH3+]", "amine")
    pipeline._prepare_ligand_variants("CC[NH3+]", "amine", max_tautomers=1)
    assert [cfg.charge_handling for cfg in configs] == [policy, policy]


@pytest.mark.parametrize("policy", [None, "preserve", "neutralize", "invalid"])
def test_headless_charge_policy_defaults_and_validation(tmp_path, policy):
    workbook = tmp_path / "input.xlsx"
    workbook.touch()
    payload = {"input_file": "input.xlsx", "output_dir": "output"}
    if policy is not None:
        payload["single"] = {"charge_handling": policy}
    path = tmp_path / "config.json"
    path.write_text(json.dumps(payload))
    if policy == "invalid":
        with pytest.raises(ValueError, match="charge_handling"):
            load_campaign_config(path, "screen")
        with pytest.raises(ValueError):
            LigandPreparationConfig(charge_handling=policy)
    else:
        assert load_campaign_config(path, "screen")["single"]["charge_handling"] == (policy or "preserve")


@pytest.mark.parametrize("policy, allowed", [("neutralize", True), ("preserve", False)])
def test_legacy_screening_resume_requires_neutralization(tmp_path, policy, allowed):
    app = object.__new__(DockMateVSApp)
    legacy = _resume_manifest(tmp_path)
    del legacy["config"]["single"]["charge_handling"]
    manifest_path = tmp_path / "run_manifest.json"
    app._write_json_atomic(manifest_path, legacy)
    current = copy.deepcopy(legacy)
    current["config"]["single"]["charge_handling"] = policy
    output = tmp_path / "docked.pdbqt"
    output.touch()
    result = _result(case_id="1ABC_LIG_sample", output_file=str(output), docking_completed=True)
    progress_path = tmp_path / "redock_progress.json"
    progress_path.write_text(json.dumps({"results": [asdict(result)]}))
    assert bool(app._load_resumable_results(manifest_path, progress_path, current)) is allowed


@pytest.mark.parametrize("policy, allowed", [("neutralize", True), ("preserve", False)])
def test_legacy_protocol_resume_requires_neutralization(tmp_path, policy, allowed):
    app = object.__new__(DockMateVSApp)
    config = _protocol_resume_config(tmp_path, tmp_path / "input.xlsx")
    actives = [{"pdb_id": "1XP1", "site_ligand": "AIH"}]
    signature = app._protocol_resume_signature(actives, config)
    del signature["single"]["charge_handling"]
    previous = {"input_file": config["input_file"], "resume_signature": signature}
    manifest_path = tmp_path / "protocol_development_manifest.json"
    manifest_path.write_text(json.dumps(previous))
    results_path = tmp_path / "protocol_development_results.csv"
    pd.DataFrame([{"pdb_id": "1XP1", "ligand_resname": "AIH"}]).to_csv(results_path, index=False)
    config["single"]["charge_handling"] = policy
    current = {"input_file": config["input_file"],
               "resume_signature": app._protocol_resume_signature(actives, config)}
    assert (app._protocol_resume_incompatibility(manifest_path, results_path, current, actives)
            is None) is allowed


@pytest.mark.parametrize("old_policy,new_policy", [(None, "preserve"),
                         ("neutralize", "preserve"), ("preserve", "neutralize")])
def test_changed_charge_policy_cannot_overwrite_saved_campaign(tmp_path, old_policy, new_policy):
    app = object.__new__(DockMateVSApp)
    app._queue = queue.Queue()
    config = {"mode": "screening", "output_dir": tmp_path,
              "single": {"charge_handling": new_policy}}
    previous = {"config": {"mode": "screening", "single": {}}}
    if old_policy:
        previous["config"]["single"]["charge_handling"] = old_policy
    manifest_path = tmp_path / "run_manifest.json"
    original = json.dumps(previous)
    manifest_path.write_text(original)
    app._run_worker([], config)
    assert app._queue.get_nowait()[0] == "preflight_failed"
    assert manifest_path.read_text() == original
    assert not (tmp_path / "redock_results.json").exists()


def test_unreadable_manifest_is_not_overwritten(tmp_path):
    app = object.__new__(DockMateVSApp)
    app._queue = queue.Queue()
    path = tmp_path / "run_manifest.json"
    path.write_text("incomplete JSON")
    app._run_worker([], {"mode": "screening", "output_dir": tmp_path})
    assert app._queue.get_nowait()[0] == "preflight_failed"
    assert path.read_text() == "incomplete JSON"
