import json
import shutil

import gemmi
import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from dockmate_vs.cli import main
from dockmate_vs.pose_export import export_cases, export_pose


def _pdbqt(mol, order, shift=0):
    lines = []
    for serial, index in enumerate(order, 1):
        atom = mol.GetAtomWithIdx(index)
        x, y, z = mol.GetConformer().GetAtomPosition(index)
        element = atom.GetSymbol()
        name = element + str(index + 1)
        lines.append(f"ATOM  {serial:5d} {name:>4s} LIG A   1    {x + shift:8.3f}{y:8.3f}{z:8.3f}  0.00  0.00     0.000 {element}\n")
    return "".join(lines)


@pytest.fixture
def saved_run(tmp_path):
    run = tmp_path / "run"
    case = run / "sample"
    variants = case / "ligand_variants"
    docked = case / "variants" / "compound_v1"
    variants.mkdir(parents=True)
    docked.mkdir(parents=True)
    mol = Chem.AddHs(Chem.MolFromSmiles("C[C@H](O)c1ccccc1"))
    assert AllChem.EmbedMolecule(mol, randomSeed=42) == 0
    with Chem.SDWriter(str(variants / "compound_v1.sdf")) as writer:
        writer.write(mol)
    # Deliberately reorder atoms, including carbons, as torsion-tree PDBQT does.
    order = list(reversed([a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() != 1]))
    (variants / "compound_v1.pdbqt").write_text(_pdbqt(mol, order))
    output = docked / "docked.pdbqt"
    output.write_text("MODEL 1\n" + _pdbqt(mol, order, 10) + "ENDMDL\nMODEL 2\n" + _pdbqt(mol, order, 20) + "ENDMDL\n")
    (case / "receptor_prepared.pdb").write_text(
        "ATOM      1  CA  ALA A   1      -8.000   1.000   1.000  1.00  0.00           C  \nEND\n"
    )
    (run / "redock_results.json").write_text(json.dumps({"results": [{
        "case_id": "sample", "output_file": str(output), "engine": "smina",
    }]}))
    return run, mol, output


def test_export_roundtrip_preserves_topology_and_requested_docked_coordinates(saved_run, tmp_path):
    run, prepared, _ = saved_run
    output = export_pose(run, "sample", pose=2, output=tmp_path / "pose.sdf")
    actual = Chem.SDMolSupplier(str(output), removeHs=False)[0]
    assert actual.GetNumAtoms() == prepared.GetNumAtoms()
    assert Chem.MolToSmiles(Chem.RemoveHs(actual)) == Chem.MolToSmiles(Chem.RemoveHs(prepared))
    heavy = [a.GetIdx() for a in prepared.GetAtoms() if a.GetAtomicNum() != 1]
    expected = prepared.GetConformer().GetPositions()[heavy] + [20, 0, 0]
    np.testing.assert_allclose(actual.GetConformer().GetPositions()[heavy], expected, atol=0.001)
    metadata = json.loads(output.with_suffix(".sdf.json").read_text())
    assert metadata["rebuilt_hydrogens"] == prepared.GetNumAtoms() - len(heavy)
    assert metadata["md_ready"] is False


def test_complex_has_unique_serials_selectable_ligand_and_companion_sdf(saved_run, tmp_path):
    run, prepared, _ = saved_run
    output = export_pose(run, "sample", format="complex-pdb", output=tmp_path / "complex.pdb", ligand_resname="DRG")
    structure = gemmi.read_structure(str(output))
    ligands = [res for chain in structure[0] for res in chain if res.name == "DRG"]
    assert len(ligands) == 1
    assert len(ligands[0]) == prepared.GetNumAtoms()
    serials = [atom.serial for chain in structure[0] for res in chain for atom in res]
    assert len(serials) == len(set(serials))
    assert output.with_suffix(".ligand.sdf").exists()


def test_export_moved_run_prefers_selected_copy_and_does_not_modify_it(saved_run, tmp_path):
    run, _, source = saved_run
    moved = tmp_path / "copied_run"
    shutil.copytree(run, moved)
    source.write_text("broken original\n")
    before = {p.relative_to(moved): p.read_bytes() for p in moved.rglob("*") if p.is_file()}
    export_pose(moved, "sample", output=tmp_path / "moved.sdf")
    after = {p.relative_to(moved): p.read_bytes() for p in moved.rglob("*") if p.is_file()}
    assert before == after


@pytest.mark.parametrize("pose", [0, -1, 3])
def test_export_rejects_invalid_pose_number(saved_run, tmp_path, pose):
    with pytest.raises(ValueError, match="positive|outside"):
        export_pose(saved_run[0], "sample", pose=pose, output=tmp_path / "out.sdf")
    assert not (tmp_path / "out.sdf").exists()


def test_export_rejects_mismatched_preparation_and_does_not_overwrite(saved_run, tmp_path):
    run, _, _ = saved_run
    output = export_pose(run, "sample", output=tmp_path / "out.sdf")
    original = output.read_bytes()
    with pytest.raises(ValueError, match="overwrite"):
        export_pose(run, "sample", output=output)
    assert output.read_bytes() == original
    prepared = run / "sample/ligand_variants/compound_v1.sdf"
    mol = Chem.SDMolSupplier(str(prepared), removeHs=False)[0]
    for index in range(mol.GetNumAtoms()):
        mol.GetConformer().SetAtomPosition(index, (0, 0, 0))
    with Chem.SDWriter(str(prepared)) as writer:
        writer.write(mol)
    with pytest.raises(ValueError, match="do not match uniquely"):
        export_pose(run, "sample", output=tmp_path / "bad.sdf")


def test_export_sdf_input_uses_prepared_bond_orders(saved_run, tmp_path):
    run, mol, original = saved_run
    output = original.with_suffix(".sd")
    order = list(reversed(range(mol.GetNumAtoms())))
    docked = Chem.RenumberAtoms(mol, order)
    with Chem.SDWriter(str(output)) as writer:
        writer.write(docked)
    results = run / "redock_results.json"
    results.write_text(json.dumps({"results": [{"case_id": "sample", "output_file": str(output)}]}))
    actual_path = export_pose(run, "sample", output=tmp_path / "sd_export.sdf")
    actual = Chem.SDMolSupplier(str(actual_path), removeHs=False)[0]
    assert actual.GetNumAtoms() == mol.GetNumAtoms()
    assert Chem.MolToSmiles(Chem.RemoveHs(actual)) == Chem.MolToSmiles(Chem.RemoveHs(mol))


def test_export_cli_lists_cases_and_exports_without_gui(saved_run, tmp_path, capsys):
    run, _, _ = saved_run
    assert main(["export", "--run", str(run), "--list-cases"]) == 0
    assert "sample\t" in capsys.readouterr().out
    assert len(export_cases(run)) == 1
    assert main(["export", "--run", str(run), "--case", "sample", "--output", str(tmp_path / "cli.sdf")]) == 0
