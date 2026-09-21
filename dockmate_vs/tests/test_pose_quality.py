from types import SimpleNamespace

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from dockmate_vs.analysis.pose_quality import (
    contact_fingerprint_from_interaction_sets,
    contact_similarity,
    summarize_posebusters_row,
)
from dockmate_vs.gui.app import DockMateVSApp


def test_posebusters_summary_reports_grouped_pass_failures():
    summary = summarize_posebusters_row({
        "mol_pred_loaded": True,
        "mol_cond_loaded": True,
        "sanitization": True,
        "inchi_convertible": True,
        "all_atoms_connected": True,
        "no_radicals": True,
        "bond_lengths": False,
        "bond_angles": True,
        "internal_steric_clash": True,
        "aromatic_ring_flatness": True,
        "internal_energy": False,
        "minimum_distance_to_protein": True,
        "volume_overlap_with_protein": True,
    })

    assert summary["posebusters_pass"] is False
    assert summary["ligand_geometry_pass"] is False
    assert summary["clash_pass"] is True
    assert summary["internal_energy_pass"] is False
    assert "bond lengths" in summary["posebusters_failed_checks"]
    assert "internal energy" in summary["posebusters_failed_checks"]


def test_plip_contact_fingerprint_and_similarity_are_type_aware():
    interactions = {
        "ZQX:L:1": SimpleNamespace(
            hydrophobic_contacts=[
                SimpleNamespace(reschain="A", restype="HIS", resnr=41)
            ],
            hbonds_ldon=[
                SimpleNamespace(reschain="A", restype="GLU", resnr=166)
            ],
            hbonds_pdon=[],
            water_bridges=[],
            saltbridge_lneg=[],
            saltbridge_pneg=[],
            pistacking=[],
            pication_laro=[],
            pication_paro=[],
            halogen_bonds=[],
            metal_complexes=[],
        )
    }

    native_contacts = contact_fingerprint_from_interaction_sets(interactions)
    pose_contacts = {"hydrophobic:A:HIS41", "hbond:A:GLY143"}
    similarity = contact_similarity(native_contacts, pose_contacts)

    assert native_contacts == {"hydrophobic:A:HIS41", "hbond:A:GLU166"}
    assert similarity["shared_contacts"] == ["hydrophobic:A:HIS41"]
    assert similarity["plip_similarity"] == pytest.approx(1 / 3)
    assert similarity["native_contact_recovery"] == pytest.approx(0.5)
    assert similarity["missing_native_contacts"] == ["hbond:A:GLU166"]
    assert similarity["new_pose_contacts"] == ["hbond:A:GLY143"]


def test_pose_quality_rows_render_in_screening_markdown():
    app = object.__new__(DockMateVSApp)
    markdown = app._summary_to_markdown({
        "total_cases": 1,
        "threshold": 2.0,
        "docking_completed": 1,
        "docking_failed": 0,
        "n_samples": 1,
        "screening_score_count": 1,
        "screening_unscored_count": 0,
        "screening_score_methods": ["vinardo (Smina score-only)"],
        "screening_score_direction": "lower",
        "screening_top_hits": [
            {
                "target_name": "Mpro",
                "pdb_id": "7ABC",
                "ligand": "LIG",
                "rank": 1,
                "compound": "hit_1",
                "score": -8.2,
                "score_source": "vinardo (Smina score-only)",
            }
        ],
        "screening_pose_quality_limit_per_structure": 5,
        "screening_pose_quality": [
            {
                "target_name": "Mpro",
                "pdb_id": "7ABC",
                "ligand": "LIG",
                "rank": 1,
                "compound": "hit_1",
                "score": -8.2,
                "posebusters_pass": False,
                "posebusters_failed_checks": ["volume overlap with protein"],
                "ligand_geometry_pass": True,
                "clash_pass": False,
                "internal_energy_pass": True,
                "plip_similarity": 0.5,
                "native_contact_recovery": 0.75,
            }
        ],
    })

    assert "## Pose Plausibility and Native-Contact Similarity" in markdown
    assert "volume overlap with protein" in markdown
    assert "| Mpro | 7ABC | LIG | 1 | hit_1 | -8.200 | Fail |" in markdown
    assert "| Pass | Fail | Pass | 0.50 | 75.0% |" in markdown


def test_posebusters_sdf_uses_prepared_topology_with_docked_coordinates(tmp_path):
    app = object.__new__(DockMateVSApp)
    case_dir = tmp_path / "case"
    prepared_dir = case_dir / "ligand_variants"
    quality_dir = case_dir / "pose_quality" / "hit_1"
    output_file = case_dir / "variants" / "ethanol_v1" / "docked.pdbqt"
    prepared_dir.mkdir(parents=True)
    quality_dir.mkdir(parents=True)
    output_file.parent.mkdir(parents=True)

    prepared = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    assert AllChem.EmbedMolecule(prepared, randomSeed=7) == 0
    prepared_sdf = prepared_dir / "ethanol_v1.sdf"
    app._write_ligand_sdf(prepared, prepared_sdf)

    docked_pose = Chem.RemoveHs(Chem.Mol(prepared), sanitize=False)
    docked_conf = docked_pose.GetConformer()
    for atom_index in range(docked_pose.GetNumAtoms()):
        position = docked_conf.GetAtomPosition(atom_index)
        position.x += 4.0
        position.y -= 2.0
        docked_conf.SetAtomPosition(atom_index, position)

    ligand_sdf, source, source_sdf, error = app._write_posebusters_ligand_sdf(
        output_file, docked_pose, quality_dir
    )

    assert source == "prepared_sdf"
    assert source_sdf == str(prepared_sdf)
    assert error is None
    assert ligand_sdf.name == "selected_pose_prepared_topology.sdf"

    written = app._load_first_sdf_mol(ligand_sdf)
    assert written is not None
    assert written.GetNumAtoms() == docked_pose.GetNumAtoms()
    assert written.GetNumBonds() == Chem.RemoveHs(prepared, sanitize=False).GetNumBonds()
    written_position = written.GetConformer().GetAtomPosition(0)
    docked_position = docked_pose.GetConformer().GetAtomPosition(0)
    assert written_position.x == pytest.approx(docked_position.x, abs=1e-4)
    assert written_position.y == pytest.approx(docked_position.y, abs=1e-4)


def test_posebusters_sdf_can_match_prepared_topology_by_heavy_atom_mcs(tmp_path):
    app = object.__new__(DockMateVSApp)
    case_dir = tmp_path / "case"
    prepared_dir = case_dir / "ligand_variants"
    quality_dir = case_dir / "pose_quality" / "hit_1"
    output_file = case_dir / "variants" / "phenol_v1" / "docked.pdbqt"
    prepared_dir.mkdir(parents=True)
    quality_dir.mkdir(parents=True)
    output_file.parent.mkdir(parents=True)

    prepared = Chem.AddHs(Chem.MolFromSmiles("Oc1ccccc1"))
    assert AllChem.EmbedMolecule(prepared, randomSeed=11) == 0
    prepared_sdf = prepared_dir / "phenol_v1.sdf"
    app._write_ligand_sdf(prepared, prepared_sdf)

    heavy = Chem.RemoveHs(prepared, sanitize=False)
    atom_order = [6, 0, 1, 2, 3, 4, 5]
    docked_pose = Chem.RenumberAtoms(heavy, atom_order)
    editable = Chem.RWMol(docked_pose)
    for atom in editable.GetAtoms():
        atom.SetIsAromatic(False)
    for bond in editable.GetBonds():
        bond.SetIsAromatic(False)
        bond.SetBondType(Chem.BondType.SINGLE)
    docked_pose = editable.GetMol()
    docked_conf = docked_pose.GetConformer()
    for atom_index in range(docked_pose.GetNumAtoms()):
        position = docked_conf.GetAtomPosition(atom_index)
        position.x += 7.0
        position.z -= 1.5
        docked_conf.SetAtomPosition(atom_index, position)

    ligand_sdf, source, source_sdf, error = app._write_posebusters_ligand_sdf(
        output_file, docked_pose, quality_dir
    )

    assert source == "prepared_sdf_mcs"
    assert source_sdf == str(prepared_sdf)
    assert error is None
    written = app._load_first_sdf_mol(ligand_sdf)
    assert written is not None
    assert written.GetNumAtoms() == docked_pose.GetNumAtoms()
    assert any(bond.GetIsAromatic() for bond in written.GetBonds())
    written_position = written.GetConformer().GetAtomPosition(0)
    docked_position = docked_pose.GetConformer().GetAtomPosition(0)
    assert written_position.x == pytest.approx(docked_position.x, abs=1e-4)
    assert written_position.z == pytest.approx(docked_position.z, abs=1e-4)
