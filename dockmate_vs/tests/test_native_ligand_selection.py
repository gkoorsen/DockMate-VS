import numpy as np
from pathlib import Path
import pytest
from Bio.PDB import PDBParser

from dockmate_vs.adaptive_docking import AdaptiveDockingPipeline
from dockmate_vs.binding_site.cocrystal import BindingSiteDefinition, BindingSiteError
from dockmate_vs.config.schema import ProteinPreparationConfig
from dockmate_vs.preparation.protein import ProteinPreparation


def test_1x8x_reference_is_exact_nonpolymer_tyr_a952(tmp_path):
    fixture = Path(__file__).parent / "fixtures/1x8x_tyr_excerpt.pdb"
    pipeline = object.__new__(AdaptiveDockingPipeline)
    pipeline.output_dir = tmp_path
    expected = [line for line in fixture.read_text().splitlines()
                if line.startswith("HETATM") and line[21:27].strip() == "A 952"]
    reference = pipeline._extract_crystal_ligand(fixture, "TYR", "A")
    assert len(expected) == 13
    assert reference.read_text().splitlines() == expected
    xyz = np.array([[float(line[s:s+8]) for s in (30, 38, 46)] for line in expected])
    site = BindingSiteDefinition(margin=4).from_cocrystal(fixture, "TYR", "A")
    np.testing.assert_allclose(site.center, xyz.mean(axis=0), atol=2e-6)
    np.testing.assert_allclose(site.size, np.ptp(xyz, axis=0)+8, atol=2e-6)


def test_1x8x_removal_preserves_all_108_protein_tyr_atom_identities_and_coordinates():
    fixture = Path(__file__).parent / "fixtures/1x8x_tyr_excerpt.pdb"
    structure = PDBParser(QUIET=True).get_structure("1X8X", fixture)

    def identities(value):
        return {(res.get_parent().id, res.id, atom.name, atom.element): atom.coord.copy()
                for res in value.get_residues() if res.id[0] == " "
                for atom in res.get_atoms()}

    before = identities(structure)
    assert len(before) == 108
    assert {key[1][1] for key in before} == {37, 127, 175, 181, 190, 254, 257, 269, 301}
    prepared = ProteinPreparation(ProteinPreparationConfig())._remove_hetero_residues(structure, ["TYR"])
    after = identities(prepared)
    assert after.keys() == before.keys()
    for key in before:
        np.testing.assert_array_equal(after[key], before[key])
    assert len(list(prepared.get_atoms())) == 108
    assert all(res.id[0] == " " for res in prepared.get_residues())


def atom(serial, name, kind="HETATM", residue="TYR", chain="A", resid=952, x=10., element="C"):
    return (f"{kind:<6}{serial:5d} {name:^4} {residue:3} {chain}{resid:4d}    "
            f"{x:8.3f}{0.:8.3f}{0.:8.3f}{1.:6.2f}{20.:6.2f}          {element:>2}\n")


def source(tmp_path):
    path = tmp_path / "tyr.pdb"
    path.write_text(atom(1, "CA", kind="ATOM", resid=37, x=-40.) +
                    atom(2, "CB", kind="ATOM", resid=37, x=-39.) +
                    atom(3, "C1") + atom(4, "O1", x=12., element="O") +
                    atom(5, "H1", x=11., element="H") +
                    atom(6, "C1", resid=953, x=60.) +
                    atom(7, "C1", chain="B", x=80.) + "END\n")
    return path


def test_reference_and_box_select_same_single_hetero_tyr(tmp_path):
    path = source(tmp_path)
    pipeline = object.__new__(AdaptiveDockingPipeline)
    pipeline.output_dir = tmp_path
    reference = pipeline._extract_crystal_ligand(path, "TYR", "A")
    lines = reference.read_text().splitlines()
    assert len(lines) == 3
    assert all(line.startswith("HETATM") and line[22:26].strip() == "952" for line in lines)
    site = BindingSiteDefinition(margin=4).from_cocrystal(path, "TYR", "A")
    np.testing.assert_allclose(site.center, [11., 0., 0.])
    np.testing.assert_allclose(site.size, [10., 8., 8.])
    assert site.ligand_atoms == 3
    assert site.ligand_heavy_atoms == 2
    assert site.to_dict()["ligand_atom_count_kind"] == "total_including_hydrogens"


def test_explicit_missing_chain_does_not_fall_back(tmp_path):
    path = source(tmp_path)
    with pytest.raises(BindingSiteError, match="Available chains for TYR: A, B"):
        BindingSiteDefinition().from_cocrystal(path, "TYR", "C")


def test_protein_tyr_is_not_a_free_native_ligand(tmp_path):
    path = tmp_path / "protein.pdb"
    path.write_text(atom(1, "CA", kind="ATOM", resid=37))
    with pytest.raises(BindingSiteError, match="protein amino-acid"):
        BindingSiteDefinition().from_cocrystal(path, "TYR", "A")


def test_auto_chain_detection_prefers_free_ligand_over_protein(tmp_path):
    from dockmate_vs.headless import HeadlessDockMateRunner

    path = tmp_path / "chain.pdb"
    path.write_text(atom(1, "CA", kind="ATOM", resid=37, chain="A") +
                    atom(2, "C1", chain="B"))
    assert HeadlessDockMateRunner()._detect_ligand_chain(path, "TYR") == "B"


def test_ligand_removal_keeps_protein_residue_of_same_name(tmp_path):
    structure = PDBParser(QUIET=True).get_structure("test", source(tmp_path))
    prep = ProteinPreparation(ProteinPreparationConfig(remove_hetero_residues=["TYR"]))
    result = prep._remove_hetero_residues(structure, ["tyr"])
    residues = list(result.get_residues())
    assert len(residues) == 1
    assert residues[0].id == (" ", 37, " ")
    assert len(list(residues[0].get_atoms())) == 2


def test_ligand_is_removed_before_fixer_can_rewrite_hetatm(monkeypatch, tmp_path):
    prep = ProteinPreparation(ProteinPreparationConfig(
        remove_hetero_residues=["TYR"], optimize_hbonds=False, validate_structure=False))
    seen = []

    def fix(path, structure):
        text = open(path).read()
        assert "HETATM" not in text
        assert "TYR A  37" in text
        seen.append(path)
        return structure

    monkeypatch.setattr(prep, "_fix_structure", fix)
    monkeypatch.setattr(prep, "_predict_pka", lambda structure: {})
    monkeypatch.setattr(prep, "_add_hydrogens", lambda structure, pka: structure)
    result = prep.prepare_receptor(str(source(tmp_path)))
    assert seen
    assert len(list(result.structure.get_residues())) == 1
