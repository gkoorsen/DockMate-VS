from types import SimpleNamespace

import dockmate_vs.preparation.protein as protein_module
from dockmate_vs.preparation.protein import (
    ProteinPreparation,
    RECEPTOR_PREPARATION_SEED,
)


def test_missing_atom_placement_uses_fixed_seed(monkeypatch, tmp_path):
    observed = {}

    class FakeFixer:
        topology = object()
        positions = object()

        def __init__(self, filename):
            observed["filename"] = filename

        def findMissingResidues(self):
            pass

        def findMissingAtoms(self):
            pass

        def addMissingAtoms(self, seed=None):
            observed["seed"] = seed

    class FakePDBFile:
        @staticmethod
        def writeFile(topology, positions, handle):
            handle.write("HEADER\n")

    class FakeParser:
        def get_structure(self, name, path):
            observed["parsed"] = (name, path)
            return "fixed-structure"

    monkeypatch.setattr(protein_module, "PDBFixer", FakeFixer)
    monkeypatch.setattr(protein_module, "PDBFile", FakePDBFile)
    preparation = object.__new__(ProteinPreparation)
    preparation.config = SimpleNamespace(
        add_missing_atoms=True,
        model_missing_loops=False,
    )
    preparation.parser = FakeParser()
    source = tmp_path / "receptor.pdb"
    source.write_text("HEADER\n")

    fixed = preparation._fix_structure(str(source), object())

    assert fixed == "fixed-structure"
    assert observed["seed"] == RECEPTOR_PREPARATION_SEED
