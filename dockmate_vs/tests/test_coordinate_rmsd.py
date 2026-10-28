from types import SimpleNamespace

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import rdFMCS

from dockmate_vs.utils.rmsd import coordinate_rmsd


def molecule(smiles, xyz):
    mol = Chem.MolFromSmiles(smiles)
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, point in enumerate(xyz):
        conf.SetAtomPosition(i, point)
    mol.AddConformer(conf)
    return mol


@pytest.mark.parametrize("symmetry", [True, False])
def test_vectorized_rmsd_matches_original_mapping_pairs_without_fitting(symmetry):
    rng = np.random.default_rng(21)
    ref = molecule("CC(C)(C)O", rng.normal(size=(5, 3)))
    probe = molecule("CC(C)(C)O", rng.normal(size=(5, 3)) + [8., 2., 0.])
    mcs = rdFMCS.FindMCS([ref, probe], ringMatchesRingOnly=False, completeRingsOnly=False,
                        matchValences=False, bondCompare=rdFMCS.BondCompare.CompareAny,
                        atomCompare=rdFMCS.AtomCompare.CompareElements)
    pattern = Chem.MolFromSmarts(mcs.smartsString)
    refs = ref.GetSubstructMatches(pattern, uniquify=not symmetry, maxMatches=1000)
    probes = probe.GetSubstructMatches(pattern, uniquify=not symmetry, maxMatches=1000)
    a = ref.GetConformer().GetPositions()
    b = probe.GetConformer().GetPositions()
    expected = min(np.sqrt(np.mean(np.sum((a[list(x)] - b[list(y)]) ** 2, axis=1)))
                   for x in refs for y in probes)
    assert coordinate_rmsd(ref, probe, use_symmetry=symmetry) == pytest.approx(expected, abs=1e-12)
    np.testing.assert_array_equal(ref.GetConformer().GetPositions(), a)
    np.testing.assert_array_equal(probe.GetConformer().GetPositions(), b)


def test_rmsd_timeout_is_explicit_and_not_a_partial_match(monkeypatch):
    mol = molecule("CO", [[0., 0., 0.], [1., 0., 0.]])
    monkeypatch.setattr(rdFMCS, "FindMCS", lambda *a, **kw: SimpleNamespace(canceled=True))
    with pytest.raises(ValueError, match="60-second limit"):
        coordinate_rmsd(mol, mol)


def test_rmsd_atom_count_mismatch_is_not_numeric():
    ref = molecule("CO", [[0., 0., 0.], [1., 0., 0.]])
    probe = molecule("C", [[0., 0., 0.]])
    with pytest.raises(ValueError, match="Heavy-atom counts differ"):
        coordinate_rmsd(ref, probe)
