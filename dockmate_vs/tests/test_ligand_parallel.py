import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

import dockmate_vs.preparation.ligand as ligand_module
from dockmate_vs.config.schema import LigandPreparationConfig
from dockmate_vs.preparation.ligand import LigandPreparation


class _InlinePool:
    """Exercise the batching path without spawning processes in the test suite."""

    def __init__(self, processes):
        self.processes = processes

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        return False

    def starmap(self, function, batches, chunksize=1):
        return [function(*batch) for batch in batches]


def test_parallel_rmsd_batches_match_serial_calculation(monkeypatch):
    mol = Chem.AddHs(Chem.MolFromSmiles("CCCC"))
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    conf_ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=11, params=params))
    assert len(conf_ids) == 11  # 55 pairs enters the parallel batching branch.

    preparation = LigandPreparation(LigandPreparationConfig(), n_cpus=2)
    serial = preparation._calculate_rmsd_matrix_serial(mol, conf_ids)

    monkeypatch.setattr(ligand_module, "Pool", _InlinePool)
    parallel = preparation._calculate_rmsd_matrix_parallel(mol, conf_ids, n_workers=2)

    # RDKit's binary conformer round-trip changes the final RMSD digits slightly.
    np.testing.assert_allclose(parallel, serial, rtol=0.0, atol=1e-6)
    np.testing.assert_allclose(parallel, parallel.T, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(np.diag(parallel), 0.0, rtol=0.0, atol=0.0)
