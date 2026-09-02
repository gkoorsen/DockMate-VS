from dockmate_vs.config.schema import LigandPreparationConfig
from dockmate_vs.preparation.ligand_cache import LigandCache


def test_cache_key_covers_complete_configuration_and_state_enumeration(tmp_path):
    cache = LigandCache(tmp_path / "ligands")
    base = LigandPreparationConfig(max_tautomers=2, max_conformers=2)
    wider_window = LigandPreparationConfig(
        max_tautomers=2,
        max_conformers=2,
        energy_window=20.0,
    )

    baseline = cache._get_cache_key("CCO", base, enumerate_states=True)

    assert cache._get_cache_key("CCO", base, enumerate_states=False) != baseline
    assert cache._get_cache_key("CCO", wider_window, enumerate_states=True) != baseline


def test_cache_key_changes_when_preparation_tool_identity_changes(tmp_path):
    cache = LigandCache(tmp_path / "ligands")
    config = LigandPreparationConfig(max_tautomers=2, max_conformers=2)
    original = cache._get_cache_key("CCO", config)

    cache._tool_signature = {"rdkit_version": "different", "obabel": None}

    assert cache._get_cache_key("CCO", config) != original
