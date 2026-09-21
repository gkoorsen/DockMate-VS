"""Pose plausibility and protein-ligand contact similarity helpers."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Mapping, Optional, Sequence, Set

import pandas as pd


POSEBUSTERS_LABELS = {
    "mol_pred_loaded": "predicted ligand loaded",
    "mol_cond_loaded": "receptor loaded",
    "sanitization": "RDKit sanitization",
    "inchi_convertible": "InChI conversion",
    "all_atoms_connected": "all atoms connected",
    "no_radicals": "no radicals",
    "molecular_formula": "molecular formula",
    "molecular_bonds": "molecular bonds",
    "double_bond_stereochemistry": "double-bond stereochemistry",
    "tetrahedral_chirality": "tetrahedral chirality",
    "bond_lengths": "bond lengths",
    "bond_angles": "bond angles",
    "internal_steric_clash": "internal steric clash",
    "aromatic_ring_flatness": "aromatic ring flatness",
    "non-aromatic_ring_non-flatness": "non-aromatic ring non-flatness",
    "double_bond_flatness": "double-bond flatness",
    "internal_energy": "internal energy",
    "protein-ligand_maximum_distance": "protein-ligand maximum distance",
    "minimum_distance_to_protein": "minimum distance to protein",
    "volume_overlap_with_protein": "volume overlap with protein",
}

LIGAND_GEOMETRY_CHECKS = (
    "sanitization",
    "inchi_convertible",
    "all_atoms_connected",
    "no_radicals",
    "bond_lengths",
    "bond_angles",
    "internal_steric_clash",
    "aromatic_ring_flatness",
    "non-aromatic_ring_non-flatness",
    "double_bond_flatness",
)
CLASH_CHECKS = (
    "internal_steric_clash",
    "minimum_distance_to_protein",
    "volume_overlap_with_protein",
)

CONTACT_ATTRIBUTES = {
    "hydrophobic": ("hydrophobic_contacts",),
    "hbond": ("hbonds_ldon", "hbonds_pdon"),
    "water_bridge": ("water_bridges",),
    "salt_bridge": ("saltbridge_lneg", "saltbridge_pneg"),
    "pi_stack": ("pistacking",),
    "pi_cation": ("pication_laro", "pication_paro"),
    "halogen": ("halogen_bonds",),
    "metal": ("metal_complexes",),
}


def _bool_or_none(value: object) -> Optional[bool]:
    if value is None:
        return None
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass
    if isinstance(value, str):
        lowered = value.strip().lower()
        if lowered in {"true", "1", "yes"}:
            return True
        if lowered in {"false", "0", "no"}:
            return False
        return None
    if isinstance(value, bool) or type(value).__name__ == "bool_":
        return bool(value)
    if isinstance(value, (int, float)):
        return bool(value)
    return None


def _present_checks(row: Mapping[str, object], checks: Sequence[str]) -> list[str]:
    return [check for check in checks if _bool_or_none(row.get(check)) is not None]


def _checks_pass(row: Mapping[str, object], checks: Sequence[str]) -> Optional[bool]:
    present = _present_checks(row, checks)
    if not present:
        return None
    return all(_bool_or_none(row.get(check)) is True for check in present)


def _posebusters_check_columns(row: Mapping[str, object]) -> list[str]:
    excluded = {"file", "molecule", "position", "mol_true_loaded"}
    columns = []
    for key, value in row.items():
        if key in excluded or str(key).startswith("rmsd_"):
            continue
        if _bool_or_none(value) is not None:
            columns.append(str(key))
    return columns


def posebusters_unavailable(error: str) -> dict:
    return {
        "posebusters_available": False,
        "posebusters_pass": None,
        "posebusters_passed_checks": 0,
        "posebusters_total_checks": 0,
        "posebusters_failed_checks": [],
        "posebusters_error": error,
        "clash_pass": None,
        "ligand_geometry_pass": None,
        "internal_energy_pass": None,
        "volume_overlap_pass": None,
    }


def summarize_posebusters_row(row: Mapping[str, object]) -> dict:
    """Summarize a PoseBusters dataframe row into stable report fields."""
    checks = _posebusters_check_columns(row)
    failed = [
        POSEBUSTERS_LABELS.get(check, check.replace("_", " "))
        for check in checks
        if _bool_or_none(row.get(check)) is False
    ]
    passed = sum(1 for check in checks if _bool_or_none(row.get(check)) is True)
    return {
        "posebusters_available": True,
        "posebusters_pass": not failed if checks else None,
        "posebusters_passed_checks": passed,
        "posebusters_total_checks": len(checks),
        "posebusters_failed_checks": failed,
        "posebusters_error": None,
        "clash_pass": _checks_pass(row, CLASH_CHECKS),
        "ligand_geometry_pass": _checks_pass(row, LIGAND_GEOMETRY_CHECKS),
        "internal_energy_pass": _checks_pass(row, ("internal_energy",)),
        "volume_overlap_pass": _checks_pass(row, ("volume_overlap_with_protein",)),
    }


def run_posebusters(ligand_sdf: Path, receptor_pdb: Path) -> dict:
    """Run PoseBusters in docking mode for one selected ligand pose."""
    try:
        from posebusters import PoseBusters  # type: ignore
    except Exception as exc:
        return posebusters_unavailable(f"PoseBusters is not installed: {exc}")

    try:
        buster = PoseBusters(config="dock", top_n=1, max_workers=0)
        frame = buster.bust([Path(ligand_sdf)], mol_cond=Path(receptor_pdb))
        if frame.empty:
            return posebusters_unavailable("PoseBusters returned no rows")
        return summarize_posebusters_row(frame.iloc[0].to_dict())
    except Exception as exc:
        return posebusters_unavailable(str(exc))


def _interaction_residue(interaction: object) -> Optional[str]:
    resnr = getattr(interaction, "resnr", None)
    if resnr is None:
        return None
    restype = str(getattr(interaction, "restype", "") or "UNK").strip().upper()
    chain = str(getattr(interaction, "reschain", "") or "").strip()
    return f"{chain}:{restype}{resnr}"


def contact_fingerprint_from_interaction_sets(interaction_sets: Mapping[object, object]) -> Set[str]:
    """Convert PLIP interaction sets into type-aware residue contact keys."""
    contacts: Set[str] = set()
    for interaction_set in interaction_sets.values():
        for contact_type, attributes in CONTACT_ATTRIBUTES.items():
            for attribute in attributes:
                for interaction in getattr(interaction_set, attribute, []) or []:
                    residue = _interaction_residue(interaction)
                    if residue:
                        contacts.add(f"{contact_type}:{residue}")
    return contacts


def plip_unavailable(error: str) -> dict:
    return {
        "plip_available": False,
        "plip_contacts": [],
        "plip_contact_count": 0,
        "plip_error": error,
    }


def plip_contact_fingerprint(complex_pdb: Path) -> dict:
    """Run PLIP for one protein-ligand complex and return contact keys."""
    try:
        from plip.structure.preparation import PDBComplex  # type: ignore
    except Exception as exc:
        return plip_unavailable(f"PLIP is not installed: {exc}")

    try:
        complex_obj = PDBComplex()
        complex_obj.load_pdb(str(complex_pdb))
        if hasattr(complex_obj, "analyze"):
            complex_obj.analyze()
        else:
            for ligand in getattr(complex_obj, "ligands", []) or []:
                complex_obj.characterize_complex(ligand)
        contacts = sorted(
            contact_fingerprint_from_interaction_sets(complex_obj.interaction_sets)
        )
        return {
            "plip_available": True,
            "plip_contacts": contacts,
            "plip_contact_count": len(contacts),
            "plip_error": None,
        }
    except Exception as exc:
        return plip_unavailable(str(exc))


def contact_similarity(reference_contacts: Iterable[str], pose_contacts: Iterable[str]) -> dict:
    reference = set(reference_contacts)
    pose = set(pose_contacts)
    shared = reference & pose
    union = reference | pose
    return {
        "plip_similarity": (len(shared) / len(union)) if union else None,
        "native_contact_recovery": (len(shared) / len(reference)) if reference else None,
        "shared_contacts": sorted(shared),
        "missing_native_contacts": sorted(reference - pose),
        "new_pose_contacts": sorted(pose - reference),
        "native_contact_count": len(reference),
        "pose_contact_count": len(pose),
        "shared_contact_count": len(shared),
    }
