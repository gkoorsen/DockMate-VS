"""Recover prepared-ligand SDF sidecars for legacy screening runs."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterable, List, Optional, Tuple, Union

from rdkit import Chem

from dockmate_vs.config.schema import LigandPreparationConfig
from dockmate_vs.preparation.ligand_cache import LigandCache


ProgressCallback = Callable[[int, int, str], None]
CancellationCallback = Callable[[], bool]


@dataclass(frozen=True)
class PreparedSdfBackfillTarget:
    """One missing SDF and the saved preparation artifact that validates it."""

    compound: str
    smiles: str
    variant_index: int
    pdbqt_path: Path
    sdf_path: Path
    case_id: str = ""


@dataclass
class PreparedSdfBackfillResult:
    """Outcome of a guarded prepared-SDF backfill."""

    written: List[Path] = field(default_factory=list)
    already_present: List[Path] = field(default_factory=list)
    errors: List[str] = field(default_factory=list)
    prepared_compounds: int = 0
    cancelled: bool = False


def _pdbqt_signature(text: str) -> Tuple[tuple, tuple]:
    """Return atom/torsion identity and coordinates from an Open Babel PDBQT."""
    topology = []
    coordinates = []
    for line in text.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            fields = line.split()
            if len(fields) < 10:
                raise ValueError("PDBQT contains a malformed atom record")
            try:
                xyz = (
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                )
            except (TypeError, ValueError):
                try:
                    xyz = tuple(float(value) for value in fields[5:8])
                except (TypeError, ValueError) as fallback_exc:
                    raise ValueError(
                        "PDBQT contains invalid atom coordinates"
                    ) from fallback_exc
            topology.append((fields[0], fields[1], fields[2], fields[-1]))
            coordinates.append(xyz)
        elif line.startswith(("ROOT", "ENDROOT", "BRANCH", "ENDBRANCH", "TORSDOF")):
            topology.append(tuple(line.split()))
    if not coordinates:
        raise ValueError("PDBQT contains no atoms")
    return tuple(topology), tuple(coordinates)


def validate_prepared_variant(
    regenerated_pdbqt: Optional[str],
    saved_pdbqt: Path,
    coordinate_tolerance: float = 0.02,
) -> Tuple[bool, Optional[str]]:
    """Verify that a regenerated variant is the exact saved preparation variant."""
    if not regenerated_pdbqt:
        return False, "regenerated variant has no PDBQT representation"
    try:
        saved_text = Path(saved_pdbqt).read_text()
        saved_topology, saved_coordinates = _pdbqt_signature(saved_text)
        new_topology, new_coordinates = _pdbqt_signature(regenerated_pdbqt)
    except (OSError, ValueError) as exc:
        return False, str(exc)

    if saved_topology != new_topology:
        return False, "atom types or rotatable-bond topology differ from the saved PDBQT"
    if len(saved_coordinates) != len(new_coordinates):
        return False, "atom count differs from the saved PDBQT"

    maximum_delta = max(
        abs(saved_value - new_value)
        for saved_xyz, new_xyz in zip(saved_coordinates, new_coordinates)
        for saved_value, new_value in zip(saved_xyz, new_xyz)
    )
    if maximum_delta > coordinate_tolerance:
        return (
            False,
            f"prepared coordinates differ from the saved PDBQT by up to "
            f"{maximum_delta:.3f} A",
        )
    return True, None


def _contains_metal(smiles: str) -> bool:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return False
    metal_symbols = {
        "LI", "NA", "K", "RB", "CS", "MG", "CA", "SR", "BA",
        "MN", "FE", "CO", "NI", "CU", "ZN", "MO", "W", "AG",
        "AU", "CD", "HG", "PT", "PD", "AL", "GA", "IN", "SN",
        "PB", "BI",
    }
    return any(atom.GetSymbol().upper() in metal_symbols for atom in molecule.GetAtoms())


def backfill_prepared_ligand_sdfs(
    targets: Iterable[PreparedSdfBackfillTarget],
    config: LigandPreparationConfig,
    *,
    n_cpus: int = 1,
    cache_dir: Union[Path, str] = "ligand_cache",
    progress: Optional[ProgressCallback] = None,
    is_cancelled: Optional[CancellationCallback] = None,
) -> PreparedSdfBackfillResult:
    """Regenerate and validate only missing SDFs, without altering docking output."""
    result = PreparedSdfBackfillResult()
    grouped = {}
    for target in targets:
        key = (target.compound, target.smiles)
        grouped.setdefault(key, []).append(target)

    cache = LigandCache(str(cache_dir))
    total = len(grouped)
    for group_number, ((compound, smiles), group_targets) in enumerate(grouped.items(), 1):
        if is_cancelled and is_cancelled():
            result.cancelled = True
            break
        if progress:
            progress(group_number - 1, total, f"Preparing {compound}")

        try:
            prepared_ligands = cache.get(
                smiles=smiles,
                mol_id=compound,
                config=config,
                n_cpus=max(1, int(n_cpus)),
                enumerate_states=not _contains_metal(smiles),
            )
            result.prepared_compounds += 1
        except Exception as exc:
            result.errors.append(f"{compound}: ligand preparation failed: {exc}")
            if progress:
                progress(group_number, total, f"Could not prepare {compound}")
            continue

        for target in group_targets:
            if target.sdf_path.exists():
                result.already_present.append(target.sdf_path)
                continue
            if target.variant_index < 1 or target.variant_index > len(prepared_ligands):
                result.errors.append(
                    f"{target.case_id or compound}: variant v{target.variant_index} is not "
                    f"among the {len(prepared_ligands)} regenerated variants"
                )
                continue

            prepared = prepared_ligands[target.variant_index - 1]
            valid, reason = validate_prepared_variant(prepared.pdbqt, target.pdbqt_path)
            if not valid:
                result.errors.append(
                    f"{target.case_id or compound} v{target.variant_index}: {reason}"
                )
                continue

            temporary = target.sdf_path.with_name(f".{target.sdf_path.name}.tmp")
            try:
                target.sdf_path.parent.mkdir(parents=True, exist_ok=True)
                prepared.to_sdf_file(str(temporary))
                supplier = Chem.SDMolSupplier(str(temporary), sanitize=False, removeHs=False)
                if not supplier or supplier[0] is None:
                    raise ValueError("regenerated SDF could not be read")
                temporary.replace(target.sdf_path)
                result.written.append(target.sdf_path)
            except Exception as exc:
                temporary.unlink(missing_ok=True)
                result.errors.append(
                    f"{target.case_id or compound} v{target.variant_index}: "
                    f"could not write SDF: {exc}"
                )

        if progress:
            progress(group_number, total, f"Prepared {compound}")

    return result
