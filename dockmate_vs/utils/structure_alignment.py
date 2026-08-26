"""
Structure superposition and binding-site grid derivation for ensemble docking.

Why this module exists
----------------------
Ensemble docking passes ONE grid box to every receptor. That is only valid if the
receptors share a coordinate frame. Structures downloaded from the PDB generally do
NOT: they sit in whatever frame their crystal form dictates. Averaging raw ligand
coordinates across such an ensemble produces a "consensus centre" in bulk solvent,
and docking then completes normally while returning scores from an empty box.

Measured on a real three-structure CK2-alpha ensemble (6RB1, 3BQC, 7ZY0) the pairwise
ligand-centroid distances were 396.5, 387.1 and 88.2 A. Every structure fell outside
a 22 A box centred on their mean.

This module provides:
  * superpose_ensemble()    - bring all members into a common frame (CA/P superposition)
  * site_center()           - per-structure binding-site centre from ONE ligand copy
  * derive_grid()           - consensus centre + box size, with a spread check
  * validate_grid()         - hard gate: does every member's site fall inside the box?

Requires gemmi (already a project dependency for structure IO).
"""

from __future__ import annotations

import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

try:
    import gemmi
except ImportError as exc:  # pragma: no cover
    raise ImportError(
        "structure_alignment requires gemmi. Install with: pip install gemmi"
    ) from exc


# Components that are crystallisation additives, cryoprotectants, ions, buffers or
# glycosylation - never the pocket ligand. Pooling these into the "ligand" centroid is
# a real failure mode: in 6TGU, ethylene glycol (EDO, 70 atoms) outnumbers the actual
# inhibitor (N92, 66 atoms) and drags the centre ~7.2 A off-site; in 1MX1 the NAG/SIA/NDG
# glycans pull it off the steroid pocket.
NON_LIGAND_COMPONENTS = frozenset({
    # water / solvent
    "HOH", "DOD", "WAT", "SOL", "TIP3", "TIP", "H2O",
    # cryo / precipitant / buffer
    "GOL", "EDO", "PEG", "PG4", "PGE", "1PE", "2PE", "P6G", "P33", "PE4", "PGO",
    "MPD", "MRD", "BU3", "DMS", "DMF", "IPA", "ACN", "TRS", "EPE", "MES", "IMD",
    "BME", "DTT", "TCE", "CAC", "BOG", "LDA", "LMT", "C8E", "OCT", "HEZ",
    "ACT", "ACY", "FMT", "EOH", "MOH", "SIN", "CIT", "FLC", "TAR", "MLA", "MLI",
    "SCN", "AZI", "NH4", "PPV", "UNX", "UNL",
    # ions
    "CL", "BR", "IOD", "I", "F", "NA", "K", "LI", "CS", "RB",
    "MG", "CA", "ZN", "MN", "FE", "FE2", "FES", "NI", "CO", "CD", "CU", "CU1",
    "HG", "AU", "AG", "PT", "PB", "BA", "SR", "AL", "GA", "YB", "SM", "EU", "TB",
    "SO4", "PO4", "PI", "NO3", "CO3", "BO3", "WO4", "MOO", "VO4", "AF3", "BEF",
    # glycans / sugars (surface glycosylation, not the pocket)
    "NAG", "NDG", "NGA", "A2G", "BMA", "MAN", "BGC", "GLC", "GAL", "GLA",
    "FUC", "FUL", "SIA", "SLB", "XYS", "XYP", "RAM", "ARA", "ARB", "RIB",
    "SUC", "TRE", "MAL", "LAT", "GCU", "ADA", "IDS", "SGN",
    # common modified residues / lipids that are not the ligand of interest
    "MSE", "PCA", "ACE", "NH2", "SEP", "TPO", "PTR", "MLY", "MLZ",
    "PLM", "MYR", "STE", "OLA", "PEE", "PC1", "CDL", "CLR", "HEM", "HEC",
})

# Cofactors: real ligands, but if the screen targets an inhibitor site these are usually
# the wrong centre to average in. Excluded by default; pass include_cofactors=True to keep.
COFACTOR_COMPONENTS = frozenset({
    "ATP", "ADP", "AMP", "ANP", "ACP", "AGS", "APC", "ADX",
    "GTP", "GDP", "GNP", "GSP", "GTG",
    "UTP", "UDP", "CTP", "CDP", "TTP", "TDP",
    "NAD", "NAI", "NAP", "NDP", "NAX", "FAD", "FDA", "FMN", "FNR",
    "SAM", "SAH", "COA", "ACO", "COO", "TPP", "TDT", "PLP", "BTN",
    "HEM", "HEC", "SF4", "F3S", "MGD",
})


@dataclass
class SiteCenter:
    """Binding-site centre extracted from one structure."""
    path: Path
    center: np.ndarray
    comp_id: str
    chain: str
    n_atoms: int
    n_copies: int
    extent: np.ndarray = field(default_factory=lambda: np.zeros(3))


@dataclass
class GridSpec:
    """Grid box derived from a co-framed ensemble."""
    center: List[float]
    size: List[float]
    max_deviation: float
    sites: List[SiteCenter]
    spread_std: List[float]

    def as_grid_params(self) -> Dict:
        """Dict in the shape the docking adapters expect."""
        return {"center": list(self.center), "size": list(self.size)}


# --------------------------------------------------------------------------- #
# Component / residue selection
# --------------------------------------------------------------------------- #

def _is_candidate_ligand(comp_id: str, include_cofactors: bool = False) -> bool:
    cid = comp_id.strip().upper()
    if not cid or cid in NON_LIGAND_COMPONENTS:
        return False
    if not include_cofactors and cid in COFACTOR_COMPONENTS:
        return False
    return True


def _comp_matches(res_name: str, wanted: str) -> bool:
    """
    Compare component ids tolerantly.

    The PDB extended CCD to five characters in 2023 (e.g. A1EIY). Legacy PDB format
    only holds three, so writing such a structure truncates the name to A1E. Anything
    that matches on residue name must allow for that or it silently finds nothing.
    """
    a = res_name.strip().upper()
    b = wanted.strip().upper()
    return a == b or a == b[:3] or b == a[:3]


def _read(path: Path) -> gemmi.Structure:
    st = gemmi.read_structure(str(path))
    st.setup_entities()
    return st


def _polymer_for_chain(st: gemmi.Structure, chain_name: Optional[str] = None):
    model = st[0]
    if chain_name is not None:
        for ch in model:
            if ch.name == chain_name:
                poly = ch.get_polymer()
                if len(poly):
                    return poly
    for ch in model:
        poly = ch.get_polymer()
        if len(poly):
            return poly
    return None


# --------------------------------------------------------------------------- #
# Site centre
# --------------------------------------------------------------------------- #

def site_center(
    path: str | Path,
    comp_id: Optional[str] = None,
    include_cofactors: bool = False,
    near: Optional[Sequence[float]] = None,
) -> Optional[SiteCenter]:
    """
    Centre of the binding-site ligand in one structure.

    Unlike averaging every HETATM, this selects a SINGLE residue copy:
      * if `comp_id` is given, that component (tolerant of 3-char truncation);
      * otherwise the largest non-additive, non-cofactor component present.
    When several copies exist (multimers), `near` picks the copy closest to a
    reference point; without `near` the first copy is used.

    Returns None when no candidate ligand is present (apo structure).
    """
    path = Path(path)
    st = _read(path)
    model = st[0]

    copies: List[Tuple[str, str, np.ndarray]] = []
    for ch in model:
        for res in ch:
            name = res.name.strip().upper()
            if comp_id is not None:
                if not _comp_matches(name, comp_id):
                    continue
            else:
                if res.het_flag != "H":
                    continue
                if not _is_candidate_ligand(name, include_cofactors):
                    continue
            pts = np.array([[a.pos.x, a.pos.y, a.pos.z] for a in res], dtype=float)
            if len(pts):
                copies.append((name, ch.name, pts))

    if not copies:
        return None

    if comp_id is None:
        # largest component wins; ties broken by first occurrence
        best_name = max({n for n, _, _ in copies},
                        key=lambda n: sum(len(p) for m, _, p in copies if m == n))
        copies = [c for c in copies if c[0] == best_name]

    if near is not None and len(copies) > 1:
        ref = np.asarray(near, dtype=float)
        chosen = min(copies, key=lambda c: float(np.linalg.norm(c[2].mean(axis=0) - ref)))
    else:
        chosen = copies[0]

    name, chain, pts = chosen
    return SiteCenter(
        path=path,
        center=pts.mean(axis=0),
        comp_id=name,
        chain=chain,
        n_atoms=len(pts),
        n_copies=len(copies),
        extent=np.ptp(pts, axis=0),
    )


# --------------------------------------------------------------------------- #
# Superposition
# --------------------------------------------------------------------------- #

def superpose_ensemble(
    structures: Sequence[str | Path],
    out_dir: str | Path,
    reference: Optional[str | Path] = None,
    comp_ids: Optional[Dict[str, str]] = None,
    include_cofactors: bool = False,
) -> Tuple[List[Path], List[Dict]]:
    """
    Superpose every member onto a reference and write aligned copies to `out_dir`.

    Superposition uses CA (protein) / P (nucleic) atoms of the chain carrying the
    binding-site ligand, so the pocket - not the crystal packing - defines the frame.
    The reference defaults to the first structure; it is copied through unchanged so
    the output set is self-consistent.

    Returns (aligned_paths, report_rows). Report rows carry rmsd, n_aligned and the
    chain used, for provenance.
    """
    structures = [Path(s) for s in structures]
    if not structures:
        raise ValueError("no structures given")

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    ref_path = Path(reference) if reference is not None else structures[0]
    comp_ids = comp_ids or {}

    def _ligand_chain(p: Path) -> Optional[str]:
        sc = site_center(p, comp_ids.get(p.stem), include_cofactors)
        return sc.chain if sc else None

    ref_st = _read(ref_path)
    ref_poly = _polymer_for_chain(ref_st, _ligand_chain(ref_path))
    if ref_poly is None or not len(ref_poly):
        raise ValueError(f"reference {ref_path.name} has no polymer chain to align on")

    aligned: List[Path] = []
    report: List[Dict] = []

    for path in structures:
        out_path = out_dir / f"{path.stem}_aligned.pdb"

        if path.resolve() == ref_path.resolve():
            st = _read(path)
            st.setup_entities()
            st.write_pdb(str(out_path))
            aligned.append(out_path)
            report.append(dict(structure=path.name, reference=ref_path.name,
                               rmsd=0.0, n_aligned=len(ref_poly),
                               chain=_ligand_chain(path) or "", is_reference=True))
            continue

        st = _read(path)
        poly = _polymer_for_chain(st, _ligand_chain(path))
        if poly is None or not len(poly):
            raise ValueError(f"{path.name}: no polymer chain to align on")

        sup = gemmi.calculate_superposition(
            ref_poly, poly, gemmi.PolymerType.PeptideL, gemmi.SupSelect.CaP
        )
        st[0].transform_pos_and_adp(sup.transform)
        st.setup_entities()
        st.write_pdb(str(out_path))

        aligned.append(out_path)
        report.append(dict(structure=path.name, reference=ref_path.name,
                           rmsd=float(sup.rmsd), n_aligned=int(sup.count),
                           chain=_ligand_chain(path) or "", is_reference=False))

    return aligned, report


# --------------------------------------------------------------------------- #
# Grid derivation and validation
# --------------------------------------------------------------------------- #

def derive_grid(
    structures: Sequence[str | Path],
    comp_ids: Optional[Dict[str, str]] = None,
    include_cofactors: bool = False,
    padding: float = 6.0,
    min_size: float = 18.0,
    max_size: float = 30.0,
    ligand_extent: Optional[float] = None,
) -> GridSpec:
    """
    Consensus grid centre and box size from a CO-FRAMED ensemble.

    Each structure contributes ONE site centre (the ligand copy nearest the reference
    site, so multimers do not skew the mean). The box is sized to cover the largest
    ligand that will be docked plus twice the residual centroid spread plus padding -
    not a fixed default, which is over-large for a compact pocket and admits spurious
    surface poses.

    Pass `ligand_extent` as the longest principal extent (A) across the screening
    library; without it the co-crystal ligand extents are used, which underestimates
    the box when screening compounds are larger.

    This function does NOT check that the ensemble is co-framed - call validate_grid()
    or superpose_ensemble() first.
    """
    structures = [Path(s) for s in structures]
    comp_ids = comp_ids or {}

    first = site_center(structures[0], comp_ids.get(structures[0].stem), include_cofactors)
    if first is None:
        raise ValueError(
            f"{structures[0].name}: no binding-site ligand found. "
            "Pass an explicit grid centre for apo structures."
        )
    ref_point = first.center

    sites: List[SiteCenter] = []
    for path in structures:
        sc = site_center(path, comp_ids.get(path.stem), include_cofactors, near=ref_point)
        if sc is None:
            raise ValueError(
                f"{path.name}: no binding-site ligand found. "
                "Pass an explicit grid centre, or exclude this structure."
            )
        sites.append(sc)

    centers = np.vstack([s.center for s in sites])
    consensus = centers.mean(axis=0)
    max_dev = float(max(np.linalg.norm(c - consensus) for c in centers))

    if ligand_extent is None:
        ligand_extent = float(np.vstack([s.extent for s in sites]).max())

    edge = float(np.ceil(ligand_extent + 2.0 * max_dev + padding))
    edge = float(np.clip(edge, min_size, max_size))

    return GridSpec(
        center=[round(float(v), 2) for v in consensus],
        size=[edge, edge, edge],
        max_deviation=round(max_dev, 2),
        sites=sites,
        spread_std=[round(float(v), 2) for v in centers.std(axis=0)],
    )


def validate_grid(
    structures: Sequence[str | Path],
    center: Sequence[float],
    size: Sequence[float],
    comp_ids: Optional[Dict[str, str]] = None,
    include_cofactors: bool = False,
    margin: float = 0.25,
) -> Tuple[bool, List[Dict]]:
    """
    Check that every member's binding site falls inside the grid box.

    This is the gate that catches an unaligned ensemble. A structure whose site centre
    lies outside the box - or within `margin` of its face - will be docked into solvent
    or against a truncated pocket, and its scores are meaningless.

    Returns (ok, rows). `ok` is False if any structure fails. Structures with no ligand
    are reported with status "no_ligand" and do not fail the check, since an explicit
    centre may have been supplied deliberately.
    """
    structures = [Path(s) for s in structures]
    comp_ids = comp_ids or {}
    center = np.asarray(center, dtype=float)
    half = np.asarray(size, dtype=float) / 2.0

    rows: List[Dict] = []
    ok = True

    for path in structures:
        sc = site_center(path, comp_ids.get(path.stem), include_cofactors, near=center)
        if sc is None:
            rows.append(dict(structure=path.name, comp_id=None, distance=None,
                             status="no_ligand"))
            continue

        delta = np.abs(sc.center - center)
        dist = float(np.linalg.norm(sc.center - center))
        if np.any(delta > half):
            status = "OUTSIDE"
            ok = False
        elif np.any(delta > half * (1.0 - margin)):
            status = "near_edge"
        else:
            status = "inside"

        rows.append(dict(structure=path.name, comp_id=sc.comp_id,
                         chain=sc.chain, n_copies=sc.n_copies,
                         distance=round(dist, 2), status=status))

    return ok, rows


def contact_residues(
    path: str | Path,
    comp_id: Optional[str] = None,
    cutoff: float = 4.5,
    include_cofactors: bool = False,
) -> List[Tuple[int, str]]:
    """
    Polymer residues within `cutoff` A of the site ligand.

    Use this to confirm a co-crystal ligand marks the intended pocket rather than an
    allosteric or surface site. A structure whose ligand contacts none of the target's
    known site residues should not be in the ensemble: it pulls the consensus centre
    away from the site being screened.
    """
    path = Path(path)
    sc = site_center(path, comp_id, include_cofactors)
    if sc is None:
        return []

    st = _read(path)
    model = st[0]

    lig_pos = []
    for ch in model:
        if ch.name != sc.chain:
            continue
        for res in ch:
            if _comp_matches(res.name, sc.comp_id):
                lig_pos = [a.pos for a in res]
                break
        if lig_pos:
            break
    if not lig_pos:
        return []

    near = set()
    for ch in model:
        for res in ch:
            if res.het_flag != "A":
                continue
            for atom in res:
                if any(atom.pos.dist(lp) <= cutoff for lp in lig_pos):
                    near.add((res.seqid.num, res.name))
                    break
    return sorted(near)
