"""Consistent selection of a single native ligand from PDB atom records."""

from Bio.PDB.Polypeptide import is_aa

NATIVE_PREPARATION_VERSION = "hetero_identity_v1"


def biopython_ligand(structure, resname, chain=None):
    """Return the same first-model residue used for reference/site selection."""
    model = next(structure.get_models())
    residues = list(model.get_residues())
    records = [("HETATM" if residue.id[0].startswith("H_") else "ATOM",
                residue.resname, residue.get_parent().id, residue.id)
               for residue in residues]
    indices = ligand_atom_indices(records, resname, chain)
    return residues[indices[0]]


def ligand_atom_indices(records, resname, chain=None):
    """Select the first matching non-polymer residue, respecting an explicit chain.

    Records are (record_type, residue_name, chain_id, residue_id) tuples in file
    order. HETATM takes precedence over ATOM; amino-acid ATOM records cannot
    serve as a free ligand. Non-amino-acid ATOM ligands remain supported.
    Residue IDs must distinguish insertion codes and separate residue copies.
    """
    name = str(resname).strip().upper()
    chain = None if chain is None else str(chain).strip()
    candidates = [i for i, (_, residue, chain_id, _) in enumerate(records)
                  if residue.strip().upper() == name
                  and (chain is None or chain_id.strip() == chain)]
    hetero = [i for i in candidates if records[i][0].strip() == "HETATM"]
    if hetero:
        candidates = hetero
    elif is_aa(name, standard=False):
        candidates = []
    if not candidates:
        available = sorted({c.strip() or "<blank>" for kind, r, c, _ in records
                            if r.strip().upper() == name and kind.strip() == "HETATM"})
        present = sorted({f"{r} {rid}" for kind, r, c, rid in records
                          if kind.strip() == "HETATM" and (chain is None or c.strip() == chain)})
        raise ValueError(
            f"No ligand found: resname={resname}, chain={chain}. "
            f"Available chains for {resname}: {', '.join(available) or 'none'}. "
            f"HETATM residues in requested chain: {', '.join(present) or 'none'}. "
            "The workbook Chain column selects the ligand copy; protein amino-acid "
            "ATOM records are not free ligands."
        )
    first = records[candidates[0]]
    return [i for i in candidates if records[i][0] == first[0]
            and records[i][2:] == first[2:]]
