"""Export docked coordinates with verified prepared-ligand chemistry.

PDBQT atom serials link the docking input to its output. The input's original
coordinates link those serials to the prepared SDF, without guessing bonds from
the docked geometry. No docking, minimisation, or force-field assignment occurs.
"""

import csv
import json
import re
from pathlib import Path

import numpy as np
from rdkit import Chem


def _results_file(run):
    root = Path(run).expanduser().resolve()
    if root.is_file():
        return root
    for name in ("redock_results.json", "redock_results.csv",
                 "protocol_development_results.csv",
                 "protocol_development/protocol_development_results.csv"):
        candidate = root / name
        if candidate.is_file():
            return candidate
    raise ValueError(f"No supported results file found in {root}")


def _local_output(raw, root):
    """Prefer this copy of a run, even if its original absolute paths exist."""
    path = Path(raw)
    for index in range(len(path.parts)):
        tail = path.parts[index:]
        if Path(*tail).is_absolute():
            candidate = Path(*tail)
        else:
            candidate = root.joinpath(*tail)
        candidate = candidate.resolve()
        if root in candidate.parents and candidate.is_file():
            return candidate
    raise ValueError(f"Docked output not found within this run: {raw}")


def _case_directory(output, root):
    for parent in output.parents:
        if parent == root:
            break
        if parent.name == "variants":
            return parent.parent
    return output.parent.parent if output.parent.name.startswith("Protocol") else output.parent


def export_cases(run):
    """Return selectable case IDs and output paths without changing the run."""
    result_file = _results_file(run)
    if result_file.suffix.lower() == ".json":
        rows = json.loads(result_file.read_text()).get("results", [])
    else:
        with result_file.open(newline="") as stream:
            rows = list(csv.DictReader(stream))
    cases = []
    seen = set()
    for row in rows:
        raw = row.get("output_file")
        if not raw:
            continue
        error = None
        try:
            output = _local_output(raw, result_file.parent)
        except ValueError as exc:
            output, error = Path(raw), str(exc)
        case_dir = _case_directory(output, result_file.parent)
        case_id = str(row.get("case_id") or case_dir.name)
        key = (case_id, output)
        if key not in seen:
            seen.add(key)
            cases.append({"case_id": case_id, "output": output, "case_dir": case_dir,
                          "row": row, "results_file": result_file, "error": error})
    return cases


def _sdf(path, index=1):
    supplier = Chem.SDMolSupplier(str(path), removeHs=False)
    if index < 1 or index > len(supplier):
        raise ValueError(f"Pose {index} is outside 1..{len(supplier)} in {path}")
    mol = supplier[index - 1]
    if mol is None or not mol.GetNumConformers():
        raise ValueError(f"Invalid molecule/coordinates at record {index} in {path}")
    if not np.isfinite(mol.GetConformer().GetPositions()).all():
        raise ValueError(f"Non-finite coordinates in {path}")
    return mol


def _pdbqt_atoms(path, pose=1):
    text = path.read_text()
    if any(line.startswith("MODEL") for line in text.splitlines()):
        models = []
        current = None
        for line in text.splitlines():
            if line.startswith("MODEL"):
                current = []
                models.append(current)
            elif line.startswith("ENDMDL"):
                current = None
            elif current is not None:
                current.append(line)
    else:
        models = [text.splitlines()]
    if pose < 1 or pose > len(models):
        raise ValueError(f"Pose {pose} is outside 1..{len(models)} in {path}")
    atoms = {}
    type_elements = {"A": "C", "NA": "N", "OA": "O", "SA": "S", "HD": "H", "HS": "H"}
    for line in models[pose - 1]:
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        serial = int(line[6:11])
        atom_type = line.split()[-1]
        element = type_elements.get(atom_type, atom_type.capitalize())
        xyz = np.array([float(line[start:start + 8]) for start in (30, 38, 46)])
        if serial in atoms or not np.isfinite(xyz).all():
            raise ValueError(f"Duplicate atom serial or invalid coordinates in {path}")
        atoms[serial] = (element, line[12:16], atom_type, xyz)
    if not atoms:
        raise ValueError(f"Pose {pose} contains no atoms in {path}")
    return atoms


def _prepared_files(case):
    output, case_dir = case["output"], case["case_dir"]
    relative = output.relative_to(case_dir)
    variant = relative.parts[1] if relative.parts[0] == "variants" else None
    candidates = []
    # rDock historically prepares its actual input inside the docking directory.
    if output.suffix.lower() in {".sd", ".sdf"} and variant:
        candidates.append(output.parent / f"{variant}.sdf")
    if variant:
        candidates.append(case_dir / "ligand_variants" / f"{variant}.sdf")
    candidates = [p for p in dict.fromkeys(candidates)
                  if p != output and p.is_file() and "reference" not in p.stem]
    if not candidates:
        candidates = [p for p in sorted(output.parent.glob("*.sdf"))
                      if p != output and "reference" not in p.stem]
        if len(candidates) > 1:
            raise ValueError("Multiple possible prepared SDFs found; restore the named variant preparation file")
    if not candidates:
        raise ValueError("Prepared ligand SDF is missing; restore the matching preparation files before export")
    topology = candidates[0]
    pdbqt = topology.with_suffix(".pdbqt")
    if output.suffix.lower() == ".pdbqt" and not pdbqt.is_file():
        raise ValueError(f"Matching prepared PDBQT is required alongside {topology}")
    return topology, pdbqt


def _pdbqt_coordinates(prepared, input_file, output_file, pose):
    original = _pdbqt_atoms(input_file)
    docked = _pdbqt_atoms(output_file, pose)
    if original.keys() != docked.keys() or any(original[k][:3] != docked[k][:3] for k in original):
        raise ValueError("Docked PDBQT atom identities differ from the prepared input")
    coordinates = prepared.GetConformer().GetPositions()
    symbols = np.array([atom.GetSymbol() for atom in prepared.GetAtoms()])
    mapping = {}
    for serial, (element, _name, _type, xyz) in original.items():
        matches = np.flatnonzero((symbols == element) & (np.linalg.norm(coordinates - xyz, axis=1) < 0.003))
        if len(matches) != 1 or int(matches[0]) in mapping:
            raise ValueError("Prepared SDF and PDBQT coordinates do not match uniquely; restore the original matching preparation pair")
        mapping[int(matches[0])] = docked[serial][3]
    if any(atom.GetIdx() not in mapping for atom in prepared.GetAtoms() if atom.GetAtomicNum() != 1):
        raise ValueError("Prepared/docked heavy-atom counts differ")
    return mapping


def _sdf_coordinates(prepared, docked):
    reference = Chem.RemoveHs(prepared)
    target = Chem.RemoveHs(docked)
    if Chem.MolToSmiles(reference) != Chem.MolToSmiles(target):
        raise ValueError("Docked SDF chemistry differs from the prepared ligand")
    match = docked.GetSubstructMatch(prepared, useChirality=True)
    if len(match) == prepared.GetNumAtoms():
        positions = docked.GetConformer().GetPositions()
        return {index: positions[other] for index, other in enumerate(match)}
    match = target.GetSubstructMatch(reference, useChirality=True)
    if len(match) != reference.GetNumAtoms():
        raise ValueError("Could not map the prepared ligand to the docked SDF")
    source_indices = [a.GetIdx() for a in prepared.GetAtoms() if a.GetAtomicNum() != 1]
    if len(source_indices) != len(match):
        raise ValueError("Nonstandard hydrogen topology requires a complete atom mapping")
    positions = target.GetConformer().GetPositions()
    return {source_indices[i]: positions[other] for i, other in enumerate(match)}


def _with_coordinates(prepared, mapping):
    result = Chem.Mol(prepared)
    conf = Chem.Conformer(result.GetNumAtoms())
    for atom in result.GetAtoms():
        atom.SetIntProp("_export_index", atom.GetIdx())
    for index, xyz in mapping.items():
        conf.SetAtomPosition(index, tuple(float(v) for v in xyz))
    result.RemoveAllConformers()
    result.AddConformer(conf)
    missing = set(range(result.GetNumAtoms())) - mapping.keys()
    if missing:
        if any(result.GetAtomWithIdx(i).GetAtomicNum() != 1 or result.GetAtomWithIdx(i).GetIsotope() for i in missing):
            raise ValueError("Missing coordinates for nonstandard atoms; cannot reconstruct safely")
        hydrogenated = Chem.AddHs(Chem.RemoveHs(result), addCoords=True)
        available = {}
        for atom in hydrogenated.GetAtoms():
            if atom.GetAtomicNum() == 1 and not atom.HasProp("_export_index"):
                parent = atom.GetNeighbors()[0].GetIntProp("_export_index")
                available.setdefault(parent, []).append(atom.GetIdx())
        for index, xyz in mapping.items():
            atom = result.GetAtomWithIdx(index)
            if atom.GetAtomicNum() != 1 or atom.GetDegree() != 1:
                continue
            choices = available.get(atom.GetNeighbors()[0].GetIdx(), [])
            if choices:
                nearest = min(choices, key=lambda i: np.linalg.norm(
                    hydrogenated.GetConformer().GetPositions()[i] - xyz))
                choices.remove(nearest)
        for index in sorted(missing):
            atom = result.GetAtomWithIdx(index)
            if atom.GetDegree() != 1:
                raise ValueError("Cannot reconstruct a nonterminal hydrogen")
            parent = atom.GetNeighbors()[0].GetIdx()
            choices = available.get(parent, [])
            if not choices:
                raise ValueError("Prepared hydrogen count is inconsistent with its valence")
            conf.SetAtomPosition(index, hydrogenated.GetConformer().GetAtomPosition(choices.pop(0)))
        # AddConformer copies the supplied conformer.
        result.RemoveAllConformers()
        result.AddConformer(conf)
    for atom in result.GetAtoms():
        atom.ClearProp("_export_index")
    Chem.SanitizeMol(result)
    return result, len(missing)


def _complex_pdb(receptor, ligand, residue_name):
    import gemmi

    structure = gemmi.read_structure(str(receptor))
    if len(structure) != 1 or not any(len(chain) for chain in structure[0]):
        raise ValueError("Prepared receptor must contain exactly one nonempty model")
    if any(res.name == residue_name for chain in structure[0] for res in chain):
        raise ValueError(f"Receptor already contains {residue_name}; choose a different --ligand-resname")
    used = {chain.name for chain in structure[0]}
    chain_name = next((c for c in "LABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789" if c not in used), None)
    if chain_name is None:
        raise ValueError("No unused one-character PDB chain ID is available")
    ligand_structure = gemmi.read_pdb_string(Chem.MolToPDBBlock(ligand))
    residue = ligand_structure[0][0][0].clone()
    residue.name = residue_name
    residue.seqid = gemmi.SeqId(1, " ")
    residue.het_flag = "H"
    chain = gemmi.Chain(chain_name)
    chain.add_residue(residue)
    structure[0].add_chain(chain)
    structure.assign_serial_numbers()
    return structure.make_pdb_string(), chain_name


def export_pose(run, case_id, pose=1, format="sdf", output=None, ligand_resname="LIG"):
    """Export a one-based file-order pose; refuse ambiguous cases and overwrites."""
    if pose < 1:
        raise ValueError("--pose is one-based and must be positive")
    if format not in {"sdf", "complex-pdb"}:
        raise ValueError("Export format must be sdf or complex-pdb")
    if not re.fullmatch(r"[A-Za-z0-9]{1,3}", ligand_resname):
        raise ValueError("--ligand-resname must contain 1-3 letters/digits")
    cases = [case for case in export_cases(run) if case["case_id"] == case_id]
    if len(cases) != 1:
        raise ValueError(f"Expected one output for case {case_id!r}, found {len(cases)}; use --list-cases")
    case = cases[0]
    if case.get("error"):
        raise ValueError(case["error"])
    topology, prepared_pdbqt = _prepared_files(case)
    prepared = _sdf(topology)
    source = case["output"]
    if source.suffix.lower() == ".pdbqt":
        mapping = _pdbqt_coordinates(prepared, prepared_pdbqt, source, pose)
    elif source.suffix.lower() in {".sd", ".sdf"}:
        mapping = _sdf_coordinates(prepared, _sdf(source, pose))
    else:
        raise ValueError(f"Unsupported docked file: {source}")
    molecule, rebuilt_hydrogens = _with_coordinates(prepared, mapping)
    molecule.SetProp("_Name", f"{case_id} pose {pose}")
    molecule.SetProp("DockMate_source", str(source))
    molecule.SetIntProp("DockMate_pose", pose)
    molecule.SetIntProp("DockMate_rebuilt_hydrogens", rebuilt_hydrogens)
    roundtrip = Chem.MolFromMolBlock(Chem.MolToMolBlock(molecule), removeHs=False)
    if (roundtrip is None or roundtrip.GetNumAtoms() != prepared.GetNumAtoms()
            or Chem.MolToSmiles(Chem.RemoveHs(roundtrip)) != Chem.MolToSmiles(Chem.RemoveHs(prepared))):
        raise ValueError("Pose geometry/topology does not preserve ligand chemistry on SDF round-trip; inspect the pose before export")
    suffix = ".sdf" if format == "sdf" else ".pdb"
    safe_id = re.sub(r"[^A-Za-z0-9_.-]+", "_", case_id)
    destination = Path(output).expanduser().resolve() if output else Path.cwd() / f"{safe_id}_pose{pose}{suffix}"
    sdf_path = destination if format == "sdf" else destination.with_suffix(".ligand.sdf")
    metadata_path = destination.with_suffix(destination.suffix + ".json")
    metadata = {"case_id": case_id, "pose": pose, "pose_order": "docking output file order (not rescored rank)",
                "source": str(source), "prepared_topology": str(topology),
                "canonical_smiles": Chem.MolToSmiles(Chem.RemoveHs(molecule)),
                "rebuilt_hydrogens": rebuilt_hydrogens,
                "hydrogen_policy": "Retain mapped coordinates; rebuild missing ordinary H with RDKit, without minimisation",
                "md_ready": False, "ligand_sdf": str(sdf_path)}
    if source.suffix.lower() == ".pdbqt":
        metadata["prepared_pdbqt"] = str(prepared_pdbqt)
    pdb_text = None
    if format == "complex-pdb":
        receptor = case["case_dir"] / "receptor_prepared.pdb"
        if not receptor.is_file():
            raise ValueError(f"Prepared receptor PDB is required: {receptor}")
        pdb_text, chain = _complex_pdb(receptor, molecule, ligand_resname.upper())
        metadata.update(receptor=str(receptor), ligand_chain=chain, ligand_resname=ligand_resname.upper())
    for path in {destination, sdf_path, metadata_path}:
        if path.exists():
            raise ValueError(f"Refusing to overwrite {path}; choose another --output")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with Chem.SDWriter(str(sdf_path)) as writer:
        writer.write(molecule)
    if pdb_text is not None:
        destination.write_text(pdb_text)
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    return destination
