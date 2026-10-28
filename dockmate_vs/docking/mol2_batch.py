"""Lossless MOL2 reading and conservative Smina batch score association.

This module does no chemical preparation or format conversion. Smina's unnamed
stdout tables may only be associated with a COMPLETE, ordered batch. An explicit
ID/mode score table can instead describe partial archived results.
"""
from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
import math
import re
from typing import Optional

import numpy as np

from .base import DockingError, DockingPose, DockingResult


@dataclass
class Mol2Record:
    name: str
    text: str
    coordinates: np.ndarray
    record_number: int


def read_mol2(path: Path) -> list[Mol2Record]:
    """Read all records, validating counts/connectivity without changing chemistry."""
    # Preserve CRLF and exact record content; bytes are never rewritten by this reader.
    text = Path(path).read_bytes().decode('utf-8')
    starts = list(re.finditer(r'^@<TRIPOS>MOLECULE\s*$', text, re.MULTILINE))
    if not starts or any(line.strip() and not line.lstrip().startswith('#')
                         for line in text[:starts[0].start()].splitlines()):
        raise DockingError(f'Invalid or empty MOL2 file: {path}')
    records = []
    for index, start in enumerate(starts):
        end = starts[index + 1].start() if index + 1 < len(starts) else len(text)
        block = text[start.start():end]
        sections = {}
        for section in block.split('@<TRIPOS>')[1:]:
            key, _, value = section.partition('\n')
            if key.strip() in sections:
                raise DockingError(f'Duplicate MOL2 section in record {index + 1}: {path}')
            sections[key.strip()] = value
        try:
            header = sections['MOLECULE'].splitlines()
            name = header[0].strip()
            n_atoms, n_bonds = map(int, header[1].split()[:2])
            atoms = [x.split() for x in sections['ATOM'].splitlines()
                     if x.strip() and not x.lstrip().startswith('#')]
            bonds = [x.split() for x in sections.get('BOND', '').splitlines()
                     if x.strip() and not x.lstrip().startswith('#')]
            if not name or n_atoms < 1 or n_bonds < 0 or (len(atoms), len(bonds)) != (n_atoms, n_bonds):
                raise ValueError('missing title or inconsistent atom/bond counts')
            atom_ids = [int(a[0]) for a in atoms]
            atom_id_set = set(atom_ids)
            bond_ids = [int(b[0]) for b in bonds]
            if len(atom_id_set) != len(atom_ids) or len(set(bond_ids)) != len(bond_ids):
                raise ValueError('duplicate atom/bond identifiers')
            if any(int(b[1]) not in atom_id_set or int(b[2]) not in atom_id_set for b in bonds):
                raise ValueError('bond refers to absent atom')
            coords = np.array([[float(v) for v in a[2:5]] for a in atoms])
            if coords.shape != (n_atoms, 3) or not np.isfinite(coords).all():
                raise ValueError('invalid coordinates')
            if any(len(a) < 6 for a in atoms) or any(len(b) < 4 for b in bonds):
                raise ValueError('missing atom or bond type')
        except (KeyError, ValueError, IndexError) as exc:
            raise DockingError(f'Invalid MOL2 record {index + 1} in {path}: {exc}') from exc
        records.append(Mol2Record(name, block, coords, index + 1))
    return records


def input_ids(path: Path) -> list[str]:
    ids = [r.name for r in read_mol2(path)]
    if len(set(ids)) != len(ids):
        raise DockingError(f'Prepared library contains duplicate compound IDs: {path}')
    return ids


def read_score_rows(path: Path, id_column='compound_id', score_column='score',
                    mode_column: Optional[str] = None, higher_is_better=False) -> dict:
    """Explicitly normalise an ID-keyed CSV to lower-is-better scores; never infer sign."""
    import csv
    with Path(path).open(newline='', encoding='utf-8-sig') as handle:
        reader = csv.DictReader(handle)
        required = {id_column, score_column} | ({mode_column} if mode_column else set())
        if not required.issubset(reader.fieldnames or []):
            raise DockingError(f'Score CSV requires columns {sorted(required)}: {path}')
        if len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise DockingError(f'Duplicate score CSV column names: {path}')
        rows = {}
        for line, row in enumerate(reader, 2):
            try:
                name = row[id_column].strip()
                mode = int(row[mode_column]) if mode_column else 1
                score = float(row[score_column])
                if not name or mode < 1 or not math.isfinite(score):
                    raise ValueError('empty ID, invalid mode or nonfinite score')
                key = (name, mode)
                if key in rows:
                    raise ValueError(f'duplicate compound/mode {key}')
                rows[key] = -score if higher_is_better else score
            except (TypeError, ValueError, AttributeError) as exc:
                raise DockingError(f'Invalid score CSV row {line}: {exc}') from exc
    return rows


def _log_tables(text: str) -> list[list[tuple]]:
    tables = []
    active = False
    for line in text.splitlines():
        if re.match(r'^\s*mode\s*\|\s*affinity\s*\|', line):
            tables.append([])
            active = True
            continue
        if not active:
            continue
        if re.match(r'^\s*\d+\s+', line):
            parts = line.split()
            try:
                if len(parts) != 4:
                    raise ValueError('expected mode, affinity and two RMSD values')
                mode = int(parts[0])
                values = tuple(float(x) for x in parts[1:])
                if mode != len(tables[-1]) + 1 or not all(math.isfinite(x) for x in values):
                    raise ValueError('invalid mode order or nonfinite score')
                if min(values[1:]) < 0:
                    raise ValueError('negative inter-pose RMSD')
                tables[-1].append((mode, *values))
            except ValueError as exc:
                raise DockingError(f'Invalid Smina score table row: {line}') from exc
        elif tables[-1] and line.strip():
            active = False
    return tables


def parse_mol2_batch(output_file: Path, expected_ids: list[str], *,
                     log_text: Optional[str] = None, score_rows: Optional[dict] = None) -> list[DockingResult]:
    """Return one result for every expected ID, refusing ambiguous associations.

    Unnamed log mapping assumes Smina's serial molecule/output loop. Callers must
    supply the complete, successful batch's log (never a timed-out process).
    CSV scores use explicit (compound_id, one-based per-compound pose mode) keys.
    """
    if not expected_ids or len(set(expected_ids)) != len(expected_ids):
        raise DockingError('Expected compound IDs must be nonempty and unique')
    if (log_text is None) == (score_rows is None):
        raise DockingError('Supply exactly one of a complete Smina log or an ID-keyed score table')
    records = read_mol2(output_file) if Path(output_file).stat().st_size else []
    groups = OrderedDict()
    previous = None
    for record in records:
        if record.name not in expected_ids:
            raise DockingError(f'Unexpected output compound ID: {record.name}')
        if record.name in groups and record.name != previous:
            raise DockingError(f'Non-contiguous/duplicated output compound: {record.name}')
        groups.setdefault(record.name, []).append(record)
        previous = record.name
    mapping = {}
    source = 'id_keyed_csv'
    if log_text is not None:
        source = 'smina_stdout_affinity'
        tables = _log_tables(log_text)
        if list(groups) != expected_ids or len(tables) != len(expected_ids):
            raise DockingError('Ambiguous unnamed Smina log: input/output IDs, order and table counts must all match; use an ID-keyed score CSV for partial archives')
        for name, table in zip(expected_ids, tables):
            if not table or len(table) != len(groups[name]):
                raise DockingError(f'Ambiguous unnamed Smina log: pose/score counts differ for {name}')
            for mode, score, lb, ub in table:
                mapping[(name, mode)] = (score, lb, ub)
    else:
        for (name, mode), score in score_rows.items():
            if name not in expected_ids or not isinstance(mode, int) or mode < 1 or not math.isfinite(score):
                raise DockingError(f'Invalid or unexpected score key: {(name, mode)}')
            if name in groups and mode > len(groups[name]):
                raise DockingError(f'Score has no matching saved pose: {(name, mode)}')
            mapping[(name, mode)] = (score, None, None)
    results = []
    for name in expected_ids:
        saved = groups.get(name, [])
        metadata = dict(output_format='mol2', output_file=str(Path(output_file).resolve()),
                        output_record_numbers=[r.record_number for r in saved], score_source=source,
                        pose_rmsd_reference='best docked mode' if log_text is not None else 'not supplied')
        missing = not saved or any((name, mode) not in mapping for mode in range(1, len(saved) + 1))
        if missing:
            reason = 'No saved pose' if not saved else 'Missing ID-keyed score for a saved pose'
            metadata['failure_category'] = 'missing_pose' if not saved else 'missing_score'
            results.append(DockingResult(name, success=False, error_message=reason, metadata=metadata))
            continue
        poses = []
        for mode, record in enumerate(saved, 1):
            score, lb, ub = mapping[(name, mode)]
            poses.append(DockingPose(score=score, rmsd_lb=lb, rmsd_ub=ub,
                                    coordinates=record.coordinates, mode_number=mode, mol2=record.text))
        results.append(DockingResult(name, poses=poses, metadata=metadata))
    return results
