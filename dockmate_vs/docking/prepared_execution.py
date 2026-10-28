"""Execute an intact prepared MOL2 library once, retaining logs and failures."""
from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import time
from typing import Optional

import numpy as np

from .base import DockingError, DockingResult
from .mol2_batch import input_ids, parse_mol2_batch, read_mol2


def checksum(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


@dataclass
class BatchExecution:
    results: list[DockingResult]
    status: str
    runtime_sec: float
    command: list[str]
    returncode: Optional[int] = None
    error_message: Optional[str] = None
    metadata: dict = field(default_factory=dict)


def dock_prepared_batch(engine, receptor_path, ligand_path, center, size,
                        output_dir, autobox_ligand=None) -> BatchExecution:
    """No preparation/splitting/retries. Timeout applies to the ENTIRE batch.

    A failed or interrupted process invalidates anonymous score association for
    the batch. Raw partial output is retained; expected IDs are all reported with
    unknown per-compound attempt status, not silently discarded or given scores.
    """
    receptor = Path(receptor_path).resolve()
    ligand = Path(ligand_path).resolve()
    engine.validate_inputs(str(receptor), str(ligand))
    ids = input_ids(ligand)
    read_mol2(receptor)
    center, size = np.asarray(center, float), np.asarray(size, float)
    if center.shape != (3,) or size.shape != (3,) or not np.isfinite([center, size]).all() or np.any(size <= 0):
        raise DockingError('Prepared batch requires finite center and positive size vectors')
    if not np.isfinite(engine.timeout_sec) or engine.timeout_sec <= 0:
        raise DockingError('Prepared batch requires a positive finite batch timeout')
    inputs = [receptor, ligand]
    extra = {}
    if autobox_ligand is not None:
        reference = Path(autobox_ligand).resolve()
        read_mol2(reference)
        inputs.append(reference)
        extra['autobox_ligand'] = str(reference)
    before = {str(p): checksum(p) for p in inputs}
    out = Path(output_dir).resolve()
    out.mkdir(parents=True, exist_ok=False)
    poses = out/'poses.mol2'
    command = engine._build_command(str(receptor), str(ligand), str(poses), center, size, **extra)
    binary = shutil.which(engine.smina_binary)
    provenance = dict(input_sha256=before, input_preparation='none', split_library=False,
                      executable_path=binary, executable_sha256=checksum(binary) if binary else None,
                      timeout_scope='entire library batch', timeout_sec=engine.timeout_sec,
                      scoring=engine.scoring_function, num_modes=engine.num_modes,
                      center=center.tolist(), size=size.tolist())
    # Persist the intended invocation even if the process cannot start.
    (out/'command.json').write_text(json.dumps(dict(command=command, **provenance), indent=2)+'\n')
    start = time.monotonic()
    returncode = None
    error = None
    status = 'completed'
    with (out/'stdout.log').open('w') as stdout, (out/'stderr.log').open('w') as stderr:
        try:
            completed = subprocess.run(command, stdout=stdout, stderr=stderr, timeout=engine.timeout_sec)
            returncode = completed.returncode
            if returncode:
                status, error = 'process_failed', f'Smina exited with code {returncode}; anonymous batch scores were not assigned'
        except subprocess.TimeoutExpired:
            status, error = 'timeout', f'Smina batch exceeded {engine.timeout_sec:g} seconds; individual attempts cannot be determined'
        except OSError as exc:
            status, error = 'launch_failed', str(exc)
    runtime = time.monotonic() - start
    after = {str(p): checksum(p) if p.exists() else None for p in inputs}
    provenance['input_sha256_after'] = after
    if before != after:
        status, error = 'input_changed', 'Prepared input checksums changed during execution'
    results = []
    if status == 'completed':
        try:
            results = parse_mol2_batch(poses, ids, log_text=(out/'stdout.log').read_text())
            if any(len(r.poses) > engine.num_modes for r in results):
                raise DockingError('Output exceeds requested num_modes')
        except (DockingError, OSError, UnicodeError) as exc:
            status, error = 'mapping_failed', str(exc)
    if status != 'completed':
        results = [DockingResult(name, success=False, error_message=error, metadata={
            'failure_category': status, 'attempt_status': 'not_started' if status == 'launch_failed' else 'unknown',
            'output_format': 'mol2', 'output_file': str(poses), 'batch_runtime_sec': runtime,
        }) for name in ids]
    for result in results:
        result.metadata['batch_runtime_sec'] = runtime
    batch = BatchExecution(results, status, runtime, command, returncode, error, provenance)
    (out/'batch.json').write_text(json.dumps(dict(
        status=status, runtime_sec=runtime, command=command, returncode=returncode,
        error_message=error, metadata=provenance, results=[r.to_dict() for r in results]
    ), indent=2, allow_nan=False)+'\n')
    return batch
