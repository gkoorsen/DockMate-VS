"""Headless screening/import of intact prepared MOL2 libraries.

This intentionally bypasses workbook parsing, ligand/receptor preparation,
variant selection and rescoring. All labels come from an explicit CSV.
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import shutil

import numpy as np
import yaml

from .docking.base import DockingError, DockingResult
from .docking.mol2_batch import input_ids, read_mol2, read_score_rows
from .docking.prepared_execution import checksum
from .docking.smina import SminaDockingEngine


def _write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def _load_labels(path, expected):
    labels = {}
    with path.open(newline='', encoding='utf-8-sig') as handle:
        reader = csv.DictReader(handle)
        if not {'compound_id', 'label'}.issubset(reader.fieldnames or []):
            raise ValueError('Labels CSV requires compound_id,label columns')
        if len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError('Duplicate labels CSV column names')
        for row in reader:
            name = (row.get('compound_id') or '').strip()
            label = (row.get('label') or '').strip().lower()
            if not name or name in labels or label not in {'active', 'decoy', 'inactive', 'unknown'}:
                raise ValueError(f'Invalid or duplicate label row: {row}')
            labels[name] = label
    if set(labels) != set(expected):
        raise ValueError('Label IDs must exactly match all prepared library IDs (including failures/unknowns)')
    return labels


def load_prepared_config(config_path, execute=True):
    path = Path(config_path).resolve()
    cfg = yaml.safe_load(path.read_text())
    allowed = {'receptor', 'reference_ligand', 'center', 'size', 'labels', 'output_dir', 'batches',
               'smina_binary', 'scoring', 'exhaustiveness', 'num_modes', 'energy_range', 'cpu',
               'seed', 'timeout_sec', 'score_columns'}
    if not isinstance(cfg, dict) or set(cfg) - allowed:
        raise ValueError('Invalid prepared-screening config or unsupported keys')
    if not isinstance(cfg.get('scoring', 'vina'), str) or not cfg.get('scoring', 'vina').strip():
        raise ValueError('scoring must be a nonempty Smina scoring-function name')
    for key in ('receptor', 'labels', 'output_dir', 'batches'):
        if key not in cfg:
            raise ValueError(f'Prepared-screening config requires {key}')
    def resolve(value):
        value = Path(value).expanduser()
        return (path.parent/value).resolve() if not value.is_absolute() else value.resolve()
    for key in ('receptor', 'labels', 'reference_ligand', 'output_dir'):
        if key in cfg:
            cfg[key] = resolve(cfg[key])
    if cfg['output_dir'].exists():
        raise ValueError('Prepared-screening output_dir must be new; existing results are never overwritten')
    if len(read_mol2(cfg['receptor'])) != 1:
        raise ValueError('Receptor must contain exactly one MOL2 record')
    if not isinstance(cfg['batches'], list) or not cfg['batches']:
        raise ValueError('At least one batch is required')
    expected = []
    for batch in cfg['batches']:
        allowed_batch = {'ligands'} if execute else {'ligands', 'poses', 'scores', 'log', 'log_complete'}
        if not isinstance(batch, dict) or set(batch) - allowed_batch or 'ligands' not in batch:
            raise ValueError('Invalid prepared batch keys')
        for key in ('ligands', 'poses', 'scores', 'log'):
            if key in batch:
                batch[key] = resolve(batch[key])
                if not batch[key].is_file():
                    raise ValueError(f'Batch file not found: {batch[key]}')
        batch['ids'] = input_ids(batch['ligands'])
        expected.extend(batch['ids'])
        if not execute:
            if 'poses' not in batch or ('scores' in batch) == ('log' in batch):
                raise ValueError('Import requires poses and exactly one of scores or log per batch')
            if 'log' in batch and batch.get('log_complete') is not True:
                raise ValueError('Import of unnamed Smina logs requires log_complete: true for a successful complete process')
    if len(set(expected)) != len(expected):
        raise ValueError('Compound IDs must be unique across all batches')
    cfg['label_map'] = _load_labels(cfg['labels'], expected)
    columns = cfg.get('score_columns', {})
    if not isinstance(columns, dict) or set(columns) - {'id_column', 'score_column', 'mode_column', 'higher_is_better'}:
        raise ValueError('Unsupported score_columns setting')
    if 'higher_is_better' in columns and not isinstance(columns['higher_is_better'], bool):
        raise ValueError('higher_is_better must be a YAML boolean')
    for key in ('id_column', 'score_column', 'mode_column'):
        if key in columns and not (key == 'mode_column' and columns[key] is None):
            if not isinstance(columns[key], str) or not columns[key].strip():
                raise ValueError(f'{key} must be a nonempty CSV column name')
    if execute:
        if 'score_columns' in cfg:
            raise ValueError('score_columns is only used for importing archived score CSVs')
        if ('center' in cfg) == ('reference_ligand' in cfg):
            raise ValueError('Choose exactly one of center or reference_ligand')
        if 'reference_ligand' in cfg:
            reference = read_mol2(cfg['reference_ligand'])
            if len(reference) != 1:
                raise ValueError('Reference ligand must contain exactly one MOL2 record')
            coords = reference[0].coordinates
            cfg['center'] = ((coords.min(0) + coords.max(0))/2).tolist()
        for key in ('center', 'size'):
            vec = np.asarray(cfg.get(key), float)
            if vec.shape != (3,) or not np.isfinite(vec).all() or (key == 'size' and np.any(vec <= 0)):
                raise ValueError(f'{key} requires three finite values (positive for size)')
        for key in ('exhaustiveness', 'num_modes', 'cpu'):
            value = cfg.get(key)
            if isinstance(value, bool) or not isinstance(value, int) or value < 1:
                raise ValueError(f'Explicit positive integer {key} is required')
        if isinstance(cfg.get('seed'), bool) or not isinstance(cfg.get('seed'), int):
            raise ValueError('An explicit integer seed is required (zero is supported)')
        for key in ('energy_range', 'timeout_sec'):
            value = cfg.get(key)
            if isinstance(value, bool) or not isinstance(value, (float, int)) or not math.isfinite(value) or value <= 0:
                raise ValueError(f'Explicit positive finite {key} is required; timeout_sec applies to each entire batch')
        binary = str(cfg.get('smina_binary', 'smina'))
        cfg['smina_binary'] = str(resolve(binary)) if '/' in binary or '\\' in binary else binary
    cfg['source_config'] = path
    return cfg


def _metrics(rows):
    """No imputation; complete-population and evaluable-only metrics are separate."""
    from sklearn.metrics import auc, average_precision_score, precision_recall_curve, roc_auc_score
    labelled = [r for r in rows if r['label'] != 'unknown']
    valid = [r for r in labelled if r['docking_completed']]
    positive = sum(r['label'] == 'active' for r in valid)
    numbers = dict(total=len(rows), labelled_total=len(labelled), unknowns=len(rows)-len(labelled),
                   labelled_actives=sum(r['label'] == 'active' for r in labelled),
                   labelled_negatives=sum(r['label'] != 'active' for r in labelled),
                   evaluable=len(valid), evaluable_actives=positive, evaluable_negatives=len(valid)-positive,
                   failed_total=sum(not r['docking_completed'] for r in rows),
                   failed_labelled=len(labelled)-len(valid))
    values = None
    if 0 < positive < len(valid):
        y = np.array([r['label'] == 'active' for r in valid])
        s = -np.array([r['best_score'] for r in valid])
        precision, recall, _ = precision_recall_curve(y, s)
        values = dict(roc_auc=float(roc_auc_score(y, s)),
                      pr_auc_trapezoidal=float(auc(recall, precision)),
                      average_precision=float(average_precision_score(y, s)))
    return dict(**numbers, all_labelled_metrics=values if len(valid) == len(labelled) else None,
                evaluable_only_metrics=values, failure_policy='no imputation; failed cases retained in table',
                molecule_ranking='minimum saved-pose score, lower is better; AUC groups equal-score ties')


def run_prepared(config_path, execute=True):
    cfg = load_prepared_config(config_path, execute)
    out = cfg['output_dir']
    out.mkdir(parents=True, exist_ok=False)
    inputs = out/'inputs'
    inputs.mkdir()
    snapshots = []
    def snapshot(source, name):
        destination = inputs/name
        shutil.copyfile(source, destination)
        source_hash = checksum(source)
        if checksum(destination) != source_hash:
            raise DockingError(f'Input changed during snapshot: {source}')
        snapshots.append(dict(source=str(source), copy=str(destination.relative_to(out)), sha256=source_hash))
        return destination
    receptor = snapshot(cfg['receptor'], 'receptor.mol2')
    snapshot(cfg['labels'], 'labels.csv')
    snapshot(cfg['source_config'], 'campaign.yml')
    reference = snapshot(cfg['reference_ligand'], 'reference.mol2') if 'reference_ligand' in cfg else None
    libraries = [snapshot(batch['ligands'], f'library_{i:03d}.mol2')
                 for i, batch in enumerate(cfg['batches'], 1)]
    for batch, library in zip(cfg['batches'], libraries):
        if input_ids(library) != batch['ids']:
            raise DockingError('Prepared library identities changed after preflight')
    if _load_labels(inputs/'labels.csv', [name for b in cfg['batches'] for name in b['ids']]) != cfg['label_map']:
        raise DockingError('Compound labels changed after preflight')
    manifest = dict(route='prepared_mol2_batch', operation='execute' if execute else 'import',
                    input_preparation='none', native_preparation_validated=False,
                    command='screen-prepared' if execute else 'import-prepared',
                    config={k: v for k, v in cfg.items() if k not in {'label_map', 'batches'}},
                    inputs=snapshots, batches=[])
    # Convert Paths only; keep scientific settings as numbers/booleans.
    manifest['config'] = {k: str(v) if isinstance(v, Path) else v for k, v in manifest['config'].items()}
    _write_json(out/'run_manifest.json', manifest)
    engine = None
    if execute:
        engine = SminaDockingEngine(smina_binary=cfg['smina_binary'],
            scoring_function=cfg.get('scoring', 'vina'), exhaustiveness=cfg['exhaustiveness'],
            num_modes=cfg['num_modes'], energy_range=cfg['energy_range'], cpu=cfg['cpu'],
            seed=cfg['seed'], timeout_sec=cfg['timeout_sec'])
        manifest['executable_version'] = engine.get_version()
    results = []
    for index, (batch, library) in enumerate(zip(cfg['batches'], libraries), 1):
        directory = out/f'batch_{index:03d}'
        if execute:
            execution = engine.dock_batch(receptor, library, cfg['center'], cfg['size'], directory,
                                          autobox_ligand=reference)
            current = execution.results
            info = dict(status=execution.status, command=execution.command,
                        runtime_sec=execution.runtime_sec, returncode=execution.returncode,
                        error_message=execution.error_message, metadata=execution.metadata)
        else:
            directory.mkdir()
            poses = directory/'poses.mol2'
            shutil.copyfile(batch['poses'], poses)
            artifacts = [dict(source=str(batch['poses']), copy=str(poses.relative_to(out)), sha256=checksum(poses))]
            kwargs = {}
            try:
                if 'scores' in batch:
                    score_file = directory/'scores.csv'
                    shutil.copyfile(batch['scores'], score_file)
                    artifacts.append(dict(source=str(batch['scores']), copy=str(score_file.relative_to(out)), sha256=checksum(score_file)))
                    all_scores = read_score_rows(score_file, **cfg.get('score_columns', {}))
                    unexpected = {name for name, mode in all_scores} - set(cfg['label_map'])
                    if unexpected:
                        raise DockingError(f'Score CSV contains unexpected compound IDs: {sorted(unexpected)}')
                    # A combined archive CSV may be shared by several libraries.
                    # Only IDs validated against the COMPLETE cohort may be selected.
                    kwargs['score_rows'] = {key: score for key, score in all_scores.items()
                                            if key[0] in batch['ids']}
                else:
                    log_file = directory/'stdout.log'
                    shutil.copyfile(batch['log'], log_file)
                    artifacts.append(dict(source=str(batch['log']), copy=str(log_file.relative_to(out)), sha256=checksum(log_file)))
                    kwargs['log_text'] = log_file.read_text()
                current = SminaDockingEngine.parse_batch_output(poses, batch['ids'], **kwargs)
                info = dict(status='imported', runtime_sec=None, artifacts=artifacts)
            except (DockingError, OSError, UnicodeError) as exc:
                current = [DockingResult(name, success=False, error_message=str(exc),
                           metadata={'failure_category': 'import_failed'}) for name in batch['ids']]
                info = dict(status='import_failed', error_message=str(exc), runtime_sec=None, artifacts=artifacts)
        for result in current:
            result.metadata['batch'] = index
            result.metadata['label'] = cfg['label_map'][result.ligand_id]
            result.metadata['scoring'] = cfg.get('scoring', 'vina')
        results.extend(current)
        manifest['batches'].append(info)
        _write_json(out/'run_manifest.json', manifest)
    rows = []
    for r in results:
        rows.append(dict(compound_id=r.ligand_id, label=cfg['label_map'][r.ligand_id],
                         engine='smina', scoring=cfg.get('scoring', 'vina'), best_score=r.best_score,
                         num_poses=len(r.poses), docking_completed=r.success,
                         runtime_sec=r.runtime, batch=r.metadata['batch'],
                         failure_category=r.metadata.get('failure_category', ''), error_message=r.error_message or '',
                         output_format='mol2', output_file=r.metadata.get('output_file', ''),
                         rmsd_status='not_applicable'))
    _write_json(out/'prepared_results.json', dict(mode='prepared_screening', results=[r.to_dict() for r in results]))
    with (out/'screening_results.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    summary = _metrics(rows)
    _write_json(out/'enrichment_summary.json', summary)
    manifest['summary'] = summary
    _write_json(out/'run_manifest.json', manifest)
    return out/'prepared_results.json'
