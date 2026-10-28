"""Prepared batch import/execution safeguards. Engine processes are mocked."""
import csv
import json
from pathlib import Path
import subprocess
from types import SimpleNamespace

import numpy as np
import pytest
import yaml

from dockmate_vs.cli import main
from dockmate_vs.docking.base import DockingError
from dockmate_vs.docking.mol2_batch import input_ids, parse_mol2_batch, read_mol2, read_score_rows
from dockmate_vs.docking.prepared_execution import checksum
from dockmate_vs.docking.smina import SminaDockingEngine
from dockmate_vs.prepared_screening import load_prepared_config, run_prepared


FIXTURE = Path(__file__).parent/'fixtures/prepared_mol2'


def synthetic_log(scores):
    """Table layout from Smina 2019 source; NOT an archived execution log."""
    return ''.join('Using random seed: 0\nmode |   affinity | dist from best mode\n'
                   '     | (kcal/mol) | rmsd l.b.| rmsd u.b.\n'
                   '-----+------------+----------+----------\n'
                   + ''.join(f'{i} {s} 0.000 0.000\n' for i, s in enumerate(batch, 1))
                   + 'Writing output ... done.\n' for batch in scores)


@pytest.fixture
def engine(monkeypatch):
    monkeypatch.setattr(SminaDockingEngine, '_verify_binary', lambda self: None)
    monkeypatch.setattr(SminaDockingEngine, 'get_version', lambda self: 'Smina test fixture')
    return SminaDockingEngine(smina_binary='smina-fixture', num_modes=1,
                             exhaustiveness=1, cpu=1, seed=0, timeout_sec=120)


def config(tmp_path, execute=False):
    cfg = dict(receptor=str(FIXTURE/'actives.mol2'), labels=str(FIXTURE/'labels.csv'),
               output_dir=str(tmp_path/'results'), batches=[
                   dict(ligands=str(FIXTURE/f'{kind}.mol2'),
                        poses=str(FIXTURE/f'{kind}-docked.mol2'), scores=str(FIXTURE/'scores.csv'))
                   for kind in ['actives', 'decoys']],
               score_columns=dict(id_column='ID', score_column='Score', higher_is_better=True))
    # A valid one-record MOL2 receptor suffices for mocked execution/import tests.
    receptor = tmp_path/'receptor.mol2'
    receptor.write_bytes(read_mol2(FIXTURE/'actives.mol2')[0].text.encode())
    cfg['receptor'] = str(receptor)
    if execute:
        cfg.pop('score_columns')
        for batch in cfg['batches']:
            batch.pop('poses'); batch.pop('scores')
        cfg.update(center=[2., 63., 67.], size=[30., 30., 30.], exhaustiveness=1,
                   num_modes=1, energy_range=3, cpu=1, seed=0, timeout_sec=120,
                   smina_binary='smina-fixture')
    path = tmp_path/'config.yml'
    path.write_text(yaml.safe_dump(cfg))
    return path, cfg


def test_author_records_keep_topology_coordinates_and_id():
    records = read_mol2(FIXTURE/'actives-docked.mol2')
    scores = read_score_rows(FIXTURE/'scores.csv', id_column='ID', score_column='Score', higher_is_better=True)
    ids = input_ids(FIXTURE/'actives.mol2')
    results = parse_mol2_batch(FIXTURE/'actives-docked.mol2', ids,
                              score_rows={k: v for k, v in scores.items() if k[0] in ids})
    assert [r.ligand_id for r in results] == ids
    for record, result in zip(records, results):
        assert result.success
        assert result.best_score == scores[(record.name, 1)]
        assert result.poses[0].mol2 == record.text
        assert '@<TRIPOS>BOND' in result.poses[0].mol2
        assert result.poses[0].pdbqt is None
        assert result.poses[0].rmsd_lb is None
        np.testing.assert_array_equal(result.poses[0].coordinates, record.coordinates)


def test_complete_log_maps_each_compound_and_multiple_modes(tmp_path):
    records = read_mol2(FIXTURE/'actives-docked.mol2')
    output = tmp_path/'multi.mol2'
    output.write_text(records[0].text*2 + records[1].text)
    results = parse_mol2_batch(output, [r.name for r in records],
                              log_text=synthetic_log([[-9., -8.], [-7.]]))
    assert [r.best_score for r in results] == [-9., -7.]
    assert [p.mode_number for p in results[0].poses] == [1, 2]
    assert results[1].metadata['output_record_numbers'] == [3]


@pytest.mark.parametrize('problem', ['missing_table', 'extra_table', 'wrong_order', 'pose_count', 'nonfinite', 'truncated'])
def test_ambiguous_log_rejected(problem):
    ids = input_ids(FIXTURE/'actives.mol2')
    log = synthetic_log([[-9.], [-8.]])
    if problem == 'missing_table': log = synthetic_log([[-9.]])
    if problem == 'extra_table': log = synthetic_log([[-9.], [-8.], [-7.]])
    if problem == 'wrong_order': ids.reverse()
    if problem == 'pose_count': log = synthetic_log([[-9., -8.], [-7.]])
    if problem == 'nonfinite': log = log.replace('-9.0', 'nan')
    if problem == 'truncated': log = log.replace('1 -9.0 0.000 0.000', '1 -9.0')
    with pytest.raises(DockingError):
        parse_mol2_batch(FIXTURE/'actives-docked.mol2', ids, log_text=log)


def test_partial_id_keyed_archive_retains_missing_compound(tmp_path):
    records = read_mol2(FIXTURE/'actives-docked.mol2')
    path = tmp_path/'partial.mol2'; path.write_text(records[1].text)
    result = parse_mol2_batch(path, [r.name for r in records],
                             score_rows={(records[0].name, 1): -9., (records[1].name, 1): -8.})
    assert result[0].success is False and result[0].best_score is None
    assert result[0].metadata['failure_category'] == 'missing_pose'
    assert result[1].best_score == -8.


def test_missing_id_keyed_score_is_not_zero_or_dropped():
    ids = input_ids(FIXTURE/'actives.mol2')
    results = parse_mol2_batch(FIXTURE/'actives-docked.mol2', ids, score_rows={(ids[1], 1): -8.})
    assert len(results) == 2 and not results[0].success
    assert results[0].best_score is None and results[1].best_score == -8.


@pytest.mark.parametrize('rows', ['A,1\nA,2\n', 'A,nan\n', ',3\n'])
def test_bad_score_csv_rejected(tmp_path, rows):
    path = tmp_path/'scores.csv'; path.write_text('compound_id,score\n'+rows)
    with pytest.raises(DockingError): read_score_rows(path)


def test_duplicate_score_column_rejected(tmp_path):
    path = tmp_path/'scores.csv'
    path.write_text('compound_id,score,score\nA,-8,-9\n')
    with pytest.raises(DockingError, match='Duplicate score CSV column'):
        read_score_rows(path)


def test_duplicate_input_ids_and_malformed_connectivity_rejected(tmp_path):
    record = read_mol2(FIXTURE/'actives.mol2')[0].text
    path = tmp_path/'input.mol2'; path.write_text(record*2)
    with pytest.raises(DockingError, match='duplicate'): input_ids(path)
    sections = record.split('@<TRIPOS>BOND')
    lines = sections[1].splitlines()
    for index, line in enumerate(lines):
        if line.strip():
            fields = line.split(); fields[1] = '999999'; lines[index] = ' '.join(fields); break
    path.write_text(sections[0]+'@<TRIPOS>BOND'+'\n'.join(lines))
    with pytest.raises(DockingError, match='absent atom'): read_mol2(path)


def test_batch_exec_uses_single_untouched_library(engine, monkeypatch, tmp_path):
    library = FIXTURE/'actives.mol2'
    before = checksum(library)
    calls = []
    def fake_run(command, stdout, stderr, timeout):
        calls.append(command)
        assert timeout == 120
        assert Path(command[command.index('--ligand')+1]).read_bytes() == library.read_bytes()
        Path(command[command.index('--out')+1]).write_bytes((FIXTURE/'actives-docked.mol2').read_bytes())
        stdout.write(synthetic_log([[-9.], [-8.]]))
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr('dockmate_vs.docking.prepared_execution.subprocess.run', fake_run)
    result = engine.dock_batch(library, library, [1, 2, 3], [30, 30, 30], tmp_path/'batch')
    assert result.status == 'completed' and len(result.results) == 2
    assert len(calls) == 1 and '--score_only' not in calls[0]
    assert calls[0][calls[0].index('--seed')+1] == '0'
    assert checksum(library) == before
    assert result.results[0].runtime is None  # Never invent per-compound runtime.
    assert (tmp_path/'batch/stdout.log').exists()
    with pytest.raises(FileExistsError):
        engine.dock_batch(library, library, [1, 2, 3], [30]*3, tmp_path/'batch')


@pytest.mark.parametrize('failure', ['timeout', 'process_failed', 'mapping_failed', 'launch_failed'])
def test_failed_batch_preserves_artifacts_and_unassigned_rows(engine, monkeypatch, tmp_path, failure):
    def fake_run(command, stdout, stderr, timeout):
        stdout.write(synthetic_log([[-9.]]))
        stderr.write('test failure diagnostic\n')
        Path(command[command.index('--out')+1]).write_text(read_mol2(FIXTURE/'actives-docked.mol2')[0].text)
        if failure == 'timeout': raise subprocess.TimeoutExpired(command, timeout)
        if failure == 'launch_failed': raise OSError('cannot start binary')
        return SimpleNamespace(returncode=2 if failure == 'process_failed' else 0)
    monkeypatch.setattr('dockmate_vs.docking.prepared_execution.subprocess.run', fake_run)
    result = engine.dock_batch(FIXTURE/'actives.mol2', FIXTURE/'actives.mol2', [1, 2, 3], [30]*3, tmp_path/'batch')
    assert result.status == failure and len(result.results) == 2
    assert all(not r.success and r.best_score is None for r in result.results)
    assert 'test failure diagnostic' in (tmp_path/'batch/stderr.log').read_text()
    assert (tmp_path/'batch/poses.mol2').stat().st_size > 0


def test_import_cli_shared_author_csv_without_subprocess(tmp_path, monkeypatch):
    path, cfg = config(tmp_path)
    def forbidden(*args, **kwargs): raise AssertionError('Import must not run any external tool')
    monkeypatch.setattr(subprocess, 'run', forbidden)
    assert main(['import-prepared', '--config', str(path)]) == 0
    out = Path(cfg['output_dir'])
    summary = json.loads((out/'enrichment_summary.json').read_text())
    assert summary['total'] == 4 and summary['failed_total'] == 0
    assert summary['all_labelled_metrics']['roc_auc'] is not None
    results = json.loads((out/'prepared_results.json').read_text())['results']
    assert all(r['best_score'] < 0 for r in results)
    assert (out/'inputs/library_001.mol2').read_bytes() == (FIXTURE/'actives.mol2').read_bytes()
    assert (out/'batch_002/poses.mol2').read_bytes() == (FIXTURE/'decoys-docked.mol2').read_bytes()
    assert main(['import-prepared', '--config', str(path)]) == 1  # No overwrite.


def test_import_partial_has_evaluable_only_metrics_and_nonzero_exit(tmp_path):
    path, cfg = config(tmp_path)
    scores = tmp_path/'scores.csv'
    original = (FIXTURE/'scores.csv').read_text().splitlines()
    scores.write_text('\n'.join(original[:-1])+'\n')
    for b in cfg['batches']: b['scores'] = str(scores)
    path.write_text(yaml.safe_dump(cfg))
    assert main(['import-prepared', '--config', str(path)]) == 1
    summary = json.loads((Path(cfg['output_dir'])/'enrichment_summary.json').read_text())
    assert summary['total'] == 4 and summary['failed_total'] == 1
    assert summary['all_labelled_metrics'] is None
    assert summary['evaluable_only_metrics'] is not None


def test_unknown_compounds_are_retained_but_excluded_from_enrichment(tmp_path):
    path, cfg = config(tmp_path)
    labels = tmp_path/'labels.csv'
    labels.write_text((FIXTURE/'labels.csv').read_text().replace('BDB9012,active', 'BDB9012,unknown'))
    cfg['labels'] = str(labels)
    path.write_text(yaml.safe_dump(cfg))
    result = run_prepared(path, execute=False)
    summary = json.loads((result.parent/'enrichment_summary.json').read_text())
    assert summary['total'] == 4 and summary['unknowns'] == 1 and summary['evaluable'] == 3
    rows = json.loads(result.read_text())['results']
    assert rows[0]['metadata']['label'] == 'unknown' and rows[0]['success']


def test_reference_box_uses_all_atom_extrema_and_preserves_zero_seed(tmp_path):
    path, cfg = config(tmp_path, execute=True)
    cfg.pop('center')
    cfg['reference_ligand'] = cfg['receptor']
    path.write_text(yaml.safe_dump(cfg))
    parsed = load_prepared_config(path)
    xyz = read_mol2(Path(cfg['reference_ligand']))[0].coordinates
    np.testing.assert_allclose(parsed['center'], (xyz.min(0) + xyz.max(0))/2)
    assert parsed['seed'] == 0
    cfg['timeout_sec'] = float('nan')
    path.write_text(yaml.safe_dump(cfg))
    with pytest.raises(ValueError, match='timeout_sec'):
        load_prepared_config(path)


@pytest.mark.parametrize('problem', ['duplicate_library', 'wrong_labels', 'unknown_key', 'unconfirmed_log', 'bad_score_column', 'null_scoring'])
def test_preflight_rejects_invalid_cohort_before_output_created(tmp_path, problem):
    path, cfg = config(tmp_path)
    if problem == 'duplicate_library': cfg['batches'][1]['ligands'] = cfg['batches'][0]['ligands']
    if problem == 'wrong_labels':
        labels = tmp_path/'labels.csv'; labels.write_text('compound_id,label\nnot-in-library,active\n')
        cfg['labels'] = str(labels)
    if problem == 'unknown_key': cfg['prepare_ligands'] = True
    if problem == 'bad_score_column': cfg['score_columns']['id_column'] = ['ID']
    if problem == 'null_scoring': cfg['scoring'] = None
    if problem == 'unconfirmed_log':
        cfg['batches'][0]['log'] = cfg['batches'][0].pop('scores')
    path.write_text(yaml.safe_dump(cfg))
    with pytest.raises(ValueError): load_prepared_config(path, execute=False)
    assert not Path(cfg['output_dir']).exists()


def test_execution_cli_retains_batching_and_failure_denominators(tmp_path, monkeypatch, engine):
    path, cfg = config(tmp_path, execute=True)
    calls = []
    def fake_run(command, stdout, stderr, timeout):
        calls.append(command)
        ligand = Path(command[command.index('--ligand')+1])
        output = Path(command[command.index('--out')+1])
        if '001' in ligand.name:
            output.write_bytes((FIXTURE/'actives-docked.mol2').read_bytes())
            stdout.write(synthetic_log([[-9.], [-8.]]))
            return SimpleNamespace(returncode=0)
        raise subprocess.TimeoutExpired(command, timeout)
    monkeypatch.setattr('dockmate_vs.docking.prepared_execution.subprocess.run', fake_run)
    assert main(['screen-prepared', '--config', str(path)]) == 1
    summary = json.loads((Path(cfg['output_dir'])/'enrichment_summary.json').read_text())
    assert len(calls) == 2  # One process per original library; never per molecule.
    assert summary['total'] == 4 and summary['failed_total'] == 2
    assert summary['all_labelled_metrics'] is None
    assert summary['evaluable_only_metrics'] is None  # Only actives survived.
