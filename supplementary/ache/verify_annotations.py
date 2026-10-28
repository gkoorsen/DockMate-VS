"""Audit full compound populations and annotation comparisons from sealed MOL2 bytes.

Standard library only. No docking, preparation, fitting or atom remapping.
Run after verify_archive.py. Output is a JSON audit report on stdout.
"""
from collections import Counter
from pathlib import Path
import csv
import io
import json
import math
import xml.etree.ElementTree as ET
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parent
PREFIX = 'ACHE_matched_live/'
DIAGNOSTIC = 'diagnostics/final_pose_check_20261026T201436/'


def rows(data):
    return list(csv.DictReader(io.StringIO(data.decode('utf-8-sig'))))


def index(records, field='compound_id'):
    result = {r[field]: r for r in records}
    assert len(result) == len(records), 'Duplicate compound IDs'
    return result


def mol2(data):
    result = {}
    for body in data.decode().split('@<TRIPOS>MOLECULE')[1:]:
        header = body.strip().splitlines()
        name = header[0].strip()
        count = int(header[1].split()[0])
        atoms = [line.split() for line in body.split('@<TRIPOS>ATOM', 1)[1]
                 .split('@<TRIPOS>', 1)[0].splitlines() if line.strip()]
        assert name not in result and len(atoms) == count
        assert len({a[0] for a in atoms}) == count
        xyz = [tuple(map(float, a[2:5])) for a in atoms]
        charges = [float(a[8]) for a in atoms]
        assert all(math.isfinite(v) for coord in xyz for v in coord)
        assert all(math.isfinite(v) for v in charges)
        result[name] = dict(coordinates=xyz, types=[a[5] for a in atoms], charges=charges)
    return result


def excel_ids(data):
    ns = {'s': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    with ZipFile(io.BytesIO(data)) as z:
        strings = []
        if 'xl/sharedStrings.xml' in z.namelist():
            strings = [''.join(el.itertext()) for el in ET.fromstring(z.read('xl/sharedStrings.xml'))]
        def value(cell):
            if cell.get('t') == 'inlineStr':
                return ''.join(t.text or '' for t in cell.findall('.//s:t', ns))
            raw = cell.findtext('s:v', default='', namespaces=ns)
            return strings[int(raw)] if cell.get('t') == 's' else raw
        # Resolve the Compounds worksheet by its workbook relationship.
        book = ET.fromstring(z.read('xl/workbook.xml'))
        sheet = next(s for s in book.findall('s:sheets/s:sheet', ns) if s.get('name') == 'Compounds')
        rid = sheet.get('{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id')
        rel = next(r for r in ET.fromstring(z.read('xl/_rels/workbook.xml.rels')) if r.get('Id') == rid)
        target = rel.get('Target')
        path = target.lstrip('/') if target.startswith('/') else 'xl/' + target
        records = ET.fromstring(z.read(path)).findall('s:sheetData/s:row', ns)
        header = records[0].findall('s:c', ns)
        column = next(c.get('r').rstrip('0123456789') for c in header if value(c) == 'compound_id')
        return [value(next(c for c in row.findall('s:c', ns)
                           if c.get('r').rstrip('0123456789') == column)) for row in records[1:]]


def main():
    compound_rows = rows((ROOT / 'COMPOUND_TABLE.csv').read_bytes())
    compounds = index(compound_rows)
    workbook_ids = excel_ids((ROOT / 'COMPOUND_TABLE.xlsx').read_bytes())
    assert workbook_ids == [r['compound_id'] for r in compound_rows]
    assert len(compounds) == len(workbook_ids) == len(set(workbook_ids)) == 1240
    assert Counter(r['label'] for r in compound_rows) == {'active': 40, 'decoy': 1200}
    diagnostics = rows((ROOT / 'annotation_diagnostics/pose_comparison.csv').read_bytes())
    diagnostic_index = index(diagnostics)
    assert set(diagnostic_index) == set(compounds)
    report = []
    with ZipFile(ROOT / 'evidence/ACHE_single_matched_live_execution.zip') as z:
        def read(name): return z.read(PREFIX + name)
        for filename in ['per_compound_comparison.csv', 'archived_pr.csv', 'archived_roc.csv',
                         'live_pr.csv', 'live_roc.csv', 'comparison.json', 'score_discrepancies.csv']:
            assert (ROOT / 'figure6_source' / filename).read_bytes() == read('analysis/' + filename)
        discrepancy_table = rows(read('analysis/score_discrepancies.csv'))
        qualifying = [r for r in compound_rows if r['status'] == 'evaluable' and
                      abs(float(r['live_affinity']) - float(r['archived_affinity'])) > 1e-12]
        assert len(discrepancy_table) == len(qualifying) == 0
        assert 'absolute_score_delta' in read('analysis/score_discrepancies.csv').decode().splitlines()[0]
        assert (ROOT / 'annotation_diagnostics/pose_comparison.csv').read_bytes() == read(DIAGNOSTIC + 'pose_comparison.csv')
        assert (ROOT / 'annotation_diagnostics/summary.json').read_bytes() == read(DIAGNOSTIC + 'summary.json')
        assert set(index(rows(read('inputs/labels.csv')))) == set(compounds)
        assert set(index(rows(read('live/screening_results.csv')))) == set(compounds)
        companion = mol2(read('diagnostics/companion_reference/ACHE-test-docked-PETS.mol2'))
        assert set(companion) == set(compounds)
        original_summary = {r['library']: r for r in json.loads(read(DIAGNOSTIC + 'summary.json'))['libraries']}
        for batch, kind, label, expected_n in [(1, 'actives', 'active', 40), (2, 'decoys', 'decoy', 1200)]:
            live_bytes = read(f'live/batch_{batch:03d}/poses.mol2')
            assert live_bytes == read(DIAGNOSTIC + f'{kind}_poses_snapshot.mol2')
            live = mol2(live_bytes)
            archived = mol2(read(f'reference/ACHE-DEKOIS-{kind}-docked.mol2'))
            expected_ids = {cid for cid, r in compounds.items() if r['label'] == label}
            assert set(live) == set(archived) == expected_ids and len(live) == expected_n
            failures = Counter()
            for cid, pose in live.items():
                record = diagnostic_index[cid]
                assert record['library'] == kind
                for key, column in [('coordinates', 'coordinates_exact_match'),
                                    ('types', 'atom_types_match'), ('charges', 'charges_match')]:
                    match = pose[key] == archived[cid][key]
                    assert record[column] == str(match), (cid, column)
                    failures[column] += not match
                for key, column in [('coordinates', 'companion_coordinates_match'),
                                    ('types', 'companion_atom_types_match'),
                                    ('charges', 'companion_charges_match')]:
                    match = pose[key] == companion[cid][key]
                    assert record[column] == str(match) and match, (cid, column)
            summary = original_summary[kind]
            assert summary['complete_readable_records'] == expected_n and not summary['unreadable_records']
            for field, column in [('coordinates_exact_matches', 'coordinates_exact_match'),
                                  ('atom_type_matches', 'atom_types_match'), ('charge_matches', 'charges_match')]:
                assert summary[field] == expected_n - failures[column]
            assert all(summary[field] == expected_n for field in
                       ['companion_coordinates_matches', 'companion_atom_type_matches', 'companion_charge_matches'])
            report.append(dict(library=kind, unique_compounds=len(live),
                               coordinate_discrepancies=failures['coordinates_exact_match'],
                               atom_type_discrepancies=failures['atom_types_match'],
                               partial_charge_discrepancies=failures['charges_match'],
                               companion_all_annotations_match=True))
        assert report[0]['atom_type_discrepancies'] == report[0]['partial_charge_discrepancies'] == 0
        assert report[1]['atom_type_discrepancies'] == 101
        assert report[1]['partial_charge_discrepancies'] == 172
        assert all(r['coordinate_discrepancies'] == 0 for r in report)
    print(json.dumps(dict(status='PASS', csv_rows=len(compounds), csv_unique_ids=len(compounds),
                         excel_rows=len(workbook_ids), excel_unique_ids=len(set(workbook_ids)),
                         diagnostic_rows=len(diagnostics), diagnostic_unique_ids=len(diagnostic_index),
                         populations_identical=True, libraries=report,
                         score_discrepancy_rows=0, score_discrepancy_filter='evaluable and abs(live_affinity - archived_affinity) > 1e-12 kcal/mol',
                         figure6_source_files_byte_identical=True,
                         method='Exact saved atom-order numerical coordinate/charge and atom-type sequence comparisons; no fitting or remapping.',
                         original_diagnostic_code=PREFIX + 'scripts/inspect_saved_poses.py'), indent=2))


if __name__ == '__main__':
    main()
