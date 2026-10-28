"""Verify the distributed supplement, frozen nested evidence and entire CSV population.
Run from any working directory: python verify_archive.py
No docking, preparation, network access, or third-party Python packages are used.
"""
from pathlib import Path
from collections import Counter
import csv, hashlib, io, json, zipfile

ROOT = Path(__file__).resolve().parent
def sha(data): return hashlib.sha256(data).hexdigest()
def rows(data): return list(csv.DictReader(io.StringIO(data.decode('utf-8-sig'))))

manifest=json.loads((ROOT/'SHA256SUMS.json').read_text())
for name,expected in manifest.items():
    data=(ROOT/name).read_bytes()
    assert sha(data)==expected['sha256'],name
    assert len(data)==expected['bytes'],name
proof=json.loads((ROOT/'TABLE_COMPLETENESS.json').read_text())
evidence=ROOT/'evidence/ACHE_single_matched_live_execution.zip'
assert sha(evidence.read_bytes())==proof['original_evidence_zip_sha256']
with zipfile.ZipFile(evidence) as z:
    prefix='ACHE_matched_live/'
    nested=json.loads(z.read(prefix+'SHA256SUMS.json'))
    assert set(z.namelist())=={prefix+p for p in nested}|{prefix+'SHA256SUMS.json'}
    for name,expected in nested.items():
        data=z.read(prefix+name)
        assert sha(data)==expected['sha256'],name
        assert len(data)==expected['bytes'],name
    for filename in ['frozen_execution_hashes.json','reference_hashes.json']:
        for name,expected in json.loads(z.read(prefix+filename)).items():
            assert sha(z.read(prefix+name))==expected,name
    compound_bytes=z.read(prefix+'analysis/per_compound_comparison.csv')
    assert (ROOT/'COMPOUND_TABLE.csv').read_bytes()==compound_bytes
    assert (ROOT/'COMPOUND_TABLE.xlsx').read_bytes()==z.read(prefix+'analysis/matched_execution_comparison.xlsx')
    compounds=rows(compound_bytes)
    native=rows(z.read(prefix+'live/screening_results.csv'))
    labels=rows(z.read(prefix+'inputs/labels.csv'))
    expected_ids=[r['compound_id'] for r in labels]
    ids=[r['compound_id'] for r in compounds]
    assert len(ids)==len(set(ids))==1240
    assert set(ids)==set(expected_ids)=={r['compound_id'] for r in native}
    assert len(native)==len(labels)==1240
    expected_labels={r['compound_id']:r['label'] for r in labels}
    assert all(r['label']==expected_labels[r['compound_id']] for r in compounds)
    assert Counter(r['label'] for r in compounds)=={'active':40,'decoy':1200}
    assert all(r['status']=='evaluable' for r in compounds)
    assert all(r['docking_completed']=='True' and r['num_poses']=='1' for r in native)
    assert ids[999]=='ZINC11499592' and ids[1000]=='ZINC20556022' and ids[1239]=='ZINC21861611'
    assert z.testzip() is None
print(json.dumps({'status':'PASS','supplement_files_verified':len(manifest),
                  'frozen_evidence_files_verified':len(nested),'compound_rows':len(ids),
                  'unique_compounds':len(set(ids)),'actives':40,'decoys':1200,
                  'rows_after_zero_based_index_999':len(ids[1000:]),
                  'last_compound':ids[-1]},indent=2))
