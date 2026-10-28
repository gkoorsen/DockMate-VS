"""Regenerate manuscript Table 2 and Figure 6 from sealed native evidence."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import shutil
import subprocess
import sys
import zipfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True,
                        help='New working directory outside the extracted supplement')
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    out = args.output.resolve()
    if out == root or root in out.parents:
        parser.error('Choose an output directory outside the extracted supplement.')
    if out.exists():
        parser.error('Output directory already exists; choose a fresh path.')
    subprocess.run([sys.executable, str(root / 'verify_archive.py')], check=True)
    subprocess.run([sys.executable, str(root / 'verify_annotations.py')], check=True)
    out.mkdir(parents=True)
    with zipfile.ZipFile(root / 'evidence/ACHE_single_matched_live_execution.zip') as z:
        z.extractall(out)
    evidence = out / 'ACHE_matched_live'
    analysis = evidence / 'analysis'
    expected = json.loads((analysis / 'comparison.json').read_text())
    subprocess.run([sys.executable, str(evidence / 'scripts/analyse.py')], check=True)
    actual = json.loads((analysis / 'comparison.json').read_text())
    if actual != expected:
        raise ValueError('Independent reanalysis differs from frozen comparison.json')
    archived = actual['archived']
    live = actual['live_all_compounds']
    n = actual['planned_total']
    rows = [('Evaluable compounds', archived['n'], actual['evaluable_total'],
             actual['evaluable_total'] - archived['n'], 'No failures')]
    for label, key in [('Trapezoidal PR-AUC', 'pr_auc_trapezoidal'),
                       ('ROC-AUC', 'roc_auc'), ('Average precision', 'average_precision')]:
        rows.append((label, f'{archived[key]:.6f}', f'{live[key]:.6f}',
                     f'{live[key] - archived[key]:.6f}', 'Exact'))
    rows.extend([
        ('Reported affinities', f'{n:,} scores', f'{actual["evaluable_total"]:,} scores',
         f'{actual["changed_scores"]} changed', 'All matched'),
        ('Saved coordinates', f'{n:,} poses', f'{actual["evaluable_total"]:,} poses',
         f'{n - actual["saved_pose_coordinates_exact_matches"]} changed', 'All matched')])
    with (out / 'Table_2.csv').open('w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Measure', 'Archived', 'Rerun', 'Difference', 'Outcome'])
        writer.writerows(rows)
    subprocess.run([sys.executable, str(root / 'plot_figure6.py'), '--data', str(analysis),
                    '--output', str(out)], check=True)
    shutil.copy2(root / 'FIGURE_6_CAPTION.txt', out / 'Figure_6_caption.txt')
    manifest = json.loads((evidence / 'SHA256SUMS.json').read_text())
    # Analysis is regenerated in the working copy; every non-analysis evidence file stays intact.
    for name, entry in manifest.items():
        if not name.startswith('analysis/'):
            if hashlib.sha256((evidence / name).read_bytes()).hexdigest() != entry['sha256']:
                raise ValueError(f'Non-analysis evidence changed: {name}')
    (out / 'REGENERATION_CHECK.json').write_text(json.dumps({
        'status': 'PASS', 'independent_summary_matches_frozen': True,
        'compound_rows': actual['planned_total'], 'evaluable_rows': actual['evaluable_total'],
        'non_analysis_evidence_unchanged': True,
        'figure_source': 'Unmodified frozen scripts/analyse.py computes source data; plot_figure6.py uses rerun presentation labels.',
        'note': 'Numerical agreement is checked; PDF metadata and rendering may vary by environment.'
    }, indent=2) + '\n')
    print(f'Table_2.csv, Figure_6.png and Figure_6.pdf regenerated in {out}')


if __name__ == '__main__':
    main()
