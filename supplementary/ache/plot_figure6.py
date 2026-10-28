"""Regenerate the final two-panel Figure 6 from complete archived/rerun data.

(a) coincident PR curves and the active-prevalence baseline;
(b) affinities aggregated by reported value within class, with count-sized markers.
Historical source column names are retained. No docking is invoked.
"""
from pathlib import Path
import argparse
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data', type=Path, default=Path(__file__).resolve().parent / 'figure6_source')
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    data = pd.read_csv(args.data / 'per_compound_comparison.csv')
    summary = json.loads((args.data / 'comparison.json').read_text())
    assert len(data) == data.compound_id.nunique() == 1240
    assert data.label.value_counts().to_dict() == {'decoy': 1200, 'active': 40}
    assert np.isfinite(data[['archived_affinity', 'live_affinity']]).all().all()
    assert np.array_equal(data.archived_affinity, data.live_affinity)
    assert np.array_equal(data.archived_affinity, data.archived_affinity.round(1))
    assert data.archived_affinity.nunique() == 101
    pr = pd.read_csv(args.data / 'archived_pr.csv')
    rerun = pd.read_csv(args.data / 'live_pr.csv')
    pd.testing.assert_frame_equal(pr, rerun, check_exact=True)
    grouped = data.groupby(['label', 'archived_affinity', 'live_affinity']).size().rename('count').reset_index()
    peak = grouped.loc[grouped['count'].idxmax()]
    assert peak['label'] == 'decoy' and peak.archived_affinity == -8.3 and peak['count'] == 35
    assert len(data[data.archived_affinity == -8.3]) == 36
    assert len(data[(data.archived_affinity == -8.3) & (data.label == 'active')]) == 1
    rho = data.archived_affinity.rank(method='average').corr(data.live_affinity.rank(method='average'))
    assert np.isclose(rho, 1, rtol=0, atol=1e-12)
    baseline = 40/1240
    area_per_compound = 18.0
    grouped['marker_area_points2'] = grouped['count'] * area_per_compound
    plt.rcParams.update({'font.size': 13, 'axes.labelsize': 15, 'legend.fontsize': 12})
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), layout='constrained', gridspec_kw={'width_ratios':[1.25,1]})
    pale, dark, orange = '#9bcbe0', '#163f73', '#e76c25'
    axes[0].plot(pr.recall, pr.precision, color=pale, lw=5, label='archived scores, recalculated')
    axes[0].plot(rerun.recall, rerun.precision, color=dark, lw=2.3, ls=(0,(4,3)), label='DockMate-VS rerun (coincident)')
    axes[0].axhline(baseline, color='#888888', lw=1.5, ls=':', label='random ranking')
    axes[0].set(xlabel='Recall', ylabel='Precision', xlim=(-.02,1.02), ylim=(-.02,1.05))
    axes[0].legend(loc='upper right', frameon=False)
    lo, hi = -14.7, -1.3
    axes[1].plot([lo,hi], [lo,hi], color='#999999', lw=1.5, label='identity', zorder=1)
    for label, marker, color, count in [('decoy','o',pale,1200),('active','D',orange,40)]:
        g = grouped[grouped.label == label]
        axes[1].scatter(g.archived_affinity, g.live_affinity, s=g.marker_area_points2,
                        marker=marker, color=color, alpha=.85, edgecolors='white', linewidths=.6,
                        label=f'{"decoys" if label == "decoy" else "actives"} (n={count})', zorder=2 if label == 'decoy' else 3)
    for count in [1,10,35]:
        axes[1].scatter([], [], s=area_per_compound*count, facecolors='none', edgecolors='#888888',
                        label=f'{count} compound' + ('' if count == 1 else 's'))
    axes[1].set(xlabel='Archived affinity (kcal/mol)', ylabel='Rerun affinity (kcal/mol)', xlim=(lo,hi), ylim=(lo,hi))
    axes[1].legend(loc='upper left', frameon=False, labelspacing=.55, handletextpad=.4)
    for ax, letter in zip(axes, ['a','b']):
        ax.spines[['top','right']].set_visible(False)
        ax.text(-.16,1.04,letter,transform=ax.transAxes,fontweight='bold',fontsize=20)
    args.output.mkdir(parents=True, exist_ok=True)
    for ext in ['png','pdf']:
        fig.savefig(args.output / ('Figure_6.'+ext), **({'dpi':200} if ext=='png' else {}))
    plt.close(fig)
    grouped.to_csv(args.output / 'Figure_6_affinity_counts.csv', index=False)
    audit = dict(status='PASS', panel_count=2, panels=['precision-recall','class-specific count-sized affinity agreement'],
                 compounds=1240, actives=40, decoys=1200, distinct_reported_affinities=101,
                 class_specific_groups=len(grouped), grouped_count_total=int(grouped['count'].sum()),
                 largest_marker=dict(label='decoy',affinity=-8.3,count=35), compounds_at_minus_8_3=36,
                 actives_at_minus_8_3=1, precision_baseline=baseline, pr_curves_pointwise_identical=True,
                 maximum_absolute_affinity_difference=float((data.live_affinity-data.archived_affinity).abs().max()),
                 average_rank_spearman=float(rho), marker_area_points2_per_compound=area_per_compound,
                 trapezoidal_pr_auc=summary['archived']['pr_auc_trapezoidal'],
                 average_precision=summary['archived']['average_precision'],
                 roc_panel=False, pixel_identity_claimed=False)
    assert round(audit['trapezoidal_pr_auc'],6)==.194743 and round(audit['average_precision'],6)==.195322
    (args.output / 'Figure_6_verification.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps(audit,indent=2))


if __name__ == '__main__':
    main()
