"""Molecular-weight-controlled diagnostics for the ESR1 screening example.

Reads the frozen screening results of the run reported in the manuscript and quantifies
how much of the reported ROC AUC survives control for molecular size.

Usage:
    python esr1_mw_controlled.py path/to/redock_results.csv

Output: esr1_mw_controlled.csv (written in the current working directory).
This script computes the numerical results supporting Figure 5 and Table 3;
it does not render a figure.
"""
import sys
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score, average_precision_score
import matplotlib.pyplot as plt

PATH = sys.argv[1] if len(sys.argv) > 1 else "redock_results.csv"
df = pd.read_csv(PATH)
y = df.control_label.to_numpy(int)
s = -df.best_score.to_numpy(float)          # higher = more favourable docking score
mw = df.molecular_weight.to_numpy(float)

auc_score = roc_auc_score(y, s)
ap = average_precision_score(y, s)
auc_mw = roc_auc_score(y, mw)
resid = s - np.polyval(np.polyfit(mw, s, 1), mw)
auc_resid = roc_auc_score(y, resid)
rs, rmw = stats.rankdata(s), stats.rankdata(mw)
resid_rank = rs - np.polyval(np.polyfit(rmw, rs, 1), rmw)
auc_resid_rank = roc_auc_score(y, resid_rank)

act, inact = np.where(y == 1)[0], np.where(y == 0)[0]


def match(seed, caliper=None):
    """Nearest-molecular-weight inactive for each active, without replacement."""
    r = np.random.default_rng(seed)
    pool, pairs = list(inact), []
    for i in r.permutation(act):
        d = np.abs(mw[pool] - mw[i])
        j = int(np.argmin(d))
        if caliper is None or d[j] <= caliper:
            pairs.append((int(i), pool.pop(j)))
    a_idx = np.array([p[0] for p in pairs]); c_idx = np.array([p[1] for p in pairs])
    idx = np.concatenate([a_idx, c_idx])
    return roc_auc_score(y[idx], s[idx]), float(np.abs(mw[a_idx] - mw[c_idx]).mean()), a_idx, c_idx


nn = [match(1000 + k) for k in range(1000)]
nn_auc = np.array([x[0] for x in nn]); nn_gap = float(np.mean([x[1] for x in nn]))
_, _, a0, c0 = nn[0]
rb = np.random.default_rng(7)
boot = np.array([roc_auc_score(y[np.concatenate([a0[k], c0[k]])], s[np.concatenate([a0[k], c0[k]])])
                 for k in (rb.integers(0, len(a0), len(a0)) for _ in range(2000))])


def strat_boot(yv, sv, reps=2000, seed=11):
    r = np.random.default_rng(seed)
    a, b = np.where(yv == 1)[0], np.where(yv == 0)[0]
    out = [roc_auc_score(yv[i], sv[i]) for i in
           (np.concatenate([r.choice(a, len(a), True), r.choice(b, len(b), True)]) for _ in range(reps))]
    return float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))


q = np.quantile(mw, [1 / 3, 2 / 3]); tert = np.digitize(mw, q)
rows = []
for name, val, sv, note in [
    ("Docking score (as reported)", auc_score, s, f"average precision {ap:.4f}"),
    ("Molecular weight alone", auc_mw, mw, "MW used directly as the ranking score"),
    ("Score, MW dependence removed (linear)", auc_resid, resid, "residual of -score on MW"),
    ("Score, MW dependence removed (rank)", auc_resid_rank, resid_rank, "residual of rank(-score) on rank(MW)"),
]:
    lo, hi = strat_boot(y, sv)
    rows.append(dict(analysis=name, auc=round(val, 4), ci_low=round(lo, 4), ci_high=round(hi, 4),
                     n=len(y), n_active=int(y.sum()), note=note))
rows.append(dict(analysis="MW-matched active/inactive pairs", auc=round(float(nn_auc.mean()), 4),
                 ci_low=round(float(np.percentile(boot, 2.5)), 4),
                 ci_high=round(float(np.percentile(boot, 97.5)), 4), n=2 * len(a0), n_active=len(a0),
                 note=f"mean |dMW| {nn_gap:.1f} Da; mean over 1000 matchings, CI from 2000 pair bootstraps"))
tert_rows, tert_ci = [], []
for k in range(3):
    m = tert == k
    lo, hi = strat_boot(y[m], s[m], seed=20 + k)
    tert_rows.append(dict(stratum=f"MW tertile {k+1}", mw_range=f"{mw[m].min():.0f}-{mw[m].max():.0f}",
                          n=int(m.sum()), n_active=int(y[m].sum()), auc=roc_auc_score(y[m], s[m])))
    tert_ci.append((lo, hi))
    rows.append(dict(analysis=f"Docking score within MW tertile {k+1} ({tert_rows[-1]['mw_range']} Da)",
                     auc=round(tert_rows[-1]["auc"], 4), ci_low=round(lo, 4), ci_high=round(hi, 4),
                     n=tert_rows[-1]["n"], n_active=tert_rows[-1]["n_active"], note="within-stratum AUC"))


def partial_spearman(a, b, c):
    ra, rb_, rc = stats.rankdata(a), stats.rankdata(b), stats.rankdata(c)
    rab, rac, rbc = np.corrcoef(ra, rb_)[0, 1], np.corrcoef(ra, rc)[0, 1], np.corrcoef(rb_, rc)[0, 1]
    return (rab - rac * rbc) / np.sqrt((1 - rac ** 2) * (1 - rbc ** 2))


rows.append(dict(analysis="Spearman rho, score vs MW", auc=round(stats.spearmanr(df.best_score, mw).statistic, 4),
                 ci_low=np.nan, ci_high=np.nan, n=len(y), n_active=int(y.sum()),
                 note="negative because lower score = more favourable"))
rows.append(dict(analysis="Partial Spearman, -score vs activity | MW",
                 auc=round(partial_spearman(s, y, mw), 4), ci_low=np.nan, ci_high=np.nan,
                 n=len(y), n_active=int(y.sum()),
                 note=f"unadjusted rho = {stats.spearmanr(s, y).statistic:.4f}"))
pd.DataFrame(rows).to_csv("esr1_mw_controlled.csv", index=False)
print(pd.DataFrame(rows).to_string(index=False, max_colwidth=46))
