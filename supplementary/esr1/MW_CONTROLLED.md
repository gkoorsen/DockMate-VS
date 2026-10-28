# Molecular-weight-controlled diagnostics (Figure 5 and Table 3)

This part of the supplement supports the size-controlled reanalysis of the ESR1 screening
subset. It is separate from the screening and property-model analysis documented in
`README.md`, whose bootstrap uses 10,000 replicates on different quantities; nothing in that
analysis verifies the values reported in Figure 5 or Table 3.

## Files

| File | Contents |
|---|---|
| `analysis/esr1_mw_controlled.py` | generates every value in Figure 5 and Table 3 from `results/redock_results.csv` |
| `publication/tables/esr1_mw_controlled.csv` | the generated table, as printed in Table 3 |
| `verification/verify_mw_controlled.py` | regenerates the table and compares it against the manuscript values |
| `verification/mw_controlled_expected.json` | the manuscript values the verifier checks against |

## Commands

From the package root, in the environment described under "Environment" in `README.md`
(`pip install -r analysis-requirements.txt`; NumPy, pandas, SciPy, scikit-learn and Matplotlib
are required):

```bash
# regenerate the table
python analysis/esr1_mw_controlled.py results/redock_results.csv

# verify it against the printed Figure 5 and Table 3 values
python verification/verify_mw_controlled.py
```

The verifier exits 0 only if every value matches at the printed four-decimal precision.
It fails on a packaged-table mismatch, non-finite values where numbers are expected,
or non-missing values where the manuscript expects a missing interval.
Running it on this package reports `rows checked: 10` and `RESULT: all values match the
manuscript`.

The analysis script writes only `esr1_mw_controlled.csv` in the current working
directory. It calculates the values supporting Figure 5 and Table 3; it does not
render a PNG. Figure 6 refers to the separate ACHE reproduction example.

## What each analysis does

| Row in Table 3 | Method |
|---|---|
| Docking score (as reported) | ROC AUC of the negated docking score against the activity label |
| Molecular weight alone | ROC AUC using molecular weight directly as the ranking score |
| Score, MW dependence removed (linear) | ROC AUC of the residual of the negated score regressed on molecular weight (least squares) |
| Score, MW dependence removed (rank) | the same on ranks: residual of rank(-score) on rank(MW) |
| MW-matched active/inactive pairs | nearest-molecular-weight inactive per active, sampled without replacement; point estimate is the mean over 1,000 matchings, interval from 2,000 pair bootstraps of one representative matching |
| Docking score within MW tertile 1-3 | ROC AUC recomputed inside each molecular-weight tertile, with the tertile range and active count reported |
| Spearman rho, score vs MW | Spearman correlation between the docking score and molecular weight |
| Partial Spearman | Spearman correlation of the negated score with activity, controlling for molecular weight, with the unadjusted value in the note |

Confidence intervals for the ROC AUC rows are stratified bootstrap percentile intervals with
2,000 replicates (`strat_boot`, seed 11), matching the Figure 5 caption. Seeds are fixed in the
script, so repeated runs reproduce the published values exactly rather than approximately.

## Scope

These diagnostics are computed from the archived per-compound results in this package. They do
not require the LIT-PCBA structures, and they do not re-run docking.
