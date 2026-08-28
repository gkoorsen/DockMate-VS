# Examples

The public examples use the DUD-E acetylcholinesterase target (ACES) and the
1E66/HUX co-crystal structure. They demonstrate protocol development followed
by labelled active/decoy screening without presenting DockMate-VS as a new
scoring function.

## ACES protocol development

`dude_aces_protocol_development_1E66.xlsx` contains the native 1E66/HUX
pose-recovery control. `campaign.protocol.yml` evaluates 24 conditions:

- Smina with Vina docking scores;
- 4 A and 6 A co-crystal box margins;
- remove-all, retain-all, and selective water handling;
- exhaustiveness 8 and 16;
- original Vina ranking and Vinardo score-only rescoring;
- seed 42, 20 modes, and at most two tautomers and conformers.

Run locally or with the core container:

```bash
dockmate-vs protocol --config examples/campaign.protocol.yml
scripts/dockmate-docker protocol examples/campaign.protocol.yml
```

The run writes to `results/dude_aces_protocol/`. Use the Summary and Charts
tabs to compare pose recovery, ranking, runtime, and protocol-factor effects.
The chart selector switches all four plots between Top-1, Top-5, and Top-10.
The summary can retain several near-equivalent pose-recovery candidates rather
than choosing a winner from negligible RMSD differences. Compare those
candidates in separate labelled enrichment runs; do not pool their scores into
one AUC.

## ACES screening benchmark

`dude_aces_screening_subset_1E66_seed42.xlsx` contains 20 clustered DUD-E ACES
actives and 200 DUD-E property-matched decoys. The compounds were selected with
a fixed random seed before docking and were not selected using DockMate-VS
scores. All compounds share the 1E66 receptor and HUX-defined site, so the GUI
recognizes this as a single-receptor assay benchmark and displays ROC,
precision-recall, score-distribution, and cumulative-recovery charts.

Run locally or with the core container:

```bash
dockmate-vs screen --config examples/campaign.screen.yml
scripts/dockmate-docker screen examples/campaign.screen.yml
```

The starter screening protocol uses Smina, Vina scoring, a 4 A co-crystal
margin, removed waters, exhaustiveness 8, 20 modes, and seed 42. It can be
amended after inspecting the protocol-development results. The run writes to
`results/dude_aces_screen/`.

## Source and interpretation

- Target page: <https://dude.docking.org/targets/aces>
- DUD-E publication: <https://doi.org/10.1021/jm300687e>
- Receptor structure: <https://www.rcsb.org/structure/1E66>
- Rebuild command: `python scripts/prepare_dude_aces_examples.py`

The workbook Metadata sheets record source checksums, selection parameters,
label semantics, and interpretation limits. DUD-E matched decoys are
physicochemically matched, topologically dissimilar presumed non-binders; they
are not experimentally confirmed inactive compounds. This compact subset is a
reproducible software-workflow example rather than a definitive assessment of
Vina, Smina, Vinardo, DUD-E, or ACES virtual-screening performance.

Scores and exact poses can vary with docking-engine version, dependency builds,
protein preparation, CPU architecture, and random seed. The generated run
manifest records the local execution environment.
