# Examples

## ESR1 native-pose protocol development

`esr1_protocol_development_1XP1.xlsx` is a minimal public example for the
Protocol Development tab. It contains the ESR1 structure 1XP1 and its native
ligand AIH. The application retrieves the public structure during preflight.

Suggested exploratory settings:

- engine: Vina;
- box: 6 A margin around AIH;
- water handling: retain all;
- exhaustiveness: 8;
- modes: 20;
- seed: 42;
- maximum tautomers: 2;
- maximum conformers: 2.

Scores and exact poses can vary with docking-engine version, dependency builds,
CPU architecture, and preparation tools. The example is intended to verify the
workflow and inspect pose-recovery outputs, not to reproduce a bit-identical
score.

The workbook does not redistribute LIT-PCBA assay compounds. Dataset downloads
must follow the source provider's terms.

## ESR1 enrichment smoke test

`esr1_enrichment_smoke_test_1XP1.xlsx` is a separate, runnable example for the
Screening tab. It docks the native 1XP1 ligand as one positive control and six
illustrative negative controls whose structures and PubChem CIDs are recorded
in the workbook. Their explicit `label` values exercise control
classification, ranking, ROC AUC, and enrichment reporting.

This deliberately small workbook is a software smoke test, not a docking-
protocol validation benchmark. The negative labels are illustrative test
labels and are not claims of biological inactivity at ESR1. Results therefore
show that the enrichment workflow operates, not that the selected protocol is
scientifically validated. The workbook contains no LIT-PCBA assay records.

## Headless campaign files

- `campaign.protocol.yml` runs a small protocol sweep on the example workbook.
- `campaign.screen.yml` runs the enrichment smoke-test workbook with a frozen
  single protocol.

Run either with the local environment (`dockmate-vs protocol|screen --config`)
or the optional core container (`scripts/dockmate-docker protocol|screen`).
