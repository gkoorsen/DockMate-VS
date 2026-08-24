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
