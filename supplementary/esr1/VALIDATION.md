# Supplement validation

Validated locally on 2026-09-10 with Python 3.9.15, NumPy 1.23.5, pandas 1.5.3,
SciPy 1.11.2, scikit-learn 1.1.1, Matplotlib 3.8.2, openpyxl 3.1.2 and
RDKit 2022.03.5 on macOS/ARM64.

- All 35 original frozen-evidence SHA-256 checksums matched before packaging.
- The separately obtained LIT-PCBA source archive matched its recorded SHA-256.
- Reconstruction reproduced all cleaning counts, including 80 actives and
  3836 inactives after removal of eight conflicting structures.
- The 80/800 subset and its exact screening order matched the source identifiers.
- The reconstructed compound CSV was byte-identical to the original CSV:
  `adbc76754325f6cbce6dd2044341c880b0abe0e0a75d2fe57da32d2fb239aa86`.
- The reconstructed screening workbook's Docking_Jobs sheet matched the CSV.
- Analysis using only the distributed structure-free result records reproduced
  the screening metrics, all bootstrap confidence intervals, descriptor
  diagnostics and all ten property-model fold AUCs to a tolerance of 1e-12.
- Six focused archive-input tests passed: exact-member reads, rejection of
  symbolic links, hard links, duplicate members, missing members and an incorrect
  source-archive checksum.

Validation outputs containing reconstructed structures were written outside the
supplement. They are not included in the ZIP. No new docking calculation or
independent predictive validation was performed. No external upload or manuscript
submission was performed.

The historical docking used earlier code and tool versions; this supplement
validates reconstruction of inputs and reanalysis of saved results, not exact
reexecution of those docking calculations with DockMate-VS 0.1.1.
