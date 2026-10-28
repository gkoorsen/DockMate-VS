ACHE REPRODUCIBILITY ANALYSIS SOURCES

The complete S2 distribution is ../packages/DockMate-VS_ACHE_Reproduction_S2.zip.
Its checksum file is ../packages/DockMate-VS_ACHE_Reproduction_S2.zip.sha256.
The ZIP includes the data, raw evidence, executed-source snapshot, report and
supporting figures; this directory contains only the maintained code.

EXTRACT THE ZIP BEFORE RUNNING THE FULL ANALYSIS
Install the requirements described in the extracted README, then run from
that extracted package directory:

  python verify_archive.py
  python verify_annotations.py
  python regenerate.py --output ../ache-reanalysis

The first command checks complete file integrity and compound accounting.
The second independently verifies all 1,240 unique compounds and both the
101 atom-type and 172 partial-charge discrepancies against sealed MOL2 data.
The final command regenerates Table 2 and the final two-panel Figure 6 using
native outputs. It does not run docking or prepare new molecules.

plot_figure6.py renders PR curves and class-specific count-sized affinity
agreement. It uses the full compound and curve data supplied in S2. ROC data
remain archived for analysis, but the final figure has no ROC panel.

The header-only score_discrepancies.csv is intentional: zero evaluable
compounds have an absolute archived-rerun score difference >1e-12 kcal/mol.
The annotation differences are recorded separately in the diagnostic CSV.

The full ZIP preserves historical paths/schema names and provenance, including
source/package-version differences, build uncertainty and host-sleep timing.
Its evidence is unchanged by the repository cleanup. Standalone manuscript
files, submission figures and packaging drafts are kept locally, outside Git.
