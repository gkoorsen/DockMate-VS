# Case-study reproducibility materials

The complete packages accompany the ESR1 and ACHE examples in the SoftwareX
manuscript. They are retained in the repository independently of whether the
journal accepts them as supplementary uploads. Both have been prepared for
submission; an actual journal upload is not asserted.

| Package | Contents | Checksum |
|---|---|---|
| [S1: ESR1](packages/DockMate-VS_ESR1_Supplement_S1.zip) | Historical campaign settings, 880-compound result data, reconstruction scripts, screening statistics and MW-controlled Figure 5/Table 3 analysis | [SHA-256](packages/DockMate-VS_ESR1_Supplement_S1.zip.sha256) |
| [S2: ACHE](packages/DockMate-VS_ACHE_Reproduction_S2.zip) | Frozen protocol, prepared inputs, executed source, executable/image identity, raw outputs, all 1,240 compounds, annotation diagnostics and independent Table 2/Figure 6 analysis | [SHA-256](packages/DockMate-VS_ACHE_Reproduction_S2.zip.sha256) |

The existing v0.1.3 software DOI does not identify these later evidence
packages. The v0.1.4 source release includes these ZIPs and checksum files
in the [Zenodo archive](https://doi.org/10.5281/zenodo.23744615). Historical software
versions and executed-source identities inside the evidence must remain intact.

## Verify the downloaded archives

From the directory containing the ZIPs and their checksum files:

```sh
shasum -a 256 -c DockMate-VS_ESR1_Supplement_S1.zip.sha256
shasum -a 256 -c DockMate-VS_ACHE_Reproduction_S2.zip.sha256
```

Extract each archive into its own working directory. Install the analysis
requirements specified by that package's README; the historical S1 and S2
analysis environments differ.

### S1: ESR1

From the extracted S1 directory:

```sh
python verify.py
python verification/verify_mw_controlled.py
python -m unittest discover -s verification -p test_verify_mw_controlled.py -v
```

The size-controlled analysis of Section 3.2 is included as
`publication/tables/esr1_mw_controlled.csv` with its generating script and a
verification command ([MW_CONTROLLED.md](esr1/MW_CONTROLLED.md)). This
2,000-replicate MW analysis is separate from the earlier 10,000-replicate
screening bootstrap. The MW script writes CSV values supporting Figure 5 and
Table 3; it does not render a figure. The verifier rejects packaged-table
mismatches, non-finite expected numbers and unexpected non-missing intervals.

S1 does not redistribute LIT-PCBA structures. Its README explains how to obtain
the source archive and reconstruct the subset, and how to reproduce the other
screening statistics without rerunning historical docking.

### S2: ACHE

From the extracted S2 directory:

```sh
python verify_archive.py
python verify_annotations.py
python regenerate.py --output ../ache-reanalysis
```

These commands verify package integrity, all 1,240 unique compound IDs and the
annotation comparisons, then independently regenerate Table 2 and the final
two-panel Figure 6 from native outputs. Figure 6 comprises precision-recall
curves and class-specific, count-sized affinity agreement; ROC data remain
available for audit but are not a panel in the final figure.

`analysis/score_discrepancies.csv` inside the sealed evidence is intentionally
header-only: none of the 1,240 evaluable compounds passes the absolute score
difference filter of greater than 1e-12 kcal/mol. Atom-type and partial-charge
differences are documented separately in the complete diagnostic table.

## Repository layout

- `packages/` contains the complete, checksum-identified distribution archives.
- `esr1/` contains the maintained S1 analysis, reconstruction and verification
  sources, plus the expected MW results needed for verifier regression tests.
- `ache/` contains the maintained S2 plotting, regeneration and audit sources.

The source directories are not extracted datasets; run the documented full
analyses from extracted packages. The package READMEs describe the exact
included evidence and limitations. The manuscript, its standalone figures,
old packaging attempts and local sign-off bundles are kept outside the
Git-tracked material. The complete supplement ZIPs intentionally retain their
supporting figures and frozen provenance.
