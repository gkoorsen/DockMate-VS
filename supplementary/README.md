# Manuscript supplements

## Download the ESR1 supplement

The complete reconstruction-based evidence package is available here:

- [DockMate-VS_ESR1_Supplement_S1.zip](packages/DockMate-VS_ESR1_Supplement_S1.zip)
- [SHA-256 checksum](packages/DockMate-VS_ESR1_Supplement_S1.zip.sha256)

Download both files, then run `shasum -a 256 -c DockMate-VS_ESR1_Supplement_S1.zip.sha256`
in their directory. Extract the ZIP and follow its README for reconstruction
and statistical analysis. The package excludes compound structures and poses.
It is hosted in this GitHub repository, not in the Zenodo 0.1.1 software archive.
For a fixed manuscript reference, use the GitHub file permalink at the commit
that added the package rather than relying on the moving `master` branch.

## Rebuild from frozen evidence

`esr1/` contains the code and documentation templates for the reconstruction-based
SoftwareX ESR1 supplement. This directory is not itself the complete evidence
package: the builder projects selected non-structural fields from the author's
frozen evidence, verifies its checksums, and creates a standalone ZIP.

```sh
python scripts/build_esr1_supplement.py \
  --frozen /path/to/frozen_esr1_seed42 \
  --output dist/esr1-supplement
```

Use a fresh output directory. The builder does not modify the frozen evidence,
rebuild software release 0.1.1, upload to Zenodo, or submit to the journal.
Do not add generated `reconstructed/` inputs to the distribution. The source
archive and original compound workbooks are deliberately excluded.

The built ZIP contains a README with reconstruction and analysis instructions,
the exact ordered subset identifiers, structure-free result records, portable
historical settings, reference statistics, and SHA-256 checksums.
