# Suggested manuscript text

Use the following availability wording after the ZIP is published on GitHub.
Replace the branch-based supplement link with its commit-specific permalink
for the submitted manuscript. The supplement is separate from the Zenodo
software archive.

## Data and code availability

DockMate-VS is available under the MIT License at
https://github.com/gkoorsen/DockMate-VS. Version 0.1.1 is archived at
https://doi.org/10.5281/zenodo.22685196; all versions are accessible through
https://doi.org/10.5281/zenodo.22283781. The reconstruction-based ESR1 supplement
is available at
https://github.com/gkoorsen/DockMate-VS/blob/master/supplementary/packages/DockMate-VS_ESR1_Supplement_S1.zip
and provides the
ordered ESR1 subset identifiers, source checksums, campaign settings,
structure-free result tables, and scripts for reconstructing the screening
inputs and reproducing the statistical analyses, including the bootstrap
intervals and property-only model. LIT-PCBA compound structures are not
redistributed: users obtain the original archive from its provider [6] and
reconstruct the inputs locally. Receptor and reference-ligand structures are
available from the Protein Data Bank under accession 1XP1 [18].

## Case-study version clarification

The ESR1 docking calculations were generated with an earlier development
version using the then-available adaptive ligand-variant selection; the saved
results were subsequently inspected with the version 0.1.1 GUI. Historical
settings and environment information are retained in the ESR1 supplement.

## Property-only diagnostic methods

A property-only logistic regression used molecular weight, logP, topological
polar surface area, rotatable-bond count, and ligand charge recorded during
preparation. Features were standardized within each training fold, and the
classifier was evaluated using 10-fold stratified, shuffled cross-validation
(seed 42). The model used L2 regularization with C=1 and the lbfgs solver; mean
ROC AUC was calculated across the ten held-out folds.
