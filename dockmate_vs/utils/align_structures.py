"""
Superpose an ensemble of structures into a common coordinate frame.

Run this BEFORE preparing receptors for ensemble docking. Ensemble docking passes one
grid box to every receptor, which is only valid if they share a coordinate frame -
structures downloaded from the PDB generally do not.

Usage
-----
    python -m dockmate_vs.utils.align_structures \\
        --structures 6RB1.pdb 3BQC.pdb 5M4F.pdb \\
        --out-dir structures_aligned

    # ligand extents from the screening library give a better box size
    python -m dockmate_vs.utils.align_structures \\
        --structures 6RB1.pdb 3BQC.pdb 5M4F.pdb \\
        --out-dir structures_aligned \\
        --reference 6RB1.pdb \\
        --comp-ids 6RB1=JWQ 3BQC=EMO 5M4F=7FC \\
        --ligand-extent 18.1

Then prepare receptors from the aligned copies and pass the printed --grid-center /
--grid-size to run_ensemble_excel_production.py.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from .structure_alignment import (
    contact_residues,
    derive_grid,
    superpose_ensemble,
    validate_grid,
)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Superpose structures into a common frame and derive a docking grid",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('--structures', nargs='+', required=True,
                        help='Structure files (PDB or mmCIF)')
    parser.add_argument('--out-dir', required=True,
                        help='Directory for aligned copies')
    parser.add_argument('--reference', default=None,
                        help='Reference structure. Defaults to the first. Prefer a '
                             'wild-type, high-resolution, single-chain entry.')
    parser.add_argument('--comp-ids', nargs='*', default=None, metavar='STEM=COMP',
                        help='Explicit ligand component per structure, e.g. 6RB1=JWQ. '
                             'Without this, the largest non-additive component is used.')
    parser.add_argument('--ligand-extent', type=float, default=None,
                        help='Longest principal extent (A) of the largest compound to be '
                             'docked. Sizes the box for the screening library rather than '
                             'the co-crystal ligands.')
    parser.add_argument('--padding', type=float, default=6.0,
                        help='Padding added to the box edge')
    parser.add_argument('--site-residues', nargs='*', type=int, default=None,
                        help='Known binding-site residue numbers. Each ligand is checked '
                             'against these; a member contacting none is flagged as '
                             'occupying a different pocket.')
    parser.add_argument('--report', default=None,
                        help='Write alignment provenance to this CSV')

    args = parser.parse_args(argv)

    comp_ids = {}
    for item in (args.comp_ids or []):
        if '=' not in item:
            parser.error(f"--comp-ids entries must be STEM=COMP, got {item!r}")
        stem, comp = item.split('=', 1)
        comp_ids[stem] = comp
        comp_ids[f"{stem}_aligned"] = comp

    structures = [Path(s) for s in args.structures]
    missing = [str(s) for s in structures if not s.exists()]
    if missing:
        parser.error(f"structures not found: {', '.join(missing)}")
    if len(structures) < 2:
        parser.error("need at least two structures to superpose")

    print(f"Superposing {len(structures)} structures "
          f"(reference: {Path(args.reference).name if args.reference else structures[0].name})")

    aligned, report = superpose_ensemble(
        structures, args.out_dir, reference=args.reference, comp_ids=comp_ids
    )

    print("\nAlignment:")
    for row in report:
        tag = "  (reference)" if row['is_reference'] else ""
        print(f"  {row['structure']:<24} RMSD {row['rmsd']:.2f} A over "
              f"{row['n_aligned']} residues, chain {row['chain']}{tag}")

    spec = derive_grid(
        aligned, comp_ids=comp_ids,
        padding=args.padding, ligand_extent=args.ligand_extent,
    )

    print("\nBinding sites in the common frame:")
    for site in spec.sites:
        copies = f"  [copy 1 of {site.n_copies}]" if site.n_copies > 1 else ""
        print(f"  {site.path.name:<24} {site.comp_id:<6} "
              f"({site.center[0]:8.2f}, {site.center[1]:8.2f}, {site.center[2]:8.2f}){copies}")

    ok, rows = validate_grid(aligned, spec.center, spec.size, comp_ids=comp_ids)

    print(f"\nConsensus centre: {spec.center[0]:.2f} {spec.center[1]:.2f} {spec.center[2]:.2f}")
    print(f"Max deviation:    {spec.max_deviation} A")
    print(f"Box size:         {spec.size[0]:.0f} A")
    print(f"Grid valid for all members: {'yes' if ok else 'NO'}")

    if not ok:
        print("\nSome sites still fall outside the box after superposition.", file=sys.stderr)
        print("This usually means one ligand occupies a different pocket "
              "(allosteric or surface) rather than the site being screened.", file=sys.stderr)
        for row in rows:
            if row['status'] == 'OUTSIDE':
                print(f"  {row['structure']}: {row['comp_id']} "
                      f"{row['distance']} A from centre", file=sys.stderr)

    if args.site_residues:
        print(f"\nContact check against site residues {args.site_residues}:")
        wanted = set(args.site_residues)
        for path in aligned:
            contacts = contact_residues(path, comp_ids.get(Path(path).stem))
            hit = sorted({num for num, _ in contacts} & wanted)
            verdict = "OK" if hit else "NO SITE CONTACT - different pocket?"
            print(f"  {Path(path).name:<24} {len(hit)}/{len(wanted)} markers {hit}  {verdict}")

    if args.report:
        merged = []
        for rep_row, site, val_row in zip(report, spec.sites, rows):
            merged.append({
                **rep_row,
                'comp_id': site.comp_id,
                'n_copies': site.n_copies,
                'center_x': round(float(site.center[0]), 3),
                'center_y': round(float(site.center[1]), 3),
                'center_z': round(float(site.center[2]), 3),
                'dist_to_consensus': val_row.get('distance'),
                'grid_status': val_row.get('status'),
            })
        with open(args.report, 'w', newline='') as fh:
            writer = csv.DictWriter(fh, fieldnames=list(merged[0].keys()))
            writer.writeheader()
            writer.writerows(merged)
        print(f"\nWrote provenance: {args.report}")

    print("\nNext steps:")
    print(f"  1. Prepare receptors from the aligned copies in {args.out_dir}/")
    print("  2. Run the experiment with:")
    print(f"       --grid-center {spec.center[0]:.2f} {spec.center[1]:.2f} {spec.center[2]:.2f} \\")
    print(f"       --grid-size {spec.size[0]:.0f} {spec.size[1]:.0f} {spec.size[2]:.0f}")

    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
