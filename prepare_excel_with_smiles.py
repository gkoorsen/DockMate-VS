#!/usr/bin/env python3
"""
Helper script to complete Excel template with SMILES.

Extracts ligand SMILES from PDB structures (for entries with PDB codes)
and creates updated Excel with SMILES column ready for validation.

Usage:
    python prepare_excel_with_smiles.py \\
        --excel anathi_template.xlsx \\
        --pdb-dir structures/ \\
        --output anathi_template_complete.xlsx
"""

import argparse
from pathlib import Path
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from loguru import logger


def extract_ligand_from_pdb(pdb_file, ligand_id):
    """
    Extract specific ligand from PDB file and convert to SMILES.
    
    Parameters
    ----------
    pdb_file : Path
        Path to PDB file
    ligand_id : str
        3-letter ligand code (e.g., 'FNO', '3WL', '4WI')
    
    Returns
    -------
    str : SMILES string or None
    """
    try:
        # Read PDB file
        with open(pdb_file, 'r') as f:
            pdb_lines = f.readlines()
        
        # Extract HETATM lines for this ligand
        ligand_lines = []
        for line in pdb_lines:
            if line.startswith('HETATM') and ligand_id in line:
                ligand_lines.append(line)
        
        if not ligand_lines:
            logger.warning(f"No HETATM records found for ligand {ligand_id} in {pdb_file.name}")
            return None
        
        # Write temporary PDB with just this ligand
        temp_pdb = pdb_file.parent / f"temp_{ligand_id}.pdb"
        with open(temp_pdb, 'w') as f:
            f.writelines(ligand_lines)
        
        # Read ligand
        mol = Chem.MolFromPDBFile(str(temp_pdb), removeHs=False)
        
        # Clean up
        temp_pdb.unlink()
        
        if mol is None:
            logger.warning(f"Could not parse ligand {ligand_id} from {pdb_file.name}")
            return None
        
        # Convert to SMILES (canonical)
        smiles = Chem.MolToSmiles(mol)
        
        return smiles
        
    except Exception as e:
        logger.error(f"Error extracting {ligand_id} from {pdb_file.name}: {e}")
        return None


def process_excel_template(excel_file, pdb_dir, output_file):
    """
    Process Excel template:
    1. Add SMILES column
    2. Extract SMILES from PDB structures for entries with PDB codes
    3. Save updated Excel
    """
    logger.info(f"Reading Excel template: {excel_file}")
    df = pd.read_excel(excel_file)
    
    # Clean column names
    df.columns = df.columns.str.strip()
    
    logger.info(f"Found {len(df)} entries")
    logger.info(f"  - {df['PDB_ID'].notna().sum()} with PDB codes")
    logger.info(f"  - {df['PDB_ID'].isna().sum()} without PDB codes (need manual SMILES)")
    
    # Add SMILES column if not present
    if 'SMILES' not in df.columns:
        df.insert(df.columns.get_loc('Target_Ligand') + 1, 'SMILES', None)
        logger.info("Added 'SMILES' column to template")
    
    # Process entries with PDB codes
    pdb_dir = Path(pdb_dir)
    extracted_count = 0
    
    logger.info("\nExtracting SMILES from PDB structures...")
    
    for idx, row in df.iterrows():
        if pd.notna(row['PDB_ID']) and pd.notna(row['Ligand']):
            pdb_id = row['PDB_ID']
            ligand_id = row['Ligand']
            
            # Find PDB file
            possible_names = [
                f"{pdb_id}_prepared.pdb",
                f"{pdb_id}_prep.pdb",
                f"{pdb_id}.pdb",
                f"{pdb_id}_clean.pdb"
            ]
            
            pdb_file = None
            for name in possible_names:
                candidate = pdb_dir / name
                if candidate.exists():
                    pdb_file = candidate
                    break
            
            if pdb_file is None:
                logger.warning(f"  PDB file not found for {pdb_id}")
                continue
            
            # Extract SMILES
            smiles = extract_ligand_from_pdb(pdb_file, ligand_id)
            
            if smiles:
                df.at[idx, 'SMILES'] = smiles
                extracted_count += 1
                logger.info(f"  ✓ {pdb_id}/{ligand_id} ({row['Target_Ligand']}): {smiles[:50]}...")
            else:
                logger.warning(f"  ✗ {pdb_id}/{ligand_id}: Could not extract")
    
    logger.info(f"\nExtracted SMILES for {extracted_count} compounds from PDB structures")
    
    # Identify entries needing manual SMILES
    needs_manual = df[df['PDB_ID'].isna() & df['SMILES'].isna()]
    
    if len(needs_manual) > 0:
        logger.info(f"\n⚠️  {len(needs_manual)} entries need SMILES added manually:")
        for idx, row in needs_manual.iterrows():
            logger.info(f"  - Row {idx}: {row['Target_Ligand']}")
        logger.info("\nPlease add SMILES for these compounds and re-save the Excel file.")
        logger.info("You can look up SMILES in:")
        logger.info("  - PubChem: https://pubchem.ncbi.nlm.nih.gov/")
        logger.info("  - ChEMBL: https://www.ebi.ac.uk/chembl/")
        logger.info("  - DrugBank: https://www.drugbank.ca/")
    
    # Save updated Excel
    output_file = Path(output_file)
    df.to_excel(output_file, index=False)
    logger.info(f"\n✅ Saved updated template to: {output_file}")
    
    # Print summary
    logger.info("\n" + "="*70)
    logger.info("SUMMARY")
    logger.info("="*70)
    logger.info(f"Total entries: {len(df)}")
    logger.info(f"SMILES extracted from PDB: {extracted_count}")
    logger.info(f"SMILES need manual entry: {len(needs_manual)}")
    logger.info(f"Decoys with SMILES: {df['decoy SMILES'].notna().sum()}")
    
    if len(needs_manual) == 0:
        logger.info("\n🎉 Excel template is COMPLETE and ready for validation!")
        logger.info(f"\nNext step:")
        logger.info(f"  python run_ensemble_from_excel_complete.py \\")
        logger.info(f"      --excel {output_file} \\")
        logger.info(f"      --pdb-dir {pdb_dir} \\")
        logger.info(f"      --output results/")
    else:
        logger.info("\n📝 Please complete the Excel file with missing SMILES, then run:")
        logger.info(f"  python run_ensemble_from_excel_complete.py \\")
        logger.info(f"      --excel {output_file} \\")
        logger.info(f"      --pdb-dir {pdb_dir} \\")
        logger.info(f"      --output results/")


def main():
    parser = argparse.ArgumentParser(
        description="Prepare Excel template with SMILES from PDB structures"
    )
    
    parser.add_argument(
        '--excel', required=True,
        help='Input Excel template'
    )
    parser.add_argument(
        '--pdb-dir', required=True,
        help='Directory with PDB structures'
    )
    parser.add_argument(
        '--output', default='anathi_template_complete.xlsx',
        help='Output Excel file with SMILES'
    )
    
    args = parser.parse_args()
    
    logger.info("="*70)
    logger.info("EXCEL TEMPLATE PREPARATION")
    logger.info("="*70)
    
    try:
        process_excel_template(args.excel, args.pdb_dir, args.output)
        return 0
    except Exception as e:
        logger.exception(f"Failed: {e}")
        return 1


if __name__ == "__main__":
    import sys
    sys.exit(main())
