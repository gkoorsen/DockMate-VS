#!/bin/bash
#
# Personalized Quick Start Script for Ensemble Validation
# Your specific paths configured
#

set -e  # Exit on error

echo "========================================================================"
echo "ENSEMBLE VALIDATION - QUICK START"
echo "========================================================================"
echo ""

# Your specific paths
PDB_DIR="/Users/gerritkoorsen/docking_platform_gui/output/redock_analysis/pdbs"
PROJECT_ROOT="/Users/gerritkoorsen/docking_platform_gui"
EXCEL_FILE="anathi_template.xlsx"  # Adjust if different location

echo "Configuration:"
echo "  PDB directory: $PDB_DIR"
echo "  Project root: $PROJECT_ROOT"
echo ""

# Check PDB directory exists
if [ ! -d "$PDB_DIR" ]; then
    echo "ERROR: PDB directory not found: $PDB_DIR"
    exit 1
fi

echo "✓ PDB directory found"
echo ""

# ============================================================================
# STEP 1: Convert PDB → PDBQT (REQUIRED!)
# ============================================================================

echo "========================================================================"
echo "STEP 1: Convert PDB → PDBQT"
echo "========================================================================"
echo ""
echo "Your docking platform requires PDBQT format, not PDB."
echo "Converting your 5 ensemble structures..."
echo ""

cd "$PDB_DIR"

# Check if prepare_receptor4.py is available
if ! command -v prepare_receptor4.py &> /dev/null; then
    echo "WARNING: prepare_receptor4.py not found in PATH"
    echo ""
    echo "You need MGLTools to convert PDB → PDBQT."
    echo "Two options:"
    echo ""
    echo "Option 1: Install MGLTools"
    echo "  Download from: http://mgltools.scripps.edu/"
    echo "  Then add to PATH or use full path to prepare_receptor4.py"
    echo ""
    echo "Option 2: Use Python/OpenBabel"
    echo "  obabel -ipdb 7EO7.pdb -opdbqt -O 7EO7_prepared.pdbqt"
    echo ""
    echo "After installing, re-run this script."
    exit 1
fi

# Convert each structure
for pdb_id in 7EO7 7XAX 7UUP 7RN1 7DPP; do
    if [ -f "${pdb_id}.pdb" ]; then
        echo "Converting ${pdb_id}.pdb → ${pdb_id}_prepared.pdbqt..."
        prepare_receptor4.py -r "${pdb_id}.pdb" -o "${pdb_id}_prepared.pdbqt" -U nphs_lps
    else
        echo "WARNING: ${pdb_id}.pdb not found in $PDB_DIR"
    fi
done

echo ""
echo "✓ PDB → PDBQT conversion complete"
echo ""

# ============================================================================
# STEP 2: Extract SMILES from PDB structures
# ============================================================================

echo "========================================================================"
echo "STEP 2: Extract SMILES from PDB structures"
echo "========================================================================"
echo ""

cd "$PROJECT_ROOT"

python scripts/prepare_excel_with_smiles.py \
    --excel "$EXCEL_FILE" \
    --pdb-dir "$PDB_DIR" \
    --output "anathi_template_complete.xlsx"

echo ""
echo "✓ SMILES extraction complete"
echo ""

# ============================================================================
# STEP 3: Manual intervention needed
# ============================================================================

echo "========================================================================"
echo "STEP 3: MANUAL ACTION REQUIRED"
echo "========================================================================"
echo ""
echo "📝 Please complete the Excel file with missing SMILES:"
echo ""
echo "   1. Open: anathi_template_complete.xlsx"
echo "   2. Add SMILES for compounds without PDB codes (9 compounds):"
echo "      - Tideglusib"
echo "      - Carmofur"
echo "      - Epigallocatechin gallate"
echo "      - GC373"
echo "      - PX-12"
echo "      - Disulfiram"
echo "      - Captan"
echo "      - GRL-1720"
echo "      - Lufotrelvir"
echo ""
echo "   3. Look up SMILES on PubChem: https://pubchem.ncbi.nlm.nih.gov/"
echo "   4. Save the Excel file"
echo ""
echo "After completing Excel, run:"
echo "   ./run_ensemble_validation.sh"
echo ""
echo "========================================================================"
