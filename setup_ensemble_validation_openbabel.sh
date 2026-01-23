#!/bin/bash
#
# Personalized Quick Start Script for Ensemble Validation
# Using OpenBabel for PDB → PDBQT conversion
#

set -e  # Exit on error

echo "========================================================================"
echo "ENSEMBLE VALIDATION - QUICK START (OpenBabel Version)"
echo "========================================================================"
echo ""

# Your specific paths
PDB_DIR="/Users/gerritkoorsen/docking_platform_gui/docking_platform_gui/preparation/output/redock_analysis/pdbs"
PROJECT_ROOT="/Users/gerritkoorsen/docking_platform_gui"
EXCEL_FILE="/Users/gerritkoorsen/docking_platform_gui/templates/anathi_template.xlsx"  # Adjust if different location

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
# STEP 1: Convert PDB → PDBQT using OpenBabel
# ============================================================================

echo "========================================================================"
echo "STEP 1: Convert PDB → PDBQT using OpenBabel"
echo "========================================================================"
echo ""

# Check if obabel is available
if ! command -v obabel &> /dev/null; then
    echo "ERROR: obabel (OpenBabel) not found in PATH"
    echo ""
    echo "Please install OpenBabel:"
    echo "  conda install -c conda-forge openbabel"
    echo ""
    echo "Or:"
    echo "  brew install open-babel  # macOS with Homebrew"
    echo "  pip install openbabel-wheel  # Python package"
    echo ""
    exit 1
fi

echo "✓ OpenBabel found: $(which obabel)"
echo ""

cd "$PDB_DIR"

echo "Converting your 5 ensemble structures..."
echo ""

# Convert each structure
# -xr flag: rigid molecule (for receptor)
# -xh flag: add hydrogens (pH 7.4)
for pdb_id in 7EO7 7XAX 7UUP 7RN1 7DPP; do
    pdb_file="${pdb_id}.pdb"
    pdbqt_file="${pdb_id}_prepared.pdbqt"
    
    if [ -f "$pdb_file" ]; then
        echo "Converting ${pdb_file} → ${pdbqt_file}..."
        
        # OpenBabel conversion with receptor settings
        obabel -ipdb "$pdb_file" \
               -opdbqt \
               -O "$pdbqt_file" \
               -xr \
               -xh \
               -p 7.4
        
        if [ -f "$pdbqt_file" ]; then
            echo "  ✓ Created: ${pdbqt_file}"
        else
            echo "  ✗ Failed to create: ${pdbqt_file}"
        fi
    else
        echo "⚠️  WARNING: ${pdb_file} not found in $PDB_DIR"
        echo "   Looking for: ${PDB_DIR}/${pdb_file}"
    fi
    echo ""
done

# Check if all conversions succeeded
echo "Checking conversions..."
all_found=true
for pdb_id in 7EO7 7XAX 7UUP 7RN1 7DPP; do
    if [ ! -f "${pdb_id}_prepared.pdbqt" ]; then
        echo "  ✗ Missing: ${pdb_id}_prepared.pdbqt"
        all_found=false
    else
        # Show file size to confirm it's not empty
        size=$(ls -lh "${pdb_id}_prepared.pdbqt" | awk '{print $5}')
        echo "  ✓ ${pdb_id}_prepared.pdbqt (${size})"
    fi
done

if [ "$all_found" = false ]; then
    echo ""
    echo "⚠️  Some PDBQT files are missing!"
    echo "Please check that all PDB files exist in: $PDB_DIR"
    exit 1
fi

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

# Check if Python script exists
if [ ! -f "scripts/prepare_excel_with_smiles.py" ]; then
    echo "⚠️  WARNING: scripts/prepare_excel_with_smiles.py not found"
    echo ""
    echo "Please ensure you've copied all scripts to:"
    echo "  $PROJECT_ROOT/scripts/"
    echo ""
    echo "Required scripts:"
    echo "  - prepare_excel_with_smiles.py"
    echo "  - run_ensemble_excel_production.py"
    echo ""
    exit 1
fi

# Check if Excel file exists
if [ ! -f "$EXCEL_FILE" ]; then
    echo "⚠️  WARNING: Excel file not found: $EXCEL_FILE"
    echo ""
    echo "Please ensure anathi_template.xlsx is in:"
    echo "  $PROJECT_ROOT/"
    echo ""
    exit 1
fi

echo "Extracting SMILES from PDB structures..."
echo ""

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
echo ""
echo "   2. Add SMILES for these 9 compounds (look up on PubChem):"
echo ""
echo "      Compound                       PubChem CID"
echo "      ─────────────────────────────────────────────"
echo "      Tideglusib                     135565641"
echo "      Carmofur                       2554"
echo "      Epigallocatechin gallate       65064"
echo "      GC373                          135565641"
echo "      PX-12                          71515"
echo "      Disulfiram                     3117"
echo "      Captan                         2346"
echo "      GRL-1720                       (search PubChem)"
echo "      Lufotrelvir                    146170876"
echo ""
echo "   3. For each compound:"
echo "      a) Go to: https://pubchem.ncbi.nlm.nih.gov/"
echo "      b) Search by name or CID"
echo "      c) Find 'Canonical SMILES' on the compound page"
echo "      d) Copy and paste into the SMILES column"
echo ""
echo "   4. Save the Excel file"
echo ""
echo "After completing, run:"
echo "   ./run_ensemble_validation.sh"
echo ""
echo "========================================================================"
echo ""
echo "Setup complete! ✅"
echo ""
