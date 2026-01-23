#!/bin/bash
#
# Run Ensemble Validation Experiment - WITH PROPER PYTHON PATH
#

set -e  # Exit on error

echo "========================================================================"
echo "ENSEMBLE VALIDATION EXPERIMENT - PRODUCTION RUN"
echo "========================================================================"
echo ""

# Your specific paths
PDB_DIR="/Users/gerritkoorsen/docking_platform_gui/output/redock_analysis/pdbs"
PROJECT_ROOT="/Users/gerritkoorsen/docking_platform_gui"
EXCEL_FILE="/data/anathi_template_complete.xlsx"
OUTPUT_DIR="results/mpro_ensemble_validation"

# SET PYTHON PATH - THIS IS CRITICAL!
export PYTHONPATH="$PROJECT_ROOT:$PYTHONPATH"

echo "Configuration:"
echo "  Excel file: $EXCEL_FILE"
echo "  PDB directory: $PDB_DIR"
echo "  Output directory: $OUTPUT_DIR"
echo "  Python path: $PYTHONPATH"
echo ""

cd "$PROJECT_ROOT"

# Check Excel file exists
if [ ! -f "$EXCEL_FILE" ]; then
    echo "ERROR: Excel file not found: $EXCEL_FILE"
    echo ""
    echo "Please complete the Excel file first:"
    echo "  ./setup_ensemble_validation_openbabel.sh"
    exit 1
fi

# Check file structure
echo "Checking file structure..."
required_files=(
    "docking_platform_gui/analysis/enrichment_analysis.py"
    "docking_platform_gui/experiments/ensemble_validation.py"
    "docking_platform_gui/experiments/docking_adapter_final.py"
    "scripts/prepare_excel_with_smiles.py"
    "scripts/run_ensemble_excel_production.py"
)

missing=0
for file in "${required_files[@]}"; do
    if [ -f "$file" ]; then
        echo "  ✓ $file"
    else
        echo "  ✗ $file MISSING"
        missing=$((missing + 1))
    fi
done

if [ $missing -gt 0 ]; then
    echo ""
    echo "ERROR: $missing required file(s) missing"
    echo ""
    echo "Please run the file structure setup first:"
    echo "  ./setup_file_structure.sh"
    echo ""
    echo "This will copy the Python modules to the correct locations."
    exit 1
fi

echo ""
echo "✓ All files present"
echo ""

# Check PDBQT files exist
echo "Checking for required PDBQT files..."
missing=0
for pdb_id in 7EO7 7XAX 7UUP 7RN1 7DPP; do
    pdbqt_file="$PDB_DIR/${pdb_id}_prepared.pdbqt"
    if [ -f "$pdbqt_file" ]; then
        echo "  ✓ $pdbqt_file"
    else
        echo "  ✗ $pdbqt_file NOT FOUND"
        missing=$((missing + 1))
    fi
done

if [ $missing -gt 0 ]; then
    echo ""
    echo "ERROR: $missing PDBQT file(s) missing"
    echo "Please run setup first:"
    echo "  ./setup_ensemble_validation_openbabel.sh"
    exit 1
fi

echo ""
echo "✓ All PDBQT files present"
echo ""

# Ask for confirmation
echo "========================================================================"
echo "EXPERIMENT SETTINGS"
echo "========================================================================"
echo ""
echo "This experiment will:"
echo "  - Dock 14 active compounds + 14 decoys"
echo "  - Test single structure (7UUP) vs ensemble (all 5)"
echo "  - Calculate AUC, EF1%, BEDROC metrics"
echo "  - Perform 1000 bootstrap iterations"
echo "  - Generate plots and recommendation"
echo ""
echo "⏰ Estimated time: 6-8 hours"
echo "💾 Output location: $OUTPUT_DIR"
echo ""

read -p "Proceed with experiment? (yes/no): " confirm

if [ "$confirm" != "yes" ]; then
    echo "Experiment cancelled"
    exit 0
fi

# Run experiment
echo ""
echo "========================================================================"
echo "STARTING EXPERIMENT"
echo "========================================================================"
echo ""
echo "Starting at: $(date)"
echo ""

python scripts/run_ensemble_excel_production.py \
    --excel "$EXCEL_FILE" \
    --pdb-dir "$PDB_DIR" \
    --output "$OUTPUT_DIR" \
    --exhaustiveness 16 \
    --cpu 4 \
    --use-cache \
    --bootstrap 1000

exit_code=$?

# Results summary
echo ""
echo "========================================================================"
if [ $exit_code -eq 0 ]; then
    echo "✅ EXPERIMENT COMPLETED SUCCESSFULLY"
    echo "========================================================================"
    echo ""
    echo "Finished at: $(date)"
    echo ""
    echo "Results saved to: $OUTPUT_DIR"
    echo ""
    echo "Key files:"
    echo "  📄 $OUTPUT_DIR/SUMMARY_REPORT.txt  ← READ THIS FIRST"
    echo "  📊 $OUTPUT_DIR/ensemble_validation_plots.png"
    echo "  📈 $OUTPUT_DIR/metrics.json"
    echo ""
    echo "To view the recommendation:"
    echo "  cat $OUTPUT_DIR/SUMMARY_REPORT.txt"
    echo ""
else
    echo "❌ EXPERIMENT FAILED"
    echo "========================================================================"
    echo ""
    echo "Check the log file for errors:"
    echo "  $OUTPUT_DIR/experiment.log"
    echo ""
fi

exit $exit_code
