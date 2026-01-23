#!/bin/bash
#
# Setup file structure for ensemble validation
# Creates necessary directories and copies files to correct locations
#

set -e

PROJECT_ROOT="/Users/gerritkoorsen/docking_platform_gui"

echo "========================================================================"
echo "SETTING UP FILE STRUCTURE"
echo "========================================================================"
echo ""

cd "$PROJECT_ROOT"

# Create necessary directories
echo "Creating directory structure..."

mkdir -p docking_platform_gui/analysis
mkdir -p docking_platform_gui/experiments
mkdir -p scripts

echo "✓ Directories created"
echo ""

# Check which files need to be copied
echo "Checking for required files in current directory..."
echo ""

files_to_copy=(
    "enrichment_analysis.py:docking_platform_gui/analysis/"
    "ensemble_validation.py:docking_platform_gui/experiments/"
    "docking_adapter_final.py:docking_platform_gui/experiments/"
    "prepare_excel_with_smiles.py:scripts/"
    "run_ensemble_excel_production.py:scripts/"
)

missing_files=()

for file_mapping in "${files_to_copy[@]}"; do
    file="${file_mapping%%:*}"
    dest="${file_mapping##*:}"
    
    if [ -f "$file" ]; then
        echo "✓ Found: $file"
        echo "  → Copying to: $dest"
        cp "$file" "$dest"
    else
        echo "✗ Missing: $file"
        missing_files+=("$file")
    fi
done

echo ""

if [ ${#missing_files[@]} -gt 0 ]; then
    echo "⚠️  WARNING: ${#missing_files[@]} file(s) missing"
    echo ""
    echo "Please ensure these files are in:"
    echo "  $PROJECT_ROOT/"
    echo ""
    echo "Missing files:"
    for file in "${missing_files[@]}"; do
        echo "  - $file"
    done
    echo ""
    echo "Download them from the outputs folder and place in:"
    echo "  $PROJECT_ROOT/"
    echo ""
    exit 1
fi

echo "✓ All files copied successfully"
echo ""

# Create __init__.py files to make packages importable
echo "Creating __init__.py files..."

touch docking_platform_gui/__init__.py
touch docking_platform_gui/analysis/__init__.py
touch docking_platform_gui/experiments/__init__.py

echo "✓ Package structure created"
echo ""

# Verify structure
echo "Verifying file structure..."
echo ""

required_files=(
    "docking_platform_gui/analysis/enrichment_analysis.py"
    "docking_platform_gui/experiments/ensemble_validation.py"
    "docking_platform_gui/experiments/docking_adapter_final.py"
    "scripts/prepare_excel_with_smiles.py"
    "scripts/run_ensemble_excel_production.py"
)

all_present=true

for file in "${required_files[@]}"; do
    if [ -f "$file" ]; then
        echo "✓ $file"
    else
        echo "✗ $file MISSING"
        all_present=false
    fi
done

echo ""

if [ "$all_present" = true ]; then
    echo "========================================================================"
    echo "✅ FILE STRUCTURE SETUP COMPLETE"
    echo "========================================================================"
    echo ""
    echo "Directory structure:"
    echo ""
    echo "  docking_platform_gui/"
    echo "  ├── analysis/"
    echo "  │   └── enrichment_analysis.py"
    echo "  ├── experiments/"
    echo "  │   ├── ensemble_validation.py"
    echo "  │   └── docking_adapter_final.py"
    echo "  └── ..."
    echo ""
    echo "  scripts/"
    echo "  ├── prepare_excel_with_smiles.py"
    echo "  └── run_ensemble_excel_production.py"
    echo ""
    echo "Next step: Run the validation experiment"
    echo "  ./run_ensemble_validation.sh"
    echo ""
else
    echo "========================================================================"
    echo "⚠️  SETUP INCOMPLETE"
    echo "========================================================================"
    echo ""
    echo "Some files are still missing. Please check above."
    exit 1
fi
