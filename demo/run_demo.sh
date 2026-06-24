#!/bin/bash
# demo/run_demo.sh - Run PhyloSOLID demo with scRNA-seq data

set -e

echo "========================================"
echo "PhyloSOLID Demo: scRNA-seq Phylogenetic Tracing"
echo "========================================"

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# ============================================
# Demo configuration - using demo/scrna/ data
# ============================================
SAMPLE="demo_scrna"
WORKDIR="$SCRIPT_DIR/scrna/output"
INPUT_DIR="$SCRIPT_DIR/scrna/input"

# Input files
MUTATION_LIST="$INPUT_DIR/identifier.txt"
BAM_FILE="$INPUT_DIR/demo_scrna.bam"
BARCODE_FILE="$INPUT_DIR/barcodes.txt"
CELLTYPE_FILE="$INPUT_DIR/celltype.txt"

# Parameters
THREADS=4
READ_LEN=120
CELL_NUM=3319

# ============================================
# Check input files
# ============================================
echo ""
echo "Checking input files..."

MISSING=0
for file in "$MUTATION_LIST" "$BAM_FILE" "$BARCODE_FILE" "$CELLTYPE_FILE"; do
    if [ ! -f "$file" ]; then
        echo "  ❌ Missing: $file"
        MISSING=1
    else
        echo "  ✓ $(basename "$file")"
    fi
done

# Check BAM index
if [ ! -f "$BAM_FILE.bai" ]; then
    echo "  ⚠️  BAM index not found (will be auto-generated)"
fi

if [ $MISSING -eq 1 ]; then
    echo ""
    echo "Error: Missing input files!"
    echo "Please ensure demo/scrna/input/ contains:"
    echo "  - identifier.txt"
    echo "  - demo_scrna.bam"
    echo "  - barcodes.txt"
    echo "  - celltype.txt"
    exit 1
fi

echo "✓ All input files found"

# ============================================
# Create output directory
# ============================================
mkdir -p "$WORKDIR"

echo ""
echo "========================================"
echo "Demo Configuration:"
echo "========================================"
echo "  Sample:          $SAMPLE"
echo "  Workdir:         $WORKDIR"
echo "  Input dir:       $INPUT_DIR"
echo "  Mutations:       $(basename "$MUTATION_LIST")"
echo "  BAM:             $(basename "$BAM_FILE")"
echo "  Barcode:         $(basename "$BARCODE_FILE")"
echo "  Cell type:       $(basename "$CELLTYPE_FILE")"
echo "  Threads:         $THREADS"
echo "  Read length:     $READ_LEN"
echo "  Cell number:     $CELL_NUM"
echo "========================================"
echo ""

# ============================================
# Run PhyloSOLID
# ============================================
echo "Running PhyloSOLID demo..."
echo "----------------------------------------"

cd "$PROJECT_DIR"

# 使用 python -m cli.main (如果 phylosolid 命令不可用)
if command -v phylosolid &> /dev/null; then
    echo "Using 'phylosolid' command"
    phylosolid --workdir "$WORKDIR" scrna \
        --sample "$SAMPLE" \
        --mutation-list "$MUTATION_LIST" \
        --bam "$BAM_FILE" \
        --barcode "$BARCODE_FILE" \
        --read-len "$READ_LEN" \
        --cellnum "$CELL_NUM" \
        --celltype-file "$CELLTYPE_FILE" \
        --threads "$THREADS" \
        --verbose
else
    echo "Using 'python -m cli.main'"
    python -m cli.main --workdir "$WORKDIR" scrna \
        --sample "$SAMPLE" \
        --mutation-list "$MUTATION_LIST" \
        --bam "$BAM_FILE" \
        --barcode "$BARCODE_FILE" \
        --read-len "$READ_LEN" \
        --cellnum "$CELL_NUM" \
        --celltype-file "$CELLTYPE_FILE" \
        --threads "$THREADS" \
        --verbose
fi

# ============================================
# Check results
# ============================================
if [ $? -eq 0 ]; then
    echo ""
    echo "----------------------------------------"
    echo "✅ Demo completed successfully!"
    echo ""
    echo "Results saved in: $WORKDIR/$SAMPLE"
    echo ""
    echo "Output structure:"
    echo "  ├── 01_features/"
    echo "  ├── 02_treeinput/"
    echo "  ├── 03_tree_building/"
    echo "  │   ├── mutation_integrator/phylo/"
    echo "  │   │   ├── final_cleaned_tree_node.json"
    echo "  │   │   ├── final_cleaned_M_*.CFMatrix"
    echo "  │   │   └── final_cleaned_I_*_for_circosPlot.txt"
    echo "  │   └── scaffold_builder/"
    echo "  │       ├── df_barcode_clones_from_phylo_tree.csv"
    echo "  │       └── demo_scrna.mutation_hierarchy.csv"
    echo "  └── pipeline_summary.yaml"
    echo ""
    echo "Key output files:"
    find "$WORKDIR/$SAMPLE" -type f \( -name "*.json" -o -name "*.CFMatrix" -o -name "*.csv" -o -name "*.yaml" \) 2>/dev/null | head -10 | sed 's/^/  /'
else
    echo ""
    echo "❌ Demo failed"
    exit 1
fi

echo ""
echo "========================================"
echo "To view results interactively, open:"
echo "  file://$WORKDIR/$SAMPLE/03_tree_building/mutation_integrator/phylo/final_cleaned_tree_node.json"
echo "========================================"