#!/bin/bash
# demo/run_demo.sh - Run PhyloSOLID demo

set -e

echo "========================================"
echo "PhyloSOLID Demo Run"
echo "========================================"

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# Set paths
SAMPLE="Org10S4D46"
MUTATION_LIST="$SCRIPT_DIR/input/02_identifier/identifier.txt"
BAM_FILE="$SCRIPT_DIR/input/01_rawdata/${SAMPLE}.bam"
BARCODE_FILE="$SCRIPT_DIR/input/01_rawdata/${SAMPLE}_CB.txt"
WORKDIR="$SCRIPT_DIR/output"

# Check if files exist
echo "Checking input files..."
for file in "$MUTATION_LIST" "$BAM_FILE" "$BARCODE_FILE"; do
    if [ ! -f "$file" ]; then
        echo "Error: $file not found!"
        exit 1
    fi
done
echo "✓ All input files found"

# Create output directory
mkdir -p "$WORKDIR"

echo ""
echo "Input files:"
echo "  Mutations: $MUTATION_LIST"
echo "  BAM:       $BAM_FILE"
echo "  Barcode:   $BARCODE_FILE"
echo "Output dir:  $WORKDIR"
echo ""

# Run PhyloSOLID
echo "Running PhyloSOLID demo..."
echo "----------------------------------------"

# 使用 python -m cli.main 而不是直接 phylosolid 命令
python -m cli.main  --workdir "$WORKDIR"  scrna \
    --sample "$SAMPLE" \
    --mutation-list "$MUTATION_LIST" \
    --bam "$BAM_FILE" \
    --barcode "$BARCODE_FILE" \
    --threads 4 \
    --read-len 91 \
    --cellnum 836 \
    --running-type benchmark \
    --verbose

echo "----------------------------------------"

# Check results
if [ $? -eq 0 ]; then
    echo ""
    echo "✅ Demo completed successfully!"
    echo "Results saved in: $WORKDIR/$SAMPLE"
    
    # Show output structure
    echo ""
    echo "Output files:"
    find "$WORKDIR/$SAMPLE" -type f -name "*.txt" -o -name "*.newick" | head -10
else
    echo ""
    echo "❌ Demo failed"
    exit 1
fi

echo ""
echo "========================================"
