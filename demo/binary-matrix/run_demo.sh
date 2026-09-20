#!/bin/bash
# demo/binary-matrix/run_demo.sh - Run PhyloSOLID binary-matrix demo

set -e

echo "========================================"
echo "PhyloSOLID Demo: binary-matrix tree building"
echo "========================================"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"

SAMPLE="demo_binput"
WORKDIR="$SCRIPT_DIR/output"
INPUT_FILE="$SCRIPT_DIR/input/demo_input.tsv"

echo ""
echo "Checking input files..."
if [ ! -f "$INPUT_FILE" ]; then
    echo "  Missing: $INPUT_FILE"
    echo "Error: Missing binary-matrix demo input"
    exit 1
fi
echo "  OK $(basename "$INPUT_FILE")"

mkdir -p "$WORKDIR"

echo ""
echo "========================================"
echo "Demo configuration"
echo "========================================"
echo "  Sample:       $SAMPLE"
echo "  Workdir:      $WORKDIR"
echo "  Input:        $INPUT_FILE"
echo "========================================"
echo ""

cd "$PROJECT_DIR"
export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/mplconfig}"
export PYTHONPATH="$PROJECT_DIR${PYTHONPATH:+:$PYTHONPATH}"

if command -v phylosolid >/dev/null 2>&1; then
    echo "Using phylosolid"
    phylosolid binary-matrix \
        --sampleid "$SAMPLE" \
        --inputfile "$INPUT_FILE" \
        --outputpath "$WORKDIR"
else
    echo "Using python -m cli.main"
    python -m cli.main binary-matrix \
        --sampleid "$SAMPLE" \
        --inputfile "$INPUT_FILE" \
        --outputpath "$WORKDIR"
fi

echo ""
echo "----------------------------------------"
echo "Drawing circos (separate visualization env)"
"$PROJECT_DIR/scripts/visualization/run_circos.sh" \
    "$WORKDIR/$SAMPLE/03_final_results/phylo" \
    "$WORKDIR/$SAMPLE/03_final_results/circos" \
    "$WORKDIR/$SAMPLE/01_scaffold_builder/df_celltype.txt"

echo ""
echo "----------------------------------------"
echo "Demo finished"
echo "Results: $WORKDIR/$SAMPLE/03_final_results/phylo/"
echo "Circos:  $WORKDIR/$SAMPLE/03_final_results/circos/"
echo "========================================"
