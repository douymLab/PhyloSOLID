#!/bin/bash
# demo/scrna/run_demo.sh - Run PhyloSOLID scRNA demo

set -e

echo "========================================"
echo "PhyloSOLID Demo: scRNA phylogenetic tracing"
echo "========================================"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"

SAMPLE="demo_scrna"
WORKDIR="$SCRIPT_DIR/output"
INPUT_DIR="$SCRIPT_DIR/input"
MUTATION_LIST="$INPUT_DIR/identifier.txt"
BAM_FILE="$INPUT_DIR/demo_scrna.bam"
BARCODE_FILE="$INPUT_DIR/barcodes.txt"
CELLTYPE_FILE="$INPUT_DIR/celltype.txt"
THREADS="${THREADS:-4}"
READ_LEN="${READ_LEN:-120}"
CELL_NUM="${CELL_NUM:-3319}"

echo ""
echo "Checking input files..."
MISSING=0
for file in "$MUTATION_LIST" "$BAM_FILE" "$BARCODE_FILE" "$CELLTYPE_FILE"; do
    if [ ! -f "$file" ]; then
        echo "  Missing: $file"
        MISSING=1
    else
        echo "  OK $(basename "$file")"
    fi
done

if [ ! -f "$BAM_FILE.bai" ] && [ ! -f "${BAM_FILE%.bam}.bai" ]; then
    echo "  BAM index not found (will be auto-generated if the pipeline supports it)"
fi

if [ "$MISSING" -eq 1 ]; then
    echo ""
    echo "Error: Missing scRNA demo inputs under demo/scrna/input/"
    echo "Need identifier.txt, demo_scrna.bam, barcodes.txt, and celltype.txt"
    exit 1
fi

mkdir -p "$WORKDIR"

echo ""
echo "========================================"
echo "Demo configuration"
echo "========================================"
echo "  Sample:       $SAMPLE"
echo "  Workdir:      $WORKDIR"
echo "  Mutations:    $(wc -l < "$MUTATION_LIST") sites"
echo "  BAM:          $(basename "$BAM_FILE")"
echo "  Barcodes:     $(wc -l < "$BARCODE_FILE")"
echo "  Threads:      $THREADS"
echo "  Read length:  $READ_LEN"
echo "========================================"
echo ""

cd "$PROJECT_DIR"
export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/mplconfig}"

RUN_CMD=(--workdir "$WORKDIR" scrna
    --sample "$SAMPLE"
    --mutation-list "$MUTATION_LIST"
    --bam "$BAM_FILE"
    --barcode "$BARCODE_FILE"
    --celltype-file "$CELLTYPE_FILE"
    --read-len "$READ_LEN"
    --cellnum "$CELL_NUM"
    --threads "$THREADS"
)

if command -v phylosolid >/dev/null 2>&1; then
    echo "Using phylosolid"
    phylosolid --verbose "${RUN_CMD[@]}"
else
    echo "Using python -m cli.main"
    python -m cli.main --verbose "${RUN_CMD[@]}"
fi

echo ""
echo "----------------------------------------"
echo "Drawing circos (separate visualization env)"
"$PROJECT_DIR/scripts/visualization/run_circos.sh" \
    "$WORKDIR/$SAMPLE/03_tree_building/05_final_results/phylo" \
    "$WORKDIR/$SAMPLE/03_tree_building/05_final_results/circos" \
    "$WORKDIR/$SAMPLE/03_tree_building/03_scaffold_builder/df_celltype.txt"

echo ""
echo "----------------------------------------"
echo "Demo finished"
echo "Results: $WORKDIR/$SAMPLE"
echo "  - 01_features/"
echo "  - 02_treeinput/"
echo "  - 03_tree_building/05_final_results/phylo/"
echo "  - 03_tree_building/05_final_results/circos/"
echo "========================================"
