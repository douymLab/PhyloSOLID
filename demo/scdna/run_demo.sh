#!/bin/bash
# demo/scdna/run_demo.sh - Run PhyloSOLID scDNA demo from target-site mini BAMs
# Bulk BAM is not required.

set -e

echo "========================================"
echo "PhyloSOLID Demo: scDNA phylogenetic tracing"
echo "========================================"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"

SAMPLE="demo_scdna"
WORKDIR="$SCRIPT_DIR/output"
INPUT_DIR="$SCRIPT_DIR/input"
BAM_DIR="$INPUT_DIR/bam"
MUTATION_LIST="$INPUT_DIR/identifier.txt"
SAMPLE_LIST="$INPUT_DIR/sample_list.txt"
THREADS="${THREADS:-4}"

echo ""
echo "Checking input files..."
MISSING=0
for file in "$MUTATION_LIST" "$SAMPLE_LIST"; do
    if [ ! -f "$file" ]; then
        echo "  Missing: $file"
        MISSING=1
    else
        echo "  OK $(basename "$file")"
    fi
done

if [ ! -d "$BAM_DIR" ]; then
    echo "  Missing BAM directory: $BAM_DIR"
    MISSING=1
else
    NBAM=$(ls "$BAM_DIR"/*.bam 2>/dev/null | wc -l)
    echo "  OK bam/ ($NBAM mini BAM files)"
    if [ "$NBAM" -eq 0 ]; then
        MISSING=1
    fi
fi

if [ "$MISSING" -eq 1 ]; then
    echo ""
    echo "Error: Missing scDNA demo inputs under demo/scdna/input/"
    exit 1
fi

echo "Indexing mini BAMs if needed..."
python - <<PY
from pathlib import Path
import pysam
bam_dir = Path("$BAM_DIR")
for bam in sorted(bam_dir.glob("*.bam")):
    bai = Path(str(bam) + ".bai")
    alt = bam.with_suffix(".bai")
    if bai.exists() or alt.exists():
        continue
    pysam.index(str(bam))
    print("  indexed", bam.name)
PY

mkdir -p "$WORKDIR"

echo ""
echo "========================================"
echo "Demo configuration"
echo "========================================"
echo "  Sample:       $SAMPLE"
echo "  Workdir:      $WORKDIR"
echo "  Mutations:    $(wc -l < "$MUTATION_LIST") sites"
echo "  BAM dir:      $BAM_DIR"
echo "  Cells:        $(wc -l < "$SAMPLE_LIST")"
echo "  Bulk BAM:     not used"
echo "  Threads:      $THREADS"
echo "  Steps:        feature_extraction tree_input tree_building"
echo "  Genotyper:    config/paths.yaml scdna.genotyper_dir or PHYLOSOLID_GENOTYPER_DIR"
echo "========================================"
echo ""

cd "$PROJECT_DIR"
export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/mplconfig}"

RUN_CMD=(--workdir "$WORKDIR" scdna
    --sample "$SAMPLE"
    --mutation-list "$MUTATION_LIST"
    --bam-dir "$BAM_DIR"
    --sample-list "$SAMPLE_LIST"
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
echo "Tree-input matrices: $WORKDIR/$SAMPLE/02_treeinput/data/"
echo "Tree-building results: $WORKDIR/$SAMPLE/03_tree_building/"
echo "  - 05_final_results/phylo/"
echo "  - 05_final_results/circos/"
echo "  - PhyloSOLID/celltree.newick (if Newick conversion succeeded)"
echo "========================================"
