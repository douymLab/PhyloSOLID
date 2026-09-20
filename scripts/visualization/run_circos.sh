#!/bin/bash
# Run PhyloSOLIDvis circos using the visualization conda env, never the tree-building env.
# Usage: run_circos.sh <phylo_dir> <circos_dir> [annotation_file]
# Skips (exit 0) if no visualization R is found.

set -e

PHYLO_DIR="${1:-}"
CIRCOS_DIR="${2:-}"
ANNOT_FILE="${3:-}"

if [ -z "$PHYLO_DIR" ] || [ -z "$CIRCOS_DIR" ]; then
    echo "Usage: $0 <phylo_dir> <circos_dir> [annotation_file]"
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
R_FILE="$SCRIPT_DIR/run_phylosolidvis.R"

find_vis_rscript() {
    if [ -n "${PHYLOSOLID_VIS_RSCRIPT:-}" ] && [ -x "$PHYLOSOLID_VIS_RSCRIPT" ]; then
        echo "$PHYLOSOLID_VIS_RSCRIPT"
        return 0
    fi

    local bases=()
    local conda_base=""
    if command -v conda >/dev/null 2>&1; then
        conda_base="$(conda info --base 2>/dev/null || true)"
        [ -n "$conda_base" ] && bases+=("$conda_base")
    fi
    if [ -n "${CONDA_PREFIX:-}" ]; then
        local walked
        walked="$(cd "$CONDA_PREFIX/.." && pwd)"
        [ "$(basename "$walked")" = "envs" ] && walked="$(cd "$walked/.." && pwd)"
        bases+=("$walked")
    fi
    bases+=("$HOME/Conda" "$HOME/anaconda3" "$HOME/miniconda3")

    local names=()
    [ -n "${PHYLOSOLID_VIS_ENV:-}" ] && names+=("$PHYLOSOLID_VIS_ENV")
    names+=(phylosolid_vis circos_env)

    local base name cand seen="|"
    for base in "${bases[@]}"; do
        [ -n "$base" ] || continue
        case "$seen" in
            *"|$base|"*) continue ;;
        esac
        seen="${seen}${base}|"
        for name in "${names[@]}"; do
            cand="${base}/envs/${name}/bin/Rscript"
            if [ -x "$cand" ]; then
                echo "$cand"
                return 0
            fi
        done
    done
    return 1
}

if [ ! -d "$PHYLO_DIR" ]; then
    echo "Circos skipped: phylo directory not found: $PHYLO_DIR"
    exit 0
fi

if ! RSCRIPT="$(find_vis_rscript)"; then
    echo "Circos skipped: visualization R not found."
    echo "  Install a separate env with: bash install_vis.sh"
    echo "  Or point to an existing env: export PHYLOSOLID_VIS_ENV=circos_env"
    echo "  Do not install PhyloSOLIDvis into pmg / phylosolid_env."
    exit 0
fi

# Refuse the tree-building envs even if someone exported them by mistake
case "$RSCRIPT" in
    */envs/pmg/bin/Rscript|*/envs/phylosolid_env/bin/Rscript)
        echo "Circos skipped: refusing to use the tree-building env ($RSCRIPT)."
        echo "  Use phylosolid_vis or circos_env instead."
        exit 0
        ;;
esac

if ! "$RSCRIPT" -e 'quit(status = if (requireNamespace("PhyloSOLIDvis", quietly=TRUE)) 0 else 2)' >/dev/null 2>&1; then
    echo "Circos skipped: $RSCRIPT does not have PhyloSOLIDvis installed."
    exit 0
fi

echo "========================================"
echo "Demo visualization: PhyloSOLIDvis circos"
echo "========================================"
echo "  Rscript:     $RSCRIPT"
echo "  inputpath:   $PHYLO_DIR"
echo "  outputpath:  $CIRCOS_DIR"
echo "  annotation:  ${ANNOT_FILE:-"(none)"}"
echo "========================================"

mkdir -p "$CIRCOS_DIR"
CMD=("$RSCRIPT" "$R_FILE" "$PHYLO_DIR" "$CIRCOS_DIR")
if [ -n "$ANNOT_FILE" ] && [ -f "$ANNOT_FILE" ]; then
    CMD+=("$ANNOT_FILE")
fi
"${CMD[@]}"
echo "Circos written to $CIRCOS_DIR"
