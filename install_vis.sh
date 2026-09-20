#!/bin/bash
# install_vis.sh - Create a separate conda env for PhyloSOLIDvis circos plots.
# Do NOT run this into the tree-building env (pmg / phylosolid_env).
# Usage: bash install_vis.sh

set -e

echo "========================================="
echo "PhyloSOLIDvis environment installation"
echo "========================================="
echo "This env is only for circos plots."
echo "Tree building stays in phylosolid_env / pmg."
echo ""

if command -v mamba &> /dev/null; then
    CONDA_CMD="mamba"
elif command -v conda &> /dev/null; then
    CONDA_CMD="conda"
else
    echo "[ERROR] Neither conda nor mamba found"
    exit 1
fi

ENV_NAME="phylosolid_vis"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ -d "$($CONDA_CMD info --base)/envs/$ENV_NAME" ]; then
    echo "[INFO] Environment already exists: $ENV_NAME"
else
    echo "[Step 1/2] Creating conda environment $ENV_NAME ..."
    $CONDA_CMD env create -f "$SCRIPT_DIR/environment_vis.yml"
fi

ENV_PREFIX=$($CONDA_CMD info --base)/envs/$ENV_NAME
if [ ! -x "$ENV_PREFIX/bin/R" ]; then
    echo "[ERROR] R not found in $ENV_PREFIX"
    exit 1
fi

echo "[Step 2/2] Installing converTree and PhyloSOLIDvis into $ENV_NAME ..."
"$ENV_PREFIX/bin/R" -e "
options(repos = c(CRAN = 'https://cloud.r-project.org'))
if (!requireNamespace('remotes', quietly=TRUE)) install.packages('remotes')
if (!requireNamespace('converTree', quietly=TRUE)) {
    remotes::install_github('xiayh17/converTree', upgrade='never')
}
if (!requireNamespace('PhyloSOLIDvis', quietly=TRUE)) {
    local_vis <- Sys.getenv('PHYLOSOLIDVIS_SRC')
    if (nzchar(local_vis) && file.exists(file.path(local_vis, 'DESCRIPTION'))) {
        message('Installing PhyloSOLIDvis from local source: ', local_vis)
        remotes::install_local(local_vis, dependencies=TRUE, upgrade='never', force=TRUE)
    } else {
        remotes::install_github('TsingYang1112/PhyloSOLIDvis', dependencies=TRUE, upgrade='never')
    }
}
library(ggplot2)
library(PhyloSOLIDvis)
message('[INFO] PhyloSOLIDvis loaded successfully')
"

echo ""
echo "========================================="
echo "Visualization environment ready"
echo "========================================="
echo "  conda activate $ENV_NAME"
echo "  Demos look for this env automatically (or circos_env)."
echo "  Override with: export PHYLOSOLID_VIS_RSCRIPT=$ENV_PREFIX/bin/Rscript"
echo "  Optional local vis source: export PHYLOSOLIDVIS_SRC=/path/to/PhyloSOLIDvis"
echo "========================================="
