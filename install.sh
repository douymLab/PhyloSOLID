#!/bin/bash
# install.sh - Install conda environment and dependencies only
# Usage: bash install.sh

set -e

echo "========================================="
echo "PhyloSOLID Environment Installation"
echo "========================================="

# Check if conda or mamba is installed
if ! command -v mamba &> /dev/null; then
    if ! command -v conda &> /dev/null; then
        echo "[ERROR] Neither conda nor mamba found"
        echo "Please install Miniconda or Mambaforge first"
        echo "Visit: https://docs.conda.io/en/latest/miniconda.html"
        exit 1
    else
        CONDA_CMD="conda"
        echo "[INFO] Using conda"
    fi
else
    CONDA_CMD="mamba"
    echo "[INFO] Using mamba"
fi

# Step 1: Create conda environment
echo ""
echo "[Step 1/2] Creating Conda environment..."
$CONDA_CMD env create -f environment.yml

# Get environment path
ENV_NAME="phylosolid_env"
ENV_PREFIX=$($CONDA_CMD info --base)/envs/$ENV_NAME

if [ ! -d "$ENV_PREFIX" ]; then
    echo "[ERROR] Environment creation failed"
    exit 1
fi

echo "[INFO] Environment created: $ENV_PREFIX"

# Step 2: Install converTree (needed for Newick conversion)
echo ""
echo "[Step 2/2] Installing converTree R package..."
$ENV_PREFIX/bin/R -e "
options(repos = c(CRAN = 'https://cloud.r-project.org'))
if (!requireNamespace('devtools', quietly=TRUE)) {
    install.packages('devtools')
}
if (!requireNamespace('remotes', quietly=TRUE)) {
    install.packages('remotes')
}
if (!requireNamespace('converTree', quietly=TRUE)) {
    message('Installing converTree...')
    remotes::install_github('xiayh17/converTree', upgrade='never', quiet=FALSE)
} else {
    message('[INFO] converTree already installed')
}
library(converTree)
message('[INFO] converTree loaded successfully')
"

echo ""
echo "========================================="
echo "Environment setup complete"
echo "========================================="
echo ""
echo "Next steps:"
echo "  1. Activate environment: conda activate $ENV_NAME"
echo "  2. Configure config/paths.yaml with your paths"
echo "  3. Install PhyloSOLID: pip install -e ."
echo "  4. Verify installation: phylosolid check-annovar --config config/paths.yaml"
echo "  5. Optional circos env (separate from this one): bash install_vis.sh"
echo "========================================="
