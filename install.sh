#!/bin/bash
# install.sh - One-click installation script for PhyloSOLID
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
    echo "[INFO] Using mamba (faster dependency resolution)"
fi

# Step 1: Create conda environment
echo ""
echo "[Step 1/3] Creating Conda environment (may take 5-10 minutes)..."
$CONDA_CMD env create -f environment.yml

# Get environment path
ENV_NAME="phylosolid_env"
ENV_PREFIX=$($CONDA_CMD info --base)/envs/$ENV_NAME

if [ ! -d "$ENV_PREFIX" ]; then
    echo "[ERROR] Environment creation failed"
    exit 1
fi

echo "[INFO] Environment created successfully: $ENV_PREFIX"

# Step 2: Install converTree R package
echo ""
echo "[Step 2/3] Installing converTree R package from GitHub..."
$ENV_PREFIX/bin/R -e "
if (!require('converTree', quietly=TRUE)) {
    message('Installing converTree...')
    if (!require('devtools', quietly=TRUE)) {
        install.packages('devtools', repos='https://cloud.r-project.org')
    }
    devtools::install_github('xiayh17/converTree', upgrade='never', quiet=FALSE)
    message('[INFO] converTree installed successfully')
} else {
    message('[INFO] converTree already installed, skipping')
}
"

# Step 3: Verify installation
echo ""
echo "[Step 3/3] Verifying critical dependencies..."
$ENV_PREFIX/bin/python -c "
import sys
packages = ['numpy', 'pandas', 'geopandas', 'vcfpy', 'scphylo', 'torch']
failed = []
for pkg in packages:
    try:
        __import__(pkg)
        print(f'[OK] {pkg}')
    except ImportError as e:
        print(f'[FAIL] {pkg}: {e}')
        failed.append(pkg)

if failed:
    print(f'\n[WARNING] Failed to import: {failed}')
    print('You may need to manually install these packages')
    sys.exit(1)
else:
    print('\n[OK] All core packages imported successfully')
"

echo ""
echo "========================================="
echo "Installation Complete"
echo "========================================="
echo ""
echo "Activate the environment:"
echo "  conda activate $ENV_NAME"
echo ""
echo "Run the demo:"
echo "  cd demo && bash run_demo.sh"
echo ""
echo "Troubleshooting:"
echo "  1. Check your internet connection"
echo "  2. Ensure you can access GitHub"
echo "  3. If using a proxy, set HTTP_PROXY and HTTPS_PROXY"
echo "========================================="