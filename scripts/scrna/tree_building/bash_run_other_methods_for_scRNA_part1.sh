#!/bin/bash
# Tree building for scRNA-seq data
# This script runs the PhyloSOLID tree building pipeline

set -e

startTime=`date +%Y%m%d-%H:%M:%S`
startTime_s=`date +%s`
echo " ----- Start at : $startTime -----"

# Get script directory and project root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"

# Load configuration for conda python path
CONFIG_FILE="${PROJECT_ROOT}/config/paths.yaml"
CONDA_PYTHON="python"  # Default to system python

if [ -f "$CONFIG_FILE" ]; then
    # Extract conda python path from config file
    CONDA_PYTHON=$(grep -A5 'conda:' "$CONFIG_FILE" | grep 'python:' | head -1 | sed -E 's/.*python: *"?([^"]*)"?.*/\1/')
    if [ -z "$CONDA_PYTHON" ]; then
        CONDA_PYTHON="python"
        echo "Warning: conda python not found in config, using system python"
    else
        echo "Using conda Python: $CONDA_PYTHON"
    fi
fi

# Add project directories to Python path for core modules
export PYTHONPATH="${PROJECT_ROOT}:${PROJECT_ROOT}/src:${PYTHONPATH}"

# Parse command line arguments
export workpath=$1
cd ${workpath}
export sample_id=$2
export cellnum=$3
export celltype_input=$4  # Original input path (may be relative)

# Convert celltype file path to absolute if it's relative
celltype_file=""
if [ -n "${celltype_input}" ] && [ "${celltype_input}" != "None" ] && [ "${celltype_input}" != "none" ]; then
    if [[ "${celltype_input}" == /* ]]; then
        # Already absolute path
        celltype_file="${celltype_input}"
    else
        # Relative path - convert to absolute using PROJECT_ROOT
        celltype_file="${PROJECT_ROOT}/${celltype_input}"
        echo "Converted relative path to absolute: ${celltype_file}"
    fi
else
    celltype_file="None"
fi

echo "===== Tree building parameters ====="
echo "workpath: ${workpath}"
echo "sample_id: ${sample_id}"
echo "cellnum: ${cellnum}"
echo "celltype_input: ${celltype_input}"
echo "celltype_file (absolute): ${celltype_file}"
echo "====================================="

echo "===== Start running: ====="

################################################################
# 1. PhyloSOLID tree building
echo "##### Run PhyloSOLID for tree"
echo "${workpath}/03_tree_building will be generated"

# Define paths
FEATURES_FILE="${workpath}/01_features/${sample_id}.benchmark_patched.feature.txt"
INPUT_PATH="${workpath}/02_treeinput/data"
OUTPUT_PATH="${workpath}/03_tree_building"

# Check if features file exists
if [ ! -f "$FEATURES_FILE" ]; then
    echo "Warning: Features file not found at $FEATURES_FILE"
    echo "Trying alternative path..."
    FEATURES_FILE="${workpath}/treeinput/features_file.txt"
    if [ ! -f "$FEATURES_FILE" ]; then
        echo "Error: Features file not found. Cannot proceed with tree building."
        exit 1
    fi
fi
echo "Using features file: $FEATURES_FILE"

# Build the base command
CMD="${CONDA_PYTHON} -m src.run_phylosilid_fullTree_scRNA \
    --sampleid ${sample_id} \
    --inputpath ${INPUT_PATH} \
    --outputpath ${OUTPUT_PATH} \
    --features_file ${FEATURES_FILE} \
    --is_predict_germ yes"

# Handle celltype file parameter
if [ "${celltype_file}" != "None" ] && [ -n "${celltype_file}" ]; then
    if [ -f "${celltype_file}" ]; then
        CMD="${CMD} --celltype_file ${celltype_file}"
        echo "✓ Using celltype file: ${celltype_file}"
        
        # Show first few lines of the file for verification
        echo "First few lines of celltype file:"
        head -3 "${celltype_file}" | sed 's/^/    /'
    else
        echo "⚠ Warning: Celltype file does not exist: ${celltype_file}"
        echo "Current directory: $(pwd)"
        echo "Parent directory contents:"
        ls -la "$(dirname "${celltype_file}")" 2>/dev/null || echo "Directory not found"
        echo "Explicitly setting celltype_file to None"
        CMD="${CMD} --celltype_file None"
    fi
else
    echo "No celltype file provided, explicitly setting to None"
    CMD="${CMD} --celltype_file None"
fi

# Print the full command for debugging
echo "Full command:"
echo "${CMD}"

# Execute the command
echo "Executing command..."
eval $CMD

# Check if command was successful
EXIT_CODE=$?
if [ $EXIT_CODE -ne 0 ]; then
    echo "Error: PhyloSOLID tree building failed with exit code $EXIT_CODE"
    exit $EXIT_CODE
fi

################################################################
# 2. Copy and convert results
echo "##### Processing results"

# Create PhyloSOLID results directory
mkdir -p ${OUTPUT_PATH}/PhyloSOLID

# Try to copy CFMatrix from scaffold_builder (most likely location)
CFMATRIX_SOURCE="${OUTPUT_PATH}/scaffold_builder/phylo_scaffold_tree/final_cleaned_M_scaffold_basedPivots.filtered_sites_inferred.CFMatrix"
CFMATRIX_DEST="${OUTPUT_PATH}/PhyloSOLID/cell_by_mut.CFMatrix"

if [ -f "$CFMATRIX_SOURCE" ]; then
    cp ${CFMATRIX_SOURCE} ${CFMATRIX_DEST}
    echo "✓ CFMatrix copied from scaffold_builder"
else
    # Try alternative location
    CFMATRIX_SOURCE="${OUTPUT_PATH}/mutation_integrator/phylo/final_cleaned_M_full_basedPivots.filtered_sites_inferred.CFMatrix"
    if [ -f "$CFMATRIX_SOURCE" ]; then
        cp ${CFMATRIX_SOURCE} ${CFMATRIX_DEST}
        echo "✓ CFMatrix copied from mutation_integrator"
    else
        echo "Warning: CFMatrix file not found in expected locations"
    fi
fi

# Convert tree format using R script
CONVERT_SCRIPT="${SCRIPT_DIR}/convert_PhyloSOLID_tree.R"
if [ -f "$CONVERT_SCRIPT" ]; then
    if [ -f "$CFMATRIX_DEST" ]; then
        Rscript ${CONVERT_SCRIPT} \
            ${CFMATRIX_DEST} \
            ${OUTPUT_PATH}/PhyloSOLID/celltree.newick
        echo "✓ Tree format converted"
    else
        echo "Warning: CFMatrix not found, skipping tree conversion"
    fi
else
    echo "Warning: convert_PhyloSOLID_tree.R not found at ${CONVERT_SCRIPT}"
fi

# Check final tree file
if [ -f "${OUTPUT_PATH}/PhyloSOLID/celltree.newick" ]; then
    echo "✓ Tree file created: ${OUTPUT_PATH}/PhyloSOLID/celltree.newick"
else
    echo "Warning: Tree file not created"
fi

echo "===== Tree building completed ====="

endTime=`date +%Y%m%d-%H:%M:%S`
endTime_s=`date +%s`
echo " ----- End at : $startTime ----> $endTime -----"
sumTime=$[ $endTime_s - $startTime_s ]
echo "Total: $sumTime seconds"