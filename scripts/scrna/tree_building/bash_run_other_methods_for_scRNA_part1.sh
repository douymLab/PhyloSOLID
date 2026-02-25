#!/bin/bash
# Tree building for scRNA-seq data

set -e

startTime=`date +%Y%m%d-%H:%M:%S`
startTime_s=`date +%s`
echo " ----- Start at : $startTime -----"

# Get script directory and project root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"

# Add src directory to Python path for core modules
export PYTHONPATH="${PROJECT_ROOT}/src:${PYTHONPATH}"

# Parse command line arguments
export workpath=$1
cd ${workpath}
export sample1=$2
export cellnum=$3
export celltype_file=$4

echo "===== Tree building parameters ====="
echo "workpath: ${workpath}"
echo "sample1: ${sample1}"
echo "cellnum: ${cellnum}"
echo "celltype_file: ${celltype_file}"
echo "====================================="

echo "===== Start running: ====="

################################################################
# 1. PhyloSOLID
echo "##### Run PhyloSOLID for tree"
echo "${workpath}/results will be generated"

# Run the main PhyloSOLID script
python -m phylosolid.run_phylosilid_fullTree_scRNA \
    --sampleid ${sample1} \
    --inputpath ${workpath}/data \
    --outputpath ${workpath}/results \
    --celltype_file ${celltype_file} \
    --features_file ${workpath}/treeinput/features_file.txt \
    --is_predict_germ yes

# Copy and convert results
cd ${workpath}/results
mkdir -p ${workpath}/results/PhyloSOLID
cd ${workpath}/results/PhyloSOLID

# Copy CFMatrix
if [ -f "${workpath}/results/scaffold_builder/phylo_scaffold_tree/final_cleaned_M_scaffold_basedPivots.filtered_sites_inferred.CFMatrix" ]; then
    cp ${workpath}/results/scaffold_builder/phylo_scaffold_tree/final_cleaned_M_scaffold_basedPivots.filtered_sites_inferred.CFMatrix \
       ${workpath}/results/PhyloSOLID/cell_by_mut.CFMatrix
else
    echo "Warning: CFMatrix file not found"
fi

# Convert tree format
if [ -f "${SCRIPT_DIR}/convert_PhyloSOLID_tree.R" ]; then
    Rscript ${SCRIPT_DIR}/convert_PhyloSOLID_tree.R \
        ${workpath}/results/PhyloSOLID/cell_by_mut.CFMatrix \
        ${workpath}/results/PhyloSOLID/celltree.newick
else
    echo "Warning: convert_PhyloSOLID_tree.R not found"
fi

endTime=`date +%Y%m%d-%H:%M:%S`
endTime_s=`date +%s`
echo " ----- End at : $startTime ----> $endTime -----"
sumTime=$[ $endTime_s - $startTime_s ]
echo "Total: $sumTime seconds"
