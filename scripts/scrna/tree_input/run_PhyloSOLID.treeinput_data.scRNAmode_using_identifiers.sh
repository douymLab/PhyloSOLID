#!/bin/bash
# scripts/scrna/tree_input/run_PhyloSOLID.treeinput_data.scRNAmode_using_identifiers.sh

startTime=`date +%Y%m%d-%H:%M:%S`
startTime_s=`date +%s`
echo " ----- Start at : $startTime -----"

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"

# Add necessary paths to PYTHONPATH
export PYTHONPATH="${SCRIPT_DIR}:${PROJECT_ROOT}/scripts/scrna/tree_input:${PYTHONPATH}"

##### Setting and loading
export workpath=$1
cd ${workpath}
export sample1=$2
export barcode1=$3
export bam_file1=$4
export mutation_list=$5
export cellnum=$6
export thread=$7

echo "===== Check the parameters ====="
echo "workpath: ${workpath}"
echo "sample1: ${sample1}"
echo "barcode1: ${barcode1}"
echo "bam_file1: ${bam_file1}"
echo "mutation_list: ${mutation_list}"
echo "cellnum: ${cellnum}"

echo "===== Start running: ====="
################################################################
# Generate tree input raw files
echo "### Step1: Generate tree input raw files ###"
mkdir -p ${workpath}/02_treeinput/

# 使用绝对路径直接运行脚本，而不是用 -m
python ${SCRIPT_DIR}/prepare_phylo_input_from_scRNA.py \
    --barcode_files ${barcode1} \
    --bams ${bam_file1} \
    --mutlist ${mutation_list} \
    --outprefix ${workpath}/02_treeinput/treeinput \
    --samples ${sample1}

################################################################
# Data preprocess before phylogeny
echo "### Step2: Data preprocess before phylogeny ###"
Rscript ${SCRIPT_DIR}/rawPreprocess_spatial.extract_identifier_sites.R \
    --inputfile ${workpath}/02_treeinput/treeinput_spot_c_${cellnum}.csv \
    --cellnum ${cellnum} \
    --outputpath ${workpath}/02_treeinput/data \
    --scid_file ${workpath}/02_treeinput/treeinput_scid_barcode.txt \
    --is_remove_cells no \
    --threshold 0.5 \
    --indid ${sample1}

endTime=`date +%Y%m%d-%H:%M:%S`
endTime_s=`date +%s`
echo " ----- End at : $endTime_s -----"
sumTime=$[ $endTime_s - $startTime_s ]
echo "$startTime ---> $endTime" "Total:$sumTime seconds"