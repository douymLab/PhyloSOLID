#!/bin/bash
# Feature extraction pipeline for scRNA-seq data
# This script calls multiple Python and R scripts to extract features from BAM files

set -e  # Exit on error

startTime=`date +%Y%m%d-%H:%M:%S`
startTime_s=`date +%s`
echo " ----- Start at : $startTime -----"

# Get script directory and project root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"

# Load configuration
CONFIG_FILE="${PROJECT_ROOT}/config/paths.yaml"

# Initialize variables
CONDA_PYTHON="python"  # Default to system python
ANNOVAR_DIR=""
HUMANDB=""
GENOME_FASTA=""
GFF3_FILE=""
MAPPABILITY_FILE=""
GNOMAD_FILE=""
RNA_EDITING_FILE=""
RUN_GET_PRIOR=""

if [ -f "$CONFIG_FILE" ]; then
    echo "Loading configuration from $CONFIG_FILE"
    
    # Extract conda python path
    CONDA_PYTHON=$(grep -A5 'conda:' "$CONFIG_FILE" | grep 'python:' | head -1 | sed -E 's/.*python: *"?([^"]*)"?.*/\1/')
    if [ -z "$CONDA_PYTHON" ]; then
        CONDA_PYTHON="python"
        echo "Warning: conda python not found in config, using system python"
    else
        echo "Using conda Python: $CONDA_PYTHON"
    fi
    
    # Extract ANNOVAR paths
    ANNOVAR_DIR=$(grep -A5 'annovar:' "$CONFIG_FILE" | grep 'script_dir:' | head -1 | sed -E 's/.*script_dir: *"?([^"]*)"?.*/\1/')
    HUMANDB=$(grep -A5 'annovar:' "$CONFIG_FILE" | grep 'humandb:' | head -1 | sed -E 's/.*humandb: *"?([^"]*)"?.*/\1/')
    
    # Extract reference paths
    GENOME_FASTA=$(grep -A10 'reference:' "$CONFIG_FILE" | grep 'genome_fasta:' | head -1 | sed -E 's/.*genome_fasta: *"?([^"]*)"?.*/\1/')
    GFF3_FILE=$(grep -A10 'reference:' "$CONFIG_FILE" | grep 'gff3_file:' | head -1 | sed -E 's/.*gff3_file: *"?([^"]*)"?.*/\1/')
    MAPPABILITY_FILE=$(grep -A10 'reference:' "$CONFIG_FILE" | grep 'mappability_file:' | head -1 | sed -E 's/.*mappability_file: *"?([^"]*)"?.*/\1/')
    GNOMAD_FILE=$(grep -A10 'reference:' "$CONFIG_FILE" | grep 'gnomad_file:' | head -1 | sed -E 's/.*gnomad_file: *"?([^"]*)"?.*/\1/')
    RNA_EDITING_FILE=$(grep -A10 'reference:' "$CONFIG_FILE" | grep 'rna_editing_file:' | head -1 | sed -E 's/.*rna_editing_file: *"?([^"]*)"?.*/\1/')
    RUN_GET_PRIOR=$(grep -A10 'reference:' "$CONFIG_FILE" | grep 'run_get_prior:' | head -1 | sed -E 's/.*run_get_prior: *"?([^"]*)"?.*/\1/')
    
    # Output read values for debugging
    echo "Read ANNOVAR_DIR: $ANNOVAR_DIR"
    echo "Read HUMANDB: $HUMANDB"
    
else
    echo "Warning: Config file not found at $CONFIG_FILE"
    echo "Using default paths from environment variables"
fi

# Remove trailing slashes
ANNOVAR_DIR=$(echo "$ANNOVAR_DIR" | sed 's/\/$//')
HUMANDB=$(echo "$HUMANDB" | sed 's/\/$//')

# Validate required paths
if [ -z "$ANNOVAR_DIR" ] || [ ! -d "$ANNOVAR_DIR" ]; then
    echo "ERROR: ANNOVAR directory not found. Please install ANNOVAR and set script_dir in config/paths.yaml"
    echo "Current ANNOVAR_DIR: $ANNOVAR_DIR"
    exit 1
fi

if [ ! -f "${ANNOVAR_DIR}/annotate_variation.pl" ]; then
    echo "ERROR: annotate_variation.pl not found in $ANNOVAR_DIR"
    ls -la "$ANNOVAR_DIR" | head -10
    exit 1
fi

if [ -z "$HUMANDB" ] || [ ! -d "$HUMANDB" ]; then
    echo "ERROR: ANNOVAR humandb not found. Please download databases and set humandb in config/paths.yaml"
    echo "Current HUMANDB: $HUMANDB"
    exit 1
fi

# Add script directory to Python path
export PYTHONPATH="${SCRIPT_DIR}:${PYTHONPATH}"

# Parse command line arguments
export sampleid=$1
export identifier_file=$2
export bamfile=$3
export readLen=$4
export out_dir_name=$5
export filtered_bam=$6
export thread=$7
export barcode_file=$8
export ase_filepath=$9
export running_type=${10}

prefix_name=${sampleid}.${running_type}_identifier
new_name="${out_dir_name}/${sampleid}.${running_type}_patched.feature.txt"

echo "===== Feature extraction parameters ====="
echo "sampleid: $sampleid"
echo "identifier_file: $identifier_file"
echo "bamfile: $bamfile"
echo "readLen: $readLen"
echo "out_dir_name: $out_dir_name"
echo "thread: $thread"
echo "barcode_file: $barcode_file"
echo "ase_filepath: $ase_filepath"
echo "running_type: $running_type"
echo "=========================================="

# Create output directory
mkdir -p ${out_dir_name}

# Get bed file of mutant sites
echo "Creating BED file from mutation list..."
awk -F"_" ' {print $1, $2, $2, $3, $4}' OFS="\t" ${identifier_file} > ${out_dir_name}/${prefix_name}.bed

################################################################
# Annovar annotation
echo "============ Annovar annotation: Start ============"
${ANNOVAR_DIR}/annotate_variation.pl -build hg38 -out ${out_dir_name}/${prefix_name}.anno \
    -dbtype wgEncodeGencodeBasicV46 ${out_dir_name}/${prefix_name}.bed ${HUMANDB}/

# Direct script path instead of -m
${CONDA_PYTHON} ${SCRIPT_DIR}/get_mutation_anno_table.py \
    -m ${identifier_file} \
    -i ${out_dir_name}/${prefix_name}.anno \
    -o ${out_dir_name}/${prefix_name}.anno_info.txt
echo "============ Annovar annotation: Completed ============"

################################################################
# Prior from gnomad
echo "============ Get prior from gnomad: Start ============"
if [ -f "$RUN_GET_PRIOR" ]; then
    ${CONDA_PYTHON} $RUN_GET_PRIOR splitprior \
        -s ${sampleid} \
        --bedfile ${out_dir_name}/${prefix_name}.bed \
        --outdir ${out_dir_name} \
        --outname .${prefix_name}.prior.out \
        --annovar $GNOMAD_FILE \
        --thread ${thread}
else
    echo "WARNING: run_get_prior.py not found, skipping prior calculation"
fi
rm -rf ${out_dir_name}/tmp
echo "============ Get prior from gnomad: Completed ============"

################################################################
# Get features
echo "============ Get features: Start ============"
# Direct script path instead of -m
${CONDA_PYTHON} ${SCRIPT_DIR}/extract_features.py \
    --mutations ${identifier_file} \
    --fasta $GENOME_FASTA \
    --gff3_file $GFF3_FILE \
    --bam ${filtered_bam} \
    --outdir ${out_dir_name} \
    --thread ${thread} \
    --outname ${prefix_name}.feature.txt \
    --mappbablity_file $MAPPABILITY_FILE \
    --annovar_annotaion_file ${out_dir_name}/${prefix_name}.anno_info.txt \
    --sample ${sampleid} \
    --readLen ${readLen} \
    --prior ${out_dir_name}/${sampleid}.${prefix_name}.prior.out
echo "============ Get features: Completed ============"

################################################################
# Get depth in mutant spots and unmutant spots
echo "============ Get depth in spots: Start ============"
# Direct script path instead of -m
${CONDA_PYTHON} ${SCRIPT_DIR}/get_mutspots_new_for_dp_and_features.py \
    --bam ${bamfile} \
    --mutation_list ${identifier_file} \
    --barcode_file ${barcode_file} \
    --run_type read \
    --features ${out_dir_name}/${prefix_name}.feature.txt \
    --outdir ${out_dir_name}/depth_in_spots \
    --threads ${thread}
echo "============ Get depth in spots: Completed ============"

################################################################
# Calculate R square between mutant allele number and expression
echo "============ Calculate R square: Start ============"
${CONDA_PYTHON} ${SCRIPT_DIR}/cal_Rsquare_as_new_feature_for_scRNA.py \
    --feature_file ${out_dir_name}/${prefix_name}.feature_depth.txt \
    --reads_filepath ${out_dir_name}/depth_in_spots \
    --output_file ${out_dir_name}/${prefix_name}.feature_depth_Rsquare.txt
echo "============ Calculate R square: Completed ============"

################################################################
# Mutation Signature Filter
echo "============ Mutation Signature Filtering: Start ============"
Rscript ${SCRIPT_DIR}/hardFilter_and_cal_mutation_sig.R \
    -i ${out_dir_name}/${prefix_name}.feature_depth_Rsquare.txt \
    -o ${out_dir_name}/${prefix_name}.feature_depth_Rsquare_sigFilter.txt
echo "============ Mutation Signature Filtering: Completed ============"

################################################################
# ASE Filter & RNA Editing
echo "============ Calculating ASE & RNA Editing Status: Start ============"

if [ "${ase_filepath}" != "no" ]; then
    echo "----- Both ASE and RNA Editing Filtering: Start -----"

    tail -n +2 ${out_dir_name}/${prefix_name}.feature_depth_Rsquare_sigFilter.txt | cut -f4 > ${out_dir_name}/${prefix_name}.mutation_list.txt

    # Direct script path instead of -m
    ${CONDA_PYTHON} ${SCRIPT_DIR}/filter_ASE.py \
        --bam ${bamfile} \
        --RNA_editing $RNA_EDITING_FILE \
        --ase_file ${ase_filepath}/${sampleid}.candidate_ASE_sites.txt \
        --mutation_list ${out_dir_name}/${prefix_name}.mutation_list.txt \
        --outdir ${out_dir_name}/${prefix_name}.ASE_RNAediting
    echo "----- Both ASE and RNA Editing Filtering: Completed -----"

    ################################################################
    # Applying ASE_filter/RNAediting_filter/UMI_consistence_filter
    echo "============ Applying filters: Start ============"
    Rscript ${SCRIPT_DIR}/hard-filter_using_ASE_RNAediting_UMI_consistence_filter.R \
        --feature_file ${out_dir_name}/${prefix_name}.feature_depth_Rsquare_sigFilter.txt \
        --ase_editing ${out_dir_name}/${prefix_name}.ASE_RNAediting/identifer_stat.txt \
        --outputfile ${out_dir_name}/${prefix_name}.feature_depth_Rsquare_sigFilter_ASE_RNAediting_UMIconsistence.txt
    echo "============ Applying filters: Completed ============"

    cp ${out_dir_name}/${prefix_name}.feature_depth_Rsquare_sigFilter_ASE_RNAediting_UMIconsistence.txt ${out_dir_name}/${prefix_name}.feature_final.txt

else
    echo "----- Only RNA Editing Filtering: Start -----"
    bedtools intersect -a <(tail -n +2 ${out_dir_name}/${prefix_name}.feature_depth_Rsquare_sigFilter.txt | awk 'BEGIN {OFS="\t"} $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/') -b <(awk 'BEGIN {OFS="\t"} $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/' $RNA_EDITING_FILE) -wa | cut -f 1-3 > ${out_dir_name}/${sampleid}.overlapping_sites_for_RNAediting.bed
    
    if [ ! -s test.bed ]; then
        echo -e "chr\tstart\tend" > ${out_dir_name}/${sampleid}.overlapping_sites_for_RNAediting.bed
    else
        sed -i '1i chr\tstart\tend' ${out_dir_name}/${sampleid}.overlapping_sites_for_RNAediting.bed
    fi
    
    awk 'BEGIN {FS=OFS="\t"} 
         FNR==NR {overlap[$1, $2, $3]=1; next}
         FNR==1 {print $0, "RNAediting_filter"; next}
         (($1, $2, $3) in overlap) {print $0, "fail"}
         !(($1, $2, $3) in overlap) {print $0, "pass"}' \
         ${out_dir_name}/${sampleid}.overlapping_sites_for_RNAediting.bed \
         ${out_dir_name}/${prefix_name}.feature_depth_Rsquare_sigFilter.txt \
         > ${out_dir_name}/${prefix_name}.feature_depth_Rsquare_sigFilter_RNAediting.txt
    echo "----- Only RNA Editing Filtering: Completed -----"

    cp ${out_dir_name}/${prefix_name}.feature_depth_Rsquare_sigFilter_RNAediting.txt ${out_dir_name}/${prefix_name}.feature_final.txt
fi

################################################################
# Extract and correct features by patch code
echo "============ Features patching: Start ============"

prior_arg=""
if [ -f "${out_dir_name}/${sampleid}.${prefix_name}.prior.out" ]; then
    prior_arg="--prior ${out_dir_name}/${sampleid}.${prefix_name}.prior.out"
fi

# Direct script path instead of -m
${CONDA_PYTHON} ${SCRIPT_DIR}/extract_feature_patch_add_r_add_pkl_refine_UMI_add_adj.py \
    --features ${out_dir_name}/${prefix_name}.feature_final.txt \
    --thread "${thread}" \
    --outdir "${out_dir_name}" \
    --outname "${new_name}" \
    --readLen "${readLen}" \
    --bam "${bamfile}" \
    ${prior_arg}

echo "============ Features patching: Completed ============"

endTime=`date +%Y%m%d-%H:%M:%S`
endTime_s=`date +%s`
echo " ----- End at : $endTime_s -----"
sumTime=$[ $endTime_s - $startTime_s ]
echo "$startTime ---> $endTime" "Total:$sumTime seconds"