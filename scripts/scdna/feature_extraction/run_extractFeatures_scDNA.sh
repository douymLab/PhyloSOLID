#!/bin/bash
#SBATCH -p dym_20240725,intel-sc3,amd-ep2
#SBATCH -q huge
#SBATCH -J batch_features
#SBATCH -c 2
#SBATCH -a 1-14
#SBATCH -o logs/batch_features_%A_%a.log
#SBATCH --mem=2G
#SBATCH --time=480:00:00


startTime=`date +%Y%m%d-%H:%M:%S`
startTime_s=`date +%s`
 
### environment
# source /storage/douyanmeiLab/zhengyunchao/miniconda3/etc/profile.d/conda.sh
# conda activate /home/douyanmeiLab/zhengyunchao/yunchao
source activate /storage/douyanmeiLab/lujinhong/tools/mambaforge/envs/breast
module load samtools/1.14


### variants
inputpath=$1
sampleName_file=$2
bamlist_path=$3
reference=$4
bw_file=$5
depth_file=$6


### inputfile
inputfile=$(ls ${inputpath}/*input.txt | sed -n "${SLURM_ARRAY_TASK_ID}p")
outputfile=$(echo $inputfile | sed 's/input\.txt$/output.bed/')

### running
python /storage/douyanmeiLab/zhengyunchao/common_mosaic_variants/scMosaicCaller/scMosaicCaller/genotyper/read_level_features_extraction.py \
        -i ${inputfile} \
        -o ${outputfile} \
        -s ${sampleName_file} \
        -b ${bamlist_path} \
        -r ${reference} \
        -u ${bw_file} \
        -d ${depth_file} \
        -n 2


endTime=`date +%Y%m%d-%H:%M:%S`
endTime_s=`date +%s`
