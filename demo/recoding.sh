
# 1. 特征提取
python -m cli.main --workdir demo/expected_output scrna \
    --sample Org4S15D63 \
    --mutation-list demo/input/Org4S15D63/02_identifier/identifier.txt \
    --bam demo/input/Org4S15D63/01_rawdata/Org4S15D63.bam \
    --barcode demo/input/Org4S15D63/01_rawdata/Org4S15D63_CB.txt \
    --threads 4 \
    --read-len 100 \
    --steps feature_extraction

# 2. Tree input
python -m cli.main --workdir demo/expected_output scrna \
    --sample Org4S15D63 \
    --mutation-list demo/input/Org4S15D63/02_identifier/identifier.txt \
    --bam demo/input/Org4S15D63/01_rawdata/Org4S15D63.bam \
    --barcode demo/input/Org4S15D63/01_rawdata/Org4S15D63_CB.txt \
    --cellnum 155 \
    --steps tree_input

# 3. 建树步骤
python -m cli.main --workdir demo/expected_output scrna \
    --sample Org4S15D63 \
    --mutation-list demo/input/Org4S15D63/02_identifier/identifier.txt \
    --bam demo/input/Org4S15D63/01_rawdata/Org4S15D63.bam \
    --barcode demo/input/Org4S15D63/01_rawdata/Org4S15D63_CB.txt \
    --celltype-file demo/input/Org4S15D63/03_celltype/celltype_file_for_Org4S15D63.txt \
    --steps tree_building


## 完整流程
python -m cli.main --workdir demo/expected_output scrna \
    --sample Org10S4D46 \
    --mutation-list demo/input/Org10S4D46/02_identifier/identifier.txt \
    --bam demo/input/Org10S4D46/01_rawdata/Org10S4D46.bam \
    --barcode demo/input/Org10S4D46/01_rawdata/Org10S4D46_CB.txt \
    --threads 4 \
    --read-len 100 \
    --cellnum 836


##### binary matrix input mode
cd /storage/douyanmeiLab/yangqing/tools/PhyloMosaicGenie/pmg/PhyloSOLID
python -m cli.main --workdir demo/expected_output binary-matrix \
    --sample UMB1465 \
    --inputfile demo/input/UMB1465/demo_input.tsv \
    --outputpath demo/expected_output/UMB1465


