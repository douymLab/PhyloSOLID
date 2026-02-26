
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
python -m cli.main --workdir demo/test_output scrna \
    --sample Org4S15D63 \
    --mutation-list demo/input/Org4S15D63/02_identifier/identifier.txt \
    --bam demo/input/Org4S15D63/01_rawdata/Org4S15D63.bam \
    --barcode demo/input/Org4S15D63/01_rawdata/Org4S15D63_CB.txt \
    --threads 4 \
    --read-len 100 \
    --cellnum 155
