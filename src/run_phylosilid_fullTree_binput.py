# Date: 2026/03/13
# Update: 2026/04/09
# Author: Qing
# Work: Test scaffold builder


##### Time #####
import time
start_time = time.perf_counter()

import sys
sys.path.append('/storage/douyanmeiLab/yangqing/tools/PhyloMosaicGenie/pmg/src')

################################################################################################
########################################## PhyloSOLID #########################################
################################################################################################
import os
import logging
import copy
import random
import pandas as pd
import numpy as np
from tqdm import tqdm
from copy import deepcopy
from pathlib import Path
import shutil

logger = logging.getLogger(__name__)

# from phylosolid import config
from phylosolid.data_loader import load_all
from phylosolid.scdna_classifier import real_time_classifier_predict
from phylosolid.germline_filter import identify_germline_variants
from phylosolid.germline_filter import *
from phylosolid.scaffold_builder import build_scaffold_tree
from phylosolid.scaffold_builder import *
from phylosolid.mutation_integrator import *

# 设置所有随机种子
RANDOM_SEED = 42  # 可以任意指定
random.seed(RANDOM_SEED)
np.random.seed(RANDOM_SEED)

# 如果有使用其他库也设置种子
try:
    import torch
    torch.manual_seed(RANDOM_SEED)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(RANDOM_SEED)
        torch.cuda.manual_seed_all(RANDOM_SEED)
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False
except ImportError:
    pass

# ------------------------------
# 配置 logging
# ------------------------------
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

# ------------------------------
# 项目参数与文件路径设置
# ------------------------------

import multiprocessing as mp
import argparse
from argparse import ArgumentParser
parser = argparse.ArgumentParser()

parser.add_argument("-s", "--sampleid", default="sampleid", type=str, help="The sampleid you can set and check.")
parser.add_argument("-o", "--outputpath", default="/storage/douyanmeiLab/yangqing/tools/PhyloMosaicGenie/pmg/src/phylosolid/samples/phylo_100k.30true/results", type=str, help="The outputpath you want to save results.")
parser.add_argument("-i", "--inputfile", default="/storage/douyanmeiLab/yangqing/tools/PhyloMosaicGenie/pmg/wenxuan/tree_vita_20260313/version3/str_info_cleaned_add_empty_new.tsv", type=str, help="The input file is formatted as a binary matrix, where rows represent cells and columns represent mutations.")

args = parser.parse_args()


# get parameters
sampleid = args.sampleid
inputfile = args.inputfile
outputpath = args.outputpath

outputpath_scaffold = os.path.join(outputpath, "processing")
outputpath_results = os.path.join(outputpath, "phylo")


# 递归创建目录，如果目录已存在也不会报错
os.makedirs(outputpath, exist_ok=True)
os.makedirs(outputpath_scaffold, exist_ok=True)
os.makedirs(outputpath_results, exist_ok=True)


# Display parameters for verification
logger.info(f"sampleid: {sampleid}")
logger.info(f"inputfile: {inputfile}")
logger.info(f"outputpath: {outputpath}")




# ------------------------------
# 参数设置
# ------------------------------

# tumor_cell_of_origin_mutations = ["chr1_46278500_C_T", "chr16_67944430_G_C", "chr17_29572953_G_A"]
# for i in tumor_cell_of_origin_mutations:
#     count_result = count_list(I_attached[i])
#     NA_prop = count_result['NA']/(count_result['NA']+count_result[0]+count_result[1])
#     print(count_result)
#     print(NA_prop)

# # {'NA': 1032, 0.0: 149, 1.0: 72}
# # {'NA': 1128, 0.0: 95, 1.0: 30}
# # {'NA': 1146, 0.0: 65, 1.0: 42}

# Default parameters matching Methods Section 3
SETTING_PARAMS = {
    "models_path": "phylosolid/models/scdna",
    
    # 1 data loader
    "p_thresh": 0.9,
    
    # 2 germline filter
    "mcf_cutoff": 0.05,
    "mcn_cutoff": 5,
    
    # Pairwise correlation criteria (Section 2.1)
    "pair_N11_min": 0,                 # minimum N11 for correlation
    "jaccard_thresh": 0.2,             # Jaccard threshold for strong correlation
    "jaccard_low": 0.1,                # lower Jaccard bound for parent-child
    "fraction_parent_child_thresh": 0.9, # fraction threshold for parent-child
    
    # 3.1 Initial filtration
    "posterior_threshold": 0.5,
    "maf_max_threshold": 0.3,
    "maf_mean_threshold": 0.1,
    
    # 3.2 Coverage-based filtration
    "na_prop_thresh_global": 0.95,      # p_NA(j) ≤ 0.9 (present in >10% cells)
    "cv_thresh": 6.0,                   # CV_j < 2
    
    # 3.3 Consensus correlation graph
    "consensus_runs": 100,             # number of randomized runs
    "consensus_clone_freq_thresh": 0.1, # minimum clone frequency threshold
    "resolution_of_graph": 0.1,
        
    # 3.4 Penalty-based placement
    "general_weight_NA": 0.001,        # NA imputation penalty parameter
    "fnfp_ratio": 1,                 # false negative/false positive ratio
    "phi": 1.0,                        # BIC complexity parameter
    
    # 4.1 Dynamic programming
    "pass_tree_cutoff": 0.9,
    "unpass_tree_cutoff": 0.1,
    
    # 4.2 fp_ratio and fn_ratio
    "fp_ratio_cutoff_across_tree": 0.2,
    "fn_ratio_cutoff_across_tree": 0.9,
    "fp_ratio_cutoff_within_subclone": 0.1,
    "fp_ratio_persite_cutoff": 0.1, 
    "fp_count_persite_cutoff": 0,
    
    "fp_ratio_per_mutation_cross_all_cells_cutoff": 0.4,
    "fp_ratio_per_cell_cross_all_muts_cutoff": 0.5,
    
    "intersection_vs_fn_flipping_ratio_cutoff": 0.2,
    "intersection_cell_count_on_mutation_cutoff": 5,
    "intersection_cell_ratio_on_mutation_cutoff": 0.2,
    "intersection_count_per_cells_cutoff": 1,
    "flipping_count_fn_per_cells_cutoff": 1,
    "flipping_to_1_count_per_cells_cutoff": 2
}

params = SETTING_PARAMS




# ------------------------------
# Step 1: Load data
# ------------------------------
logger.info("===== Step1: Loading data ...")

I_raw = pd.read_csv(inputfile, sep='\t', encoding='utf-8', index_col=0)
I_raw.columns = I_raw.columns.str.replace(':', '_', regex=False)
I = I_raw.copy()
logger.info(f"Loaded data: {len(I_raw)} cells, {len(I_raw.columns)} mutations")


# Generate P_somatic, V_somatic, A_somatic, C_somatic, I_somatic, df_reads_somatic, df_features_new 几个数据框
I_somatic = I_raw.copy()
I_somatic_withNA3 = I_somatic.replace({np.nan: 3}).astype(int)

P_somatic = I_raw.copy()
V_somatic = I_raw.copy()
A_somatic = I_raw.replace({0: 0, 1: 10}).fillna(0).astype(int)
C_somatic = I_raw.replace({0: 10, 1: 10}).fillna(0).astype(int)
df_reads_somatic = pd.concat([
    pd.DataFrame([[f"{(I_raw[col]==1).sum()*10}/{I_raw[col].count()*10}" for col in I_raw.columns]], 
                 index=['bulk'], columns=I_raw.columns),
    I_raw.applymap(lambda x: "0/10" if x == 0 else ("10/10" if x == 1 else np.nan))
])
df_features = pd.DataFrame([
    (I_raw == 1).sum().astype(int),
    (I_raw == 1).sum() / len(I_raw)
], index=['mutant_cellnum', 'mutant_cell_fraction'], columns=I_raw.columns)
df_features_new, empty_mutations = update_features_matrix(I_somatic, df_reads_somatic, df_features, params["mcf_cutoff"])
df_features_new = add_mutation_proportions_to_features(df_features_new, I_somatic)

somatic_mutations = list(I_raw.columns)




# ------------------------------
# Step 2: Scaffold builder
# ------------------------------
logger.info("===== Step2: Construct scaffold tree ...")

##### Load Celltype Data
barcodes = df_reads_somatic.index.tolist()  # 获取所有条形码
df_celltype = pd.DataFrame({
    "barcode": barcodes,
    "cell_type": ["default_type"] * len(barcodes)
})

df_celltype.to_csv(os.path.join(outputpath_scaffold, "df_celltype.txt"), sep="\t")
logger.info(f"Celltype data loaded: {df_celltype.shape[0]} cells")


##### Scaffold builder
logging.info("Running scaffold building ...")
immune_mutations = ['chr14_105707803_C_A', 'chr14_106276422_A_G', 'chr14_106276454_C_G', 'chr22_22881588_G_A', 'chr22_22895482_G_A', 'chr22_22895504_C_A', 'chr22_22901226_C_G', 'chr22_22901169_G_C', 'chr22_22895709_C_G']
results_of_scaffold = build_scaffold_tree(
    P_somatic = P_somatic, 
    V_somatic = V_somatic, 
    A_somatic = A_somatic, 
    C_somatic = C_somatic, 
    I_somatic = I_somatic,
    df_reads_somatic = df_reads_somatic,
    df_features_new = df_features_new,
    params = params,
    is_filter_quality = "no",
    outputpath = outputpath_scaffold,
    sampleid = sampleid,
    immune_mutations = immune_mutations,
    df_celltype = df_celltype
)


T_scaffold, M_scaffold, df_flipping_spots, df_total_flipping_count, final_cleaned_I_selected_withNA3, final_cleaned_M_scaffold, backbone_mutations, mutation_group, spots_to_split, group_mutations, remained_mutations, high_cv_mutations = results_of_scaffold
# scaffold_mutations = [i for i in initial_scaffold_mutations if i not in remained_mutations_by_scaffold_building]
scaffold_mutations = list(M_scaffold.columns)
non_scaffold_mutations = [i for i in somatic_mutations if i not in scaffold_mutations]


print_tree(T_scaffold)


logging.info(f"Identified {len(scaffold_mutations)} scaffold variants")

# logging.info(f"The number of all_mutations is: {len(all_mutations)}")
# logging.info(f"The number of predicted_germline_mutations is: {len(predicted_germline_mutations)}")
logging.info(f"The number of somatic_mutations is: {len(somatic_mutations)}")
logging.info(f"The number of scaffold_mutations is: {len(scaffold_mutations)}")
logging.info(f"The number of non_scaffold_mutations is: {len(non_scaffold_mutations)}")



# ------------------------------
# Move files to destination
# ------------------------------
source = Path(outputpath_scaffold) / "phylo_scaffold_tree"
destination = Path(outputpath_results)

for item in source.iterdir():
    shutil.move(str(item), str(destination / item.name))

source.rmdir()



# ------------------------------
# End of Process
# ------------------------------
end_time = time.perf_counter


##### Time #####
finish_time = time.perf_counter()
print("Program finished in {:.4f} seconds".format(finish_time-start_time))




