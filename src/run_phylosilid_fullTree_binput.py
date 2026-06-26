#!/usr/bin/env python3
import os
os.environ['PYTHONHASHSEED'] = '42'


# Date: 2026/03/13
# Update: 2026/06/18
# Author: Qing
# Work: Binary matrix mode


##### Time #####
import time
start_time = time.perf_counter()


################################################################################################
########################################## PhyloSOLID #########################################
################################################################################################
import logging
import copy
import random
import pandas as pd
import numpy as np
from tqdm import tqdm
from copy import deepcopy
from pathlib import Path
import shutil
import json

logger = logging.getLogger(__name__)

# Import from src (project structure)
from src.scaffold_builder import build_scaffold_tree
from src.scaffold_builder import *
from src.mutation_integrator import *


# ------------------------------
# Logging configuration
# ------------------------------
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

# ------------------------------
# Parse command line arguments
# ------------------------------
import multiprocessing as mp
import argparse
from argparse import ArgumentParser
parser = argparse.ArgumentParser()

parser.add_argument("-s", "--sampleid", default="", type=str, help="Sample ID")
parser.add_argument("-o", "--outputpath", default="./output", type=str, help="Output path for results")
parser.add_argument("-i", "--inputfile", default="", type=str, help="Input binary matrix file (rows=cells, columns=mutations)")
parser.add_argument("--seed", default=42, type=int, help="Random seed for reproducibility")

args = parser.parse_args()

# 设置所有随机种子（统一管理）
from src.reproducibility import set_seed, deterministic_choice
set_seed(args.seed)


# Get parameters
sampleid = args.sampleid
inputfile = args.inputfile
outputpath = args.outputpath + "/" + sampleid

outputpath_scaffold = os.path.join(outputpath, "processing")
outputpath_results = os.path.join(outputpath, "phylo")

# Create directories
os.makedirs(outputpath, exist_ok=True)
os.makedirs(outputpath_scaffold, exist_ok=True)
os.makedirs(outputpath_results, exist_ok=True)

# Display parameters
logger.info(f"sampleid: {sampleid}")
logger.info(f"inputfile: {inputfile}")
logger.info(f"outputpath: {outputpath}")

# ------------------------------
# Parameter settings (matching Methods Section 3)
# ------------------------------
SETTING_PARAMS = {
    "models_path": "phylosolid/models/scdna",
    
    # 1 data loader
    "p_thresh": 0.9,
    
    # 2 germline filter
    "mcf_cutoff": 0.05,
    "mcn_cutoff": 5,
    
    # Pairwise correlation criteria (Section 2.1)
    "pair_N11_min": 0,
    "jaccard_thresh": 0.2,
    "jaccard_low": 0.1,
    "fraction_parent_child_thresh": 0.9,
    
    # 3.1 Initial filtration
    "posterior_threshold": 0.5,
    "maf_max_threshold": 0.3,
    "maf_mean_threshold": 0.1,
    
    # 3.2 Coverage-based filtration
    "na_prop_thresh_global": 0.95,
    "cv_thresh": 6.0,
    
    # 3.3 Consensus correlation graph
    "consensus_runs": 100,
    "consensus_clone_freq_thresh": 0.1,
    "resolution_of_graph": 0.1,
    
    # 3.4 Penalty-based placement
    "general_weight_NA": 0.001,
    "fnfp_ratio": 1,
    "phi": 1.0,
    
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
    "fp_count_per_mutation_cross_all_cells_cutoff": 0,
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
# Helper functions
# ------------------------------
def update_features_matrix(I, df_reads, df_features, mcf_cutoff):
    """Update features matrix with mutation statistics"""
    # Update mutant_cell_fraction_detected
    df_features.loc['mutant_cell_fraction_detected'] = 0
    
    for col in I.columns:
        total_cells = len(I[col].dropna())
        if total_cells > 0:
            detected_fraction = (I[col] == 1).sum() / total_cells
            df_features.loc['mutant_cell_fraction_detected', col] = detected_fraction
        else:
            df_features.loc['mutant_cell_fraction_detected', col] = 0
    
    # Identify empty mutations
    empty_mutations = [col for col in I.columns if df_features.loc['mutant_cell_fraction_detected', col] == 0]
    
    return df_features, empty_mutations

def add_mutation_proportions_to_features(df_features, I):
    """Add mutation proportions to features matrix"""
    for col in I.columns:
        total_cells = len(I[col].dropna())
        if total_cells > 0:
            df_features.loc['mutant_cell_fraction', col] = (I[col] == 1).sum() / total_cells
        else:
            df_features.loc['mutant_cell_fraction', col] = 0
    return df_features

def print_tree(T):
    """Print tree structure"""
    if hasattr(T, 'print'):
        T.print()
    else:
        print("Tree object cannot be printed directly")

# ------------------------------
# Step 1: Load data
# ------------------------------
logger.info("===== Step1: Loading data ...")

I_raw = pd.read_csv(inputfile, sep='\t', encoding='utf-8', index_col=0)
I_raw.columns = I_raw.columns.str.replace(':', '_', regex=False)
I = I_raw.copy()
logger.info(f"Loaded data: {len(I_raw)} cells, {len(I_raw.columns)} mutations")

# Generate P_somatic, V_somatic, A_somatic, C_somatic, I_somatic, df_reads_somatic, df_features_new
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

somatic_mutations = list((reorder_columns_by_mutant_stats(I_raw, df_features_new)[0]).columns)

# ------------------------------
# Step 2: Scaffold builder
# ------------------------------
logger.info("===== Step2: Construct scaffold tree ...")

# Load Celltype Data
barcodes = df_reads_somatic.index.tolist()
df_celltype = pd.DataFrame({
    "barcode": barcodes,
    "cell_type": ["default_type"] * len(barcodes)
})

df_celltype.to_csv(os.path.join(outputpath_scaffold, "df_celltype.txt"), sep="\t")
logger.info(f"Celltype data loaded: {df_celltype.shape[0]} cells")

# Run scaffold builder
logger.info("Running scaffold building ...")
immune_mutations = []

results_of_scaffold = build_scaffold_tree(
    P_somatic=P_somatic,
    V_somatic=V_somatic,
    A_somatic=A_somatic,
    C_somatic=C_somatic,
    I_somatic=I_somatic,
    df_reads_somatic=df_reads_somatic,
    df_features_new=df_features_new,
    params=params,
    is_filter_quality="no",
    outputpath=outputpath_scaffold,
    sampleid=sampleid,
    immune_mutations=immune_mutations,
    df_celltype=df_celltype
)

# Unpack results
T_scaffold, M_scaffold, df_flipping_spots, df_total_flipping_count, final_cleaned_I_selected_withNA3, \
final_cleaned_M_scaffold, backbone_mutations, mutation_group, spots_to_split, group_mutations, \
remained_mutations, high_cv_mutations = results_of_scaffold

scaffold_mutations = list(M_scaffold.columns)
non_scaffold_mutations = [i for i in somatic_mutations if i not in scaffold_mutations]

print_tree(T_scaffold)

logger.info(f"Identified {len(scaffold_mutations)} scaffold variants")
logger.info(f"Number of somatic_mutations: {len(somatic_mutations)}")
logger.info(f"Number of scaffold_mutations: {len(scaffold_mutations)}")
logger.info(f"Number of non_scaffold_mutations: {len(non_scaffold_mutations)}")

# ------------------------------
# Move files to destination
# ------------------------------
source = Path(outputpath_scaffold) / "phylo_scaffold_tree"
destination = Path(outputpath_results)

for item in source.iterdir():
    shutil.move(str(item), str(destination / item.name))

try:
    source.rmdir()
except:
    pass

# ------------------------------
# End of Process
# ------------------------------
end_time = time.perf_counter()

print("Program finished in {:.4f} seconds".format(end_time - start_time))