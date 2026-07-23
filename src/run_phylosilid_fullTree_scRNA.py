#!/usr/bin/env python3
import os
os.environ['PYTHONHASHSEED'] = '42'

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)


# Date: 2025/09/16
# Update: 2025/10/13
# Latest: 2026/07/22
# Author: Qing


##### Time #####
import time
start_time = time.perf_counter()


################################################################################################
########################################## PhyloSOLID ##########################################
################################################################################################
import logging
import copy
import random
import pandas as pd
import numpy as np
from tqdm import tqdm
from copy import deepcopy

logger = logging.getLogger(__name__)

from src.data_loader import load_all
from src.scrna_classifier import real_time_classifier_predict
from src.germline_filter import identify_germline_variants
from src.germline_filter import *
from src.scaffold_builder import build_scaffold_tree
from src.scaffold_builder import *
from src.mutation_integrator import *


# ------------------------------
# Configure logging
# ------------------------------
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

# ------------------------------
# Project parameters and file paths
# ------------------------------
import multiprocessing as mp
import argparse
from argparse import ArgumentParser
parser = argparse.ArgumentParser()

parser.add_argument("-s", "--sampleid", default="", type=str, help="The sampleid you can set and check.")
parser.add_argument("-i", "--inputpath", default="yourpath/data", type=str, help="The inputpath contains the preprocessing results from raw posterior-reads data.")
parser.add_argument("-o", "--outputpath", default="yourpath/results", type=str, help="The outputpath you want to save results.")
parser.add_argument("-c", "--celltype_file", default=None, type=str, help="The celltype_file you should provide. If you can't generate this file, please set 'None'.")
parser.add_argument("--features_file", default=None, type=str, help="The features_file you should provide. If you can't generate this file, please set 'None'.")
parser.add_argument("--is_predict_germ", default="no", choices=["yes", "no"], type=str, help="Select 'yes' or 'no' to determine whether to predict germline mutations.")
parser.add_argument("--is_detect_passtree_by_dp", default="no", choices=["yes", "no"], type=str, help="Select 'yes' or 'no' to determine whether to run Dynamic programing step.")
parser.add_argument("--is_filter_quality", default="yes", choices=["yes", "no"], type=str, help="Select 'yes' or 'no' to determine whether to filter mutations in scaffold steps by coverage quality.")
parser.add_argument("--cv_rank_thresh", default=0.3, type=float, help="CV rank threshold for coverage-based filtration (top percentage of lowest CV)")
parser.add_argument("--remove_artifact_mutations", default="yes", choices=["yes", "no"], type=str, help="Select 'yes' or 'no' to determine whether to permanently remove artifact mutations.")
parser.add_argument("--seed", default=42, type=int, help="Random seed for reproducibility")

args = parser.parse_args()

# Set all random seeds
from src.reproducibility import set_seed, deterministic_choice
set_seed(args.seed)


# get parameters
sampleid = args.sampleid
inputpath = args.inputpath
outputpath = args.outputpath
celltype_file = args.celltype_file
features_file = args.features_file
is_predict_germ = args.is_predict_germ
is_detect_passtree_by_dp = args.is_detect_passtree_by_dp
is_filter_quality = args.is_filter_quality
cv_rank_thresh = args.cv_rank_thresh
remove_artifact_mutations = args.remove_artifact_mutations


outputpath_classifier = os.path.join(outputpath, "classifier_filter")
outputpath_germline = os.path.join(outputpath, "germline_filter")
outputpath_scaffold = os.path.join(outputpath, "scaffold_builder")
outputpath_full = os.path.join(outputpath, "mutation_integrator")
os.makedirs(outputpath_classifier, exist_ok=True)
os.makedirs(outputpath_germline, exist_ok=True)
os.makedirs(outputpath_scaffold, exist_ok=True)
os.makedirs(outputpath_full, exist_ok=True)


# Display parameters for verification
logger.info(f"sampleid: {sampleid}")
logger.info(f"inputpath: {inputpath}")
logger.info(f"outputpath: {outputpath}")
logger.info(f"celltype_file: {celltype_file}")
logger.info(f"features_file: {features_file}")
logger.info(f"is_predict_germ: {is_predict_germ}")
logger.info(f"is_detect_passtree_by_dp: {is_detect_passtree_by_dp}")
logger.info(f"is_filter_quality: {is_filter_quality}")
logger.info(f"cv_rank_thresh: {cv_rank_thresh}")
logger.info(f"remove_artifact_mutations: {remove_artifact_mutations}")


# ------------------------------
# Parameter settings
# ------------------------------

SETTING_PARAMS = {
    "models_path": "phylosolid/models/scdna",
    
    # 1 data loader
    "p_thresh": 0.5,
    
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
    "cv_rank_thresh": cv_rank_thresh,
    
    # 3.3 Consensus correlation graph
    "consensus_runs": 100,
    "consensus_clone_freq_thresh": 0.1,
    "resolution_of_graph": 1,
        
    # 3.4 Penalty-based placement
    "general_weight_NA": 0.001,
    "fnfp_ratio": 0.1,
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
    
    "fp_ratio_per_mutation_cross_all_cells_cutoff": 0.2,
    "fp_count_per_mutation_cross_all_cells_cutoff": 5,
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
# Load raw sequencing data and build initial binary genotype matrix.
data = load_all(inputpath)
P_raw, V_raw, C_raw, A_raw = data["P"], data["V"], data["C"], data["A"]
df_features = data['features']
df_reads_raw = data['df_reads']
I_raw = build_binary_I(P_raw, V_raw, C_raw, params["p_thresh"])

logger.info(f"Loaded data: {len(P_raw)} cells, {len(I_raw.columns)} mutations")


##### Output binary matrix with NA=0
I_raw_withNA0 = I_raw.replace({np.nan: 0}).astype(int)
I_raw_withNA0.to_csv(os.path.join(outputpath, "I_raw_withNA0.txt"), sep="\t")
I_raw_withNA0_T = I_raw_withNA0.T
I_raw_withNA0_T.to_csv(os.path.join(outputpath, "I_raw_withNA0_T.txt"), sep="\t")

I_raw_withNA3 = I_raw.replace({np.nan: 3}).astype(int)
I_raw_withNA3.to_csv(os.path.join(outputpath, "I_raw_withNA3.txt"), sep="\t")
I_raw_withNA3_T = I_raw_withNA3.T
I_raw_withNA3_T.to_csv(os.path.join(outputpath, "I_raw_withNA3_T.txt"), sep="\t")

##### Output for BSCITE bulk input
df_corrected = df_reads_raw.copy()

def split_and_sum(series):
    before_sum = 0
    after_sum = 0
    
    for value in series:
        if pd.notna(value) and '/' in str(value):
            parts = str(value).split('/')
            if len(parts) == 2:
                try:
                    before_sum += int(parts[0])
                    after_sum += int(parts[1])
                except ValueError:
                    continue
                    
    return f"{before_sum}/{after_sum}"

for col in df_corrected.columns:
    df_corrected.loc['pseudo_bulk', col] = split_and_sum(df_corrected[col])

output_bulk_data = []

sample1 = 'pseudo_bulk'
sample2 = 'bulk'

remaining_cells = [idx for idx in df_corrected.index if idx not in ['pseudo_bulk', 'bulk']]
sample3 = deterministic_choice(remaining_cells, salt="scRNA_bulk_sample")

print(f"Using samples: {sample1}, {sample2}, {sample3}")

for i, snp_name in enumerate(df_corrected.columns):
    sample1_value = df_corrected.loc[sample1, snp_name]
    sample2_value = df_corrected.loc[sample2, snp_name] 
    sample3_value = df_corrected.loc[sample3, snp_name]
    
    samples_mutant = []
    samples_reference = []
    
    for value in [sample1_value, sample2_value, sample3_value]:
        if pd.isna(value):
            samples_mutant.append(0)
            samples_reference.append(0)
        elif '/' in str(value):
            mutant, reference = str(value).split('/')
            samples_mutant.append(int(mutant))
            samples_reference.append(int(reference))
        else:
            samples_mutant.append(0)
            samples_reference.append(0)
    
    mutant_counts = f"{samples_mutant[0]};{samples_mutant[1]};{samples_mutant[2]}"
    reference_counts = f"{samples_reference[0]};{samples_reference[1]};{samples_reference[2]}"
    
    chromosome, position = get_random_chromosome_position(snp_name)
    
    output_bulk_data.append({
        'ID': f'mut{i}',
        'Chromosome': chromosome,
        'Position': position,
        'MutantCount': mutant_counts,
        'ReferenceCount': reference_counts
    })

output_bulk_df = pd.DataFrame(output_bulk_data)
output_bulk_df.to_csv(os.path.join(outputpath, "input_BULK.txt"), sep='\t', index=False)


##### Remove cells with no mutations (all zeros)
I_filtered = I_raw[I_raw.eq(1).any(axis=1)]
I = reorder_columns_by_mutant_stats(I_filtered, df_features)[0]
all_mutations = list(I.columns)

cells_in_I = I.index
cols_in_I = I.columns
P = P_raw.loc[cells_in_I, cols_in_I]
V = V_raw.loc[cells_in_I, cols_in_I]
C = C_raw.loc[cells_in_I, cols_in_I]
A = A_raw.loc[cells_in_I, cols_in_I]
bulk_row = df_reads_raw.loc[['bulk']]
df_reads_cells = df_reads_raw.loc[cells_in_I, cols_in_I]
df_reads = pd.concat([bulk_row, df_reads_cells])

df_features_new, empty_mutations = update_features_matrix(I, df_reads, df_features, params["mcf_cutoff"])


# ------------------------------
# Save parameters to JSON
# ------------------------------
import json

params_file = os.path.join(outputpath, "params_used.json")
with open(params_file, 'w') as f:
    json.dump({
        "sampleid": sampleid,
        "inputpath": inputpath,
        "outputpath": outputpath,
        "seed": args.seed,
        "is_predict_germ": is_predict_germ,
        "is_detect_passtree_by_dp": is_detect_passtree_by_dp,
        "is_filter_quality": is_filter_quality,
        "cv_rank_thresh": cv_rank_thresh,
        "remove_artifact_mutations": remove_artifact_mutations,
        "params": params
    }, f, indent=4, default=str)

logger.info(f"Parameters saved to: {params_file}")


# ------------------------------
# Step 2: Classifier for identifying mosaic mutations
# ------------------------------
logger.info("===== Step2: Classifier ...")

models_path = params['models_path']

if features_file is None or features_file == "None" or features_file == "no":
    candidate_mutations = list(I_raw.columns)
else:
    df_for_classifier = pd.read_csv(features_file, sep="\t")
    if 'uniqueid' in df_for_classifier.columns:
        df_for_classifier['mutation_id'] = df_for_classifier['uniqueid'].str.split('.').str[0]
    else:
        df_for_classifier['uniqueid'] = df_for_classifier['identifier']+"."+sampleid
        df_for_classifier['mutation_id'] = df_for_classifier['uniqueid'].str.split('.').str[0]
    
    results_of_classifier = real_time_classifier_predict(df_for_classifier, sampleid, outputpath_classifier)
    
    df_results_of_classifier = results_of_classifier.copy()
    df_candidate_mosaic = df_results_of_classifier.loc[df_results_of_classifier['predicted_label']=='mosaic',:]
    candidate_mutations = [i for i in list(set(list(df_candidate_mosaic['mutation_id']))) if i in list(I_raw.columns)]

logging.info(f"Predicted {len(candidate_mutations)} candidate mosaic mutations")


P_candidate = P[candidate_mutations].copy()
V_candidate = V[candidate_mutations].copy()
A_candidate = A[candidate_mutations].copy()
C_candidate = C[candidate_mutations].copy()
I_candidate = I[candidate_mutations].copy()
df_reads_candidate = df_reads[candidate_mutations].copy()


# ------------------------------
# Step 3. Germline filtering
# ------------------------------
logger.info("===== Step3: Predict germline mutations ...")

if is_predict_germ=="yes":
    logging.info("Running germline filtering ...")
    stats_df, germline_mutations = identify_germline_variants(
        P=P, V=V, C=C, df_reads=df_reads, df_features_new=df_features_new, 
        p_thresh=params["p_thresh"],
        mcf_cutoff=params["mcf_cutoff"],
        mcn_cutoff=params["mcn_cutoff"],
        outputpath=outputpath_germline,
        sampleid=sampleid
    )
else:
    germline_mutations = set()

predicted_germline_mutations = list(germline_mutations)
rescued_germline_mutations = [mut for mut in predicted_germline_mutations 
                             if mut in df_features_new.columns 
                             and df_features_new.loc['mutant_cell_fraction_detected', mut] < 0.5]

removed_germline_mutations = [i for i in predicted_germline_mutations if i not in rescued_germline_mutations]


logging.info(f"Identified {len(predicted_germline_mutations)} germline variants")


plot_heatmap_with_germline_mutations(I, predicted_germline_mutations, os.path.join(outputpath_germline, sampleid+".heatmap_with_predicted_germline_mutations_and_histograms.pdf"))
removed_artifact_mutations = []
somatic_mutations_init = [i for i in all_mutations if i not in removed_germline_mutations and i not in removed_artifact_mutations]
somatic_mutations = list((reorder_columns_by_mutant_stats(I[somatic_mutations_init], df_features_new)[0]).columns)
P_somatic = P[somatic_mutations].copy()
V_somatic = V[somatic_mutations].copy()
A_somatic = A[somatic_mutations].copy()
C_somatic = C[somatic_mutations].copy()
I_somatic = I[somatic_mutations].copy()
df_reads_somatic = df_reads[somatic_mutations].copy()

I_somatic_withNA3 = I_somatic.replace({np.nan: 3}).astype(int)
I_somatic_withNA3.to_csv(os.path.join(outputpath_germline, "I_somatic_withNA3.txt"), sep="\t")
df_features_new = add_mutation_proportions_to_features(df_features_new, I_somatic)


# ------------------------------
# Step 4: Scaffold builder
# ------------------------------
logger.info("===== Step4: Construct scaffold tree ...")

if celltype_file is None or celltype_file == "None":
    barcodes = df_reads_somatic.index.tolist()
    df_celltype = pd.DataFrame({
        "barcode": barcodes,
        "cell_type": ["default_type"] * len(barcodes)
    })
else:
    df_celltype = pd.read_csv(celltype_file, sep="\t")

df_celltype.to_csv(os.path.join(outputpath_scaffold, "df_celltype.txt"), sep="\t")
logger.info(f"Celltype data loaded: {df_celltype.shape[0]} cells")


##### Scaffold builder
logging.info("Running scaffold building ...")
immune_mutations = []
results_of_scaffold = build_scaffold_tree(
    P_somatic = P_somatic, 
    V_somatic = V_somatic, 
    A_somatic = A_somatic, 
    C_somatic = C_somatic, 
    I_somatic = I_somatic,
    df_reads_somatic = df_reads_somatic,
    df_features_new = df_features_new,
    params = params,
    is_filter_quality = is_filter_quality,
    outputpath_scaffold = outputpath_scaffold,
    sampleid = sampleid,
    df_celltype = df_celltype,
    immune_mutations = immune_mutations
)


T_scaffold, M_scaffold, df_flipping_spots, df_total_flipping_count, final_cleaned_I_selected_withNA3, final_cleaned_M_scaffold, backbone_mutations, mutation_group, spots_to_split, group_mutations, remained_mutations, conflict_mutations, high_cv_mutations = results_of_scaffold
scaffold_mutations = list(M_scaffold.columns)
non_scaffold_mutations = [i for i in somatic_mutations if i not in scaffold_mutations]


mutation_clones_for_scaffold = get_mutation_clone_and_backbone_mut_as_keys_by_first_level_with_frequency(T_scaffold, I_somatic)
df_barcode_clones_for_scaffold = assign_clone_labels(M_scaffold, mutation_clones_for_scaffold)

df_barcode_clones_for_scaffold.to_csv(os.path.join(outputpath_scaffold, "df_barcode_clones_from_phylo_tree.csv"), sep=',', index=False)


logging.info(f"Identified {len(scaffold_mutations)} scaffold variants")


logging.info(f"The number of all_mutations is: {len(all_mutations)}")
logging.info(f"The number of predicted_germline_mutations is: {len(predicted_germline_mutations)}")
logging.info(f"The number of somatic_mutations is: {len(somatic_mutations)}")
logging.info(f"The number of scaffold_mutations is: {len(scaffold_mutations)}")
logging.info(f"The number of non_scaffold_mutations is: {len(non_scaffold_mutations)}")


# ------------------------------
# Step 5: DP pass tree & prepare data
# ------------------------------
logger.info("===== Step5: Run DP for passtree ...")

if is_detect_passtree_by_dp=="yes":
    logging.info("Running dynamic programming ...")
    
    pass_tree_cutoff=params['pass_tree_cutoff']
    unpass_tree_cutoff=params['unpass_tree_cutoff']
    
    df_DP_results, passtree_mutations, onecell_mutations = run_dp_pass_tree(
        data = data, 
        df_features_new = df_features_new, 
        M_scaffold = M_scaffold, 
        outputpath_full = outputpath_full, 
        scaffold_mutations = scaffold_mutations,
        p_thresh=p_thresh,
        pass_tree_cutoff=pass_tree_cutoff,
        unpass_tree_cutoff=unpass_tree_cutoff,
        is_log_value_for_likelihoods=True
    )
    uppasstree_mutations = [i for i in all_mutations if i not in passtree_mutations]
else:
    logging.info("Skip running dynamic programming ...")
    passtree_mutations = all_mutations


##### Prepare data for full-resolved tree building

logging.info("Prepare the data for full-resolved tree building ...")

attached_mutations = [i for i in passtree_mutations if i not in scaffold_mutations and i not in removed_germline_mutations and i not in removed_artifact_mutations]
logging.info(f"The number of attached_mutations is: {len(attached_mutations)}")

I_attached_selected = I[scaffold_mutations + attached_mutations]
I_attached_selected_sorted = I_attached_selected[I_attached_selected.apply(lambda col: (col == 1).sum(), axis=0).sort_values(ascending=False).index]
I_attached_sorted_non_empty = I_attached_selected_sorted[I_attached_selected_sorted.eq(1).any(axis=1)]

P_attached_sorted_non_empty = P.loc[
    I_attached_sorted_non_empty.index,
    I_attached_sorted_non_empty.columns
]

I_attached_split, P_attached_split = split_spots_by_immune_mutations(spots_to_split, [i for i in immune_mutations if i in I_attached_sorted_non_empty.columns], I_attached_sorted_non_empty, P_attached_sorted_non_empty)
I_attached, sorting_stats_of_I_attached = reorder_columns_by_mutant_stats(
    I_attached_split, 
    df_features_new,
    min_cell_threshold=30,
    bin_size=5,
    descending=True
)
P_attached = P_attached_split[I_attached.columns]

IRank_mutations = I_attached.columns.tolist()
IRank_mutations_reversed = I_attached.columns[::-1].tolist()
all_conflict_mutations = conflict_mutations.copy()


# ------------------------------
# Step 6: Initial full-resolved tree (Steps 6.1-6.4)
# ------------------------------

# ============================================================================
# This section places accessory mutations (M_accessory) onto the scaffold tree
# (T_scaffold) using the Bayesian placement procedure.
#
# Input:  T_scaffold, M_scaffold, M_accessory
# Output: Initial fully-resolved tree (T_current, M_current)
# ============================================================================

logger.info("=" * 80)
logger.info("===== Step6: Initial full-resolved tree ...")
logger.info("PART I: INITIAL SCAFFOLD INTEGRATION")
logger.info("=" * 80)
logger.info("")

##### Prepare for mutation placement
logger.info("Calculating penalties and refining placements ...")
T_current = copy.deepcopy(T_scaffold)

new_rows_for_current = I_attached.index.difference(M_scaffold.index)
new_data_for_current = pd.DataFrame(0, index=new_rows_for_current, columns=M_scaffold.columns)
M_current_each_mut = pd.concat([M_scaffold, new_data_for_current])
all_nodes_in_T_scaffold = T_scaffold.all_names_no_root()
M_current = merge_mutations(M_current_each_mut, all_nodes_in_T_scaffold)
M_current.insert(0, 'ROOT', 1)

# Penalty parameters
omega_NA = params['general_weight_NA'] if params['general_weight_NA'] else 0.001
fnfp_ratio = params['fnfp_ratio']
phi = params['phi']

root_mutations = []


# ----------------------------------------------------------------------------
# Step 6.1: First-pass integration of accessory mutations
# ----------------------------------------------------------------------------
logger.info("STEP 6.1: First-pass integration of accessory mutations")
logger.info("-" * 80)
logger.info("Integrating accessory mutations (M_accessory) onto the scaffold")
logger.info("tree (T_scaffold) using the Bayesian placement procedure.")
logger.info("")
logger.info("  M_scaffold: scaffold/backbone mutations (retained)")
logger.info("  M_accessory: accessory mutations to be placed")
logger.info("  M_extraneous: mutations with no clear intersection (flagged)")
logger.info("-" * 80)

sorted_attached_mutations = [i for i in I_attached.columns if i in attached_mutations]

if len(sorted_attached_mutations) > 0:
    external_mutations_of_attached_on_scaffold, conflict_mutations_temp, T_current, M_current, root_mutations = attach_mutations_to_current_tree(
        sorted_attached_mutations=sorted_attached_mutations,
        T_current=T_current,
        M_current=M_current,
        I_attached=I_attached,
        P_attached=P_attached,
        ω_NA=omega_NA,
        fnfp_ratio=fnfp_ratio,
        φ=phi,
        logger=logger,
        root_mutations=root_mutations
    )
    all_conflict_mutations.extend(conflict_mutations_temp)
else:
    external_mutations_of_attached_on_scaffold = []

logger.info(f"  |M_scaffold|: {len(scaffold_mutations)}")
logger.info(f"  |M_accessory|: {len(sorted_attached_mutations)}")
logger.info(f"  |M_accessory successfully placed|: {len(sorted_attached_mutations) - len(external_mutations_of_attached_on_scaffold)}")
logger.info(f"  |M_extraneous|: {len(external_mutations_of_attached_on_scaffold)}")
logger.info("")


# ----------------------------------------------------------------------------
# Step 6.2: Second-pass integration of extraneous mutations
# ----------------------------------------------------------------------------
logger.info("STEP 6.2: Second-pass integration of extraneous mutations")
logger.info("-" * 80)
logger.info("Re-attempting integration of extraneous mutations (M_extraneous)")
logger.info("that failed to find clear intersections in the first pass.")
logger.info("-" * 80)

if len(external_mutations_of_attached_on_scaffold) > 0:
    sorted_external_mutations = [i for i in I_attached.columns if i in external_mutations_of_attached_on_scaffold]
    final_external_mutations_of_attached_on_scaffold, conflict_mutations_temp, T_current, M_current, root_mutations = attach_mutations_to_current_tree(
        sorted_attached_mutations=sorted_external_mutations,
        T_current=T_current,
        M_current=M_current,
        I_attached=I_attached,
        P_attached=P_attached,
        ω_NA=omega_NA,
        fnfp_ratio=fnfp_ratio,
        φ=phi,
        logger=logger,
        root_mutations=root_mutations
    )
    all_conflict_mutations.extend(conflict_mutations_temp)
else:
    final_external_mutations_of_attached_on_scaffold = []

logger.info(f"  |M_extraneous remaining after second pass|: {len(final_external_mutations_of_attached_on_scaffold)}")
logger.info("")


# ----------------------------------------------------------------------------
# Step 6.3: Cluster-based integration of remaining extraneous mutations
# ----------------------------------------------------------------------------
logger.info("STEP 6.3: Cluster-based integration of remaining extraneous mutations")
logger.info("-" * 80)
logger.info("Grouping remaining M_extraneous by intersection patterns and")
logger.info("integrating them as coherent phylogenetic clusters.")
logger.info("-" * 80)

external_mutations = final_external_mutations_of_attached_on_scaffold
logging.info(f"The number of external_mutations is: {len(external_mutations)}")

final_external_mutations = []
if len(external_mutations) > 0:
    sorted_external_mutations = [i for i in I_attached.columns if i in external_mutations]
    final_external_mutations, conflict_mutations_temp, T_current, M_current, root_mutations = attach_mutations_to_current_tree(
        sorted_attached_mutations=sorted_external_mutations,
        T_current=T_current,
        M_current=M_current,
        I_attached=I_attached,
        P_attached=P_attached,
        ω_NA=omega_NA,
        fnfp_ratio=fnfp_ratio,
        φ=phi,
        logger=logger,
        root_mutations=root_mutations
    )
    all_conflict_mutations.extend(conflict_mutations_temp)

# Mutations still not placed -> cluster-based integration
logging.info(f"The number of final_external_mutations is: {len(final_external_mutations)}")

remained_mutations = []
if len(final_external_mutations) > 0:
    subtree_groups = cluster_external_mutations_by_intersection(I_attached, final_external_mutations)
    
    logger.info("Processing remaining external mutations by building subtrees")
    
    remained_mutations, conflict_mutations_temp, T_current, M_current, root_mutations = process_external_mutations_by_subtree_groups(
        subtree_groups=subtree_groups,
        T_current=T_current,
        M_current=M_current,
        I_attached=I_attached,
        P_attached=P_attached,
        ω_NA=omega_NA,
        fnfp_ratio=fnfp_ratio,
        φ=phi,
        logger=logger,
        root_mutations=root_mutations
    )
    all_conflict_mutations.extend(conflict_mutations_temp)

logger.info(f"  |M_extraneous resolved via clustering|: {len(final_external_mutations)}")
logger.info(f"  |M_remained|: {len(remained_mutations)}")
logger.info("")


# ----------------------------------------------------------------------------
# Step 6.4: Final integration - root assignment of remained mutations
# ----------------------------------------------------------------------------
logger.info("STEP 6.4: Final integration - root assignment of remained mutations")
logger.info("-" * 80)
logger.info("Performing final integration attempt for M_remained mutations.")
logger.info("Mutations that remain unresolved are assigned to ROOT as M_root.")
logger.info("-" * 80)

final_remained_mutations = []
if len(remained_mutations) > 0:
    sorted_remained_mutations = [i for i in I_attached.columns if i in remained_mutations]
    final_remained_mutations, conflict_mutations_temp, T_current, M_current, root_mutations = attach_mutations_to_current_tree(
        sorted_attached_mutations=sorted_remained_mutations,
        T_current=T_current,
        M_current=M_current,
        I_attached=I_attached,
        P_attached=P_attached,
        ω_NA=omega_NA,
        fnfp_ratio=fnfp_ratio,
        φ=phi,
        logger=logger,
        root_mutations=root_mutations
    )
    all_conflict_mutations.extend(conflict_mutations_temp)

logger.info(f"  |M_root|: {len(root_mutations)}")
logger.info(f"  |M_ambiguous|: {len(all_conflict_mutations)}")
logger.info("")


# Step 6 Summary
# ----------------------------------------------------------------------------
logger.info("=" * 80)
logger.info("Step 6 COMPLETED: Initial fully-resolved tree constructed")
logger.info("-" * 80)

# Print tree structure
logger.info("Current tree structure:")
logger.info("-" * 40)
print_tree_logger(T_current)
logger.info("-" * 40)

# Count mutations (exclude ROOT)
mutations_on_tree = [m for m in M_current.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist() if m != 'ROOT']

logger.info(f"  Mutations on tree (excluding ROOT): {len(mutations_on_tree)}")
logger.info(f"  Cells: {M_current.shape[0]}")
logger.info(f"  Merged columns (including ROOT): {M_current.shape[1]}")
logger.info("=" * 80)
logger.info("")


# ------------------------------
# Step 7: Prune-and-regraft (PRG) refinement (Steps 7.1-7.4)
# ------------------------------

# ============================================================================
# This section identifies mutations with elevated discordance rates using
# four complementary metrics, prunes them from the tree, and re-integrates
# them using the Bayesian placement procedure.
#
# Metrics:
#   7.1: r_FP^(clone) - Intra-clonal false positive discordance rate
#   7.2: r_FP^(global), r_FN^(global) - Global discordance rates
#   7.3: r_FP^(locus) - Locus-specific false positive discordance rate
#   7.4: eta - Ancestral retention fraction
#
# Input:  Initial fully-resolved tree (T_current, M_current)
# Output: Refined tree (T_current, M_current)
# ============================================================================

logger.info("=" * 80)
logger.info("Step 7: Prune-and-regraft (PRG) refinement")
logger.info("=" * 80)
logger.info("")


# ----------------------------------------------------------------------------
# Step 7.1: Intra-clonal FP discordance pruning (M_FP^(clone))
# ----------------------------------------------------------------------------
logger.info("STEP 7.1: Intra-clonal FP discordance pruning")
logger.info("-" * 80)
logger.info("  r_FP^(clone)(j) = delta_FP(j) / |C_j|")
logger.info("  where C_j is the set of cells in the clone where mutation j")
logger.info("  is phylogenetically placed.")
logger.info("")
logger.info("  M_FP^(clone): mutations flagged by intra-clonal FP discordance")
logger.info("  M_FP-daughter^(clone): daughter mutations co-flagged")
logger.info("  G_FP^(clone): clone-associated group containing M_FP^(clone)")
logger.info("-" * 80)

# Recalculate backbone nodes (T_current may have changed since Step 6)
current_backbone_nodes = get_first_level_backbone_nodes(T_current)
expanded_mutations_of_current_backbone_nodes = [mutation for node in current_backbone_nodes for mutation in node.split('|')]

mutation_clones_for_subclone = get_mutation_clone_and_backbone_mut_as_keys_by_first_level_with_frequency(T_current, I_attached)

# ---- Compute intra-clonal FP discordance rates ----
T_checkpoint_fpratio_within_subclone = copy.deepcopy(T_current)
M_checkpoint_fpratio_within_subclone = M_current.copy()

M_for_fp_ratio_and_fn_ratio_fpratio_within_subclone = M_checkpoint_fpratio_within_subclone.drop(columns=['ROOT'], errors='ignore')
mutations_on_T_current_fpratio_within_subclone = M_for_fp_ratio_and_fn_ratio_fpratio_within_subclone.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()
M_for_fp_ratio_and_fn_ratio_fpratio_within_subclone = split_merged_columns(M_for_fp_ratio_and_fn_ratio_fpratio_within_subclone, mutations_on_T_current_fpratio_within_subclone)

df_fp_ratio_fpratio_within_subclone, fp_mutations_dict_for_out_subclone_muts_fpratio_within_subclone, fp_mutations_dict_for_in_subclone_muts_fpratio_within_subclone = calculate_fp_ratios_within_subclone(
    M_for_fp_ratio_and_fn_ratio_fpratio_within_subclone, I_attached, mutation_clones_for_subclone
)

fp_ratio_cutoff_within_subclone = params['fp_ratio_cutoff_within_subclone']

# ---- Identify M_FP^(clone) ----
# Mutations with intra-clonal FP discordance rate >= threshold
rehanged_mutations_by_fpratio_within_subclone = df_fp_ratio_fpratio_within_subclone[
    df_fp_ratio_fpratio_within_subclone['fp_ratio_within_subclone_for_in_subclone_muts'] >= fp_ratio_cutoff_within_subclone
]['identifier'].tolist()

# Exclude backbone mutations from pruning
rehanged_mutations_by_fpratio_within_subclone_but_backbone = [
    i for i in rehanged_mutations_by_fpratio_within_subclone 
    if i not in list(set(expanded_mutations_of_current_backbone_nodes + scaffold_mutations))
]

logger.info(f"  |M_FP^(clone)|: {len(rehanged_mutations_by_fpratio_within_subclone)}")
logger.info(f"  |M_FP^(clone) \ M_backbone|: {len(rehanged_mutations_by_fpratio_within_subclone_but_backbone)}")

# ---- Identify M_FP-daughter^(clone) ----
# Get ordered branch groups for FP-flagged mutations
ordered_branch_groups_for_rehanged_mutations_by_fpratio_within_subclone_but_backbone = find_ordered_branch_groups_for_rehanged_mutations_with_keys_as_earlist(
    T_current, rehanged_mutations_by_fpratio_within_subclone_but_backbone
)

# Map each flagged mutation to its FP-inducing mutations
filtered_fp_mutations_dict_by_fpratio_within_subclone = {
    mut: other_muts 
    for mut, other_muts in fp_mutations_dict_for_in_subclone_muts_fpratio_within_subclone.items() 
    if mut in rehanged_mutations_by_fpratio_within_subclone_but_backbone
}

# Map FP-flagged mutations to their node names in the tree
nodes_rehanged_mutations_by_fpratio_within_subclone_but_backbone = list(set([
    find_mutation_column(mutation, M_current.columns) 
    for mutation in rehanged_mutations_by_fpratio_within_subclone_but_backbone
]))

# Get all daughter mutations for each FP-flagged node
node_dict = {node.name: node for node in T_current.traverse()}
daughters_nodes_dict_by_fpratio_within_subclone = {
    mutation: get_all_daughter_mutations(node_dict[mutation])
    for mutation in nodes_rehanged_mutations_by_fpratio_within_subclone_but_backbone
    if mutation in node_dict
}
daughters_mutations_dict_by_fpratio_within_subclone = list(set([
    daughter for daughters_list in daughters_nodes_dict_by_fpratio_within_subclone.values() 
    for daughter in daughters_list
    if daughter not in rehanged_mutations_by_fpratio_within_subclone_but_backbone
]))

logger.info(f"  |M_FP-daughter^(clone)|: {len(daughters_mutations_dict_by_fpratio_within_subclone)}")
logger.info(f"  |G_FP^(clone)|: {len(ordered_branch_groups_for_rehanged_mutations_by_fpratio_within_subclone_but_backbone)} groups")

# ---- Pool mutations for PRG ----
daughters_to_leaf_mutations_fpratio_within_subclone = []
fp_mutations_fpratio_within_subclone = []

for branch_mut in ordered_branch_groups_for_rehanged_mutations_by_fpratio_within_subclone_but_backbone:
    daughter_list = ordered_branch_groups_for_rehanged_mutations_by_fpratio_within_subclone_but_backbone[branch_mut]
    
    daughters_to_leaf_mutations_fpratio_within_subclone = list({
        item 
        for key in daughter_list
        if key in daughters_mutations_dict_by_fpratio_within_subclone
        for item in [key] + daughters_mutations_dict_by_fpratio_within_subclone[key]
    })
    
    fp_mutations_fpratio_within_subclone_init = list({
        item 
        for key in daughter_list
        if key in filtered_fp_mutations_dict_by_fpratio_within_subclone
        for item in [key] + filtered_fp_mutations_dict_by_fpratio_within_subclone[key]
    })
    fp_mutations_fpratio_within_subclone = [
        i for i in fp_mutations_fpratio_within_subclone_init 
        if i not in daughter_list 
        and i not in daughters_to_leaf_mutations_fpratio_within_subclone 
        and i not in expanded_mutations_of_current_backbone_nodes
    ]

# QC filter: require mutant cell number > 5
daughters_to_leaf_mutations_fpratio_within_subclone_qc = [
    mut for mut in daughters_to_leaf_mutations_fpratio_within_subclone 
    if df_features_new[mut]['mutant_cellnum'] > 5
]
fp_mutations_fpratio_within_subclone_qc = set(
    fp_mutations_fpratio_within_subclone + 
    [i for i in daughters_to_leaf_mutations_fpratio_within_subclone 
     if i not in daughters_to_leaf_mutations_fpratio_within_subclone_qc]
)

# Sort mutations for deterministic processing
sorted_fp_mutations_fpratio_within_subclone = [
    i for i in I_attached.columns 
    if i in fp_mutations_fpratio_within_subclone_qc
]
sorted_daughters_to_leaf_mutations_fpratio_within_subclone = [
    i for i in I_attached.columns 
    if i in daughters_to_leaf_mutations_fpratio_within_subclone_qc
]
sorted_rehanged_mutations_all_fpratio_within_subclone = (
    sorted_fp_mutations_fpratio_within_subclone + 
    sorted_daughters_to_leaf_mutations_fpratio_within_subclone
)

# ---- Prune and regraft M_FP^(clone) and M_FP-daughter^(clone) ----
external_mutations_fpratio_within_subclone_by_sorted_fp_mutations_fpratio_within_subclone = []
external_mutations_fpratio_within_subclone_by_sorted_daughters_to_leaf_mutations_fpratio_within_subclone = []

if len(sorted_rehanged_mutations_all_fpratio_within_subclone) > 0:
    T_removed_fpratio_within_subclone, M_removed_fpratio_within_subclone = remove_mutations_from_tree_and_matrix(
        T_checkpoint_fpratio_within_subclone, M_checkpoint_fpratio_within_subclone, 
        sorted_rehanged_mutations_all_fpratio_within_subclone
    )
    logger.info(f"  Tree after pruning: {M_removed_fpratio_within_subclone.shape[0]} cells, {M_removed_fpratio_within_subclone.shape[1]} mutations")
    
    T_current = copy.deepcopy(T_removed_fpratio_within_subclone)
    M_current = M_removed_fpratio_within_subclone.copy()

if len(sorted_fp_mutations_fpratio_within_subclone) > 0:
    external_mutations_fpratio_within_subclone_by_sorted_fp_mutations_fpratio_within_subclone, conflict_mutations_temp, T_current, M_current, root_mutations = attach_mutations_to_current_tree(
        sorted_attached_mutations=sorted_fp_mutations_fpratio_within_subclone,
        T_current=T_current,
        M_current=M_current,
        I_attached=I_attached,
        P_attached=P_attached,
        ω_NA=omega_NA,
        fnfp_ratio=fnfp_ratio,
        φ=phi,
        logger=logger,
        root_mutations=root_mutations
    )
    all_conflict_mutations.extend(conflict_mutations_temp)

if len(sorted_daughters_to_leaf_mutations_fpratio_within_subclone) > 0:
    external_mutations_fpratio_within_subclone_by_sorted_daughters_to_leaf_mutations_fpratio_within_subclone, conflict_mutations_temp, T_current, M_current, root_mutations = attach_mutations_to_current_tree(
        sorted_attached_mutations=sorted_daughters_to_leaf_mutations_fpratio_within_subclone,
        T_current=T_current,
        M_current=M_current,
        I_attached=I_attached,
        P_attached=P_attached,
        ω_NA=omega_NA,
        fnfp_ratio=fnfp_ratio,
        φ=phi,
        logger=logger,
        root_mutations=root_mutations
    )
    all_conflict_mutations.extend(conflict_mutations_temp)

# ---- Re-attempt any remaining external mutations ----
external_mutations_fpratio_within_subclone = list(set(
    external_mutations_fpratio_within_subclone_by_sorted_fp_mutations_fpratio_within_subclone + 
    external_mutations_fpratio_within_subclone_by_sorted_daughters_to_leaf_mutations_fpratio_within_subclone
))

if len(external_mutations_fpratio_within_subclone) > 0:
    logger.info(f"  Re-attempting {len(external_mutations_fpratio_within_subclone)} external mutations...")
    sorted_external_mutations_fpratio_within_subclone = [
        i for i in I_attached.columns 
        if i in external_mutations_fpratio_within_subclone
    ]
    final_external_mutations_fpratio_within_subclone, conflict_mutations_temp, T_current, M_current, root_mutations = attach_mutations_to_current_tree(
        sorted_attached_mutations=sorted_external_mutations_fpratio_within_subclone,
        T_current=T_current,
        M_current=M_current,
        I_attached=I_attached,
        P_attached=P_attached,
        ω_NA=omega_NA,
        fnfp_ratio=fnfp_ratio,
        φ=phi,
        logger=logger,
        root_mutations=root_mutations
    )
    all_conflict_mutations.extend(conflict_mutations_temp)
    logger.info(f"  Final external mutations: {len(final_external_mutations_fpratio_within_subclone)}")

logger.info("")


# ---- Recompute metrics after pruning ----
mutation_clones_for_subclone_v2 = get_mutation_clone_and_backbone_mut_as_keys_by_first_level_with_frequency(T_current, I_attached)

T_checkpoint_fpratio_within_subclone_v2 = copy.deepcopy(T_current)
M_checkpoint_fpratio_within_subclone_v2 = M_current.copy()

M_for_fp_ratio_and_fn_ratio_fpratio_within_subclone_v2 = M_checkpoint_fpratio_within_subclone_v2.drop(columns=['ROOT'], errors='ignore')
mutations_on_T_current_fpratio_within_subclone_v2 = M_for_fp_ratio_and_fn_ratio_fpratio_within_subclone_v2.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()
M_for_fp_ratio_and_fn_ratio_fpratio_within_subclone_v2 = split_merged_columns(M_for_fp_ratio_and_fn_ratio_fpratio_within_subclone_v2, mutations_on_T_current_fpratio_within_subclone_v2)

df_fp_ratio_fpratio_within_subclone_v2, fp_mutations_dict_for_out_subclone_muts_fpratio_within_subclone_v2, fp_mutations_dict_for_in_subclone_muts_fpratio_within_subclone_v2 = calculate_fp_ratios_within_subclone(
    M_for_fp_ratio_and_fn_ratio_fpratio_within_subclone_v2, I_attached, mutation_clones_for_subclone_v2
)

# Merge pre- and post-pruning metrics for comparison
df_fp_ratio_fpratio_within_subclone_final = pd.merge(
    df_fp_ratio_fpratio_within_subclone, 
    df_fp_ratio_fpratio_within_subclone_v2, 
    on='identifier', 
    suffixes=('.1', '.2')
)


# ----------------------------------------------------------------------------
# Step 7.2: Global FP/FN discordance pruning (M_FP^(global), M_FN^(global))
# ----------------------------------------------------------------------------
logger.info("STEP 7.2: Global FP/FN discordance pruning")
logger.info("-" * 80)
logger.info("  r_FP^(global)(j) = delta_FP(j) / N_total")
logger.info("  r_FN^(global)(j) = delta_FN(j) / N_total")
logger.info("")
logger.info("  M_FP^(global): mutations with r_FP^(global) >= threshold")
logger.info("  M_FN^(global): mutations with r_FN^(global) >= threshold")
logger.info("  M_FP-daughter^(global): daughter mutations co-flagged")
logger.info("  G_FP^(global): global FP-associated group")
logger.info("  G_FN^(global): global FN-associated group")
logger.info("-" * 80)

# Recalculate backbone nodes (T_current may have changed in Step 7.1)
current_backbone_nodes = get_first_level_backbone_nodes(T_current)
expanded_mutations_of_current_backbone_nodes = [mutation for node in current_backbone_nodes for mutation in node.split('|')]

# ---- Compute global FP/FN discordance rates ----
T_checkpoint_fpfnratio_across_tree = copy.deepcopy(T_current)
M_checkpoint_fpfnratio_across_tree = M_current.copy()

M_for_fp_ratio_and_fn_ratio_fpfnratio_across_tree = M_checkpoint_fpfnratio_across_tree.drop(columns=['ROOT'], errors='ignore')
mutations_on_T_current_fpfnratio_across_tree = M_for_fp_ratio_and_fn_ratio_fpfnratio_across_tree.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()
M_for_fp_ratio_and_fn_ratio_fpfnratio_across_tree = split_merged_columns(M_for_fp_ratio_and_fn_ratio_fpfnratio_across_tree, mutations_on_T_current_fpfnratio_across_tree)

df_fp_ratio_and_fn_ratio_fpfnratio_across_tree, fp_mutations_dict_fpfnratio_across_tree = calculate_fp_fn_ratios_across_tree(
    M_for_fp_ratio_and_fn_ratio_fpfnratio_across_tree, I_attached
)

fp_ratio_cutoff_across_tree = params['fp_ratio_cutoff_across_tree']
fn_ratio_cutoff_across_tree = params['fn_ratio_cutoff_across_tree']

# ---- Identify M_FP^(global) ----
# Mutations with global FP discordance rate >= threshold
rehanged_fp_mutations_by_fpfnratio_across_tree = df_fp_ratio_and_fn_ratio_fpfnratio_across_tree[
    df_fp_ratio_and_fn_ratio_fpfnratio_across_tree['fp_ratio'] >= fp_ratio_cutoff_across_tree
]['identifier'].tolist()

# Exclude backbone mutations from pruning
rehanged_fp_mutations_by_fpfnratio_across_tree_but_backbone = [
    i for i in rehanged_fp_mutations_by_fpfnratio_across_tree 
    if i not in list(set(expanded_mutations_of_current_backbone_nodes + scaffold_mutations))
]

logger.info(f"  |M_FP^(global)|: {len(rehanged_fp_mutations_by_fpfnratio_across_tree)}")
logger.info(f"  |M_FP^(global) \ M_backbone|: {len(rehanged_fp_mutations_by_fpfnratio_across_tree_but_backbone)}")

# ---- Identify G_FP^(global) ----
# Get ordered branch groups for FP-flagged mutations
ordered_branch_groups_for_rehanged_fp_mutations_by_fpfnratio_across_tree_but_backbone = find_ordered_branch_groups_for_rehanged_mutations_with_keys_as_earlist(
    T_current, rehanged_fp_mutations_by_fpfnratio_across_tree_but_backbone
)

# Map each FP-flagged mutation to its FP-inducing mutations
filtered_fp_mutations_dict_by_fpfnratio_across_tree = {
    mut: other_muts 
    for mut, other_muts in fp_mutations_dict_fpfnratio_across_tree.items() 
    if mut in rehanged_fp_mutations_by_fpfnratio_across_tree_but_backbone
}

# Map FP-flagged mutations to their node names in the tree
nodes_rehanged_fp_mutations_by_fpfnratio_across_tree_but_backbone = list(set([
    find_mutation_column(mutation, M_current.columns) 
    for mutation in rehanged_fp_mutations_by_fpfnratio_across_tree_but_backbone
]))

# Get all daughter mutations for each FP-flagged node
node_dict = {node.name: node for node in T_current.traverse()}
daughters_nodes_dict_by_fpfnratio_across_tree = {
    mutation: get_all_daughter_mutations(node_dict[mutation])
    for mutation in nodes_rehanged_fp_mutations_by_fpfnratio_across_tree_but_backbone
    if mutation in node_dict
}

# ---- Identify M_FP-daughter^(global) ----
daughters_mutations_dict_by_fpfnratio_across_tree = list(set([
    daughter for daughters_list in daughters_nodes_dict_by_fpfnratio_across_tree.values() 
    for daughter in daughters_list
    if daughter not in rehanged_fp_mutations_by_fpfnratio_across_tree_but_backbone
]))

logger.info(f"  |M_FP-daughter^(global)|: {len(daughters_mutations_dict_by_fpfnratio_across_tree)}")
logger.info(f"  |G_FP^(global)|: {len(ordered_branch_groups_for_rehanged_fp_mutations_by_fpfnratio_across_tree_but_backbone)} groups")

# ---- Pool mutations for PRG ----
daughters_to_leaf_mutations_fpfnratio_across_tree = []
fp_mutations_fpfnratio_across_tree = []

for branch_mut in ordered_branch_groups_for_rehanged_fp_mutations_by_fpfnratio_across_tree_but_backbone:
    daughter_list = ordered_branch_groups_for_rehanged_fp_mutations_by_fpfnratio_across_tree_but_backbone[branch_mut]
    
    daughters_to_leaf_mutations_fpfnratio_across_tree = list({
        item 
        for key in daughter_list
        if key in daughters_mutations_dict_by_fpfnratio_across_tree
        for item in [key] + daughters_mutations_dict_by_fpfnratio_across_tree[key]
    })
    
    fp_mutations_fpfnratio_across_tree_init = list({
        item 
        for key in daughter_list
        if key in filtered_fp_mutations_dict_by_fpfnratio_across_tree
        for item in [key] + filtered_fp_mutations_dict_by_fpfnratio_across_tree[key]
    })
    fp_mutations_fpfnratio_across_tree = [
        i for i in fp_mutations_fpfnratio_across_tree_init 
        if i not in daughter_list 
        and i not in daughters_to_leaf_mutations_fpfnratio_across_tree 
        and i not in list(set(expanded_mutations_of_current_backbone_nodes + scaffold_mutations))
    ]

# QC filter: require mutant cell number > 5
daughters_to_leaf_mutations_fpfnratio_across_tree_qc = [
    mut for mut in daughters_to_leaf_mutations_fpfnratio_across_tree 
    if df_features_new[mut]['mutant_cellnum'] > 5
]
fp_mutations_fpfnratio_across_tree_qc = set(
    fp_mutations_fpfnratio_across_tree + 
    [i for i in daughters_to_leaf_mutations_fpfnratio_across_tree 
     if i not in daughters_to_leaf_mutations_fpfnratio_across_tree_qc]
)

# Sort mutations for deterministic processing
sorted_fp_mutations_fpfnratio_across_tree = [
    i for i in I_attached.columns 
    if i in fp_mutations_fpfnratio_across_tree_qc
]
sorted_daughters_to_leaf_mutations_fpfnratio_across_tree = [
    i for i in I_attached.columns 
    if i in daughters_to_leaf_mutations_fpfnratio_across_tree_qc
]
sorted_rehanged_mutations_all_fpfnratio_across_tree = (
    sorted_fp_mutations_fpfnratio_across_tree + 
    sorted_daughters_to_leaf_mutations_fpfnratio_across_tree
)

# ---- Prune and regraft M_FP^(global) and M_FP-daughter^(global) ----
external_mutations_fpfnratio_across_tree_by_sorted_fp_mutations_fpfnratio_across_tree = []
external_mutations_fpfnratio_across_tree_by_sorted_daughters_to_leaf_mutations_fpfnratio_across_tree = []

if len(sorted_rehanged_mutations_all_fpfnratio_across_tree) > 0:
    T_removed_fpfnratio_across_tree, M_removed_fpfnratio_across_tree = remove_mutations_from_tree_and_matrix(
        T_checkpoint_fpfnratio_across_tree, M_checkpoint_fpfnratio_across_tree, 
        sorted_rehanged_mutations_all_fpfnratio_across_tree
    )
    logger.info(f"  Tree after pruning: {M_removed_fpfnratio_across_tree.shape[0]} cells, {M_removed_fpfnratio_across_tree.shape[1]} mutations")
    
    T_current = copy.deepcopy(T_removed_fpfnratio_across_tree)
    M_current = M_removed_fpfnratio_across_tree.copy()

if len(sorted_fp_mutations_fpfnratio_across_tree) > 0:
    external_mutations_fpfnratio_across_tree_by_sorted_fp_mutations_fpfnratio_across_tree, conflict_mutations_temp, T_current, M_current, root_mutations = attach_mutations_to_current_tree(
        sorted_attached_mutations=sorted_fp_mutations_fpfnratio_across_tree,
        T_current=T_current,
        M_current=M_current,
        I_attached=I_attached,
        P_attached=P_attached,
        ω_NA=omega_NA,
        fnfp_ratio=fnfp_ratio,
        φ=phi,
        logger=logger,
        root_mutations=root_mutations
    )
    all_conflict_mutations.extend(conflict_mutations_temp)

if len(sorted_daughters_to_leaf_mutations_fpfnratio_across_tree) > 0:
    external_mutations_fpfnratio_across_tree_by_sorted_daughters_to_leaf_mutations_fpfnratio_across_tree, conflict_mutations_temp, T_current, M_current, root_mutations = attach_mutations_to_current_tree(
        sorted_attached_mutations=sorted_daughters_to_leaf_mutations_fpfnratio_across_tree,
        T_current=T_current,
        M_current=M_current,
        I_attached=I_attached,
        P_attached=P_attached,
        ω_NA=omega_NA,
        fnfp_ratio=fnfp_ratio,
        φ=phi,
        logger=logger,
        root_mutations=root_mutations
    )
    all_conflict_mutations.extend(conflict_mutations_temp)

# ---- Identify M_FN^(global) ----
# Mutations with global FN discordance rate >= threshold
rehanged_fn_mutations_by_fpfnratio_across_tree = df_fp_ratio_and_fn_ratio_fpfnratio_across_tree[
    df_fp_ratio_and_fn_ratio_fpfnratio_across_tree['fn_ratio'] >= fn_ratio_cutoff_across_tree
]['identifier'].tolist()

sorted_fn_mutations_fpfnratio_across_tree = [
    i for i in I_attached.columns 
    if i in rehanged_fn_mutations_by_fpfnratio_across_tree
]

logger.info(f"  |M_FN^(global)|: {len(rehanged_fn_mutations_by_fpfnratio_across_tree)}")

# ---- Prune and regraft M_FN^(global) ----
external_mutations_fpfnratio_across_tree_by_sorted_fn_mutations_fpfnratio_across_tree = []

if len(rehanged_fn_mutations_by_fpfnratio_across_tree) > 0:
    T_removed_fpfnratio_across_tree, M_removed_fpfnratio_across_tree = remove_mutations_from_tree_and_matrix(
        T_checkpoint_fpfnratio_across_tree, M_checkpoint_fpfnratio_across_tree, 
        rehanged_fn_mutations_by_fpfnratio_across_tree
    )
    logger.info(f"  Tree after pruning: {M_removed_fpfnratio_across_tree.shape[0]} cells, {M_removed_fpfnratio_across_tree.shape[1]} mutations")
    
    T_current = copy.deepcopy(T_removed_fpfnratio_across_tree)
    M_current = M_removed_fpfnratio_across_tree.copy()

if len(sorted_fn_mutations_fpfnratio_across_tree) > 0:
    external_mutations_fpfnratio_across_tree_by_sorted_fn_mutations_fpfnratio_across_tree, conflict_mutations_temp, T_current, M_current, root_mutations = attach_mutations_to_current_tree(
        sorted_attached_mutations=sorted_fn_mutations_fpfnratio_across_tree,
        T_current=T_current,
        M_current=M_current,
        I_attached=I_attached,
        P_attached=P_attached,
        ω_NA=omega_NA,
        fnfp_ratio=fnfp_ratio,
        φ=phi,
        logger=logger,
        root_mutations=root_mutations
    )
    all_conflict_mutations.extend(conflict_mutations_temp)

# ---- Re-attempt any remaining external mutations ----
external_mutations_fpfnratio_across_tree = list(set(
    external_mutations_fpfnratio_across_tree_by_sorted_fp_mutations_fpfnratio_across_tree +
    external_mutations_fpfnratio_across_tree_by_sorted_daughters_to_leaf_mutations_fpfnratio_across_tree +
    external_mutations_fpfnratio_across_tree_by_sorted_fn_mutations_fpfnratio_across_tree
))

if len(external_mutations_fpfnratio_across_tree) > 0:
    logger.info(f"  Re-attempting {len(external_mutations_fpfnratio_across_tree)} external mutations...")
    sorted_external_mutations_fpfnratio_across_tree = [
        i for i in I_attached.columns 
        if i in external_mutations_fpfnratio_across_tree
    ]
    final_external_mutations_fpfnratio_across_tree, conflict_mutations_temp, T_current, M_current, root_mutations = attach_mutations_to_current_tree(
        sorted_attached_mutations=sorted_external_mutations_fpfnratio_across_tree,
        T_current=T_current,
        M_current=M_current,
        I_attached=I_attached,
        P_attached=P_attached,
        ω_NA=omega_NA,
        fnfp_ratio=fnfp_ratio,
        φ=phi,
        logger=logger,
        root_mutations=root_mutations
    )
    all_conflict_mutations.extend(conflict_mutations_temp)
    logger.info(f"  Final external mutations: {len(final_external_mutations_fpfnratio_across_tree)}")

logger.info("")


# ---- Recompute metrics after pruning ----
T_checkpoint_fpfnratio_across_tree_v2 = copy.deepcopy(T_current)
M_checkpoint_fpfnratio_across_tree_v2 = M_current.copy()

M_for_fp_ratio_and_fn_ratio_fpfnratio_across_tree_v2 = M_checkpoint_fpfnratio_across_tree_v2.drop(columns=['ROOT'], errors='ignore')
mutations_on_T_current_fpfnratio_across_tree_v2 = M_for_fp_ratio_and_fn_ratio_fpfnratio_across_tree_v2.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()
M_for_fp_ratio_and_fn_ratio_fpfnratio_across_tree_v2 = split_merged_columns(M_for_fp_ratio_and_fn_ratio_fpfnratio_across_tree_v2, mutations_on_T_current_fpfnratio_across_tree_v2)

df_fp_ratio_and_fn_ratio_fpfnratio_across_tree_v2, fp_mutations_dict_fpfnratio_across_tree_v2 = calculate_fp_fn_ratios_across_tree(
    M_for_fp_ratio_and_fn_ratio_fpfnratio_across_tree_v2, I_attached
)

# Merge pre- and post-pruning metrics for comparison
df_fp_ratio_and_fn_ratio_fpfnratio_across_tree_final = pd.merge(
    df_fp_ratio_and_fn_ratio_fpfnratio_across_tree, 
    df_fp_ratio_and_fn_ratio_fpfnratio_across_tree_v2, 
    on='identifier', 
    suffixes=('.1', '.2')
)

# Merge with Step 7.1 results
combined_df_fp_ratios_within_subclone_and_fpfn_ratios_across_tree = pd.merge(
    df_fp_ratio_fpratio_within_subclone_final,
    df_fp_ratio_and_fn_ratio_fpfnratio_across_tree_final,
    on='identifier',
    how='left'
)


# ----------------------------------------------------------------------------
# Step 7.3: Locus-specific FP discordance pruning (M_FP^(locus))
# ----------------------------------------------------------------------------
logger.info("STEP 7.3: Locus-specific FP discordance pruning")
logger.info("-" * 80)
logger.info("  r_FP^(locus)(j) = delta_FP(j) / coverage_j")
logger.info("")
logger.info("  M_FP^(locus): mutations flagged by locus-specific FP discordance")
logger.info("  - These mutations exhibit site-specific bias")
logger.info("  - May indicate sequencing errors or mapping artifacts")
logger.info("-" * 80)

# Recalculate backbone nodes (T_current may have changed in Step 7.2)
current_backbone_nodes = get_first_level_backbone_nodes(T_current)
expanded_mutations_of_current_backbone_nodes = [mutation for node in current_backbone_nodes for mutation in node.split('|')]

mutation_clones_for_persitefp = get_mutation_clone_and_backbone_mut_as_keys_by_first_level_with_frequency(T_current, I_attached)

# ---- Compute locus-specific FP discordance rates ----
T_checkpoint_fp_ratio_persitefp = copy.deepcopy(T_current)
M_checkpoint_fp_ratio_persitefp = M_current.copy()

M_for_fp_ratio_persitefp = M_checkpoint_fp_ratio_persitefp.drop(columns=['ROOT'], errors='ignore')
mutations_on_T_current_persitefp = M_for_fp_ratio_persitefp.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()
M_for_fp_ratio_persitefp = split_merged_columns(M_for_fp_ratio_persitefp, mutations_on_T_current_persitefp)

df_fp_ratio_persitefp = calculate_fp_ratios_persite_within_subclone(M_for_fp_ratio_persitefp, I_attached, mutation_clones_for_persitefp)

fp_ratio_persite_cutoff = params['fp_ratio_persite_cutoff']

# ---- Identify M_FP^(locus) ----
# Mutations with locus-specific FP discordance rate >= threshold
rehanged_mutations_by_persitefp = df_fp_ratio_persitefp[
    df_fp_ratio_persitefp['fp_ratio_persite'] >= fp_ratio_persite_cutoff
]['identifier'].tolist()

# Exclude backbone mutations from pruning
rehanged_mutations_by_persitefp_but_backbone = [
    i for i in rehanged_mutations_by_persitefp 
    if i not in list(set(expanded_mutations_of_current_backbone_nodes + scaffold_mutations))
]

logger.info(f"  |M_FP^(locus)|: {len(rehanged_mutations_by_persitefp)}")
logger.info(f"  |M_FP^(locus) \ M_backbone|: {len(rehanged_mutations_by_persitefp_but_backbone)}")

sorted_rehanged_mutations_by_persitefp_but_backbone = [
    i for i in I_attached.columns 
    if i in rehanged_mutations_by_persitefp_but_backbone
]

# ---- Prune and regraft M_FP^(locus) ----
external_mutations_by_sorted_rehanged_mutations_by_persitefp_but_backbone = []

if len(sorted_rehanged_mutations_by_persitefp_but_backbone) > 0:
    T_removed_fp_ratio_persitefp, M_removed_fp_ratio_persitefp = remove_mutations_from_tree_and_matrix(
        T_checkpoint_fp_ratio_persitefp, M_checkpoint_fp_ratio_persitefp, 
        sorted_rehanged_mutations_by_persitefp_but_backbone
    )
    logger.info(f"  Tree after pruning: {M_removed_fp_ratio_persitefp.shape[0]} cells, {M_removed_fp_ratio_persitefp.shape[1]} mutations")
    
    T_current = copy.deepcopy(T_removed_fp_ratio_persitefp)
    M_current = M_removed_fp_ratio_persitefp.copy()

if len(sorted_rehanged_mutations_by_persitefp_but_backbone) > 0:
    external_mutations_by_sorted_rehanged_mutations_by_persitefp_but_backbone, conflict_mutations_temp, T_current, M_current, root_mutations = attach_mutations_to_current_tree(
        sorted_attached_mutations=sorted_rehanged_mutations_by_persitefp_but_backbone,
        T_current=T_current,
        M_current=M_current,
        I_attached=I_attached,
        P_attached=P_attached,
        ω_NA=omega_NA,
        fnfp_ratio=fnfp_ratio,
        φ=phi,
        logger=logger,
        root_mutations=root_mutations
    )
    all_conflict_mutations.extend(conflict_mutations_temp)

# ---- Re-attempt any remaining external mutations ----
if len(external_mutations_by_sorted_rehanged_mutations_by_persitefp_but_backbone) > 0:
    logger.info(f"  Re-attempting {len(external_mutations_by_sorted_rehanged_mutations_by_persitefp_but_backbone)} external mutations...")
    sorted_external_mutations_by_sorted_rehanged_mutations_by_persitefp_but_backbone = [
        i for i in I_attached.columns 
        if i in external_mutations_by_sorted_rehanged_mutations_by_persitefp_but_backbone
    ]
    final_external_mutations_fp_ratio_persitefp, conflict_mutations_temp, T_current, M_current, root_mutations = attach_mutations_to_current_tree(
        sorted_attached_mutations=sorted_external_mutations_by_sorted_rehanged_mutations_by_persitefp_but_backbone,
        T_current=T_current,
        M_current=M_current,
        I_attached=I_attached,
        P_attached=P_attached,
        ω_NA=omega_NA,
        fnfp_ratio=fnfp_ratio,
        φ=phi,
        logger=logger,
        root_mutations=root_mutations
    )
    all_conflict_mutations.extend(conflict_mutations_temp)
    logger.info(f"  Final external mutations: {len(final_external_mutations_fp_ratio_persitefp)}")

logger.info("")


# ---- Recompute metrics after pruning ----
mutation_clones_for_persitefp_v2 = get_mutation_clone_and_backbone_mut_as_keys_by_first_level_with_frequency(T_current, I_attached)

T_checkpoint_fp_ratio_persitefp_v2 = copy.deepcopy(T_current)
M_checkpoint_fp_ratio_persitefp_v2 = M_current.copy()

M_for_fp_ratio_persitefp_v2 = M_checkpoint_fp_ratio_persitefp_v2.drop(columns=['ROOT'], errors='ignore')
mutations_on_T_current_persitefp = M_for_fp_ratio_persitefp_v2.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()
M_for_fp_ratio_persitefp_v2 = split_merged_columns(M_for_fp_ratio_persitefp_v2, mutations_on_T_current_persitefp)

df_fp_ratio_persitefp_v2 = calculate_fp_ratios_persite_within_subclone(M_for_fp_ratio_persitefp_v2, I_attached, mutation_clones_for_persitefp_v2)

# Merge pre- and post-pruning metrics for comparison
df_fp_ratio_persitefp_final = pd.merge(
    df_fp_ratio_persitefp, 
    df_fp_ratio_persitefp_v2, 
    on='identifier', 
    suffixes=('.1', '.2')
)

# Merge with Step 7.1 and 7.2 results
final_combined_df_fp_ratios_within_subclone_and_fpfn_ratios_across_tree_and_persite_fp_ratio = pd.merge(
    combined_df_fp_ratios_within_subclone_and_fpfn_ratios_across_tree,
    df_fp_ratio_persitefp_final,
    on='identifier',
    how='left'
)

final_combined_df_fp_ratios_within_subclone_and_fpfn_ratios_across_tree_and_persite_fp_ratio.to_csv(
    os.path.join(outputpath_full, "final_combined_df_fp_ratios_within_subclone_and_fpfn_ratios_across_tree_and_persite_fp_ratio.csv"), 
    sep=","
)


# ----------------------------------------------------------------------------
# Step 7.4: Ancestral retention-based orphaned mutation identification
# ----------------------------------------------------------------------------
logger.info("STEP 7.4: Ancestral retention-based orphaned mutation identification")
logger.info("-" * 80)
logger.info("  eta(j, p(j)) = N_intersection / N_mutant(j)")
logger.info("")
logger.info("  M_orphan^(retention): mutations with eta < threshold")
logger.info("  M_orphan-progeny: daughter mutations of M_orphan^(retention)")
logger.info("")
logger.info("These mutations lack sufficient clonal parent support and will be")
logger.info("relocated to phylogenetically independent branches.")
logger.info("-" * 80)

# Recalculate backbone nodes (T_current may have changed in Step 7.3)
current_backbone_nodes = get_first_level_backbone_nodes(T_current)
expanded_mutations_of_current_backbone_nodes = [mutation for node in current_backbone_nodes for mutation in node.split('|')]

# ---- Compute ancestral retention fraction ----
T_checkpoint_outgroup = copy.deepcopy(T_current)
M_checkpoint_outgroup = M_current.copy()

M_for_fp_ratio_and_fn_ratio_outgroup = M_checkpoint_outgroup.drop(columns=['ROOT'], errors='ignore')
mutations_on_T_current_outgroup = M_for_fp_ratio_and_fn_ratio_outgroup.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()
M_for_fp_ratio_and_fn_ratio_outgroup = split_merged_columns(M_for_fp_ratio_and_fn_ratio_outgroup, mutations_on_T_current_outgroup)

df_intersection_and_inter_vs_fn_flipping_ratio_per_mutation = calculate_intersection_and_inter_vs_fn_flipping_ratio_per_mutation(
    T_checkpoint_outgroup, M_checkpoint_outgroup, I_attached
)

df_intersection_and_inter_vs_fn_flipping_ratio_per_mutation.to_csv(
    os.path.join(outputpath_full, "df_intersection_and_inter_vs_fn_flipping_ratio_per_mutation.csv"), 
    sep=","
)

intersection_vs_fn_flipping_ratio_cutoff = params['intersection_vs_fn_flipping_ratio_cutoff']
intersection_cell_count_on_mutation_cutoff = params['intersection_cell_count_on_mutation_cutoff']
intersection_cell_ratio_on_mutation_cutoff = params['intersection_cell_ratio_on_mutation_cutoff']

# ---- Identify M_orphan^(retention) ----
outgroup_mutations = df_intersection_and_inter_vs_fn_flipping_ratio_per_mutation[(
    (df_intersection_and_inter_vs_fn_flipping_ratio_per_mutation['parent_retention_ratio'] <= intersection_vs_fn_flipping_ratio_cutoff) & 
    ((df_intersection_and_inter_vs_fn_flipping_ratio_per_mutation['intersection_cell_count_on_mutation'] <= intersection_cell_count_on_mutation_cutoff) | 
     (df_intersection_and_inter_vs_fn_flipping_ratio_per_mutation['intersection_cell_ratio_on_mutation'] <= intersection_cell_ratio_on_mutation_cutoff))
)]['mutation'].tolist()

# Exclude backbone mutations
outgroup_mutations_but_backbone = [
    i for i in outgroup_mutations 
    if i not in list(set(expanded_mutations_of_current_backbone_nodes))
]

logger.info(f"  |M_orphan^(retention)|: {len(outgroup_mutations)}")
logger.info(f"  |M_orphan^(retention) \ M_backbone|: {len(outgroup_mutations_but_backbone)}")

# ---- Identify M_orphan-progeny ----
sorted_rehanged_mutations_all_outgroup = []

if outgroup_mutations_but_backbone:
    nodes_outgroup_mutations_but_backbone = list(set([
        find_mutation_column(mutation, M_checkpoint_outgroup.columns) 
        for mutation in outgroup_mutations_but_backbone
    ]))
    
    node_dict = {node.name: node for node in T_checkpoint_outgroup.traverse()}
    daughter_nodes_of_outgroup_mutations_but_backbone = {
        mutation: get_all_daughter_mutations(node_dict[mutation])
        for mutation in nodes_outgroup_mutations_but_backbone
        if mutation in node_dict
    }
    daughter_mutations_of_outgroup_mutations_but_backbone = list(set([
        daughter for daughters_list in daughter_nodes_of_outgroup_mutations_but_backbone.values() 
        for daughter in daughters_list
        if daughter not in outgroup_mutations_but_backbone
    ]))
    
    sorted_outgroup_mutations_but_backbone = [
        i for i in I_attached.columns 
        if i in outgroup_mutations_but_backbone
    ]
    sorted_daughter_mutations_of_outgroup_mutations_but_backbone = [
        i for i in I_attached.columns 
        if i in daughter_mutations_of_outgroup_mutations_but_backbone
    ]
    sorted_rehanged_mutations_all_outgroup = (
        sorted_outgroup_mutations_but_backbone + 
        sorted_daughter_mutations_of_outgroup_mutations_but_backbone
    )
    
    logger.info(f"  |M_orphan-progeny|: {len(daughter_mutations_of_outgroup_mutations_but_backbone)}")

# ---- Prune and regraft M_orphan^(retention) and M_orphan-progeny ----
external_mutations_by_intersection_vs_fn = []

if len(sorted_rehanged_mutations_all_outgroup) > 0:
    logger.info(f"  Pruning and regrafting {len(sorted_rehanged_mutations_all_outgroup)} orphaned mutations...")
    
    for remove_mut_by_once in sorted_rehanged_mutations_all_outgroup:
        T_removed_outgroup, M_removed_outgroup = remove_mutations_from_tree_and_matrix(
            T_checkpoint_outgroup, M_checkpoint_outgroup, [remove_mut_by_once]
        )
        logger.info(f"    Tree after pruning: {M_removed_outgroup.shape[0]} cells, {M_removed_outgroup.shape[1]} mutations")
        
        T_current = copy.deepcopy(T_removed_outgroup)
        M_current = M_removed_outgroup.copy()
        
        external_mutations_by_intersection_vs_fn_temp, conflict_mutations_temp, T_current, M_current, root_mutations = attach_mutations_to_current_tree(
            sorted_attached_mutations=[remove_mut_by_once],
            T_current=T_current,
            M_current=M_current,
            I_attached=I_attached,
            P_attached=P_attached,
            ω_NA=omega_NA,
            fnfp_ratio=fnfp_ratio,
            φ=phi,
            logger=logger,
            root_mutations=root_mutations
        )
        all_conflict_mutations.extend(conflict_mutations_temp)
        external_mutations_by_intersection_vs_fn.extend(external_mutations_by_intersection_vs_fn_temp)

logger.info("")


# ----------------------------------------------------------------------------
# Step 7 Summary
# ----------------------------------------------------------------------------
logger.info("=" * 80)
logger.info("Step 7 COMPLETED: Tree refined by discordance-guided pruning")
logger.info("-" * 80)

# Print tree structure after PRG refinement
logger.info("Refined tree structure:")
logger.info("-" * 40)
print_tree_logger(T_current)
logger.info("-" * 40)

# Count mutations (exclude ROOT)
mutations_on_tree = [m for m in M_current.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist() if m != 'ROOT']

logger.info(f"  Mutations on tree (excluding ROOT): {len(mutations_on_tree)}")
logger.info(f"  Cells: {M_current.shape[0]}")
logger.info(f"  Merged columns (including ROOT): {M_current.shape[1]}")
logger.info(f"  |M_ambiguous|: {len(all_conflict_mutations)}")
logger.info("=" * 80)
logger.info("")


# ============================================================================
# CHECKPOINT: Weighted Discordance Index (Omega) Computation
# ============================================================================
# This checkpoint computes the weighted discordance index (Omega) before quality
# control, serving as a baseline for assessing QC effectiveness.
# ============================================================================

logger.info("=" * 80)
logger.info("CHECKPOINT: Weighted Discordance Index (Omega) computation")
logger.info("-" * 80)
logger.info("Computing the weighted discordance index (Omega) to quantify the")
logger.info("aggregate phylogenetic inconsistency prior to quality control.")
logger.info("")
logger.info("  Omega = N_deltaFP + lambda * N_deltaFN,  where lambda = 0.1")
logger.info("")
logger.info("  - N_deltaFP: total false positive discordance count")
logger.info("  - N_deltaFN: total false negative discordance count")
logger.info("  - lambda: empirical FN/FP discordance weight")
logger.info("-" * 80)

M_for_omega = M_current.drop(columns=['ROOT'], errors='ignore')
mutations_on_tree_for_omega = M_for_omega.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()
M_for_omega = split_merged_columns(M_for_omega, mutations_on_tree_for_omega)

M_for_omega_clean = M_for_omega.loc[:, (M_for_omega != 0).any(axis=0)]
M_for_omega_clean = M_for_omega_clean.loc[(M_for_omega_clean != 0).any(axis=1)]

I_for_omega = I_attached.loc[M_for_omega_clean.index, M_for_omega_clean.columns].replace({np.nan: 3}).astype(int)

N_deltaFP = ((I_for_omega == 1) & (M_for_omega_clean == 0)).sum().sum()
N_deltaFN = ((I_for_omega == 0) & (M_for_omega_clean == 1)).sum().sum()
N_NAto0 = ((I_for_omega == 3) & (M_for_omega_clean == 0)).sum().sum()
N_NAto1 = ((I_for_omega == 3) & (M_for_omega_clean == 1)).sum().sum()

omega_before_qc = N_deltaFP + params['fnfp_ratio'] * N_deltaFN

logger.info("")
logger.info("  ┌─────────────────────────────────────────────────────────────────────┐")
logger.info("  │              WEIGHTED DISCORDANCE INDEX (PRE-QC)                   │")
logger.info("  ├─────────────────────────────────────────────────────────────────────┤")
logger.info(f"  │  Weighted Discordance Index (Omega)        : {omega_before_qc:>10.2f}      │")
logger.info(f"  │    - delta_FP discordance                  : {N_deltaFP:>10}          │")
logger.info(f"  │    - delta_FN discordance                  : {N_deltaFN:>10}          │")
logger.info(f"  │    - NA->0 imputations                     : {N_NAto0:>10}          │")
logger.info(f"  │    - NA->1 imputations                     : {N_NAto1:>10}          │")
logger.info(f"  │    - FN/FP weight (lambda)                 : {params['fnfp_ratio']:>10.1f}      │")
logger.info("  ├─────────────────────────────────────────────────────────────────────┤")
logger.info(f"  │  Cells                                   : {M_for_omega_clean.shape[0]:>10}│")
logger.info(f"  │  Mutations                               : {M_for_omega_clean.shape[1]:>10}│")
logger.info("  └─────────────────────────────────────────────────────────────────────┘")

omega_checkpoint_file = os.path.join(outputpath_full, "weighted_discordance_index_pre_qc.txt")
with open(omega_checkpoint_file, 'w') as f:
    f.write("=" * 60 + "\n")
    f.write("WEIGHTED DISCORDANCE INDEX (PRE-QC)\n")
    f.write("=" * 60 + "\n\n")
    f.write(f"Omega (Weighted Discordance Index): {omega_before_qc:.2f}\n")
    f.write(f"  N_deltaFP: {N_deltaFP}\n")
    f.write(f"  N_deltaFN: {N_deltaFN}\n")
    f.write(f"  lambda (FN/FP weight): {params['fnfp_ratio']}\n\n")
    f.write(f"Cells: {M_for_omega_clean.shape[0]}\n")
    f.write(f"Mutations: {M_for_omega_clean.shape[1]}\n")
    f.write("=" * 60 + "\n")

logger.info(f"  Checkpoint saved to: {omega_checkpoint_file}")
logger.info("")


# ------------------------------
# Step 8: Quality control & filtration (Steps 8.1-8.4)
# ------------------------------

# ============================================================================
# Step 8 identifies and removes low-quality cells and mutations:
#   8.1: C_orphan - Phylogenetically orphaned cells (lack ancestral support)
#   8.2: M_artifact - Artifactual loci (persistently high global FP discordance)
#   8.3: C_chimeric - Chimeric cells (doublets, excessive FP discordance)
#   8.4: Final tree status summary (pre-output)
#
# Input:  Refined tree (T_current, M_current)
# Output: Final high-confidence tree (T_current, M_current)
# ============================================================================

logger.info("=" * 80)
logger.info("Step 8: Quality control & filtration")
logger.info("=" * 80)
logger.info("")

# ---- Get QC control parameters ----
logger.info(f"  Remove artifact mutations: {remove_artifact_mutations}")
logger.info("")


# ----------------------------------------------------------------------------
# Step 8.1: Identification of phylogenetically orphaned cells (C_orphan)
# ----------------------------------------------------------------------------
logger.info("STEP 8.1: Identification of phylogenetically orphaned cells (C_orphan)")
logger.info("-" * 80)
logger.info("Identifying phylogenetically orphaned cells that lack sufficient")
logger.info("ancestral mutation support.")
logger.info("")
logger.info("  Criteria for C_orphan classification:")
logger.info("    - Intersection count <= 1 (minimal clonal support)")
logger.info("    - FN discordance >= 1 (persistent false negatives)")
logger.info("    - NA->1 imputations >= 2 (unreliable imputation)")
logger.info("-" * 80)

# Recalculate backbone nodes (T_current may have changed in Step 7)
current_backbone_nodes = get_first_level_backbone_nodes(T_current)
expanded_mutations_of_current_backbone_nodes = [mutation for node in current_backbone_nodes for mutation in node.split('|')]

# ---- Compute cell-level metrics ----
T_checkpoint_wireless_cells = copy.deepcopy(T_current)
M_checkpoint_wireless_cells = M_current.copy()

M_for_fp_ratio_and_fn_ratio_wireless_cells = M_checkpoint_wireless_cells.drop(columns=['ROOT'], errors='ignore')
mutations_on_T_current_wireless_cells = M_for_fp_ratio_and_fn_ratio_wireless_cells.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()
M_for_fp_ratio_and_fn_ratio_wireless_cells = split_merged_columns(M_for_fp_ratio_and_fn_ratio_wireless_cells, mutations_on_T_current_wireless_cells)

df_intersection_and_flipping_to_1_count_per_cell = calculate_intersection_and_flipping_to_1_count_per_cell(
    M_for_fp_ratio_and_fn_ratio_wireless_cells, I_attached
)

df_intersection_and_flipping_to_1_count_per_cell.to_csv(
    os.path.join(outputpath_full, "df_intersection_and_flipping_to_1_count_per_cell.csv"), 
    sep=","
)

# ---- Identify C_orphan ----
intersection_count_per_cells_cutoff = params['intersection_count_per_cells_cutoff']
flipping_count_fn_per_cells_cutoff = params['flipping_count_fn_per_cells_cutoff']
flipping_to_1_count_per_cells_cutoff = params['flipping_to_1_count_per_cells_cutoff']

df_wireless_cells_filter = df_intersection_and_flipping_to_1_count_per_cell.loc[
    ((df_intersection_and_flipping_to_1_count_per_cell['intersection_count'] == intersection_count_per_cells_cutoff) & 
     (df_intersection_and_flipping_to_1_count_per_cell['flipping_count_fn'] >= flipping_count_fn_per_cells_cutoff) & 
     (df_intersection_and_flipping_to_1_count_per_cell['flipping_to_1_count'] >= flipping_to_1_count_per_cells_cutoff))
]

identified_wireless_cells = list(df_wireless_cells_filter.index)

# ---- Identify doublet candidates from orphaned mutations ----
conflicting_cells_as_doublets_from_parents_format_nested_list = df_intersection_and_inter_vs_fn_flipping_ratio_per_mutation.loc[
    df_intersection_and_inter_vs_fn_flipping_ratio_per_mutation['mutation'].isin(outgroup_mutations), 
    'intersection_cells_on_mutation_parents'
]
conflicting_cells_as_doublets_from_parents_format_flat_list = sum(conflicting_cells_as_doublets_from_parents_format_nested_list.tolist(), [])

conflicting_cells_as_doublets_from_children_format_nested_list = df_intersection_and_inter_vs_fn_flipping_ratio_per_mutation.loc[
    df_intersection_and_inter_vs_fn_flipping_ratio_per_mutation['mutation'].isin(outgroup_mutations), 
    'intersection_cells_on_mutation_children'
]
conflicting_cells_as_doublets_from_children_format_flat_list = sum(conflicting_cells_as_doublets_from_children_format_nested_list.tolist(), [])

# ---- Remove identified cells ----
to_be_removed_cells = list(set(
    identified_wireless_cells + 
    conflicting_cells_as_doublets_from_parents_format_flat_list + 
    conflicting_cells_as_doublets_from_children_format_flat_list
))

logger.info(f"  |C_orphan| (wireless cells): {len(identified_wireless_cells)}")
logger.info(f"  |C_doublet_candidates|: {len(conflicting_cells_as_doublets_from_parents_format_flat_list)}")
logger.info(f"  |Total cells to remove|: {len(to_be_removed_cells)}")

if len(to_be_removed_cells) > 0:
    with open(os.path.join(outputpath_full, "likely_doublet_cells_removed_from_tree_by_fn_flipping.csv"), 'w') as f:
        for cell in to_be_removed_cells:
            f.write(cell + '\n')
    
    M_current = M_current.drop(to_be_removed_cells, errors='ignore')
    logger.info(f"  Removed {len(to_be_removed_cells)} cells from the tree")

logger.info("")


# ----------------------------------------------------------------------------
# Step 8.2: Artifactual loci (M_artifact) and chimeric cells (C_chimeric)
# ----------------------------------------------------------------------------
logger.info("STEP 8.2: Artifactual loci and chimeric cell detection")
logger.info("-" * 80)
logger.info("Identifying and permanently removing:")
logger.info("")
logger.info("  (1) Artifactual loci (M_artifact): mutations with persistently")
logger.info("      high global FP discordance (r_FP^(global) >= threshold)")
logger.info("  (2) Chimeric cells (C_chimeric): cells exhibiting excessive")
logger.info("      cross-mutation FP discordance (>50%)")
logger.info("")
logger.info("These represent the most severe data quality issues and are")
logger.info("permanently excluded from the final phylogeny.")
logger.info("-" * 80)

# Recalculate backbone nodes (T_current may have changed in Step 8.1)
current_backbone_nodes = get_first_level_backbone_nodes(T_current)
expanded_mutations_of_current_backbone_nodes = [mutation for node in current_backbone_nodes for mutation in node.split('|')]

# ---- Initialize variables (ensure they exist regardless of remove_artifact_mutations) ----
to_be_removed_mutations_by_fp_mutations_cross_all_cells = []
external_mutations_cross_all_cells_by_sorted_rehanged_mutations_all_by_fp_mutations_cross_all_cells = []

# ---- Compute comprehensive FP metrics ----
T_checkpoint_artifact_and_doublet = copy.deepcopy(T_current)
M_checkpoint_artifact_and_doublet = M_current.copy()

M_for_artifact_and_doublet = M_checkpoint_artifact_and_doublet.drop(columns=['ROOT'], errors='ignore')
mutations_on_T_current_artifact_and_doublet = M_for_artifact_and_doublet.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()
M_for_artifact_and_doublet = split_merged_columns(M_for_artifact_and_doublet, mutations_on_T_current_artifact_and_doublet)

df_fp_ratio_per_mutation_cross_all_cells, df_fp_ratio_per_cell_cross_all_muts, overall_metrics, fp_mutations_dict_cross_all_cells = calculate_comprehensive_fp_metrics(
    M_for_artifact_and_doublet, 
    I_attached.loc[M_checkpoint_artifact_and_doublet.index, M_for_artifact_and_doublet.columns]
)

df_fp_ratio_per_mutation_cross_all_cells.to_csv(
    os.path.join(outputpath_full, "df_fp_ratio_per_mutation_cross_all_cells.csv"), 
    sep=","
)
df_fp_ratio_per_cell_cross_all_muts.to_csv(
    os.path.join(outputpath_full, "df_fp_ratio_per_cell_cross_all_muts.csv"), 
    sep=","
)

fp_ratio_per_mutation_cross_all_cells_cutoff = params['fp_ratio_per_mutation_cross_all_cells_cutoff']
fp_count_per_mutation_cross_all_cells_cutoff = params['fp_count_per_mutation_cross_all_cells_cutoff']
fp_ratio_per_cell_cross_all_muts_cutoff = params['fp_ratio_per_cell_cross_all_muts_cutoff']

# ---- Identify M_artifact (permanently removed) ----
rehanged_fp_mutations_cross_all_cells = df_fp_ratio_per_mutation_cross_all_cells[
    (df_fp_ratio_per_mutation_cross_all_cells['fp_cells_ratio_per_mutation'] >= fp_ratio_per_mutation_cross_all_cells_cutoff) & 
    (df_fp_ratio_per_mutation_cross_all_cells['fp_cells_count'] >= fp_count_per_mutation_cross_all_cells_cutoff)
]['identifier'].tolist()

# Exclude backbone mutations from artifact identification
rehanged_fp_mutations_cross_all_cells_but_backbone = [
    i for i in rehanged_fp_mutations_cross_all_cells 
    if i not in list(set(expanded_mutations_of_current_backbone_nodes + scaffold_mutations))
]

# These are candidates for permanent removal
artifact_candidates = rehanged_fp_mutations_cross_all_cells_but_backbone

# ---- Apply removal control ----
if remove_artifact_mutations == "yes":
    to_be_removed_mutations_by_fp_mutations_cross_all_cells = artifact_candidates
    
    logger.warning("=" * 80)
    logger.warning(f"🔴 PERMANENTLY REMOVED ARTIFACTUAL LOCI (M_artifact)")
    logger.warning("-" * 80)
    logger.warning(f"  |M_artifact|: {len(to_be_removed_mutations_by_fp_mutations_cross_all_cells)} total")
    for mut in to_be_removed_mutations_by_fp_mutations_cross_all_cells:
        logger.warning(f"    - {mut}")
    logger.warning("=" * 80)
    
    pd.Series(to_be_removed_mutations_by_fp_mutations_cross_all_cells).to_csv(
        os.path.join(outputpath_full, "artifact_mutations_permanently_removed.csv"),
        index=False, header=['mutation']
    )
    
    # Remove artifact mutations from the tree
    M_current = M_current.drop(columns=to_be_removed_mutations_by_fp_mutations_cross_all_cells, errors='ignore')
    logger.info(f"  Removed {len(to_be_removed_mutations_by_fp_mutations_cross_all_cells)} artifact mutations from the tree")
    
    # ---- Get branch groups and daughter mutations for reattachment ----
    # (Only if there are artifact candidates to process)
    if len(artifact_candidates) > 0:
        ordered_branch_groups_for_rehanged_fp_mutations_cross_all_cells_but_backbone = find_ordered_branch_groups_for_rehanged_mutations_with_keys_as_earlist(
            T_current, artifact_candidates
        )
        
        filtered_fp_mutations_dict_cross_all_cells = {
            mut: other_muts 
            for mut, other_muts in fp_mutations_dict_cross_all_cells.items() 
            if mut in artifact_candidates
        }
        
        nodes_rehanged_fp_mutations_cross_all_cells_but_backbone = list(set([
            find_mutation_column(mutation, M_current.columns) 
            for mutation in artifact_candidates
        ]))
        
        node_dict = {node.name: node for node in T_current.traverse()}
        daughters_nodes_dict_cross_all_cells = {
            mutation: get_all_daughter_mutations(node_dict[mutation])
            for mutation in nodes_rehanged_fp_mutations_cross_all_cells_but_backbone
            if mutation in node_dict
        }
        daughters_mutations_dict_cross_all_cells = list(set([
            daughter for daughters_list in daughters_nodes_dict_cross_all_cells.values() 
            for daughter in daughters_list
            if daughter not in artifact_candidates
        ]))
        
        daughters_to_leaf_mutations_cross_all_cells = []
        fp_mutations_cross_all_cells = []
        
        for branch_mut in ordered_branch_groups_for_rehanged_fp_mutations_cross_all_cells_but_backbone:
            daughter_list = ordered_branch_groups_for_rehanged_fp_mutations_cross_all_cells_but_backbone[branch_mut]
            
            daughters_to_leaf_mutations_cross_all_cells = list({
                item 
                for key in daughter_list
                if key in daughters_mutations_dict_cross_all_cells
                for item in [key] + daughters_mutations_dict_cross_all_cells[key]
            })
            
            fp_mutations_cross_all_cells_init = list({
                item 
                for key in daughter_list
                if key in filtered_fp_mutations_dict_cross_all_cells
                for item in [key] + filtered_fp_mutations_dict_cross_all_cells[key]
            })
            fp_mutations_cross_all_cells = [
                i for i in fp_mutations_cross_all_cells_init 
                if i not in daughter_list 
                and i not in daughters_to_leaf_mutations_cross_all_cells 
                and i not in list(set(expanded_mutations_of_current_backbone_nodes + scaffold_mutations))
            ]
        
        daughters_to_leaf_mutations_cross_all_cells_qc = [
            mut for mut in daughters_to_leaf_mutations_cross_all_cells 
            if df_features_new[mut]['mutant_cellnum'] > 5
        ]
        fp_mutations_cross_all_cells_qc = set(
            fp_mutations_cross_all_cells + 
            [i for i in daughters_to_leaf_mutations_cross_all_cells 
             if i not in daughters_to_leaf_mutations_cross_all_cells_qc]
        )
        
        sorted_fp_mutations_cross_all_cells = [
            i for i in I_attached.columns 
            if i in fp_mutations_cross_all_cells_qc
        ]
        sorted_daughters_to_leaf_mutations_cross_all_cells = [
            i for i in I_attached.columns 
            if i in daughters_to_leaf_mutations_cross_all_cells_qc
        ]
        sorted_rehanged_mutations_all_cross_all_cells = (
            sorted_fp_mutations_cross_all_cells + 
            sorted_daughters_to_leaf_mutations_cross_all_cells
        )
        
        # Mutations to reattach (excluding permanently removed ones)
        sorted_rehanged_mutations_all_by_fp_mutations_cross_all_cells = [
            i for i in sorted_rehanged_mutations_all_cross_all_cells 
            if i not in artifact_candidates
        ]
        remove_mutations_for_rebuild = list(set(
            artifact_candidates + 
            sorted_rehanged_mutations_all_by_fp_mutations_cross_all_cells
        ))
        
        # ---- Prune and reattach ----
        if len(remove_mutations_for_rebuild) > 0:
            T_removed_cross_all_cells, M_removed_cross_all_cells = remove_mutations_from_tree_and_matrix(
                T_checkpoint_artifact_and_doublet, M_checkpoint_artifact_and_doublet, 
                remove_mutations_for_rebuild
            )
            logger.info(f"  Tree after pruning: {M_removed_cross_all_cells.shape[0]} cells, {M_removed_cross_all_cells.shape[1]} mutations")
            
            T_current = copy.deepcopy(T_removed_cross_all_cells)
            M_current = M_removed_cross_all_cells.copy()
        
        if len(sorted_rehanged_mutations_all_by_fp_mutations_cross_all_cells) > 0:
            external_mutations_cross_all_cells_by_sorted_rehanged_mutations_all_by_fp_mutations_cross_all_cells, conflict_mutations_temp, T_current, M_current, root_mutations = attach_mutations_to_current_tree(
                sorted_attached_mutations=sorted_rehanged_mutations_all_by_fp_mutations_cross_all_cells,
                T_current=T_current,
                M_current=M_current,
                I_attached=I_attached,
                P_attached=P_attached,
                ω_NA=omega_NA,
                fnfp_ratio=fnfp_ratio,
                φ=phi,
                logger=logger,
                root_mutations=root_mutations
            )
            all_conflict_mutations.extend(conflict_mutations_temp)

else:
    # remove_artifact_mutations == "no"
    logger.info(f"  ⚠️ Artifact mutation removal is DISABLED (remove_artifact_mutations=no)")
    logger.info(f"  |M_artifact| candidates identified but NOT removed: {len(artifact_candidates)}")
    
    # Still save the list for reference
    pd.Series(artifact_candidates).to_csv(
        os.path.join(outputpath_full, "artifact_mutations_identified_but_not_removed.csv"),
        index=False, header=['mutation']
    )
    # to_be_removed_mutations_by_fp_mutations_cross_all_cells is already initialized as []
    # external_mutations_cross_all_cells_by_sorted_rehanged_mutations_all_by_fp_mutations_cross_all_cells is already initialized as []

# ---- Identify C_chimeric (chimeric cells / doublets) ----
identified_doublet_cells = df_fp_ratio_per_cell_cross_all_muts.loc[
    df_fp_ratio_per_cell_cross_all_muts['fp_muts_ratio_per_cell'] >= fp_ratio_per_cell_cross_all_muts_cutoff
]['cell_id'].tolist()

logger.info("")
logger.info(f"  |C_chimeric| (doublet cells): {len(identified_doublet_cells)}")
logger.info("  - These cells exhibit excessive cross-mutation FP discordance")
logger.info("  - Likely represent doublets or technical artifacts")

if len(identified_doublet_cells) > 0:
    with open(os.path.join(outputpath_full, "likely_doublet_cells_removed_from_tree_by_fp_ratio.csv"), 'w') as f:
        for cell in identified_doublet_cells:
            f.write(cell + '\n')
    
    M_current = M_current.drop(identified_doublet_cells, errors='ignore')
    logger.info(f"  Removed {len(identified_doublet_cells)} chimeric cells from the tree")

logger.info("")


# ----------------------------------------------------------------------------
# Step 8.3: Attach conflict mutations to ROOT
# ----------------------------------------------------------------------------
logger.info("STEP 8.3: Attach conflict mutations to ROOT")
logger.info("-" * 80)
logger.info("Mutations with ambiguous placements (M_ambiguous) are attached")
logger.info("to the ROOT node as a final resolution step.")
logger.info("-" * 80)

logger.info(f"  |M_ambiguous| before resolution: {len(all_conflict_mutations)}")

final_remained_mutations = []
final_conflict_mutations = []

if len(all_conflict_mutations) > 0:
    subtree_groups = cluster_external_mutations_by_intersection(I_attached, all_conflict_mutations)
    logger.info(f"  Processing {len(subtree_groups)} conflict clusters")
    
    final_remained_mutations, final_conflict_mutations, T_current, M_current, root_mutations = process_external_mutations_by_subtree_groups(
        subtree_groups=subtree_groups,
        T_current=T_current,
        M_current=M_current,
        I_attached=I_attached,
        P_attached=P_attached,
        ω_NA=omega_NA,
        fnfp_ratio=fnfp_ratio,
        φ=phi,
        logger=logger,
        root_mutations=root_mutations
    )

logger.info(f"  |M_ambiguous| after resolution: {len(final_conflict_mutations)}")
logger.info("")


# ----------------------------------------------------------------------------
# Step 8.4: Final tree status summary (pre-output)
# ----------------------------------------------------------------------------
logger.info("=" * 80)
logger.info("STEP 8.4: Final tree status summary (pre-output)")
logger.info("=" * 80)

# ---- Final tree status ----
M_final = M_current.drop(columns=['ROOT'], errors='ignore')
mutations_on_T_final = M_final.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()
M_final_split = split_merged_columns(M_final, mutations_on_T_final)
final_cleaned_M = M_final_split.loc[(M_final_split != 0).any(axis=1)]

logger.info("Final tree status (before output):")
logger.info(f"  Cells: {final_cleaned_M.shape[0]}")
logger.info(f"  Mutations: {final_cleaned_M.shape[1]}")
logger.info(f"  Shape: {final_cleaned_M.shape}")
logger.info("-" * 80)

# ---- Mutations not on tree ----
# Initialize summary_data with all keys
summary_data = {
    'final_remained_mutations': len(final_remained_mutations),
    'final_conflict_mutations': len(final_conflict_mutations),
    'root_mutations': len(root_mutations),
    'artifact_mutations_permanently_removed': len(to_be_removed_mutations_by_fp_mutations_cross_all_cells),
}

# Try to get external mutations count (may not exist)
try:
    external_count = len(set(
        external_mutations_by_intersection_vs_fn +
        external_mutations_cross_all_cells_by_sorted_rehanged_mutations_all_by_fp_mutations_cross_all_cells
    ))
    summary_data['all_external_mutations'] = external_count
except NameError:
    summary_data['all_external_mutations'] = 0

logger.info("Mutations NOT on tree:")
total = 0
for name, count in summary_data.items():
    logger.info(f"  {name}: {count}")
    total += count
logger.info(f"  Total: {total}")
logger.info("=" * 80)

# ---- Save to file ----
summary_file = os.path.join(outputpath_full, "final_tree_status_summary_pre_output.txt")
with open(summary_file, 'w') as f:
    f.write("=" * 60 + "\n")
    f.write("FINAL TREE STATUS SUMMARY (PRE-OUTPUT)\n")
    f.write("=" * 60 + "\n\n")
    f.write(f"Cells: {final_cleaned_M.shape[0]}\n")
    f.write(f"Mutations: {final_cleaned_M.shape[1]}\n")
    f.write(f"Shape: {final_cleaned_M.shape}\n\n")
    f.write("-" * 60 + "\n")
    f.write("Mutations NOT on tree:\n")
    for name, count in summary_data.items():
        f.write(f"  {name}: {count}\n")
    f.write(f"  Total: {total}\n")
    f.write("=" * 60 + "\n")

logger.info(f"Saved to: {summary_file}")
logger.info("")


# ----------------------------------------------------------------------------
# Step 8 Summary
# ----------------------------------------------------------------------------
logger.info("=" * 80)
logger.info("Step 8 COMPLETED: Quality control and filtration finished")
logger.info("-" * 80)

# Print tree structure after QC
logger.info("Final tree structure after QC:")
logger.info("-" * 40)
print_tree_logger(T_current)
logger.info("-" * 40)

# Count mutations (exclude ROOT)
mutations_on_tree = [m for m in M_current.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist() if m != 'ROOT']

logger.info(f"  Mutations on tree (excluding ROOT): {len(mutations_on_tree)}")
logger.info(f"  Cells: {M_current.shape[0]}")
logger.info(f"  Merged columns (including ROOT): {M_current.shape[1]}")
logger.info(f"  |M_artifact| (permanently removed): {len(to_be_removed_mutations_by_fp_mutations_cross_all_cells)}")
logger.info(f"  |C_chimeric| (removed): {len(identified_doublet_cells)}")
logger.info(f"  |C_orphan| (removed): {len(to_be_removed_cells)}")
logger.info("=" * 80)
logger.info("")


# ------------------------------
# Step 9: Post-processing & output (Steps 9.1-9.5)
# ------------------------------

# ============================================================================
# Step 9 formats the final tree and matrices, computes flipping statistics,
# and exports all results to files.
#
# Sub-steps:
#   9.1: Prepare final data - Expand merged columns and format tree
#   9.2: Export matrices and tree files
#   9.3: Identify flipping spots (delta_FP, delta_FN, NA->0, NA->1)
#   9.4: Calculate total flipping counts and final Omega
#   9.5: Export tree formats (JSON, TXT) and clone assignments
#
# Input:  Final tree (T_current, M_current)
# Output: All result files in output directory
# ============================================================================

logger.info("=" * 80)
logger.info("Step 9: Post-processing & output")
logger.info("=" * 80)
logger.info("")


# ----------------------------------------------------------------------------
# Step 9.1: Prepare final data
# ----------------------------------------------------------------------------
logger.info("STEP 9.1: Prepare final data")
logger.info("-" * 80)

# ---- Drop ROOT column and add root mutations back ----
M_current_filtered = M_current.drop(columns=['ROOT'], errors='ignore')

for mut_on_root in root_mutations:
    M_current_filtered.insert(0, mut_on_root, 1)

# ---- Expand merged columns ----
mutations_on_T_current = M_current_filtered.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()

T_full = copy.deepcopy(T_current)
M_full = split_merged_columns(M_current_filtered, mutations_on_T_current)

logger.info("Final full-resolved tree:")
print_tree_logger(T_full)

logger.info(f"  Final tree cells: {M_full.shape[0]}")
logger.info(f"  Final tree mutations: {M_full.shape[1]}")
logger.info("")


# ----------------------------------------------------------------------------
# Step 9.2: Output results
# ----------------------------------------------------------------------------
logger.info("STEP 9.2: Output result files")
logger.info("-" * 80)

# ---- Create output directory ----
phylo_dir = os.path.join(outputpath_full, "phylo")
os.makedirs(phylo_dir, exist_ok=True)

# ---- Export I matrix with NA=3 ----
I_full_withNA3 = I_attached.replace({np.nan: 3}).astype(int)
I_full_withNA3.to_csv(os.path.join(phylo_dir, "I_full_withNA3.txt"), sep="\t")

# ---- Export M matrix ----
WriteTfile(os.path.join(phylo_dir, "M_full_basedPivots.filtered_sites_inferred"), M_full, M_full.index.tolist(), M_full.columns.tolist(), judge="yes")

# ---- Clean and export final matrices ----
final_cleaned_M_full = M_full.loc[:, (M_full != 0).any(axis=0)]
final_cleaned_M_full = final_cleaned_M_full.loc[(final_cleaned_M_full != 0).any(axis=1)]

kept_rows = final_cleaned_M_full.index
kept_cols = final_cleaned_M_full.columns

final_cleaned_I_full_withNA3 = I_full_withNA3.loc[kept_rows, kept_cols]

WriteTfile(os.path.join(phylo_dir, "final_cleaned_M_full_basedPivots.filtered_sites_inferred"), 
           final_cleaned_M_full, final_cleaned_M_full.index.tolist(), final_cleaned_M_full.columns.tolist(), judge="yes")
final_cleaned_I_full_withNA3.to_csv(os.path.join(phylo_dir, "final_cleaned_I_full_withNA3_for_circosPlot.txt"), sep="\t")

logger.info(f"  Output directory: {phylo_dir}")
logger.info("")


# ----------------------------------------------------------------------------
# Step 9.3: Identify flipping spots
# ----------------------------------------------------------------------------
logger.info("STEP 9.3: Identify flipping spots")
logger.info("-" * 80)
logger.info("Computing discordance spots for each mutation:")
logger.info("  - delta_FN: tree predicts 1, observed 0 (false negative)")
logger.info("  - delta_FP: tree predicts 0, observed 1 (false positive)")
logger.info("  - NA->1: missing data imputed as mutant")
logger.info("  - NA->0: missing data imputed as reference")
logger.info("-" * 80)

df_bin_withNA3_for_flipping = final_cleaned_I_full_withNA3.copy()
df_phylogeny = final_cleaned_M_full.copy()

# ---- Compute flipping spots ----
false_negative_flipping_spots = df_bin_withNA3_for_flipping.apply(
    lambda col: find_flipping_spots(col, df_phylogeny[col.name], condition_in_bin=0, condition_phylogeny=1)
)
false_positive_flipping_spots = df_bin_withNA3_for_flipping.apply(
    lambda col: find_flipping_spots(col, df_phylogeny[col.name], condition_in_bin=1, condition_phylogeny=0)
)
NAto1_flipping_spots = df_bin_withNA3_for_flipping.apply(
    lambda col: find_flipping_spots(col, df_phylogeny[col.name], condition_in_bin=3, condition_phylogeny=1)
)
NAto0_flipping_spots = df_bin_withNA3_for_flipping.apply(
    lambda col: find_flipping_spots(col, df_phylogeny[col.name], condition_in_bin=3, condition_phylogeny=0)
)

# ---- Handle empty results ----
if false_negative_flipping_spots.empty:
    false_negative_flipping_spots = {col: [] for col in df_bin_withNA3_for_flipping.columns}

if false_positive_flipping_spots.empty:
    false_positive_flipping_spots = {col: [] for col in df_bin_withNA3_for_flipping.columns}

if NAto1_flipping_spots.empty:
    NAto1_flipping_spots = {col: [] for col in df_bin_withNA3_for_flipping.columns}

if NAto0_flipping_spots.empty:
    NAto0_flipping_spots = {col: [] for col in df_bin_withNA3_for_flipping.columns}

# ---- Build flipping spots dataframe ----
df_flipping_spots = pd.DataFrame({
    'Mutation': df_bin_withNA3_for_flipping.columns,
    'delta_FN_spots': [', '.join(false_negative_flipping_spots.get(col, [])) for col in df_bin_withNA3_for_flipping.columns],
    'delta_FP_spots': [', '.join(false_positive_flipping_spots.get(col, [])) for col in df_bin_withNA3_for_flipping.columns],
    'NA_to_1_spots': [', '.join(NAto1_flipping_spots.get(col, [])) for col in df_bin_withNA3_for_flipping.columns],
    'NA_to_0_spots': [', '.join(NAto0_flipping_spots.get(col, [])) for col in df_bin_withNA3_for_flipping.columns]
})
df_flipping_spots.to_csv(os.path.join(phylo_dir, "df_flipping_spots.txt"), sep="\t", index=False)

logger.info("")


# ----------------------------------------------------------------------------
# Step 9.4: Calculate total flipping counts
# ----------------------------------------------------------------------------
logger.info("STEP 9.4: Calculate total flipping counts")
logger.info("-" * 80)
logger.info("Computing aggregate discordance statistics for the final tree.")
logger.info("")
logger.info("  Omega = N_deltaFP + lambda * N_deltaFN,  where lambda = 0.1")
logger.info("-" * 80)

# ---- Compute total discordance counts ----
total_FN_flipping = ((df_bin_withNA3_for_flipping == 0) & (df_phylogeny == 1)).sum().sum()
total_FP_flipping = ((df_bin_withNA3_for_flipping == 1) & (df_phylogeny == 0)).sum().sum()
total_NAto0 = ((df_bin_withNA3_for_flipping == 3) & (df_phylogeny == 0)).sum().sum()
total_NAto1 = ((df_bin_withNA3_for_flipping == 3) & (df_phylogeny == 1)).sum().sum()

omega_final = total_FP_flipping + params['fnfp_ratio'] * total_FN_flipping

logger.info("")
logger.info("  ┌─────────────────────────────────────────────────────────────────────┐")
logger.info("  │              WEIGHTED DISCORDANCE INDEX (FINAL)                    │")
logger.info("  ├─────────────────────────────────────────────────────────────────────┤")
logger.info(f"  │  Weighted Discordance Index (Omega)        : {omega_final:>10.2f}      │")
logger.info(f"  │    - delta_FP discordance                  : {total_FP_flipping:>10}          │")
logger.info(f"  │    - delta_FN discordance                  : {total_FN_flipping:>10}          │")
logger.info(f"  │    - NA->0 imputations                     : {total_NAto0:>10}       │")
logger.info(f"  │    - NA->1 imputations                     : {total_NAto1:>10}       │")
logger.info(f"  │    - FN/FP weight (lambda)                 : {params['fnfp_ratio']:>10.1f}      │")
logger.info("  └─────────────────────────────────────────────────────────────────────┘")

# ---- Save total flipping counts ----
df_total_flipping_count = pd.DataFrame({
    'total_delta_FP': [total_FP_flipping],
    'total_delta_FN': [total_FN_flipping],
    'total_NA_to_0': [total_NAto0],
    'total_NA_to_1': [total_NAto1],
    'weighted_discordance_index_Omega': [omega_final]
})
df_total_flipping_count.to_csv(os.path.join(phylo_dir, "df_total_flipping_count.txt"), sep="\t", index=False)

# ---- Save per-site flipping counts ----
df_flip_counts_tree = calculate_flip_counts_per_site(df_bin_withNA3_for_flipping, df_phylogeny)
df_flip_counts_tree.to_csv(os.path.join(phylo_dir, "df_flipping_count_for_each_mut.txt"), sep="\t", index=True)

print(df_total_flipping_count.iloc[0,:])

logger.info(f"The shape of final_cleaned_M_full.shape: {final_cleaned_M_full.shape}")
print(final_cleaned_M_full.shape)
logger.info("")


# ----------------------------------------------------------------------------
# Step 9.5: Tree format and clone information
# ----------------------------------------------------------------------------
logger.info("STEP 9.5: Tree format and clone information")
logger.info("-" * 80)
logger.info("Exporting tree in multiple formats and assigning clone labels.")
logger.info("-" * 80)

# ---- Export tree as JSON ----
tree_dict = tree_to_dict(T_full)

with open(os.path.join(phylo_dir, 'final_cleaned_tree_node.json'), 'w') as f:
    json.dump(tree_dict, f, indent=4)

# ---- Export tree as text ----
T_full.save_to_file(os.path.join(phylo_dir, 'final_cleaned_tree_node.txt'))

# ---- Assign clone labels to cells ----
mutation_clones = get_mutation_clone_and_backbone_mut_as_keys_by_first_level_with_frequency(T_full, I_attached)
df_barcode_clones = assign_clone_labels(M_full, mutation_clones)

df_barcode_clones.to_csv(os.path.join(phylo_dir, "df_barcode_clones_from_phylo_tree.csv"), sep=',', index=False)

# ---- Export tree in Newick format ----
try:
    newick_str = tree_to_newick(T_full)
    with open(os.path.join(phylo_dir, 'final_cleaned_tree.newick'), 'w') as f:
        f.write(newick_str + ';')
    logger.info("  Exported Newick format: final_cleaned_tree.newick")
except NameError:
    logger.warning("  tree_to_newick function not available, skipping Newick export")

logger.info("")


# ----------------------------------------------------------------------------
# Step 9.6: Save run summary
# ----------------------------------------------------------------------------
logger.info("STEP 9.6: Save run summary")
logger.info("-" * 80)

run_summary_file = os.path.join(outputpath_full, "run_summary.txt")

# Get the number of final mutations (exclude ROOT)
final_mutations_count = len([m for m in M_full.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist() if m != 'ROOT'])

with open(run_summary_file, 'w') as f:
    f.write("=" * 80 + "\n")
    f.write("PhyloSOLID RUN SUMMARY\n")
    f.write("=" * 80 + "\n\n")
    f.write(f"Sample ID: {sampleid}\n")
    f.write(f"Run date: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    f.write(f"Seed: {args.seed}\n")
    f.write(f"Output path: {outputpath_full}\n\n")
    
    f.write("-" * 80 + "\n")
    f.write("PARAMETERS\n")
    f.write("-" * 80 + "\n")
    f.write(f"  is_predict_germ: {is_predict_germ}\n")
    f.write(f"  is_detect_passtree_by_dp: {is_detect_passtree_by_dp}\n")
    f.write(f"  is_filter_quality: {is_filter_quality}\n")
    f.write(f"  cv_rank_thresh: {cv_rank_thresh}\n")
    f.write(f"  remove_artifact_mutations: {remove_artifact_mutations}\n")
    f.write(f"  fnfp_ratio: {params['fnfp_ratio']}\n")
    f.write(f"  fp_ratio_cutoff_within_subclone: {params['fp_ratio_cutoff_within_subclone']}\n")
    f.write(f"  fp_ratio_cutoff_across_tree: {params['fp_ratio_cutoff_across_tree']}\n")
    f.write(f"  fn_ratio_cutoff_across_tree: {params['fn_ratio_cutoff_across_tree']}\n")
    f.write(f"  fp_ratio_persite_cutoff: {params['fp_ratio_persite_cutoff']}\n\n")
    
    f.write("-" * 80 + "\n")
    f.write("INPUT STATISTICS\n")
    f.write("-" * 80 + "\n")
    f.write(f"  Total mutations: {len(all_mutations)}\n")
    f.write(f"  Germline removed: {len(predicted_germline_mutations)}\n")
    f.write(f"  Somatic mutations: {len(somatic_mutations)}\n")
    f.write(f"  Scaffold mutations: {len(scaffold_mutations)}\n")
    f.write(f"  Accessory mutations: {len(attached_mutations)}\n")
    f.write(f"  Cells: {M_current.shape[0]}\n\n")
    
    f.write("-" * 80 + "\n")
    f.write("MUTATION CLASSIFICATION\n")
    f.write("-" * 80 + "\n")
    f.write(f"  M_scaffold: {len(scaffold_mutations)}\n")
    f.write(f"  M_accessory (integrated): {len(attached_mutations)}\n")
    f.write(f"  M_artifact (REMOVED): {len(to_be_removed_mutations_by_fp_mutations_cross_all_cells)}\n")
    f.write(f"  M_root (root-assigned): {len(root_mutations)}\n")
    f.write(f"  M_ambiguous (conflict): {len(all_conflict_mutations)}\n\n")
    
    f.write("-" * 80 + "\n")
    f.write("CELL CLASSIFICATION\n")
    f.write("-" * 80 + "\n")
    f.write(f"  C_resolved (final): {final_cleaned_M_full.shape[0]}\n")
    f.write(f"  C_orphan (REMOVED): {len(to_be_removed_cells)}\n")
    f.write(f"  C_chimeric (REMOVED): {len(identified_doublet_cells)}\n\n")
    
    f.write("-" * 80 + "\n")
    f.write("DISCORDANCE METRICS\n")
    f.write("-" * 80 + "\n")
    f.write(f"  Ω (pre-QC): {omega_before_qc:.2f}\n")
    f.write(f"  Ω (final): {omega_final:.2f}\n")
    f.write(f"  Ω reduction: {omega_before_qc - omega_final:.2f}\n\n")
    
    f.write("-" * 80 + "\n")
    f.write("OUTPUT FILES\n")
    f.write("-" * 80 + "\n")
    f.write(f"  Main results: {phylo_dir}/\n")
    f.write("  Files exported:\n")
    f.write("    - I_full_withNA3.txt\n")
    f.write("    - M_full_basedPivots.filtered_sites_inferred\n")
    f.write("    - final_cleaned_M_full_basedPivots.filtered_sites_inferred\n")
    f.write("    - final_cleaned_I_full_withNA3_for_circosPlot.txt\n")
    f.write("    - df_flipping_spots.txt\n")
    f.write("    - df_total_flipping_count.txt\n")
    f.write("    - df_flipping_count_for_each_mut.txt\n")
    f.write("    - final_cleaned_tree_node.json\n")
    f.write("    - final_cleaned_tree_node.txt\n")
    f.write("    - df_barcode_clones_from_phylo_tree.csv\n")
    f.write("=" * 80 + "\n")

logger.info(f"  Run summary saved to: {run_summary_file}")
logger.info("")


# ----------------------------------------------------------------------------
# Step 9 Summary
# ----------------------------------------------------------------------------
logger.info("=" * 80)
logger.info("Step 9 COMPLETED: All results exported successfully")
logger.info("-" * 80)
logger.info(f"  Output directory: {phylo_dir}")
logger.info("  Files exported:")
logger.info("    - I_full_withNA3.txt")
logger.info("    - M_full_basedPivots.filtered_sites_inferred")
logger.info("    - final_cleaned_M_full_basedPivots.filtered_sites_inferred")
logger.info("    - final_cleaned_I_full_withNA3_for_circosPlot.txt")
logger.info("    - df_flipping_spots.txt")
logger.info("    - df_total_flipping_count.txt")
logger.info("    - df_flipping_count_for_each_mut.txt")
logger.info("    - final_cleaned_tree_node.json")
logger.info("    - final_cleaned_tree_node.txt")
logger.info("    - df_barcode_clones_from_phylo_tree.csv")
logger.info("=" * 80)


# ============================================================================
# FINAL SUMMARY
# ============================================================================
logger.info("")
logger.info("=" * 80)
logger.info("PHYLOSOLID: PHYLOGENETIC RECONSTRUCTION COMPLETED")
logger.info("=" * 80)
logger.info("")
logger.info("  ┌─────────────────────────────────────────────────────────────────────┐")
logger.info("  │  MUTATION CLASSIFICATION SUMMARY                                   │")
logger.info("  ├─────────────────────────────────────────────────────────────────────┤")
logger.info(f"  │  M_scaffold                     : {len(scaffold_mutations):>10}│")
logger.info(f"  │  M_accessory (integrated)       : {len(attached_mutations):>10}│")
logger.info(f"  │  M_artifact (REMOVED)           : {len(to_be_removed_mutations_by_fp_mutations_cross_all_cells):>10}│")
logger.info(f"  │  M_root (root-assigned)         : {len(root_mutations):>10}│")
logger.info(f"  │  M_ambiguous (conflict)         : {len(all_conflict_mutations):>10}│")
logger.info("  ├─────────────────────────────────────────────────────────────────────┤")
logger.info("  │  CELL CLASSIFICATION SUMMARY                                       │")
logger.info("  ├─────────────────────────────────────────────────────────────────────┤")
logger.info(f"  │  C_resolved (final)             : {final_cleaned_M_full.shape[0]:>10}│")
logger.info(f"  │  C_orphan (REMOVED)             : {len(to_be_removed_cells):>10}│")
logger.info(f"  │  C_chimeric (REMOVED)           : {len(identified_doublet_cells):>10}│")
logger.info("  ├─────────────────────────────────────────────────────────────────────┤")
logger.info("  │  DISCORDANCE METRICS                                               │")
logger.info("  ├─────────────────────────────────────────────────────────────────────┤")
logger.info(f"  │  Omega (pre-QC)                 : {omega_before_qc:>10.2f}      │")
logger.info(f"  │  Omega (final)                  : {omega_final:>10.2f}      │")
logger.info(f"  │  Omega reduction                : {omega_before_qc - omega_final:>10.2f}      │")
logger.info("  └─────────────────────────────────────────────────────────────────────┘")
logger.info("")
logger.info("=" * 80)
logger.info("PhyloSOLID completed successfully!")
logger.info("=" * 80)


# ------------------------------
# End of Process
# ------------------------------
finish_time = time.perf_counter()
print("Program finished in {:.4f} seconds".format(finish_time - start_time))

