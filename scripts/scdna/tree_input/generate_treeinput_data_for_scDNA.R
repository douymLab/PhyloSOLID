# Date: 2023/06/01
# Update: 2024/09/20
# Author: Qing Yang
# Work: Filtering take high-confidence mutations and uses them to build low-resolution but high-confidence phylogenetic tree, then calculate somatic posterior per based on tree, and finally extract and test features.


##### Time #####
import time
start_time = time.perf_counter()


##### Input parameters
import multiprocessing as mp
import argparse
from argparse import ArgumentParser
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--inputpath", default="/storage/douyanmeiLab/yangqing/tools/PhyloMosaicGenie/pmg/data/", type=str, help="The inputpath contains the preprocessing results from raw posterior-reads data.")
parser.add_argument("-o", "--outputpath", default="/storage/douyanmeiLab/yangqing/tools/PhyloMosaicGenie/pmg/features/", type=str, help="The outputpath you want to save the features informatin and low-resolution tree.")
parser.add_argument("-p", "--phylo_or_valid", default="valid", choices=["phylo", "valid"], type=str, help="Select 'phylo' or 'valid' as this parameter to perform the tree inference or somatic mutation verification task.")
parser.add_argument("-n", "--pseudo_node", default="no", choices=["yes", "no"], type=str, help="Select 'yes' or 'no' to identify if there is a pseudo-node that not includes any mutations when you construct a phylogenetic tree.")
parser.add_argument("-w", "--general_weight_na", default=0.001, type=float, help="Set a small value as the weight of the general NA flipping.")
parser.add_argument("-s", "--fnfp_ratio_subclone", default=0.5, type=float, help="The ratio of the false-negative rate to the false-positive rate is used to convert the weight to flipping.")
parser.add_argument("-r", "--fnfp_ratio_late", default=1e-5, type=float, help="The ratio of the false-negative rate to the false-positive rate is used to convert the weight to flipping.") ## cancer or relable: 1e-5; normal:1e-2~0.5
parser.add_argument("-g", "--grnd", default="/storage/douyanmeiLab/yangqing/tools/PhyloMosaicGenie/realData/grnd/tree1465_grnd.CFMatrix", type=str, help="The ground-truth tree file (cellIDxmutID) if the '-p' parameter you selected is 'valid'.")
parser.add_argument("-a", "--alpha_beta", default="/storage/douyanmeiLab/yangqing/tools/PhyloMosaicGenie/pmg/input/UMB1465_alpha_beta_sorted.txt", type=str, help="The alpha_beta values per cell file.")
parser.add_argument("-l", "--labelfile", default="", type=str, help="The 'label.txt' if you have created. Otherwise, please ignore the parameter.")
parser.add_argument("-t", "--threshold", default=0.9, type=float, help="The threshold of somatic posterior per cell used to determine whether a unit is initially identified as a mutation.")
parser.add_argument("-m", "--mut_filter", default="yes", choices=["yes", "no"], type=str, help="Select 'yes' or 'no' to determine whether to filter somatic mutations.")
parser.add_argument("-u", "--is_update_somatic_posterior", default="no", choices=["yes", "no"], type=str, help="Select 'yes' or 'no' to determine whether to update somatic posterior per site based on tree.")
parser.add_argument("-j", "--pass_tree_cutoff", default=0.9, type=float, help="The cutoff of somatic posterior per site based on tree whether a site pass phylogeny.")
parser.add_argument("-k", "--unpass_tree_cutoff", default=0.1, type=float, help="The cutoff of somatic posterior per site only in one cell based on tree whether a site pass phylogeny.")
parser.add_argument("-c", "--correct_fp_manually", default="no", type=str, help="The file contains cellID and mutID you want to manually false positive flip to 0. If this parameter is not set, we will default to 'no' to avoid manual correction")
parser.add_argument("-e", "--extra_NoneCluster_cutoff", default="no", type=str, help="A coverage cutoff of an additional cluster with high coverage but no mutant alleles can be generated. If this parameter is not set, the default is set to 'no', then no such additional cluster is set")
parser.add_argument("-v", "--is_log_value_for_likelihoods", default="yes", choices=["yes", "no"], type=str, help="The value format of your likelihoods.")
# parser.add_argument("-t", "--thread", default=mp.cpu_count(), type=int, help="Cpu count.")
args = parser.parse_args()

# Display parameters for verification
print("Parameters:")
print("Input path: ", args.inputpath)
print("Output path: ", args.outputpath)
print("Phylo or valid: ", args.phylo_or_valid)
print("Pseudo node: ", args.pseudo_node)
print("The General weight_na ratio for pivot mutations: ", args.general_weight_na)
print("The FN/FP ratio for pivot mutations: ", args.fnfp_ratio_subclone)
print("The FN/FP ratio for late mutations: ", args.fnfp_ratio_late)
print("Ground truth file: ", args.grnd)
print("Alpha-beta file: ", args.alpha_beta)
print("Label file: ", args.labelfile)
print("Threshold: ", args.threshold)
print("Mutation filter: ", args.mut_filter)
print("Pass tree cutoff: ", args.pass_tree_cutoff)
print("Unpass tree cutoff: ", args.unpass_tree_cutoff)
# Uncomment if needed
# print("Thread count: ", args.thread)


##### Load libraies
import os
import math
import random
import numpy as np
import pandas as pd
from tqdm import tqdm
import networkx as nx
import scipy.stats as stats
from scipy.stats import mode
from scipy.special import logsumexp, comb, perm, beta
from statsmodels.stats.proportion import proportion_confint
import itertools
import ast
import re
from itertools import chain, combinations
from collections import Counter
from collections import defaultdict
from multiprocessing import Process, Queue
import matplotlib.pyplot as plt
import seaborn as sns
import scphylo as scp
from scphylo.pl._helper import (
    _add_barplot,
    _add_chromplot,
    _clonal_cell_mutation_list,
    _get_tree,
    _newick_info2_mutation_list,
)
import matplotlib
matplotlib.use('Agg')  # Use a non-interactive backend
import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning, message='divide by zero encountered in log')
warnings.filterwarnings("ignore", category=FutureWarning)


##### Parameters
inputpath = args.inputpath
outputpath = args.outputpath
if not os.path.exists(outputpath):
    os.makedirs(outputpath)

fnfp_ratio_subclone = args.fnfp_ratio_subclone
fnfp_ratio_late = args.fnfp_ratio_late
cutoff = args.threshold
general_weight_na = args.general_weight_na

##### Read files
in_posterior = pd.read_csv(inputpath+"/data.posterior_matrix.txt", sep='\t', index_col=0).T
in_llmut = pd.read_csv(inputpath+"/data.likelihood_mut_matrix.txt", sep='\t', index_col=0).T
in_llunmut = pd.read_csv(inputpath+"/data.likelihood_unmut_matrix.txt", sep='\t', index_col=0).T
in_reads = pd.read_csv(inputpath+"/data.allele_count.txt", sep='\t', index_col=0).T
in_baseq = pd.read_csv(inputpath+"/data.baseq_dict.txt", sep='\t', index_col=0).T
in_features = pd.read_csv(inputpath+"/features.preprocess_items.txt", sep='\t', index_col=0).T

in_alpha_beta = pd.read_csv(args.alpha_beta, sep='\t')

# input data display
in_posterior
in_llmut
in_llunmut
in_reads
in_baseq
in_features
in_alpha_beta


##### Functions
def log_dataframe_info(df, df_name):
    # Function to log basic information about the dataframe
    print(f" ===== {df_name} Information ===== ")
    print(f"Shape: {df.shape}")
    # print(f"Columns: {df.columns.tolist()}")
    print(f"Missing values: {df.isnull().sum().sum()}")
    print("\n")

def divide_columns(row):
    if pd.notna(row['total']) and row['total'] != 0:
        return int(row['alt']) / int(row['total'])
    else:
        return np.nan

def reads2df(reads):
    df = pd.DataFrame(reads, columns=['data'])
    df['alt'] = df['data'].str.split('/').str[0]
    df['total'] = df['data'].str.split('/').str[1]
    df['af'] = df.apply(divide_columns, axis=1)
    df['is_mut'] = pd.to_numeric(df['alt']).apply(lambda x: 1 if x >= 1 else 0)
    df['is_mut'] = df['is_mut'].fillna(0)
    return df    

def get_binaryVec(vector, df_reads_persite, cutoff):
    more_cutoff = np.where(vector >= cutoff, 1, 0)
    mut_info = df_reads_persite['is_mut'].to_numpy()
    subclone = np.logical_and(more_cutoff, mut_info).astype(int)
    return subclone

def get_ternaryVec_NAto3(vector, df_reads_persite, cutoff):
    more_cutoff = np.where((vector >= cutoff) & (vector <= 1), 1, np.where(vector == 3, 3, 0))
    mut_info = df_reads_persite['is_mut'].to_numpy()
    subclone = np.where((more_cutoff == 1) & (mut_info == 1), 1,
        np.where((more_cutoff == 0) & (mut_info == 1), 0,
            np.where((more_cutoff == 1) & (mut_info == 0), 0,
                np.where((more_cutoff == 3) | (mut_info == 3), 3, 0))))
    return subclone

def get_binaryVec_by_cutoff(vector, cutoff):
    more_cutoff = np.where(vector >= cutoff, 1, 0)
    return more_cutoff

def get_binaryVec_by_mutAllele(df_reads):
    mut_info = df_reads['is_mut'].to_numpy()
    return mut_info

def posterior2bin(df_posterior, df_reads, cutoff):
    M_posterior = df_posterior.values
    M_reads = df_reads.iloc[1:,].values
    # Get the primarily binary martrix
    M_bin = []
    for i, v in enumerate(M_posterior.T):
        M_bin.append(get_binaryVec(v, reads2df(M_reads[:,i]), cutoff))
    M_bin = np.array(M_bin).T
    df_bin = pd.DataFrame(M_bin)
    df_bin.columns = df_posterior.columns
    df_bin.index = df_posterior.index
    df_bin.index.name = "cellIDxmutID"
    return [M_posterior, M_reads, M_bin, df_bin]

def posterior2ter_NAto3(df_posterior, df_reads, cutoff):
    M_posterior = df_posterior.values
    M_reads = df_reads.iloc[1:,].values
    # Get the primarily binary martrix
    M_bin = M_posterior.copy()
    M_bin[np.isnan(M_bin)] = 3
    for i, v in enumerate(M_bin.T):
        M_bin[:,i] = get_ternaryVec_NAto3(v, reads2df(M_reads[:,i]), cutoff)
    df_bin = pd.DataFrame(M_bin).astype(int)
    df_bin.columns = df_posterior.columns
    df_bin.index = df_posterior.index
    df_bin.index.name = "cellIDxmutID"
    return [M_posterior, M_reads, M_bin, df_bin]

def posterior2bin_bothPosteriorMutallele(df_posterior, df_reads, cutoff):
    M_posterior = df_posterior.values
    M_reads = df_reads.iloc[1:,].values
    # Get the primarily binary martrix
    M_bin = []
    for i, v in enumerate(M_posterior.T):
        M_bin.append(get_binaryVec(v, reads2df(M_reads[:,i]), cutoff))
    M_bin = np.array(M_bin).T
    df_bin = pd.DataFrame(M_bin)
    df_bin.columns = df_posterior.columns
    df_bin.index = df_posterior.index
    df_bin.index.name = "cellIDxmutID"
    return df_bin

def posterior2ter_NAto3_bothPosteriorMutallele(df_posterior, df_reads, cutoff):
    M_posterior = df_posterior.values
    M_reads = df_reads.iloc[1:,].values
    # Get the primarily binary martrix
    M_bin = M_posterior.copy()
    M_bin[np.isnan(M_bin)] = 3
    for i, v in enumerate(M_bin.T):
        M_bin[:,i] = get_ternaryVec_NAto3(v, reads2df(M_reads[:,i]), cutoff)
    df_bin = pd.DataFrame(M_bin).astype(int)
    df_bin.columns = df_posterior.columns
    df_bin.index = df_posterior.index
    df_bin.index.name = "cellIDxmutID"
    return df_bin

def posterior2bin_onlyPosterior(df_posterior, cutoff):
    M_posterior = df_posterior.values
    # Get the primarily binary martrix
    M_bin = []
    for i, v in enumerate(M_posterior.T):
        M_bin.append((v >= cutoff).astype(int))
    M_bin = np.array(M_bin).T
    df_bin = pd.DataFrame(M_bin)
    df_bin.columns = df_posterior.columns
    df_bin.index = df_posterior.index
    df_bin.index.name = "cellIDxmutID"
    return df_bin

def posterior2ter_NAto3_onlyPosterior(df_posterior, cutoff):
    M_posterior = df_posterior.values
    # Get the primarily binary martrix
    M_bin = M_posterior.copy()
    M_bin[np.isnan(M_bin)] = 3
    for i, v in enumerate(M_bin.T):
        M_bin[:,i] = np.where((v >= cutoff) & (v <= 1), 1, np.where(v == 3, 3, 0))
    df_bin = pd.DataFrame(M_bin).astype(int)
    df_bin.columns = df_posterior.columns
    df_bin.index = df_posterior.index
    df_bin.index.name = "cellIDxmutID"
    return df_bin

def filter_preliminary(value):
    if pd.isna(value):
        return False
    numerator, denominator = map(int, value.split('/'))
    return numerator / denominator >= 0.2 and numerator >= 3

def data_cleaning(in_posterior, in_reads, in_features, cutoff, is_filter_mut):
    ##### Step1: Pre-filtering for mutations (that appear in all or none of the cells at the same time)
    [M_posterior_raw, M_reads_raw, M_bin_raw, df_bin_raw] = posterior2bin(in_posterior, in_reads, cutoff)
    df_bin_raw.to_csv(outputpath+"/df_bin_raw.txt", sep="\t")
    # df_bin_raw.to_csv(outputpath+"/df_bin_raw.txt", sep="\t")
    ### Filter1: Delete the columns that all values equal 0 or 1
    all_zeros_indices = np.where((M_bin_raw == 0).all(axis=0))[0]
    all_zeros_mutid = list(in_posterior.columns[all_zeros_indices])
    all_ones_indices = np.where((M_bin_raw == 1).all(axis=0))[0]
    all_ones_mutid = list(in_posterior.columns[all_ones_indices])
    remain_mutid = [mutid for mutid in in_posterior.columns if mutid not in (all_zeros_mutid+all_ones_mutid)]
    no_zeros_ones_posterior = in_posterior.loc[:, remain_mutid]
    no_zeros_ones_reads = in_reads.loc[:, remain_mutid]
    no_zeros_ones_features = in_features.loc[:, remain_mutid]
    # no_zeros_ones_posterior.to_csv(outputpath+"/no_zeros_ones_posterior.txt", sep="\t")
    # no_zeros_ones_reads.to_csv(outputpath+"/no_zeros_ones_reads.txt", sep="\t")
    # no_zeros_ones_features.to_csv(outputpath+"/no_zeros_ones_features.txt", sep="\t")
    [M_posterior_no_zeros_ones, M_reads_no_zeros_ones, M_bin_no_zeros_ones, df_bin_no_zeros_ones] = posterior2bin(no_zeros_ones_posterior, no_zeros_ones_reads, cutoff)
    # df_bin_no_zeros_ones.to_csv(outputpath+"/df_bin_no_zeros_ones.txt", sep="\t")
    print("The number of candidate sites exclude all zeros or all ones for cells is: "+str(len(remain_mutid)))
    ##### Step2: Filtering for mutations
    if is_filter_mut=="yes":    
        ### Filter2: Exclude low-confidence mutations by permissive conditions
        ## clean1: At least one cell
        preliminarily_filtered = no_zeros_ones_reads.iloc[1:, :].apply(lambda x: x.map(filter_preliminary)).any()
        print("The number of sites that passed the preliminarily filter for low confidence mutations is: "+str(preliminarily_filtered.sum()))
        # The number of sites that passed the preliminarily filter for low confidence mutations is: 123
        ## clen2: The proportion of cells screened as NA
        na_proportion_persite = no_zeros_ones_posterior.apply(lambda col: col.isna().mean())
        # na_proportion_persite = no_zeros_ones_reads[no_zeros_ones_reads.index != 'bulk'].apply(lambda x: pd.to_numeric(x.str.split('/').str[1], errors='coerce')).isna().mean()
        # na_proportion_persite = no_zeros_ones_posterior.apply(lambda col: col.isna().sum() / len(col))
        no_zeros_ones_features.loc['na_proportion_persite'] = na_proportion_persite
        mean_na_proportion = np.mean(na_proportion_persite)
        std_na_proportion = np.std(na_proportion_persite)
        # This takes the mean plus twice the standard deviation as the upper threshold
        na_upper_threshold = mean_na_proportion + 1.5 * std_na_proportion
        na_filtered = no_zeros_ones_features.T["na_proportion_persite"]<=na_upper_threshold
        print("The number of sites that passed the NA proportion filter for low confidence mutations is: "+str(na_filtered.sum()))
        # The number of sites that passed the NA proportion filter for low confidence mutations is: 363
        ## clean3: filter the outliers of avg_cov_percell
        cov_lower_threshold = np.percentile(no_zeros_ones_features.loc['avg_cov_percell'], 50)
        cov_filtered = no_zeros_ones_features.T["avg_cov_percell"]>=cov_lower_threshold
        print("The number of sites that passed the avg_cov_percell filter for low confidence mutations is: "+str(cov_filtered.sum()))
        # The number of sites that passed the avg_cov_percell filter for low confidence mutations is: 289
        ## clean4: As a rule of thumb, avg_mutAF_percell>=0.1 && max_mutAF_percell>=0.3
        mutAF_filtered = ((no_zeros_ones_features.T["avg_mutAF_percell"]>=0.3) & (no_zeros_ones_features.T["max_mutAF_percell"]>=0.5))
        # Calculate the outlier of avg_cov_percell
        print("The number of sites that passed the further filter for low confidence mutations is: "+str(mutAF_filtered.sum()))
        # The number of sites that passed the further filter for low confidence mutations is: 163
        ##### Mutations filtering results
        filtered_sites = no_zeros_ones_features.T[preliminarily_filtered & na_filtered & cov_filtered & mutAF_filtered].index.tolist()
        print("The number of sites that passed the two steps filter for low confidence mutations is: "+str(len(filtered_sites)))
        # The number of sites that passed the two steps filter for low confidence mutations is: 85
        if len(filtered_sites)==0:
            print("FilterError: All of your mutations are filtered, please check your input or change your filtering parameters such as setting 'mut_filter' as 'no', or something else.")
        else:
            filtered_posterior = no_zeros_ones_posterior[filtered_sites]
            filtered_reads = no_zeros_ones_reads[filtered_sites]
            filtered_features = no_zeros_ones_features[filtered_sites]
            [filtered_posterior_input, filtered_reads_input, filtered_bin, df_filtered_bin] = posterior2bin(filtered_posterior, filtered_reads, cutoff)
    elif is_filter_mut=="no":
        ##### Don't filter any mutations
        filtered_sites = no_zeros_ones_features.columns.tolist()
        filtered_posterior = no_zeros_ones_posterior
        filtered_reads = no_zeros_ones_reads
        filtered_features = no_zeros_ones_features
        [filtered_posterior_input, filtered_reads_input, filtered_bin, df_filtered_bin] = posterior2bin(filtered_posterior, filtered_reads, cutoff)
    # Save the intermediate file
    if len(filtered_sites)!=0:
        filtered_posterior.to_csv(outputpath+"/filtered_posterior.txt", sep="\t")
        filtered_reads.to_csv(outputpath+"/filtered_reads.txt", sep="\t")
        filtered_features.to_csv(outputpath+"/filtered_features.txt", sep="\t")
        df_filtered_bin.to_csv(outputpath+"/filtered_bin.txt", sep="\t")
        print("The number of candidate sites for building pivot tree is: "+str(len(filtered_sites)))
        ##### Step3: Filtering for cells
        ### Additional subclone named as "no_mut_subclone = "
        no_mut_cellidx = np.where(np.all(filtered_bin == 0, axis=1))[0]
        no_mut_cellname = filtered_posterior.index[no_mut_cellidx]
        # no_mut_subclone = pd.DataFrame(np.zeros((len(no_mut_cellname), len(filtered_sites)+1)), index=no_mut_cellname, columns=filtered_sites+['no_mut_subclone'])
        # no_mut_subclone.iloc[:, -1] = 1
        ### Pending processing matrix consists of cells that contain at least one mutation
        pending_cellidx = np.where(np.any(filtered_bin == 1, axis=1))[0]
        pending_cellname = filtered_posterior.index[pending_cellidx]
        pending_posterior = filtered_posterior.loc[pending_cellname,:]
        pending_reads = filtered_reads.loc[np.insert(pending_cellname, 0, 'bulk'),:]
        pending_features = filtered_features
        [pending_posterior_input, pending_reads_input, pending_bin, df_pending_bin] = posterior2bin(pending_posterior, pending_reads, cutoff)
        # [pending_posterior_input, pending_reads_input, df_pending_bin] = posterior2ter_NAto3(pending_posterior, pending_reads, cutoff)
        # Check to see if there are cells do not have any mutations
        return [filtered_sites, pending_posterior, pending_reads, pending_features, pending_bin, no_zeros_ones_posterior, no_zeros_ones_reads, no_zeros_ones_features, no_mut_cellname]

def intersect_is_self(v1, v0):
    # If True, v1 is a subset of v0
    v0_indices = [i for i, val in enumerate(v0) if val == 1]
    v1_indices = [i for i, val in enumerate(v1) if val == 1]
    v0_set = set(v0_indices)
    v1_set = set(v1_indices)
    if v1_set == v1_set.intersection(v0_set):
        return True
    else:
        return False

def intersect_is_empty(v1, v0):
    v0_indices = [i for i, val in enumerate(v0) if val == 1]
    v1_indices = [i for i, val in enumerate(v1) if val == 1]
    v0_set = set(v0_indices)
    v1_set = set(v1_indices)
    if len(v1_set.intersection(v0_set)) == 0:
        return True
    else:
        return False

def intersect_is_conflict(v1, v0):
    v0_indices = [i for i, val in enumerate(v0) if val == 1]
    v1_indices = [i for i, val in enumerate(v1) if val == 1]
    v0_set = set(v0_indices)
    v1_set = set(v1_indices)
    if len(v1_set.intersection(v0_set)) != 0 and len(v1_set.difference(v0_set)) != 0:
        if intersect_is_self(v0, v1):
            return False
        else:
            return True
    else:
        return False

def find_conflicting_mutations(target_mut, in_bin_withoutNA):
    conflict_muts = []
    for mut in in_bin_withoutNA.columns:
        if mut == target_mut:
            continue
        if intersect_is_conflict(in_bin_withoutNA[target_mut].values, in_bin_withoutNA[mut].values):
            conflict_muts.append(mut)
    return conflict_muts, len(conflict_muts)

def find_all_conflicting_mutations(in_bin_withoutNA):
    conflict_dict = {}
    for target_mut in in_bin_withoutNA.columns:
        conflict_muts, conflict_count = find_conflicting_mutations(target_mut, in_bin_withoutNA)
        conflict_dict[target_mut] = {"conflicts": conflict_muts, "count": conflict_count}
    return conflict_dict

def convert_conflict_results_to_dataframe(conflict_results):
    df = pd.DataFrame.from_dict(conflict_results, orient='index').reset_index()
    df.columns = ['target_mut', 'conflicts', 'count']
    df['conflicts'] = df['conflicts'].apply(lambda x: ', '.join(x))
    return df

def find_intersected_vectors(subclone_vectors, bin_vec):
    intersected_vectors = []
    intersected_indices = []
    for idx, sub_vec in enumerate(subclone_vectors):
        intersection = np.logical_and(sub_vec, bin_vec)
        if any(intersection):
            intersected_vectors.append(sub_vec)
            intersected_indices.append(idx)
    return [intersected_vectors, intersected_indices]

def find_max_sum_subvector(subclone_mutlist, subclone_binV_list):
    max_sum = -np.inf
    max_sum_index = None
    for i, subvector in enumerate(subclone_binV_list):
        subvector_sum = sum(subvector)
        if subvector_sum > max_sum:
            max_sum = subvector_sum
            max_sum_index = i
    return [subclone_binV_list[max_sum_index], subclone_mutlist[max_sum_index]]

def find_intersecting_subclone_and_leader(df_subclone, bin_vec):
    max_mutid = None
    max_sum = -np.inf
    intersecting_vector = None
    intersect_info_list = []
    intersect_all_mutid = []
    intersect_all_binV = []
    intersect_df_idxs = []
    for idx, row in df_subclone.iterrows():
        subclone_vector = np.array(row['subclone_vector'])
        subclone_mutlist = row['subclone_mutlist']
        subclone_binV_list = row['subclone_binV_list']
        if sum(subclone_vector * bin_vec)>0:
            [max_binV, max_mutid] = find_max_sum_subvector(subclone_mutlist, subclone_binV_list)
            intersect_info = [subclone_vector, max_binV, max_mutid]
            intersect_info_list.append(intersect_info)
            intersect_all_mutid = intersect_all_mutid+subclone_mutlist
            intersect_all_binV = intersect_all_binV+subclone_binV_list
            intersect_df_idxs = intersect_df_idxs+[idx]
    return [intersect_info_list, intersect_all_mutid, intersect_all_binV, intersect_df_idxs]

def find_intersecting_vector_in_existing(existing_vectors, bin_vec, already_processed_indices):
    intersect_info_list = []
    for idx, existing_vector in enumerate(existing_vectors):
        intersection_num = sum(existing_vector * bin_vec)
        if intersection_num > 0:
            mutid_in_matrix = already_processed_indices[idx]
            intersect_info = [existing_vector, mutid_in_matrix, intersection_num]
            intersect_info_list.append(intersect_info)
    return intersect_info_list

def get_conflict_entries_from_all_intersetct_entries(intersect_entries, pivot_bin):
    conflict_entries = []
    for considered_entry in intersect_entries:
        if intersect_is_conflict(considered_entry[0], pivot_bin):
            conflict_entries.append(considered_entry)
    return conflict_entries

def calculate_difference_sum_withNA(vec_posterior, vec_binary, cutoff, fnfp_ratio, general_weight_na):
    # 计算 weight_na_to_1 和 weight_na_to_0
    # general_weight_na = 1/len(vec_posterior)
    # general_weight_na = 0.01
    vec_fixed = np.where(vec_posterior >= cutoff, 1, np.where(vec_posterior < cutoff, 0, np.nan))
    non_nan_values = vec_fixed[~np.isnan(vec_fixed)]
    count_1 = np.sum(non_nan_values == 1)
    count_0 = np.sum(non_nan_values == 0)
    total_count = count_1 + count_0
    weight_na_to_1 = count_1 / total_count if total_count > 0 else 0
    weight_na_to_1 = general_weight_na * weight_na_to_1
    weight_na_to_0 = count_0 / total_count if total_count > 0 else 0
    weight_na_to_0 = general_weight_na * weight_na_to_0
    # 计算 penalty
    fixed_bin = (vec_posterior >= cutoff).astype(int)
    non_na_positions = ~np.isnan(vec_posterior)
    value_fn_positions = np.where((fixed_bin == 0) & (vec_binary == 1) & non_na_positions)[0]
    value_fp_positions = np.where((fixed_bin == 1) & (vec_binary == 0) & non_na_positions)[0]
    na_positions = np.where(np.isnan(vec_posterior))[0]
    # 计算各类 penalty
    na_to_1_penalty = weight_na_to_1 * np.sum(vec_binary[na_positions] == 1)
    na_to_0_penalty = weight_na_to_0 * np.sum(vec_binary[na_positions] == 0)
    fp_penalty = np.sum(np.abs(vec_posterior[value_fp_positions] - vec_binary[value_fp_positions]))
    fn_penalty = fnfp_ratio * np.sum(np.abs(vec_posterior[value_fn_positions] - vec_binary[value_fn_positions]))
    # 总 penalty
    penalty_score = na_to_1_penalty + na_to_0_penalty + fp_penalty + fn_penalty
    return penalty_score

# calculate_difference_sum_withNA(np.array([np.nan,0.1,1,0.9,0.9,0]), np.array([1,0,1,1,0,0]), cutoff=0.5, fnfp_ratio=1e-5)
# 1.5

def calculate_difference_sum_withNA_multiV(intersect_mutid_list, df_posterior, vec_binary, cutoff, fnfp_ratio, general_weight_na):
    penalty_multiV = 0
    for intersect_mutid in intersect_mutid_list:
        penalty_eachV = calculate_difference_sum_withNA(df_posterior.iloc[:,intersect_mutid].values, vec_binary, cutoff, fnfp_ratio, general_weight_na)
        penalty_multiV = penalty_multiV + penalty_eachV
    return penalty_multiV

def calculate_difference_sum_withNA_multiV_under_intersect_subset(intersect_mutid_list, df_posterior, vec_binary, cutoff, fnfp_ratio, M_copy, general_weight_na):
    penalty_multiV = 0
    for intersect_mutid in intersect_mutid_list:
        new_vector = ((M_copy[:,intersect_mutid]*vec_binary)>0).astype(int)
        penalty_eachV = calculate_difference_sum_withNA(df_posterior.iloc[:,intersect_mutid].values, new_vector, cutoff, fnfp_ratio, general_weight_na)
        penalty_multiV = penalty_multiV + penalty_eachV
    return penalty_multiV

def calculate_difference_sum_withNA_multiV_under_conflict_fp(intersect_mutid_list, df_posterior, vec_binary, cutoff, fnfp_ratio, M_copy, general_weight_na):
    penalty_multiV = 0
    for intersect_mutid in intersect_mutid_list:
        new_vector = ((M_copy[:,intersect_mutid]-vec_binary)>0).astype(int)
        penalty_eachV = calculate_difference_sum_withNA(df_posterior.iloc[:,intersect_mutid].values, new_vector, cutoff, fnfp_ratio, general_weight_na)
        penalty_multiV = penalty_multiV + penalty_eachV
    return penalty_multiV

def write_subclone_persite(M_bin, mut_index, bin_vec, df_posterior, df, cutoff, fnfp_ratio_subclone, general_weight_na):
    is_germ = False
    M_copy = M_bin.copy()
    my_df = df.copy()
    new_binV = bin_vec
    if my_df.empty:
        mutlist = [mut_index]
        binV_list = [bin_vec]
        my_df.loc[len(my_df)] = [bin_vec, mutlist, binV_list, [mut_index]]
    else:
        subclone_vectors = list(my_df['subclone_vector'].values)
        if all(np.apply_along_axis(lambda v: intersect_is_empty(bin_vec, v), axis=1, arr=subclone_vectors)):
            # Judgment: The intersection of the add vector with all the current child clones is empty.
            mutlist = [mut_index]
            binV_list = [bin_vec]
            my_df.loc[len(my_df)] = [bin_vec, mutlist, binV_list, [mut_index]]
        elif any(np.apply_along_axis(lambda v: intersect_is_self(bin_vec, v), axis=1, arr=subclone_vectors)):
            # Judgment: The adding vector is subset of one of present subclones.
            vector_to_find = [v for v in subclone_vectors if intersect_is_self(bin_vec, v)][0]
            df_i = [df_i for df_i,row in my_df.iterrows() if np.array_equal(vector_to_find, row['subclone_vector'])][0]
            my_df.loc[df_i]['subclone_mutlist'].append(mut_index)
            my_df.loc[df_i]['subclone_binV_list'].append(bin_vec)
            # => subclone_leader wouldn't change
        elif any(np.apply_along_axis(lambda v: intersect_is_self(v, bin_vec), axis=1, arr=subclone_vectors)):
            # Judgment: The adding vector is superset of one of present subclones.
            df_i_find = []
            for vector_to_find in [v for v in subclone_vectors if intersect_is_self(v, bin_vec)]:
                df_i = [df_i for df_i,row in my_df.iterrows() if np.array_equal(vector_to_find, row['subclone_vector'])][0]
                df_i_find.append(df_i)
            # Replace entries higher in the index
            find_i = min(df_i_find)
            my_df.at[find_i, 'subclone_vector'] = bin_vec
            my_df.loc[find_i]['subclone_mutlist'].append(mut_index)
            my_df.loc[find_i]['subclone_binV_list'].append(bin_vec)
            my_df.loc[find_i]['subclone_leader'] = [mut_index]
            df_i_find.remove(find_i)
            for remove_i in df_i_find:
                for mv in my_df.loc[remove_i]['subclone_mutlist']:
                    my_df.loc[find_i]['subclone_mutlist'].append(mv)
                    my_df.loc[find_i]['subclone_binV_list'].append(bin_vec)
                my_df = my_df.drop(remove_i)
        else:
            # Judgment: The addition vector conflicts with one of the current child clones => flipping
            penalty_score = np.inf
            refine_subclone = None
            [intersect_entries, intersect_all_mutid, intersect_all_binV, intersect_df_idxs] = find_intersecting_subclone_and_leader(my_df, bin_vec)
            intersect_mutid_list = [entry[2] for entry in intersect_entries]
            ### Case1: remain previours vectors as subclone
            for considered_entry in intersect_entries:
                refine_vec = considered_entry[0]
                # np.array_equal(considered_entry[0], considered_entry[1])
                new_penalty = calculate_difference_sum_withNA(df_posterior.iloc[:, considered_entry[2]].values, refine_vec*considered_entry[1], cutoff, fnfp_ratio_subclone, general_weight_na) + calculate_difference_sum_withNA(df_posterior.iloc[:, mut_index].values, refine_vec*bin_vec, cutoff, fnfp_ratio_subclone, general_weight_na)
                if new_penalty < penalty_score:
                    penalty_score = new_penalty
                    refine_subclone = refine_vec
                    refine_case = "case1: previous"
            #### Case2: adopt current vector as subclone
            refine_vec = bin_vec
            new_penalty = calculate_difference_sum_withNA(df_posterior.iloc[:, mut_index].values, refine_vec, cutoff, fnfp_ratio_subclone, general_weight_na) + calculate_difference_sum_withNA_multiV(intersect_mutid_list, df_posterior, refine_vec, cutoff, fnfp_ratio_subclone, general_weight_na)
            if new_penalty < penalty_score:
                penalty_score = new_penalty
                refine_subclone = refine_vec
                refine_case = "case2: current"
            ### Case3: adopt their union as subclone
            refine_vec = (bin_vec + np.sum([entry[0] for entry in intersect_entries], axis=0)>0).astype(int)
            new_penalty = calculate_difference_sum_withNA(df_posterior.iloc[:, mut_index].values, refine_vec*bin_vec, cutoff, fnfp_ratio_subclone, general_weight_na) + calculate_difference_sum_withNA_multiV(intersect_mutid_list, df_posterior, refine_vec, cutoff, fnfp_ratio_subclone, general_weight_na)
            if new_penalty < penalty_score:
                penalty_score = new_penalty
                refine_subclone = refine_vec
                refine_case = "case3: extending"
            ### Refine df_subclone
            if penalty_score > theta:
                M_copy[:,mut_index] = np.ones(len(bin_vec), dtype=int)
                is_germ = True
            else:
                if refine_case=="case1: previous":
                    M_copy[:,mut_index] = refine_subclone * bin_vec
                    add_index = [df_i for df_i,row in my_df.iterrows() if np.array_equal(refine_subclone, row['subclone_vector'])][0]
                    my_df.loc[add_index]['subclone_mutlist'].append(mut_index)
                    my_df.loc[add_index]['subclone_binV_list'].append(refine_subclone * bin_vec)
                elif refine_case=="case2: current":
                    for i in intersect_all_mutid:
                        M_copy[:,intersect_entries[i][2]] = intersect_entries[i][0]*refine_subclone
                    my_df.drop(intersect_df_idxs, inplace = True)
                    mutlist = [mut_index]
                    binV_list = [bin_vec]
                    my_df.loc[len(my_df)] = [refine_subclone, mutlist, binV_list, [mut_index]]
                elif refine_case=="case3: extending":
                    M_copy[:,mut_index] = refine_subclone * bin_vec
                    mutlist = intersect_all_mutid + [mut_index]
                    binV_list = intersect_all_binV + [bin_vec]
                    my_df.drop(intersect_df_idxs, inplace = True)
                    my_df.reset_index(drop=True, inplace=True)
                    max_pre_index = [intersect_entries[np.argmax([np.sum(e[0]) for e in intersect_entries])][2]]
                    max_pre_sum = sum(intersect_entries[np.argmax([np.sum(e[0]) for e in intersect_entries])][0])
                    rather_leader = [mut_index] if sum(bin_vec) > max_pre_sum else max_pre_index
                    my_df.loc[len(my_df)] = [refine_subclone, mutlist, binV_list, rather_leader]
    my_df.reset_index(drop=True, inplace=True)
    return [my_df, M_copy, is_germ]

def check_subclone(my_dict):
    # Check whether all subclones cover all cells, and can identify cells that is not be coverd so that we can get the sites that need to be re-process.
    check_vector = np.zeros(len(next(iter(my_dict.values()))), dtype=int)
    for idx, v in my_dict.items():
        check_vector = np.logical_or(check_vector, v)
    check_vector = check_vector.astype(int)
    return check_vector

def sortIndex(M_bin, M_posterior, cutoff):
    sum_binary_0 = np.sum((M_posterior < cutoff) & ~np.isnan(M_posterior), axis=0)
    sum_binary_1 = np.sum((M_posterior >= cutoff) & ~np.isnan(M_posterior), axis=0)
    sum_binary = np.where((sum_binary_1+sum_binary_0) == 0, 0, sum_binary_1 / (sum_binary_1+sum_binary_0))
    sum_posterior = np.nansum(M_posterior, axis=0)
    ISet=[]
    for sum_unique in sorted(list(set(sum_binary)), reverse=True):
        indices = [index for index, value in enumerate(sum_binary) if value == sum_unique]
        sub_ISet = [indices[i] for i in np.argsort(-np.argsort(sum_posterior[indices]))]
        ISet = ISet + sub_ISet
    return ISet

def get_early_mutations(df_features):
    # The number in the top 90% of all possible cell numbers.
    cutoff_early_mut = np.percentile(list(set(df_features.loc['mutant_cellnum',:].astype(int))), 90)
    early_muts_columns = df_features.loc['mutant_cellnum'] >= cutoff_early_mut
    early_muts = df_features.columns[early_muts_columns].tolist()
    early_muts_index = [index for index, element in enumerate(df_features.columns) if element in early_muts]
    return [early_muts, early_muts_index]

def get_subclone_with_flip(M_bin, df_posterior, ISet, theta, cutoff, fnfp_ratio_subclone, general_weight_na):
    M_posterior = df_posterior.values
    M_copy = M_bin.astype(int).copy()
    df_subclone = pd.DataFrame(columns=["subclone_vector", "subclone_mutlist", "subclone_binV_list", 'subclone_leader'])
    germ_sites = []
    for pivot_index in ISet:
        print(str(pivot_index)+" : "+str(df_posterior.columns.tolist()[pivot_index]))
        pivot_vector=M_posterior[:,pivot_index]
        pivot_bin=M_bin[:,pivot_index]
        [df_subclone, M_copy, is_germ] = write_subclone_persite(M_copy, pivot_index, pivot_bin.astype(int), df_posterior, df_subclone, cutoff, fnfp_ratio_subclone, general_weight_na)
        if is_germ:
            germ_sites.append(pivot_index)
    # Add other infomartion
    M_copy = M_copy.astype(int)
    df_subclone['subclone_key'] = ['subclone_' + str(i) for i in range(1,df_subclone.shape[0]+1)]
    df_subclone['subclone_leader_mutid'] = df_posterior.columns[[i[0] for i in df_subclone['subclone_leader']]]
    ### Check whether the df_subclone contains all cells
    subclone_vectors = list(df_subclone['subclone_vector'].values)
    sums_per_position = np.sum(subclone_vectors, axis=0)
    non_one_indices = list(np.where(sums_per_position != 1)[0])
    if np.all(sums_per_position < 2):
        print("=> df_subclone is conflict-free!")
        if len(non_one_indices)==0:
            print("=> df_subclone contains all cells!")
        else:
            print("=> df_subclone does not contain all cells!")
    else:
        print("=> OMG! df_subclone is conflict!")
    # Convert df_subclone['subclone_vector'] from array to list
    df_subclone['subclone_vector'] = [s.tolist() for s in df_subclone['subclone_vector']]
    df_subclone.reset_index(drop=True, inplace=True)
    return [df_subclone, M_copy, germ_sites]

def find_conflict_pairs(M_input):
    conflict_pairs = []
    for i in range(M_input.shape[1]):
        v1 = M_input[:, i]
        for j in range(i+1, M_input.shape[1]):
            v0 = M_input[:, j]
            if intersect_is_conflict(v1, v0):
                conflict_pairs.append((i, j))
    # print("相互之间存在冲突的列：", conflict_pairs)
    # unique_conflict_pairs = sorted(list(set([item for sublist in conflict_pairs for item in sublist])))
    return conflict_pairs

def split_and_retain_larger_group_adaptive(data):
    # 找到非零元素的索引
    non_zero_indices = np.where(data != 0)[0]
    # 提取非零元素
    non_zero_values = data[non_zero_indices]
    # 计算平均值和标准差
    mean_value = np.mean(non_zero_values)
    std_dev = np.std(non_zero_values)
    # 设定一个更加灵活的阈值
    cutoff_value = mean_value - 0.5 * std_dev  # 使用均值减去 0.5 倍的标准差作为阈值
    # 将较大的一组元素的索引保留为1，较小的组保留为0
    larger_group_indices = non_zero_indices[non_zero_values >= cutoff_value]
    # 创建一个与原始数据相同大小的全零数组
    result = np.zeros_like(data)
    # 将较大组的索引对应的位置设置为1
    result[larger_group_indices] = 1
    return result

def get_conflict_to_one_clone_mutidx_list(one_clone, Sremaining, M_copy):
    conflict_to_one_clone_mutidx_list = []
    for idx in Sremaining:
        if intersect_is_conflict(one_clone, M_copy[:,idx]):
            conflict_to_one_clone_mutidx_list.append(idx)
    return conflict_to_one_clone_mutidx_list

def get_subsets_to_one_clone_mutidx_list(one_clone, Sremaining, M_copy):
    conflict_to_one_clone_mutidx_list = []
    for idx in Sremaining:
        if intersect_is_self( M_copy[:,idx], one_clone):
            conflict_to_one_clone_mutidx_list.append(idx)
    return conflict_to_one_clone_mutidx_list

def check_and_process_conflict(M_out, df_pos, cutoff, fnfp_ratio_late, ISet_process, subclone_vectors, pivot_sites, general_weight_na):
    M_final = M_out.copy()
    ISet_remain = [i for i in ISet_process if i not in pivot_sites]
    while not scp.ul.is_conflict_free_gusfield(pd.DataFrame(M_final)):
        ### Step1: Find conflict sites
        conflict_pairs = find_conflict_pairs(M_final)
        ### Step2: Process conflict
        ISet_conflict_free = []
        for pivot_index in ISet_remain:
            pivot_bin = M_final[:,pivot_index]
            if any(np.apply_along_axis(lambda v: intersect_is_self(pivot_bin, v), axis=1, arr=subclone_vectors)):
                ISet_conflict_free.append(pivot_index)
        print("The number of conflict-free mutations is : " + str(len(ISet_conflict_free)))
        ISet_conflict_free_copy = ISet_conflict_free.copy()
        ISet_conflict = [i for i in ISet_remain if i not in ISet_conflict_free]
        print("After the initial inference, the number of mutations that are still conflict is as follows: " + str(len(ISet_conflict)))
        all_elements = [item for sublist in conflict_pairs for item in sublist]
        element_counts = Counter(all_elements)
        sorted_counts = sorted(element_counts.items(), key=lambda x: x[1], reverse=True)
        process_con = sorted_counts[0]
        process_index = process_con[0]
        conflict_list = []
        for pair in conflict_pairs:
            if process_index in pair:
                other_element = pair[0] if pair[1] == process_index else pair[1]
                conflict_list.append(other_element)
        # process (false negative)
        # row_sums = np.sum(M_final[:, [process_index] + conflict_list], axis=1)
        if len(conflict_list)==1:
            min_score = float('inf')
            # case1: conflict_list follows process_index
            case_bin = M_final[:, process_index]
            conflict_index = conflict_list[0]
            flipping_score = calculate_difference_sum_withNA(df_pos.iloc[:, conflict_index], case_bin, cutoff, fnfp_ratio_late, general_weight_na) + calculate_difference_sum_withNA(df_pos.iloc[:, process_index], case_bin, cutoff, fnfp_ratio_late, general_weight_na)
            if flipping_score<min_score:
                process_case = "case1"
                process_case_bin = case_bin
            # case2: process_index follow conflict_list
            case_bin = (np.sum(M_final[:, conflict_list], axis=1) > 0).astype(int)
            conflict_index = conflict_list[0]
            flipping_score = calculate_difference_sum_withNA(df_pos.iloc[:, conflict_index], case_bin, cutoff, fnfp_ratio_late, general_weight_na) + calculate_difference_sum_withNA(df_pos.iloc[:, process_index], case_bin, cutoff, fnfp_ratio_late, general_weight_na)
            if flipping_score<min_score:
                process_case = "case2"
                process_case_bin = case_bin
            # case3: process_index + conflict_list ()
            # final determin
            if process_case=="case1":
                M_final[:, conflict_index] = M_final[:, conflict_index]*M_final[:, process_index]
            elif process_case=="case2":
                M_final[:, process_index] = M_final[:, conflict_index]*M_final[:, process_index]
        else:
            # if scp.ul.is_conflict_free_gusfield(pd.DataFrame(M_final[:,conflict_list])):
            min_score = float('inf')
            for conflict_index in conflict_list:
                flipping_score = calculate_difference_sum_withNA(df_pos.iloc[:, conflict_index], M_final[:, conflict_index], cutoff, fnfp_ratio_late, general_weight_na) + calculate_difference_sum_withNA(df_pos.iloc[:, process_index], M_final[:, conflict_index], cutoff, fnfp_ratio_late, general_weight_na)
                if flipping_score<min_score:
                    process_case_bin = M_final[:, conflict_index]
            M_final[:,process_index] = process_case_bin*M_final[:, process_index]
    return M_final

def check_conflict(M_out, pivot_index, already_processed_indices):
    """
    检查 M_out[:, pivot_index] 是否与 already_processed_indices 中的每一个向量存在冲突。
    参数:
        M_out: 二维数组，包含要检查的所有向量。
        pivot_index: 整数，指定当前待检查向量的索引。
        already_processed_indices: 列表，包含要与之比较的向量的索引。
    返回:
        bool: 如果存在冲突，返回 True；否则返回 False。
    """
    pivot_vector = M_out[:, pivot_index]
    for idx, existing_vector in enumerate(M_out[:, already_processed_indices].T):
        intersection_num = sum(pivot_vector * existing_vector)
        if intersect_is_conflict(pivot_vector, existing_vector):
            # print(f"{already_processed_indices[idx]}: {intersection_num}: conflict!!!")
            return True
        # else:
        #     print(f"{already_processed_indices[idx]}: {intersection_num}")
    return False

def lowPtreeFNonly(M_bin, df_pos, subclone_vectors): # very greedy algorithm that constructs a ptree matrix from M_inputs by adding 1's. Assumes that there are only false negatives.
    M_pos = df_pos.values
    ISet = sortIndex(M_input, M_pos)
    ### Step2: Process mutations that are conflict-free with any subclone
    M_copy = M_input.copy()
    ISet_process = ISet.copy()
    while len(ISet_process)>0:
        pivot_index=ISet_process[0]          # index of the pivot vector
        pivot_vector=M_copy[:,pivot_index]   # vector used for pivoting the current iteration
        Sremaining=ISet_process.copy()
        cum_vector=np.copy(pivot_vector)     # holds the union vector
        while_cont=1
        while while_cont==1:                 # main loop
            while_cont=0       
            for j in Sremaining:
                cap_j=cum_vector*M_copy[:,j] # intersection of the pivot and the jth column
                if np.any(cap_j):            # continue as long as there is a column having non-empty intersection
                    cum_vector=cum_vector+M_copy[:,j]
                    while_cont=1
                    Sremaining.remove(j)
        M_copy[:,pivot_index]=cum_vector
        ISet_process.remove(pivot_index)
    return M_copy

# def greedyPtreeFNFP_withNA(M_input, df_pos, cutoff, fnfp_ratio_late, df_subclone, ISet_process, pivot_sites, general_weight_na): # Modified Greedy algorithm for dealing with matrices with both NA and false positive values.
#     M_pos = df_pos.values
#     M_copy = M_input.copy()
#     ### Step1: Sort by posterior and binary values
#     ISet_remain = ISet_process.copy()
#     ### Step2: Process mutations that are conflict-free with any subclone
#     while len(ISet_remain)>0:
#         pivot_index = ISet_remain[0]
#         Sremaining = ISet_remain.copy()
#         pivot_bin = M_copy[:,pivot_index]
#         pivot_pos = M_pos[:,pivot_index]
#         if pivot_index in pivot_sites:
#             M_copy[:,pivot_index] = pivot_bin.copy()
#         else:
#             cum_vector = np.copy(pivot_bin)
#             while_cont = 1
#             # The histogram for the union
#             cum_hist = np.zeros(M_input.shape[0])
#             cum_hist_int = np.zeros(M_input.shape[0])
#             cum_flipping_hist = np.zeros(M_input.shape[0])
#             # Continue uniting vectors until no candidate remains
#             while while_cont == 1:
#                 while_cont = 0
#                 for j in Sremaining:
#                     cap_i = pivot_bin*M_copy[:,j]
#                     min_vec_size = np.min([np.sum(M_copy[:,j]),np.sum(pivot_bin)])
#                     cap_size = np.sum(cap_i)
#                     # Check if the columns have a meaningful overlap, any intersection is counted in the histogram (paramter lambda set as 0)
#                     if cap_size/min_vec_size > 0:
#                         cum_hist = cum_hist+np.nan_to_num(M_pos[:,j])
#                         cum_hist_int = cum_hist_int+M_copy[:,j].astype(int)
#                         cum_flipping_hist = cum_flipping_hist+calculate_difference_sum_withNA(M_pos[:,j], pivot_bin, cutoff, fnfp_ratio_late, general_weight_na)
#                         # Found an overlapping vector so keep on going
#                         while_cont = 1
#                         # United vector is removed from the search set
#                         Sremaining.remove(j)
#             # The elements that repeated few times are considered to be false positives (paramter miu set as)
#             # cum_vector = (cum_hist>1) & (cum_hist_int>0)
#             cum_vector = cum_hist_int>0
#             ncap = np.sum(cum_vector*pivot_bin)
#             # clean up the false positives wrt. the established pivot
#             for j in ISet_remain:
#                 # Intersection of union with column j
#                 capj = cum_vector*M_copy[:,j]
#                 capj_score = calculate_difference_sum_withNA(M_pos[:,j], capj, cutoff, fnfp_ratio_late, general_weight_na)
#                 # Difference of column j from the union
#                 difj = M_copy[:,j]>capj
#                 difj_score = calculate_difference_sum_withNA(M_pos[:,j], difj, cutoff, fnfp_ratio_late, general_weight_na)
#                 if capj_score < difj_score:
#                     M_copy[:,j] = capj
#                 else:
#                     M_copy[:,j] = difj
#             # Correct the false negatives in the pivot
#             M_copy[:,pivot_index] = cum_vector.astype(int)
#         # removing the pivot from the search space
#         ISet_remain.remove(pivot_index)
#     ### Step3: Process mutations that are conflict with exsiting subclone (mainly false positive flip)
#     # check and process conflict
#     if scp.ul.is_conflict_free_gusfield(pd.DataFrame(M_copy)):
#         M_out = M_copy.copy()
#     else:
#         subclone_vectors = list(df_subclone['subclone_vector'].values)
#         M_out = check_and_process_conflict(M_copy, df_pos, cutoff, fnfp_ratio_late, ISet_process, subclone_vectors, pivot_index, general_weight_na)
#     return M_copy

def get_selected_combinations(conflict_entries):
    # Single element combination
    single_element_combinations = [[entry] for entry in conflict_entries]
    # Missing a combination of elements
    missing_one_combinations = [list(combo) for combo in combinations(conflict_entries, len(conflict_entries) - 1)]
    # Total element combination
    all_elements_combination = [conflict_entries]
    # Merge all combinations
    selected_combinations = single_element_combinations + missing_one_combinations + all_elements_combination
    return selected_combinations

def greedyPtreeFNFP_withNA(M_input, df_pos, cutoff, fnfp_ratio_late, df_subclone, ISet_process, pivot_sites, general_weight_na): 
    # Modified Greedy algorithm for dealing with matrices with both NA and false positive values.
    refine_case_list = []
    M_pos = df_pos.values
    ##### Step1: Sort by posterior and binary values
    ##### Step2: Process mutations that are conflict-free with any subclone
    M_copy = M_input.copy()
    M_out = M_input.copy()
    # conflict_num_threshold = len(ISet_process)/len(df_subclone)
    conflict_num_threshold = np.inf
    print("====================================================================")
    print("The number of mutations processed below is: " + str(len(ISet_process)))
    already_processed_indices = []
    max_iterations = 10  # Set the maximum number of iterations to prevent dead loops
    for pivot_index in ISet_process:
        print(str(pivot_index) + " : " + str(df_pos.columns.tolist()[pivot_index]))
        pivot_bin = M_out[:, pivot_index]
        pivot_pos = M_pos[:, pivot_index]
        if pivot_index in pivot_sites:
            M_out[:, pivot_index] = pivot_bin.copy()
            refine_case = "case0: pivot"
            print(refine_case)
        else:
            existing_vectors = M_out[:, already_processed_indices].T
            if existing_vectors.shape[0] == 0:
                M_out[:, pivot_index] = pivot_bin.copy()
                refine_case = "case0: continue"
                print(refine_case)
            elif all(np.apply_along_axis(lambda v: not intersect_is_conflict(pivot_bin, v), axis=1, arr=existing_vectors)):
                M_out[:, pivot_index] = pivot_bin.copy()
                refine_case = "case0: merge"
                print(refine_case)
            else:
                has_conflict = check_conflict(M_out, pivot_index, already_processed_indices)
                print(has_conflict)
                iterations = 0  # Initializes the iteration counter
                last_conflict_state = has_conflict  # Record the last conflict status
                while has_conflict and iterations < max_iterations:
                    iterations += 1
                    pivot_bin = M_out[:, pivot_index]
                    penalty_score = np.inf
                    refine_subclone = None
                    intersect_entries = find_intersecting_vector_in_existing(existing_vectors, pivot_bin, already_processed_indices)
                    intersect_mutid_list = [entry[1] for entry in intersect_entries]
                    conflict_entries = get_conflict_entries_from_all_intersetct_entries(intersect_entries, pivot_bin)
                    conflict_mutid_list = [entry[1] for entry in conflict_entries]
                    ##### Handle multiple conflict situations
                    ### Case0: move to top conflict because of too many conflicts with previous muts
                    if len(conflict_entries) >= conflict_num_threshold:
                        M_out[:, pivot_index] = np.array(df_subclone.loc[df_subclone['subclone_mutlist'].apply(lambda x: pivot_index in x), 'subclone_vector'].iloc[0])
                        refine_case = "case0: move2sub"
                        print(refine_case)
                        break
                    else:
                        ### Case1: remain previours vectors as fixed_vector
                        for considered_entry in conflict_entries:
                            refine_vec = considered_entry[0]
                            new_penalty = calculate_difference_sum_withNA(df_pos.iloc[:, considered_entry[1]].values, refine_vec * considered_entry[0], cutoff, fnfp_ratio_late, general_weight_na) + calculate_difference_sum_withNA(df_pos.iloc[:, pivot_index].values, refine_vec * pivot_bin, cutoff, fnfp_ratio_late, general_weight_na)
                            print(new_penalty)
                            if new_penalty < penalty_score:
                                penalty_score = new_penalty
                                refine_subclone = refine_vec
                                refine_case = "case1: previous"
                        #### Case2: adopt current vector as fixed_vector
                        refine_vec = pivot_bin
                        new_penalty = calculate_difference_sum_withNA(df_pos.iloc[:, pivot_index].values, refine_vec, cutoff, fnfp_ratio_late, general_weight_na) + calculate_difference_sum_withNA_multiV(conflict_mutid_list, df_pos, refine_vec, cutoff, fnfp_ratio_late, general_weight_na)
                        print(new_penalty)
                        if new_penalty < penalty_score:
                            penalty_score = new_penalty
                            refine_subclone = refine_vec
                            refine_case = "case2: current"
                        ### Case3: adopt their union as fixed_vector
                        # # all combination conditions by all conflicts
                        # all_combinations = [list(combo) for r in range(1, len(conflict_entries) + 1) for combo in itertools.combinations(conflict_entries, r)]
                        # # select_combination is that includes single element combinations, all element combinations, and combinations that are missing an element
                        # selected_combinations = [
                        #     list(combo) 
                        #     for r in [1, len(conflict_entries) - 1, len(conflict_entries)] 
                        #     for combo in combinations(conflict_entries, r)
                        # ]
                        selected_combinations = get_selected_combinations(conflict_entries)
                        # calculate case and penalty 
                        sub_penalty_each_conflict = np.inf
                        for combo_entries in selected_combinations:
                            each_refine_vec = (pivot_bin + np.sum([entry[0] for entry in combo_entries], axis=0)>0).astype(int)
                            each_extend_conflict_mutid = [i[1] for i in combo_entries]
                            each_fp_conflict_mutid = [i for i in conflict_mutid_list if i not in each_extend_conflict_mutid]
                            for case_under_3 in ['case3-1: exrtend_pivot', 'case3-2: extend_mother']:
                                if case_under_3=="case3-1: exrtend_pivot":
                                    ## Case-1: extend pivot_bin to pivot_bin ∩ conflict_each
                                    each_conflict_penalty_score = calculate_difference_sum_withNA(df_pos.iloc[:, pivot_index].values, each_refine_vec, cutoff, fnfp_ratio_late, general_weight_na) + calculate_difference_sum_withNA_multiV_under_intersect_subset(each_extend_conflict_mutid, df_pos, each_refine_vec, cutoff, fnfp_ratio_late, M_out, general_weight_na) + calculate_difference_sum_withNA_multiV_under_conflict_fp(each_fp_conflict_mutid, df_pos, each_refine_vec, cutoff, fnfp_ratio_late, M_out, general_weight_na)
                                    print(each_conflict_penalty_score)
                                    print(case_under_3)
                                if case_under_3=="case3-2: extend_mother":
                                ## Case-2: extend conflict_each to pivot_bin ∩ conflict_each
                                    each_conflict_penalty_score = calculate_difference_sum_withNA(df_pos.iloc[:, pivot_index].values, ((pivot_bin*each_refine_vec)>0).astype(int), cutoff, fnfp_ratio_late, general_weight_na) + calculate_difference_sum_withNA_multiV(each_extend_conflict_mutid, df_pos, each_refine_vec, cutoff, fnfp_ratio_late, general_weight_na) + calculate_difference_sum_withNA_multiV_under_conflict_fp(each_fp_conflict_mutid, df_pos, each_refine_vec, cutoff, fnfp_ratio_late, M_out, general_weight_na)
                                    print(each_conflict_penalty_score)
                                    print(case_under_3)
                                if each_conflict_penalty_score < sub_penalty_each_conflict:
                                    sub_penalty_each_conflict = each_conflict_penalty_score
                                    sub_refine_vec = each_refine_vec
                                    sub_refine_case = case_under_3
                                    sub_extend_conflict_mutid = each_extend_conflict_mutid
                                    sub_fp_conflict_mutid = each_fp_conflict_mutid
                        print(sub_penalty_each_conflict)
                        print(sub_refine_case)
                        # minimum sub case
                        new_penalty = sub_penalty_each_conflict
                        refine_vec = sub_refine_vec
                        refine_case3 = sub_refine_case
                        extend_conflict_mutid = sub_extend_conflict_mutid
                        fp_conflict_mutid = sub_fp_conflict_mutid
                        if new_penalty < penalty_score:
                            penalty_score = new_penalty
                            refine_subclone = sub_refine_vec
                            refine_case = refine_case3
                        ### Refine df_subclone
                        print(refine_case)
                        print(penalty_score)
                        if penalty_score > theta:
                            M_out[:,pivot_index] = np.ones(len(pivot_bin), dtype=int)
                            is_germ = True
                        else:
                            if refine_case=="case1: previous":
                                M_out[:,pivot_index] = refine_subclone * pivot_bin
                            elif refine_case=="case2: current": ## 这种情况是最不可能的，罚分会明显高于其他情况；若最后是这个 case 则大概率会 move to germline
                                for flip_mutid in conflict_mutid_list:
                                    M_out[:,flip_mutid] = M_out[:,flip_mutid]*refine_subclone
                            elif refine_case=="case3-1: exrtend_pivot":
                                M_out[:,pivot_index] = refine_subclone
                                for fn_mutid in extend_conflict_mutid:
                                    M_out[:,fn_mutid] = ((M_out[:,fn_mutid]*refine_subclone)>0).astype(int)
                                for fp_mutid in fp_conflict_mutid:
                                    M_out[:,fp_mutid] = ((M_out[:,fp_mutid]-refine_subclone)>0).astype(int)
                            elif refine_case=="case3-2: extend_mother":
                                M_out[:,pivot_index] = ((pivot_bin*refine_subclone)>0).astype(int)
                                for fn_mutid in extend_conflict_mutid:
                                    M_out[:,fn_mutid] = refine_subclone
                                for fp_mutid in fp_conflict_mutid:
                                    M_out[:,fp_mutid] = ((M_out[:,fp_mutid]-refine_subclone)>0).astype(int)
                    # Recheck the conflict status
                    has_conflict = check_conflict(M_out, pivot_index, already_processed_indices)
                    # If the conflict status does not change, exit the loop
                    if has_conflict == last_conflict_state:
                        print("The conflict status does not change, and the loop exits...")
                        break
                    last_conflict_state = has_conflict  # Update the status of the last conflict
                    print(refine_case)
                # Check whether the maximum number of iterations has been exceeded
                if iterations >= max_iterations:
                    print(f"If the maximum number of iterations is exceeded, the default sub-clone vector is set: {pivot_index}")
                    M_out[:, pivot_index] = np.array(df_subclone.loc[df_subclone['subclone_mutlist'].apply(lambda x: pivot_index in x), 'subclone_vector'].iloc[0])
                    refine_case = "case4: out_iter"
        already_processed_indices.append(pivot_index)
        refine_case_list.append(refine_case)
    if scp.ul.is_conflict_free_gusfield(pd.DataFrame(M_out)):
        M_final = M_out.copy()
    else:
        subclone_vectors = list(df_subclone['subclone_vector'].values)
        M_final = check_and_process_conflict(M_out, df_pos, cutoff, fnfp_ratio_late, ISet_process, subclone_vectors, pivot_sites, general_weight_na)
    return M_final

def greedyPtreeFNFP_underSubclone(M_input, df_pos, cutoff, fnfp_ratio_late, df_subclone, ISet_process, pivot_sites, general_weight_na): # Modified Greedy algorithm for dealing with matrices with both NA and false positive values.
    # hist_cutoff = math.ceil(M_input.shape[1]/M_input.shape[0])
    subclone_vectors = list(df_subclone['subclone_vector'].values)
    M_pos = df_pos.values
    ### Step1: Sort by posterior and binary values
    ### Step2: Process mutations that are conflict-free with any subclone
    M_subclone = M_input.copy()
    Sremaining = ISet_process.copy()
    print("====================================================================")
    print("The number of mutations processed below is: " + str(len(ISet_process)))
    for idx,subclone_row in df_subclone.iterrows():
        # print(idx)
        # print(subclone_row)
        [subclone_vector, subclone_mutlist, subclone_binV_list, subclone_leader, subclone_key, subclone_leader_mutid, raw_subclone_vector] = subclone_row
        subclone_vector = np.array(subclone_vector)
        for belongs_i in range(len(subclone_mutlist)):
            pivot_index = subclone_mutlist[belongs_i]
            print(str(pivot_index)+" : "+str(df_pos.columns.tolist()[pivot_index]))
            if pivot_index in pivot_sites:
                M_subclone[:, pivot_index] = subclone_vector
            else:
                M_subclone[:, pivot_index] = subclone_vector * M_input[:,pivot_index]
    M_free = greedyPtreeFNFP_withNA(M_subclone, df_pos, cutoff, fnfp_ratio_late, df_subclone, ISet_process, pivot_sites, general_weight_na)
    return M_free

def findPivots(M_bin, df_posterior, theta, early, cutoff, fnfp_ratio_subclone, fnfp_ratio_late, general_weight_na):
    M_posterior = df_posterior.values
    ISet_init = sortIndex(M_bin, M_posterior, cutoff)
    # Find pivots and construct low-resolution tree
    ISet_init_copy = ISet_init.copy()
    ISet_init_copy.remove(early)
    ISet = [early] + ISet_init_copy
    # Extract subclone dataframe with the logic of building tree
    [df_subclone, M_copy, germ_sites] = get_subclone_with_flip(M_bin, df_posterior, ISet, theta, cutoff, fnfp_ratio_subclone, general_weight_na)
    ### Rule out germline mutations as a sublcone
    no_germ_ISet = ISet.copy()
    M_germ = np.ones((M_bin.shape[0], len(germ_sites)), dtype=int)
    while any(all(elem == 1 for elem in sub_list) for sub_list in list(df_subclone['subclone_vector'])):
        # print("Refine df_subclone: Identify one or several mutations that are believed to cover all cells during the process of determining subclones. These mutations are considered germline mutations. Reevaluate the remaining mutations to define df_subclone.")
        print("Refine df_subclone: one or several mutations may be germline.")
        index_of_all_ones = next(index for index, sublist in enumerate(list(df_subclone['subclone_vector'])) if all(elem == 1 for elem in sublist))
        new_germ_site = list(df_subclone['subclone_leader'])[index_of_all_ones]
        germ_sites = list(set(germ_sites + new_germ_site))
        no_germ_ISet = [i for i in no_germ_ISet if i not in germ_sites]
        # Construct subtree M_germ
        M_germ_list = np.array(list(df_subclone['subclone_vector'])[index_of_all_ones])
        if len(M_germ) == 0:
            # If M_germ is empty, initialize the two-dimensional list with each element starting as a list
            M_germ = np.array([[val] for val in M_germ_list])
        else:
            # Otherwise, add the value of M_germ_list column by column to the existing M_germ column
            M_germ = np.column_stack((M_germ, M_germ_list))
        # Find pivot and extract subclone.
        [df_subclone, M_copy_no_germ, new_germ_sites] = get_subclone_with_flip(M_bin, df_posterior, no_germ_ISet, theta, cutoff, fnfp_ratio_subclone, general_weight_na)
        germ_sites = list(set(germ_sites + new_germ_sites))
        no_germ_ISet = [i for i in no_germ_ISet if i not in germ_sites]
        # Construct subtree M_germ
        M_germ_list = np.array(list(df_subclone['subclone_vector'])[index_of_all_ones])
        if len(M_germ) == 0:
            # If M_germ is empty, initialize the two-dimensional list with each element starting as a list
            M_germ = np.array([[val] for val in M_germ_list])
        # else:
        #     # Otherwise, add the value of M_germ_list column by column to the existing M_germ column
        #     M_germ = np.column_stack((M_germ, M_germ_list))
    ### Process mutations, (unless germ_sites + standby_sites), (including pivot_sites + phy_sites)
    # Pivot mutations
    pivot_sites = [item for sublist in df_subclone['subclone_leader'] for item in sublist]
    # Sites that need to be futherly processed
    phy_sites = [s for s in ISet if s not in pivot_sites and s not in germ_sites]
    # ISet_process
    ISet_process = [i for i in ISet if i in pivot_sites+phy_sites]
    ### The cells contained in df_subclone are screened for further processing
    filtered_useful_cells_vector = [sum(x) for x in zip(*df_subclone['subclone_vector'])]
    filtered_useful_cells_indices = [i for i, x in enumerate(filtered_useful_cells_vector) if x > 0]
    df_subclone['raw_subclone_vector'] = df_subclone['subclone_vector']
    df_subclone['subclone_vector'] = df_subclone['raw_subclone_vector'].apply(lambda x: [x[i] for i in filtered_useful_cells_indices])
    ### Tree input: filter cells
    filtered_M_input = M_bin[filtered_useful_cells_indices,:]
    filtered_df_pos = df_posterior.iloc[filtered_useful_cells_indices]
    filtered_M_germ = M_germ[filtered_useful_cells_indices,:]
    ### Tree input: reorder muts
    for idx,pivot_line in enumerate(df_subclone.iterrows()):
        target_index = pivot_line[1]['subclone_leader'][0]
        filtered_M_input[:,target_index] = np.array(pivot_line[1]['subclone_vector'])
    # filtered_M_input[:, germ_sites] = 1
    for i in range(len(germ_sites)):
        germ_onesite = germ_sites[i]
        M_germ_onesite = filtered_M_germ[:,i]
        filtered_M_input[:,germ_onesite] = M_germ_onesite
    ### Build pivot tree
    M_init_all = greedyPtreeFNFP_underSubclone(filtered_M_input, filtered_df_pos, cutoff, fnfp_ratio_late, df_subclone, ISet_process, pivot_sites, general_weight_na)
    M_pivot_phy = M_init_all[:, ISet_process]
    ### Update M_tree including all sites (new order which is differ from filtered_df_pos)
    tree_sites_withOrder = germ_sites + ISet_process
    tree_mutid_withOrder = list(filtered_df_pos.columns[tree_sites_withOrder])
    if (len(M_germ)>0):
        M_tree_withOrder = np.hstack((M_germ[filtered_useful_cells_indices,:], M_pivot_phy))
    else:
        M_tree_withOrder = M_pivot_phy
    # Determine whether the tree is conflict-free
    if scp.ul.is_conflict_free_gusfield(pd.DataFrame(M_tree_withOrder)):
        print("The current tree is conflict-free!")
        # Calculate flipping score
        flipping_penalty = []
        for site in pivot_sites:
            flipping_penalty_i = calculate_difference_sum_withNA(filtered_df_pos.values[:,site], M_tree_withOrder[:,site], cutoff, fnfp_ratio_subclone, general_weight_na)
            flipping_penalty.append(flipping_penalty_i)
        for site in [s for s in tree_sites_withOrder if s not in pivot_sites]:
            flipping_penalty_i = calculate_difference_sum_withNA(filtered_df_pos.values[:,site], M_tree_withOrder[:,site], cutoff, fnfp_ratio_late, general_weight_na)
            flipping_penalty.append(flipping_penalty_i)
        flipping_fp_num = []
        flipping_fn_num = []
        flipping_na_to_1_num = []
        flipping_na_to_0_num = []
        for site in tree_sites_withOrder:
            flipping_fp_num_i = np.sum((filtered_df_pos.values[:, site] >= cutoff) & (M_tree_withOrder[:, site] == 0))
            flipping_fp_num.append(flipping_fp_num_i)
            flipping_fn_num_i = np.sum((filtered_df_pos.values[:, site] < cutoff) & (M_tree_withOrder[:, site] == 1))
            flipping_fn_num.append(flipping_fn_num_i)
            flipping_na_to_1_num_i = np.sum(np.isnan(filtered_df_pos.values[:, site]) & (M_tree_withOrder[:, site] == 1))
            flipping_na_to_1_num.append(flipping_na_to_1_num_i)
            flipping_na_to_0_num_i = np.sum(np.isnan(filtered_df_pos.values[:, site]) & (M_tree_withOrder[:, site] == 0))
            flipping_na_to_0_num.append(flipping_na_to_0_num_i)
        flipping_penalty_total = np.nansum(flipping_penalty)
        flipping_fp_num_total = np.nansum(flipping_fp_num)
        flipping_fn_num_total = np.nansum(flipping_fn_num)
        flipping_na_to_1_num_total = np.nansum(flipping_na_to_1_num)
        flipping_na_to_0_num_total = np.nansum(flipping_na_to_0_num)
        flipping_num_total = flipping_fp_num_total + flipping_fn_num_total + flipping_na_to_1_num_total + flipping_na_to_0_num_total
        print("Flipping score: "+str(flipping_penalty_total))
    else:
        print("The current tree is conflict! May something is wrong!")
        flipping_penalty_total = np.inf
    germ_sites_mutid = filtered_df_pos.columns[germ_sites].tolist()
    return [M_tree_withOrder, df_subclone, flipping_penalty_total, pivot_sites, tree_mutid_withOrder, germ_sites_mutid, filtered_useful_cells_indices, flipping_fp_num_total, flipping_fn_num_total, flipping_na_to_1_num_total, flipping_na_to_0_num_total, flipping_num_total]

def minFlippingPtree(M_bin, df_posterior, theta, df_features, cutoff, fnfp_ratio_subclone, fnfp_ratio_late, general_weight_na):
    ##### Get significant non-nan mutations
    # Calculate significant lower cutoff using binomial-test
    mean_nonna_count_in_allmuts = (df_posterior.count()).mean()
    # 300.57894736842104
    nonna_prop_lower_threshold = proportion_confint(mean_nonna_count_in_allmuts, df_posterior.shape[0], alpha=0.05, method='normal')[0]
    # nonna_prop_upper_threshold = proportion_confint(mean_nonna_count_in_allmuts, df_posterior.shape[0], alpha=0.05, method='normal')[1]
    df_nonna_prop = 1-(df_posterior.isna()).mean()
    sig_nonna_mutid_list = df_nonna_prop[df_nonna_prop >= nonna_prop_lower_threshold].index
    ##### Idenify early mutations and observe the distribution of cell fraction
    [early_muts, early_muts_index] = get_early_mutations(df_features)
    sig_early_muts = sig_nonna_mutid_list.intersection(early_muts).tolist()
    sig_early_muts_index = [early_muts_index[i] for i in range(len(early_muts)) if early_muts[i] in sig_early_muts]
    ##### Initialization and circulation
    min_score = np.inf
    using_early = []
    df_traversal_early_recoder = pd.DataFrame(columns=["early_mut_index", "early_mutid", "flipping_penalty_total",  "df_subclone_leader_mutid","pivot_sites", 'length_filtered_sites_order', 'germ_sites_mutid', "flipping_fp_num_total", "flipping_fn_num_total", "flipping_na_to_1_num_total", "flipping_na_to_0_num_total", "flipping_num_total"])
    for i in range(len(sig_early_muts_index)):
        early = sig_early_muts_index[i]
        # Early mutations are the few mutations that are most likely to be pivot, and the executing loop determines which one has the smallest total flip penalty as the first pivot, and the germline mutations were already ruled out first.
        print("=====> Stepwise-" + str(i+1) + ":")
        M_bin_copy = M_bin.copy()
        df_posterior_copy = df_posterior.copy()
        [M_all, df_subclone, flipping_penalty_total, pivot_sites, filtered_sites_order, germ_sites_mutid, filtered_useful_cells_indices, flipping_fp_num_total, flipping_fn_num_total, flipping_na_to_1_num_total, flipping_na_to_0_num_total, flipping_num_total] = findPivots(M_bin_copy, df_posterior_copy, theta, early, cutoff, fnfp_ratio_subclone, fnfp_ratio_late, general_weight_na)
        each_early_mut_index = early
        each_early_mutid = early_muts[i]
        each_flipping_penalty_total = flipping_penalty_total
        each_df_subclone_leader_mutid = [s for s in df_subclone['subclone_leader_mutid']]
        each_pivot_sites = pivot_sites
        each_length_filtered_sites_order = len(filtered_sites_order)
        each_germ_sites_mutid = germ_sites_mutid
        each_flipping_fp_num_total = flipping_fp_num_total
        each_flipping_fn_num_total = flipping_fn_num_total
        each_flipping_na_to_1_num_total = flipping_na_to_1_num_total
        each_flipping_na_to_0_num_total = flipping_na_to_0_num_total
        each_flipping_num_total = flipping_num_total
        # Adds the current result to df_traversal_early_recoder
        new_row = {
            "early_mut_index": each_early_mut_index,
            "early_mutid": each_early_mutid,
            "flipping_penalty_total": each_flipping_penalty_total,
            "df_subclone_leader_mutid": each_df_subclone_leader_mutid,
            "pivot_sites": each_pivot_sites,
            "length_filtered_sites_order": each_length_filtered_sites_order,
            "germ_sites_mutid": each_germ_sites_mutid,
            "flipping_fp_num_total": each_flipping_fp_num_total,
            "flipping_fn_num_total": each_flipping_fn_num_total,
            "flipping_na_to_1_num_total": each_flipping_na_to_1_num_total,
            "flipping_na_to_0_num_total": each_flipping_na_to_0_num_total,
            "flipping_num_total": each_flipping_num_total,
        }
        new_row_df = pd.DataFrame([new_row])
        df_traversal_early_recoder = pd.concat([df_traversal_early_recoder, new_row_df], ignore_index=True)
        if flipping_penalty_total < min_score:
            min_score = flipping_penalty_total
            [M_all_min, df_subclone_min, flipping_penalty_total_min, pivot_sites_min, filtered_sites_order_min, filtered_useful_cells_indices_min] = [M_all, df_subclone, flipping_penalty_total, pivot_sites, filtered_sites_order, filtered_useful_cells_indices]
            using_early = early
    print("=====================================================")
    print("===== Low-resolution tree is completely constructed! ")
    print("=====================================================")
    print("The output M_tree shape: "+str(M_all_min.shape))
    print("The df_subclone shape: "+str(df_subclone.shape))
    print("The flipping_penalty_total_min: "+str(flipping_penalty_total_min))
    print("The using early mutations is: "+str(df_features.columns[using_early]))
    return [M_all_min.astype(int), df_subclone_min, flipping_penalty_total_min, pivot_sites_min, filtered_sites_order_min, df_traversal_early_recoder, filtered_useful_cells_indices_min]

def has_duplicates(lst):
    seen = set()
    for item in lst:
        if tuple(item) in seen:
            return True
        seen.add(tuple(item))
    return False

def get_allBranchSet_as_dict(phylogeny):
    df = phylogeny.apply(lambda col: ''.join(map(str, col)), axis=0)
    unique_panels = df.unique()
    result_dict = {panel: 0 for panel in unique_panels}
    return result_dict

def str2array(string_value):
    return np.array([int(char) for char in string_value])

def get_allBranchSet(M_tree):
    columns_to_delete = np.all(M_tree == 1, axis=0)
    M_tree_noRoot = M_tree[:, ~columns_to_delete]
    clusters_allBranchSet = [str2array(s) for s in get_allBranchSet_as_dict(pd.DataFrame(M_tree_noRoot)).keys()]
    return clusters_allBranchSet

def write_subclone_persite_Binary_noConflict(mut_index, vector, df):
    my_df = df.copy()
    waiting_index = []
    if my_df.empty:
        mutlist = [mut_index]
        my_df.loc[len(my_df)] = [vector, mutlist]
        subclone_num = 1
    else:
        subclone_vectors = list(my_df['subclone_vector'].values)
        if all(np.apply_along_axis(lambda v: intersect_is_empty(vector, v), axis=1, arr=subclone_vectors)):
            # Judgment: The intersection of the add vector with all the current child clones is empty.
            mutlist = [mut_index]
            my_df.loc[len(my_df)] = [vector, mutlist]
            subclone_num = 1
        elif any(np.apply_along_axis(lambda v: intersect_is_self(vector, v), axis=1, arr=subclone_vectors)):
            # Judgment: The adding vector is subset of one of present subclones.
            vector_to_find = [v for v in subclone_vectors if intersect_is_self(vector, v)][0]
            df_i = [df_i for df_i,row in my_df.iterrows() if np.array_equal(vector_to_find, row['subclone_vector'])][0]
            my_df.loc[df_i]['subclone_mutlist'].append(mut_index)
            subclone_num = 1
        elif any(np.apply_along_axis(lambda v: intersect_is_self(v, vector), axis=1, arr=subclone_vectors)):
            # Judgment: The adding vector is superset of one of present subclones.
            df_i_find = []
            for vector_to_find in [v for v in subclone_vectors if intersect_is_self(v, vector)]:
                df_i = [df_i for df_i,row in my_df.iterrows() if np.array_equal(vector_to_find, row['subclone_vector'])][0]
                df_i_find.append(df_i)
            # Replace entries higher in the index
            find_i = min(df_i_find)
            my_df.at[find_i, 'subclone_vector'] = vector
            my_df.loc[find_i]['subclone_mutlist'].append(mut_index)
            df_i_find.remove(find_i)
            for remove_i in df_i_find:
                for mv in my_df.loc[remove_i]['subclone_mutlist']:
                    my_df.loc[find_i]['subclone_mutlist'].append(mv)
                my_df = my_df.drop(remove_i)
            subclone_num = 1
        else:
            # Judgment: The addition vector conflicts with one of the current child clones
            subclone_num = 1
            df_i_find = []
            sum_to_find = []
            for df_i,subclone_vector in enumerate(subclone_vectors):
                # There may be multiple existing subclones for different summation results. First find all subclones that are confilct with current process vector.
                if intersect_is_conflict(vector, subclone_vector):
                    df_i_find.append(df_i)
                    sum_to_find.append(sum(subclone_vector))
            if all(np.sum(vector)>sum_to_find):
                waiting_index = []
                for iter_i in df_i_find:
                    waiting_index = waiting_index + my_df.loc[iter_i]['subclone_mutlist']
                find_i = min(df_i_find)
                my_df.at[find_i, 'subclone_vector'] = vector
                my_df.at[find_i, 'subclone_mutlist'] = [mut_index]
            else:
                subclone_num = 0
                waiting_index = [mut_index]
    return[my_df, subclone_num, waiting_index]

def get_1stBranchSet(M_tree):
    columns_to_delete = np.all(M_tree == 1, axis=0)
    M_tree_noRoot = M_tree[:, ~columns_to_delete]
    ISet=np.argsort(sum(M_tree_noRoot))[::-1]
    df_1stBranchSet = pd.DataFrame(columns=["subclone_vector", "subclone_mutlist"])
    for pivot_index in ISet:
        pivot_bin=M_tree_noRoot[:,pivot_index]
        [df_1stBranchSet, pivot_subclone_num, waiting_index] = write_subclone_persite_Binary_noConflict(pivot_index, pivot_bin.astype(int), df_1stBranchSet)
    clusters_1stBranchSet = [np.array(s) for s in list(df_1stBranchSet['subclone_vector'].values)]
    return [df_1stBranchSet, clusters_1stBranchSet]

def get_earlyBranchSet(no_pivot_clusters, pivot_clusters):
    early_clusters = []
    for pivot_cluster in pivot_clusters:
        early_clusters.append(pivot_cluster)
        subset_clusters = [cluster for cluster in no_pivot_clusters if intersect_is_self(cluster, pivot_cluster)]
        if len(subset_clusters)==0:
            continue
        else:
            M_subset_clusters = np.vstack(subset_clusters).T
            secondary_clusters = get_1stBranchSet(M_subset_clusters)[1]
            early_clusters = early_clusters + secondary_clusters
    return early_clusters

def get_leafBranchSet(selected_clusters, pivot_clusters):
    leaf_clusters = []
    for pivot_cluster in pivot_clusters:
        selected_clusters_under_pivot = [cluster for cluster in selected_clusters if intersect_is_self(cluster, pivot_cluster)]
        if len(selected_clusters_under_pivot)==1:
            leaf_clusters = leaf_clusters + selected_clusters_under_pivot
        else:
            filtered_clusters = [cluster for cluster in selected_clusters_under_pivot if not np.array_equal(cluster, pivot_cluster)]
            leaf_clusters = leaf_clusters + filtered_clusters
    return leaf_clusters

def is_subset(cluster, node_cluster):
    # Gets the position of 1 in the cluster
    cluster_one_positions = np.where(cluster == 1)[0]
    # Determine if the node_cluster has at least one, but not all, 1s at these locations
    one_positions = np.sum(node_cluster[cluster_one_positions])
    # At least one 1 but not all 1
    return 0 < one_positions < len(cluster_one_positions)

def transform_Q_to_E(phred_score):
    """
    Transform phred quality score (Q) to the error probability (E)
    """
    return pow(10,-1*int(phred_score)/10)

import ast
def get_qvalues_list(data_str, ref_allele, alt_allele):
    try:
        data_dict = ast.literal_eval(data_str)
    except ValueError:
        return [], []
    ref_qvalues = []
    alt_qvalues = []
    if 'no_reads' in data_str:
        return [], []
    else:
        # traverse data_dict dictionary
        for key, count in data_dict.items():
            try:
                allele, baseq = key.split('_')
                baseq = int(baseq)
            except ValueError:
                continue
            if allele == ref_allele:
                ref_qvalues.extend([baseq] * count)
            elif allele == alt_allele:
                alt_qvalues.extend([baseq] * count)
        return [ref_qvalues, alt_qvalues]

# def calculate_likelihood_het_by_prod(q_values):
#     likelihood_het = 1.0
#     for q_value in q_values:
#         E = transform_Q_to_E(q_value)
#         likelihood_het *= 0.5 * (1 - E) + 0.5 * E / 3
#     return likelihood_het

def calculate_likelihood_het_by_BB(BB_paras):
    [ref_count, mut_count, alpha_value, beta_value] = BB_paras
    likelihood_het = beta((ref_count+alpha_value), (mut_count+beta_value))/beta(alpha_value, beta_value)
    return likelihood_het

def string_to_tuple(input_value):
    if isinstance(input_value, tuple):
        return input_value  # 如果是元组，直接返回
    elif isinstance(input_value, str):
        # 去掉括号并处理字符串
        input_value = input_value.strip('()')
        elements = [elem.strip() for elem in input_value.split(',')]
        return tuple(elements)
    else:
        raise ValueError("Input must be a string or a tuple")

def string_to_dict(input_value):
    """
    Convert a string representation of a dictionary to an actual dictionary.
    If the input is already a dictionary, return it as is.
    The keys in the input string should be valid Python identifiers.
    If the value is '*', it will be preserved as a string.
    """
    # 如果输入已经是字典，直接返回
    if isinstance(input_value, dict):
        return input_value
    # 检查输入是否是字符串
    if isinstance(input_value, str):
        # 特殊情况处理：检查是否是 '{no_reads: *}' 格式
        if input_value.strip() == "{no_reads: *}":
            return {'no_reads': '*'}
        # 使用正则表达式替换所有的键，加上引号，同时保留 '*'
        input_value = re.sub(r'(\*?\w+)', r"'\1'", input_value)  # 将键用单引号包裹
        # 使用 ast.literal_eval 将字符串转换为字典
        try:
            result_dict = ast.literal_eval(input_value)
        except (ValueError, SyntaxError):
            print("Invalid dictionary format.")
            return {}
        # 将值转换为数值类型，保留 '*' 字符串
        for key in result_dict:
            if isinstance(result_dict[key], str) and result_dict[key] == '*':
                continue  # 保持 '*' 字符串不变
            if isinstance(result_dict[key], (int, float, str)):
                try:
                    result_dict[key] = float(result_dict[key]) if '.' in str(result_dict[key]) else int(result_dict[key])
                except ValueError:
                    continue  # 处理非数值类型的值
        return result_dict
    # 如果输入既不是字典也不是字符串，返回空字典
    print("Input must be a string or a dictionary.")
    return {}

def calculate_likelihood_het_by_prod(mut_type, sc_reads_dict, mu_condition_unphased):
    """
    Use product of error rates to calculate single cell genotype likelihood when unphased.
    e.g., mu_condition_unphased=(j_allele,k_allele,not_j_allele[i],k_allele)
    """
    mu_condition_unphased = string_to_tuple(mu_condition_unphased)
    sc_reads_dict = string_to_dict(sc_reads_dict)
    # 处理特殊键
    filtered_dict = {}
    special_keys = [key for key in sc_reads_dict if re.match(r'^\*_\d+$', key)]
    if len(special_keys) == 1 and len(sc_reads_dict) == 1:
        # 如果只有一个特殊键
        filtered_dict = {'no_reads': '*'}
    else:
        # 过滤掉特殊键
        filtered_dict = {key: value for key, value in sc_reads_dict.items() if key not in special_keys}
    if mut_type=="mutated":
        j_allele=mu_condition_unphased[2]
        k_allele=mu_condition_unphased[3]
    elif mut_type=="unmutated":
        j_allele=mu_condition_unphased[0]
        k_allele=mu_condition_unphased[1]
    log_sc_gt_likelihood = np.log(1)
    if j_allele == k_allele:
        for key, value in sc_reads_dict.items():
            try:
                read_block = re.split('_', key)
                q_value = int(read_block[1])
            except (ValueError, IndexError):
                continue
            # print(f"Key: {key}, Value: {value}, Type of Value: {type(value)}")
            if read_block[0] == j_allele:
                log_sc_gt_likelihood += value * np.log(0.5 * (1 - transform_Q_to_E(q_value)) + 0.5 * (1 - transform_Q_to_E(q_value)))
            elif read_block[0] != j_allele and read_block[0] != '*':
                log_sc_gt_likelihood += value * np.log(0.5 * transform_Q_to_E(q_value) / 3 + 0.5 * transform_Q_to_E(q_value) / 3)
    else:
        for key, value in sc_reads_dict.items():
            try:
                read_block = re.split('_', key)
                q_value = int(read_block[1])
            except (ValueError, IndexError):
                continue
            if read_block[0] == j_allele or read_block[0] == k_allele:
                log_sc_gt_likelihood += value * np.log(0.5 * (1 - transform_Q_to_E(q_value)) + 0.5 * transform_Q_to_E(q_value) / 3)
            elif read_block[0] != j_allele and read_block[0] != k_allele and read_block[0] != '*':
                log_sc_gt_likelihood += value * np.log(0.5 * transform_Q_to_E(q_value) / 3 + 0.5 * transform_Q_to_E(q_value) / 3)
    return log_sc_gt_likelihood

def convert_to_log(df, epsilon=1e-10):
    """
    Convert the raw probability data to log format.
    Avoid the log-zero problem by taking the log(max(val, epsilon)).
    """
    # Use epsilon to prevent log(0) problems
    return np.log(np.maximum(df, epsilon))

def exp_normalize(log_prob):
    max_log_prob = log_prob.max()
    y = np.exp(log_prob - max_log_prob)
    return y / y.sum()

def normalize_columns(llmut_col, llunmut_col):
    # 转换为 log 格式
    log_llmut_col = llmut_col
    log_llunmut_col = llunmut_col
    paired_likelihoods_prod = np.vstack([log_llmut_col, log_llunmut_col]).T
    norm_llmut_prod, norm_llunmut_prod = np.log(np.array([exp_normalize(pair) for pair in paired_likelihoods_prod]).T)
    return pd.Series({"norm_llmut": norm_llmut_prod, "norm_llunmut": norm_llunmut_prod})

def apply_normalization(df, is_log="yes"):
    if is_log=="yes":
        process_df = df
    else:
        # Converts the original probability value to log format
        process_df = convert_to_log(df)
    norm_results = pd.DataFrame(index=df.index)
    num_columns = df.shape[1] // 2
    for i in range(num_columns):
        llmut_col = process_df.iloc[:, i]
        llunmut_col = process_df.iloc[:, num_columns + i]
        normalized = normalize_columns(llmut_col, llunmut_col)
        norm_results[f'norm_llmut_{i}'] = normalized['norm_llmut']
        norm_results[f'norm_llunmut_{i}'] = normalized['norm_llunmut']
    return norm_results

def calculate_likelihood_per_site(new_mutid, in_posterior, orderd_index, in_reads_noBulk, in_features, in_baseq, df_alpha_beta):
    # print(new_mutid)
    ##### get posterior, reads and likelihood information 
    df_posterior_likelihoods = pd.concat([
        pd.DataFrame(list(in_posterior[new_mutid]), columns=['posterior'])
    ], axis=1)
    df_posterior_likelihoods.index = list(orderd_index)
    ### split reads count and baseq information
    site_reads = in_reads_noBulk.loc[orderd_index, new_mutid]
    site_reads = site_reads.fillna('0/0')
    # get the reads count of each type(ref, mut, total)
    split_reads_info_2cols_ref_total = (site_reads.str.split('/', expand=True)).astype(int)
    ref_count = (split_reads_info_2cols_ref_total[1]).astype(int)-(split_reads_info_2cols_ref_total[0]).astype(int)
    df_site_reads = pd.concat([site_reads, ref_count, split_reads_info_2cols_ref_total], axis=1)
    df_site_reads.columns = ['reads_info', 'ref_count', 'mut_count', 'total_count']
    ##### get baseq information 
    # ref and alt allele
    ref_allele = in_features[new_mutid]['ref']
    alt_allele = in_features[new_mutid]['mut']
    site_baseq = in_baseq.loc[orderd_index, new_mutid]
    # # get_qvalues_list(site_baseq, ref_allele, alt_allele)
    # all_qvalues = site_baseq.apply(lambda x: (np.concatenate((get_qvalues_list(str(x), ref_allele, alt_allele)[0], get_qvalues_list(str(x), ref_allele, alt_allele)[1]))).astype(int))
    # depth = all_qvalues.apply(len)
    # df_site_baseq
    df_site_baseq = pd.concat([
        site_baseq
    ], axis=1)
    df_site_baseq.columns = ['baseq_info']
    # df_site_baseq = pd.concat([
    #     site_baseq,
    #     pd.DataFrame({'all_qvalues':all_qvalues}),
    #     pd.DataFrame({'depth':depth})
    # ], axis=1)
    # df_site_baseq.columns = ['baseq_info', 'baseq_list', 'depth']
    ##### combine all information into df_site_combined
    df_site_combined = pd.concat([df_posterior_likelihoods, df_site_reads, df_site_baseq, df_alpha_beta], axis=1)
    ################################################################
    ##### Recalculate likelihood
    site_genotype = in_features[new_mutid]['genotype']
    ### new_site_llmut_prod = df_site_combined.apply(lambda row: cal_product_mut(row['baseq_info'], mu_condition_unphased), axis=1)
    likelihoods_unmutant = df_site_combined.apply(lambda row: calculate_likelihood_het_by_prod("unmutated", row['baseq_info'], site_genotype), axis=1)
    likelihoods_het_by_prod = df_site_combined.apply(lambda row: calculate_likelihood_het_by_prod("mutated", row['baseq_info'], site_genotype), axis=1)
    # calculate_likelihood_het_by_BB(np.array(df_site_combined.iloc[0,:][['ref_count', 'mut_count', 'alpha', 'beta']]))
    likelihoods_het_by_BB = np.log(df_site_combined.apply(lambda row: calculate_likelihood_het_by_BB(row[['ref_count', 'mut_count', 'alpha', 'beta']]), axis=1))
    ##### Normalization   
    ### Product
    paired_likelihoods_prod = np.vstack((likelihoods_het_by_prod, likelihoods_unmutant)).T
    norm_llmut_prod, norm_llunmut_prod = np.log(np.array([exp_normalize(pair) for pair in paired_likelihoods_prod]).T)
    ### Beta-binomial
    paired_likelihoods_BB = np.vstack((likelihoods_het_by_BB, likelihoods_unmutant)).T
    norm_llmut_BB, norm_llunmut_BB = np.log(np.array([exp_normalize(pair) for pair in paired_likelihoods_BB]).T)
    ##### combine data into df_site_combined
    df_site_combined['likelihoods_mutant_prod'] = norm_llmut_prod
    df_site_combined['likelihoods_unmutant_prod'] = norm_llunmut_prod
    df_site_combined['likelihoods_mutant_BB'] = norm_llmut_BB
    df_site_combined['likelihoods_unmutant_BB'] = norm_llunmut_BB
    return [norm_llmut_prod, norm_llunmut_prod, norm_llmut_BB, norm_llunmut_BB]

def log_sum_exp(log_probs):
    max_log_prob = np.max(log_probs)
    return max_log_prob + np.log(np.sum(np.exp(log_probs - max_log_prob)))

def logdiffexp(a, b):
       return logsumexp([a, b], b=[1, -1])

# def logdiffexp(x, y):
#     """ 计算 log(exp(x) - exp(y))，用于计算排除某些情况后的概率 """
#     return x + np.log(1 - np.exp(y - x))

# def logdiffexp_2para(x, y):
#     # 计算 log(exp(x) - exp(y) - exp(z))
#     max_value = np.max([x, y])
#     return max_value + np.log1p(-np.exp(x - max_value) - np.exp(y - max_value))

# def logdiffexp_3para(x, y, z):
#     # 计算 log(exp(x) - exp(y) - exp(z))
#     max_value = np.max([x, y, z])
#     return max_value + np.log1p(-np.exp(x - max_value) - np.exp(y - max_value) - np.exp(z - max_value))

def DP_calSomaticPosterior_withoutTree(site_llmut, site_llunmut, cell_num):
    log_sc_gt_likelihood_dict = {}
    log_sc_gt_likelihood_given_l_list = []
    log_sc_gt_posterior_given_l_list = []
    log_sc_gt_likelihood_dict[(0,0)]=np.log(1)
    for i in range(cell_num):
        log_sc_gt_likelihood_dict[(-1,i)] = -np.inf
    for i in range(1,int(cell_num)+1):
        log_sc_gt_likelihood_dict[(i,0)] = -np.inf
    for l in range(int(cell_num)+1):
        for n in range(1,int(cell_num)+1):
            # M_l_n = M_l_n-1 * P(Dn|Gn_is_not_mutate) + M_l-1_n-1 * P(Dn|Gn_is_mutated)
            log_sc_gt_likelihood_dict[(l,n)] = logsumexp([log_sc_gt_likelihood_dict[(l,n-1)] + site_llunmut[n-1], log_sc_gt_likelihood_dict[(l-1,n-1)] + site_llmut[n-1]])
        log_sc_gt_likelihood_given_l_list.append(log_sc_gt_likelihood_dict[(l,cell_num)])
    posterior_l_is_0_dp_path=exp_normalize(np.array(log_sc_gt_likelihood_given_l_list))[0]
    somatic_posterior_per_site = float(1 - posterior_l_is_0_dp_path)
    return somatic_posterior_per_site

def get_optimal_cell_combination(path_dict, cell_num, max_l):
    # Retrace the path to find mutant combinations
    l, n = max_l, cell_num
    mut_cell_indices = [0] * cell_num
    while (l, n) in path_dict and path_dict[(l, n)] is not None:
        prev_l, prev_n = path_dict[(l, n)]
        if prev_l == l - 1:  # mutation
            mut_cell_indices[n-1] = 1  # n-1 cell mutation
        l, n = prev_l, prev_n
    return mut_cell_indices

def DP_calSomaticPosterior_late(site_llmut, site_llunmut, leaf_cluster):
    ### Late: t least 1 cell
    non_existing = np.sum(site_llunmut[leaf_cluster==0])
    probs_llmut = site_llmut[leaf_cluster==1]
    probs_llunmut = site_llunmut[leaf_cluster==1]
    cell_num = len(probs_llmut)
    # The initializing state
    likelihood_dict = {}
    likelihood_list = []
    path_dict = {}  # To track the path
    # Initialize DP tables
    likelihood_dict[(0,0)] = np.log(1)
    path_dict[(0,0)] = None  # No previous state
    for i in range(cell_num):
        likelihood_dict[(-1,i)] = -np.inf
        path_dict[(-1,i)] = None
    for i in range(1, int(cell_num)+1):
        likelihood_dict[(i,0)] = -np.inf
        path_dict[(i,0)] = None
    # DP: n cells, l mutant cells
    for l in range(int(cell_num)+1):
        for n in range(1, int(cell_num)+1):
            # M_l_n = M_l_n-1 * P(Dn|Gn_is_not_mutate) + M_l-1_n-1 * P(Dn|Gn_is_mutated)
            if n == 1:  # Apply non_existing only at the first step
                option1 = likelihood_dict[(l, n-1)] + probs_llunmut[n-1] + non_existing
                option2 = likelihood_dict[(l-1, n-1)] + probs_llmut[n-1] + non_existing
            else:
                option1 = likelihood_dict[(l, n-1)] + probs_llunmut[n-1]
                option2 = likelihood_dict[(l-1, n-1)] + probs_llmut[n-1]
            likelihood_dict[(l, n)] = logsumexp([option1, option2])
            if option1 > option2:
                path_dict[(l,n)] = (l, n-1)
            else:
                path_dict[(l,n)] = (l-1, n-1)
        likelihood_list.append(likelihood_dict[(l,cell_num)])
    all_prob = logsumexp(likelihood_list)
    at_least_one_mutcell_prob = np.log(np.exp(all_prob) - np.exp(likelihood_list[0]))
    one_mutcell_prob = likelihood_list[1]
    # Find the combination of mutations with the greatest probability
    max_l = np.argmax([likelihood_dict[(l, cell_num)] for l in range(cell_num + 1)])
    mut_cells = get_optimal_cell_combination(path_dict, cell_num, max_l)
    return [at_least_one_mutcell_prob, likelihood_dict[(max_l, cell_num)], mut_cells, one_mutcell_prob]

def DP_calSomaticPosterior_early(site_llmut, site_llunmut, pivot_clusters):
    ### Early: at least 2 clusters but not all clusters
    cluster_llmut = np.array([np.sum(site_llmut[c == 1]) for c in pivot_clusters])
    cluster_llunmut = np.array([np.sum(site_llunmut[c == 1]) for c in pivot_clusters])
    cluster_num = len(cluster_llmut)
    # The initializing state
    likelihood_dict = {}
    likelihood_list = []
    path_dict = {}
    # Initialize DP tables
    likelihood_dict[(0,0)] = np.log(1)
    path_dict[(0,0)] = None  # No previous state
    for i in range(cluster_num):
        likelihood_dict[(-1, i)] = -np.inf
        path_dict[(-1,i)] = None
    for i in range(1, int(cluster_num)+1):
        likelihood_dict[(i, 0)] = -np.inf
        path_dict[(i,0)] = None
    # DP: n cells, l mutant cells
    for l in range(int(cluster_num)+1):
        for n in range(1, int(cluster_num)+1):
            # Calculate the two possible values for logsumexp
            option1 = likelihood_dict[(l,n-1)] + cluster_llunmut[n-1]
            option2 = likelihood_dict[(l-1,n-1)] + cluster_llmut[n-1]
            likelihood_dict[(l,n)] = logsumexp([option1, option2])
            # Choose the maximum of the two options
            if option1 > option2:
                path_dict[(l,n)] = (l, n-1)
            else:
                path_dict[(l,n)] = (l-1, n-1)
        likelihood_list.append(likelihood_dict[(l,cluster_num)])
    all_prob = logsumexp(likelihood_list)
    at_least_two_mutcluster_but_not_allcluster_prob = logdiffexp(all_prob, logsumexp([likelihood_list[0], likelihood_list[1], likelihood_list[-1]]))
    # Find the combination of mutations with the greatest probability
    max_l = range(2, cluster_num)[np.argmax([likelihood_dict[(l, cluster_num)] for l in range(2, cluster_num)])]  # at least 2 clusters but all clusters
    mut_clusters = get_optimal_cell_combination(path_dict, cluster_num, max_l)
    return [at_least_two_mutcluster_but_not_allcluster_prob, likelihood_dict[(max_l, cluster_num)], mut_clusters]

def DP_calSomaticPosterior_internal(site_llmut, site_llunmut, internal_subset):
    ### Internal: at least 2 clusters
    non_existing = np.sum(site_llunmut[sum(internal_subset)==0])
    cluster_llmut = np.array([np.sum(site_llmut[c==1]) for c in internal_subset])
    cluster_llunmut = np.array([np.sum(site_llunmut[c==1]) for c in internal_subset])
    cluster_num = len(cluster_llmut)
    # The initializing state
    likelihood_dict = {}
    path_dict = {}  # To track the path
    likelihood_list = [] 
    # Initial state
    likelihood_dict[(0,0)] = np.log(1)
    path_dict[(0,0)] = None
    for i in range(cluster_num):
        likelihood_dict[(-1,i)] = -np.inf
        path_dict[(-1,i)] = None
    for i in range(1, int(cluster_num)+1):
        likelihood_dict[(i,0)] = -np.inf
        path_dict[(i,0)] = None
    # DP: n cells, l mutant cells
    for l in range(int(cluster_num)+1):
        for n in range(1, int(cluster_num)+1):
            option1 = likelihood_dict[(l, n-1)] + cluster_llunmut[n-1] + non_existing
            option2 = likelihood_dict[(l-1, n-1)] + cluster_llmut[n-1] + non_existing
            likelihood_dict[(l,n)] = logsumexp([option1, option2])
            if option1 > option2:
                path_dict[(l, n)] = (l, n-1)
            else:
                path_dict[(l, n)] = (l-1, n-1)
        likelihood_list.append(likelihood_dict[(l,cluster_num)])
    all_prob = logsumexp(likelihood_list)
    at_least_two_mutcluster_prob = logdiffexp(all_prob, logsumexp([likelihood_list[0],likelihood_list[1]]))
    # Find the combination of mutations with the greatest probability
    max_l = range(2, cluster_num+1)[np.argmax([likelihood_dict[(l, cluster_num)] for l in range(2, cluster_num+1)])]  # at least 2 clusters but all clusters
    mut_clusters = get_optimal_cell_combination(path_dict, cluster_num, max_l)
    return [at_least_two_mutcluster_prob, likelihood_dict[(max_l, cluster_num)], mut_clusters]

def get_newSomaticPosterior(site_llmut, site_llunmut, node_clusters_dict):
    # node/cluster
    leaf_clusters = node_clusters_dict['leaf_clusters']
    internal_clusters = node_clusters_dict['internal_clusters']
    pivot_clusters = node_clusters_dict['pivot_clusters']
    node_clusters = leaf_clusters+internal_clusters+pivot_clusters
    node_num = np.sum(leaf_clusters) + len(internal_clusters)
    ### calculate posterior under conditions
    node_position_list = []
    branch_state_list = []
    cluster_list = []
    mut_posterior_list = []
    max_mut_posterior_list = []
    onecell_mut_posterior_list = []
    # Condition: late branches
    for idx,leaf_cluster in enumerate(leaf_clusters):
        node_position_list.append("leaf_clusters")
        branch_state_list.append("late")
        # dp
        late_results = DP_calSomaticPosterior_late(site_llmut, site_llunmut, leaf_cluster)
        # cluster
        late_cluster_predicted = np.zeros_like(leaf_cluster)
        indices_with_1 = np.where(leaf_cluster == 1)[0]
        late_cluster_predicted[indices_with_1] = late_results[2]
        cluster_list.append(late_cluster_predicted)
        # posterior
        mut_posterior_list.append(late_results[0])
        max_mut_posterior_list.append(late_results[1])
        onecell_mut_posterior_list.append(late_results[3])
    # Condition: internal branches
    for idx,internal_cluster in enumerate(internal_clusters):
        node_position_list.append("internal_clusters")
        branch_state_list.append("internal")
        # dp
        internal_subset = [node_cluster for node_cluster in node_clusters if is_subset(internal_cluster, node_cluster)]
        if len(internal_subset)>1:
            internal_results = DP_calSomaticPosterior_internal(site_llmut, site_llunmut, internal_subset)
            # cluster
            indices_with_1 = np.where(np.array(internal_results[2]) == 1)[0]
            internal_cluster_predicted = np.sum(np.array(internal_subset)[indices_with_1], axis=0)
            cluster_list.append(internal_cluster_predicted)
            # posterior
            mut_posterior_list.append(internal_results[0])
            max_mut_posterior_list.append(internal_results[1])
        elif len(internal_subset)==1:
            # cluster
            cluster_list.append(internal_subset[0])
            # posterior
            mut_posterior_list.append(np.sum(np.concatenate([site_llmut[internal_cluster==1], site_llunmut[internal_cluster==0]])))
            max_mut_posterior_list.append(np.sum(np.concatenate([site_llmut[internal_cluster==1], site_llunmut[internal_cluster==0]])))
    # Condition: early branches
    if len(pivot_clusters)>2:
        node_position_list.append("pivot_clusters")
        branch_state_list.append("early")
        # dp
        early_results = DP_calSomaticPosterior_early(site_llmut, site_llunmut, pivot_clusters)
        # cluster
        indices_with_1 = np.where(np.array(early_results[2]) == 1)[0]
        early_cluster_predicted = np.sum(np.array(pivot_clusters)[indices_with_1], axis=0)
        cluster_list.append(early_cluster_predicted)
        # posterior
        mut_posterior_list.append(early_results[0])
        max_mut_posterior_list.append(early_results[1])
    else:
        node_position_list.append("pivot_clusters")
        branch_state_list.append("none")
        cluster_list.append(np.zeros(len(site_llmut), dtype=int))
        mut_posterior_list.append(-np.inf)
        max_mut_posterior_list.append(-np.inf)
    # Collect final results
    df_allSP = pd.DataFrame({'node_position': node_position_list, 'cluster': cluster_list, 'branch_state': branch_state_list, 'somatic_posterior_conditional': max_mut_posterior_list})
    max_row = df_allSP.loc[df_allSP['somatic_posterior_conditional'].idxmax()]
    # Calculate posterior
    wrong_posterior = np.exp(np.sum(site_llmut))
    unmut_posterior = np.exp(np.sum(site_llunmut))
    mut_posterior = np.exp(logsumexp(mut_posterior_list))*(1/node_num)
    onecell_posterior = np.exp(logsumexp(onecell_mut_posterior_list))*(1/node_num)
    ### normalized SP
    total_prob = wrong_posterior + unmut_posterior + mut_posterior
    normalized_mut = mut_posterior / total_prob
    normalized_onecell = onecell_posterior / total_prob
    normalized_unmut = unmut_posterior / total_prob
    normalized_wrong = wrong_posterior / total_prob
    normalized_mut_max = np.exp(max_row['somatic_posterior_conditional'])*(1/node_num) / total_prob
    ### out frame
    df_outSP = max_row.copy()
    df_outSP = df_outSP.drop(['somatic_posterior_conditional'])
    df_outSP['somatic_posterior_per_site_max'] = normalized_mut_max
    df_outSP['somatic_posterior_per_site_onecell'] = normalized_onecell
    df_outSP['somatic_posterior_per_site'] = normalized_mut
    df_outSP['artifact_posterior_per_site'] = normalized_unmut
    df_outSP['germline_posterior_per_site'] = normalized_wrong
    return [df_allSP, df_outSP]

def all_newSomaticPosterior(df_llmut_prod, df_llunmut_prod, df_llmut_BB, df_llunmut_BB, M_lowR):
    ### Step1: Get all nodes/clusters we need
    all_clusters = get_allBranchSet(M_lowR)
    [df_subclone, pivot_clusters] = get_1stBranchSet(M_lowR)
    no_pivot_clusters = [cluster for cluster in all_clusters if not any(np.array_equal(cluster, pivot_clusters) for pivot_clusters in pivot_clusters)]
    selected_clusters = get_earlyBranchSet(no_pivot_clusters, pivot_clusters)
    leaf_clusters = get_leafBranchSet(selected_clusters, pivot_clusters)
    internal_clusters = [cluster for cluster in selected_clusters if not any(np.array_equal(cluster, leaf_cluster) for leaf_cluster in leaf_clusters)]
    # Fornmated low resolution tree 
    M_pivot = np.array(selected_clusters).T
    # Generate important clusters dictionary
    node_clusters_dict = {
        "leaf_clusters": leaf_clusters, 
        "internal_clusters": internal_clusters, 
        "pivot_clusters": pivot_clusters
    }
    ### Step2: initiate 2 dataframe
    # Out somatic_posterior_per_site based on tree (normalization)
    df_newSP_out_init = pd.DataFrame(columns=['node_position', 'cluster', 'branch_state', 'somatic_posterior_per_site_max', 'somatic_posterior_per_site_onecell', 'somatic_posterior_per_site', 'artifact_posterior_per_site', 'germline_posterior_per_site'])
    df_newSP_out = df_newSP_out_init.dropna(axis=1, how='all')
    # All conditions somatic_posterior_per_site based on tree (list all)
    df_list = []
    for node_position, clusters in node_clusters_dict.items():
        if node_position == 'pivot_clusters':
            # 对于 pivot_clusters，只生成一行
            df_list.append(pd.DataFrame({
                'node_position': [node_position],
                'cluster': ['uncertain'],
                'branch_state': ['early']
            }))
        else:
            for cluster in clusters:
                branch_state = 'internal' if node_position == 'internal_clusters' else 'leaf'
                df_list.append(pd.DataFrame({
                    'node_position': [node_position],
                    'cluster': [cluster],
                    'branch_state': [branch_state]
                }))
    df_newSP_allNodes = pd.concat(df_list, ignore_index=True)
    ### Step3: Calculation SP for each site
    M_llmut_prod = df_llmut_prod.values
    M_llunmut_prod = df_llunmut_prod.values
    M_llmut_BB = df_llmut_BB.values
    M_llunmut_BB = df_llunmut_BB.values
    withoutTree_posterior_prod = []
    withoutTree_posterior_BB = []
    df_newSP_allNodes_prod = df_newSP_allNodes.copy()
    df_newSP_allNodes_BB = df_newSP_allNodes.copy()
    df_newSP_prod_out = df_newSP_out.copy()
    df_newSP_BB_out = df_newSP_out.copy()
    for i in tqdm(range(0, M_llunmut_prod.shape[1])):
        new_mutid = df_llunmut_prod.columns[i]
        site_llmut_prod = M_llmut_prod[:,i]
        site_llunmut_prod = M_llunmut_prod[:,i]
        site_llmut_BB = M_llmut_BB[:,i]
        site_llunmut_BB = M_llunmut_BB[:,i]
        ## prod results:
        each_newSP_prod = get_newSomaticPosterior(site_llmut_prod, site_llunmut_prod, node_clusters_dict)
        # df_allNodes
        new_col_in_df_newSP_allNodes_prod = pd.DataFrame({new_mutid: each_newSP_prod[0]['somatic_posterior_conditional']})
        df_newSP_allNodes_prod = pd.concat([df_newSP_allNodes_prod, new_col_in_df_newSP_allNodes_prod], axis=1, ignore_index=True)
        # out posterior (mut, unmut, wrong)
        df_newSP_prod_out = pd.concat([df_newSP_prod_out, pd.DataFrame([each_newSP_prod[1]])], ignore_index=True)
        withoutTree_posterior_prod.append(DP_calSomaticPosterior_withoutTree(site_llmut_prod, site_llunmut_prod, len(site_llunmut_prod)))
        ## BB results: 
        each_newSP_BB = get_newSomaticPosterior(site_llmut_BB, site_llunmut_BB, node_clusters_dict)
        # df_allNodes
        new_col_in_df_newSP_allNodes_BB = pd.DataFrame({new_mutid: each_newSP_BB[0]['somatic_posterior_conditional']})
        df_newSP_allNodes_BB = pd.concat([df_newSP_allNodes_BB, new_col_in_df_newSP_allNodes_BB], axis=1, ignore_index=True)
        # out posterior (mut, unmut, wrong)
        df_newSP_BB_out = pd.concat([df_newSP_BB_out, pd.DataFrame([each_newSP_BB[1]])], ignore_index=True)
        withoutTree_posterior_BB.append(DP_calSomaticPosterior_withoutTree(site_llmut_BB, site_llunmut_BB, len(site_llunmut_BB)))
    df_newSP_allNodes_prod.columns = df_newSP_allNodes_BB.columns = ['node_position', 'cluster', 'branch_state']+list(df_llunmut_prod.columns)
    df_newSP_prod_out.index = df_newSP_BB_out.index = df_llunmut_prod.columns
    return [df_subclone, M_pivot, df_newSP_allNodes_prod, df_newSP_prod_out, df_newSP_allNodes_BB, df_newSP_BB_out, withoutTree_posterior_prod, withoutTree_posterior_BB]

def calculate_ratio(reads_info_item):
    if isinstance(reads_info_item, str) and '/' in reads_info_item:
        parts = reads_info_item.split('/')
        numerator = float(parts[0])
        denominator = float(parts[1])
        if denominator != 0:
            return numerator / denominator
        else:
            return np.nan
    else:
        return reads_info_item

def test_features(df_features, tp_mutid, fp_mutid, kind):
    # Test features
    features_num = df_features.shape[1]
    ncols_num = 3
    nrows_num = int(np.ceil(features_num/ncols_num))
    fig, axes = plt.subplots(nrows=nrows_num, ncols=ncols_num, figsize=((ncols_num*4), (nrows_num*3)))
    feature_test = []
    for idx,feature in enumerate(df_features.columns):
        # Gets the current subgraph (Divide by ncols)
        row = idx // ncols_num  # Row in which the subgraph is located
        col = idx % ncols_num  # Column of the subgraph
        ax = axes[row, col]
        # data by group
        TP=list(df_features.loc[tp_mutid,][feature])
        FP=list(df_features.loc[fp_mutid,][feature])
        # 执行Wilcoxon秩和检验
        statistic, p_value = stats.ranksums(TP, FP)
        # 输出检验结果
        # print("Wilcoxon Rank Sum Test Statistic:", statistic)
        print(feature, ":", p_value)
        feature_test.append(p_value)
        # boxplot
        if kind=="boxplot":
            sns.boxplot(data=[TP, FP], ax=ax)
        elif kind=="violinplot":
            sns.violinplot(data=[TP, FP], ax=ax)
        # ax.set_ylim(0, 1)
        ax.set_title(feature)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(['TP', 'FP'])
        ax.text(0.5, 0.9, f"p-value: {p_value:.2e}", transform=ax.transAxes, ha="center")
    # show or save
    plt.tight_layout()
    # plt.show()
    plt.savefig(outputpath+"/compared_features.filtered_sites."+kind+".pdf")  # 图片保存的文件名和格式
    # feature test dataframe
    df_feature_test = pd.DataFrame({
        'feature_name': df_features.columns,
        'p_value': feature_test,
        })
    df_feature_test.to_csv(outputpath+"/out_feature_test.filtered_sites.txt", sep="\t")
    return df_feature_test

def CalBranchSupports(num_bootstraps, out_phylogeny, fp_pool, num_columns_to_select, branch_set):
    ### Calculate branch stability by bootstrap.
    # num_columns_to_select = 19 (20% fpr)
    # Branch confidence: bootstrap
    for iterator in range(num_bootstraps):
        fp_column_names_to_select = random.sample(fp_pool.columns.tolist(), num_columns_to_select)
        df_select = pd.concat([out_phylogeny, fp_pool[fp_column_names_to_select]], axis=1)
        M_bin = np.array(df_select)
        [bret, pret, M_out] = greedyPtreeNA(M_bin,sum(M_bin),oc=0.1,hc=25)
        bootstrap_phylogeny = pd.DataFrame(M_out)
        bootstrap_phylogeny.index = out_phylogeny.index
        bootstrap_phylogeny.columns = list(out_phylogeny.columns)+fp_column_names_to_select
        bootstrap_phylogeny.index.name = "cellIDxmutID"
        # Intermediate tree visualization
        bootstrap_tree = scp.ul.to_tree(bootstrap_phylogeny)
        scp.pl.clonal_tree(bootstrap_tree, output_file=outputpath+"/tree_scphylo.allsites.bootstrap_"+str(num_bootstraps)+"_iterate"+str(iterator)+".pdf")
        # unique columns for bootstrap tree
        bootstrap_branch = list(get_allBranchSet_as_dict(bootstrap_phylogeny).keys())
        for branch in bootstrap_branch:
            if branch in branch_set:
                branch_set[branch] += 1
    # Calculate the bootstrap support value
    branch_supports = branch_set.copy()
    for key in branch_supports:
        branch_supports[key] /= num_bootstraps
        branch_supports[key] = "{:.1%}".format(branch_supports[key])
        # branch_supports[key] = "{:.2f}".format(branch_supports[key])
    return branch_supports

# Get >=？ dp but no mutAllele spots
def check_valid_spot(value, cutoff):
    if pd.isna(value):
        return False  # 如果是 NaN，则不符合条件
    try:
        alt, total = map(int, value.split('/'))  # 分割 'alt' 和 'total'
        return alt == 0 and total >= cutoff  # 判断是否符合条件
    except ValueError:
        return False  # 如果分割出错，则不符合条件

def get_valid_spot_ids(df, cutoff):
    # 去掉 'bulk' 行，只对其他行进行处理
    # filtered_df = df.drop('bulk')
    filtered_df = df.copy()
    # 对每一列使用 apply 和 lambda 函数来筛选出符合条件的行名
    spot_ids = filtered_df.apply(lambda col: col.index[col.apply(lambda x: check_valid_spot(x, cutoff))].tolist())
    return spot_ids

# Function to determine phylogeny_label based on the criteria
def determine_phylogeny_label(row, pass_tree_cutoff, unpass_tree_cutoff):
    if row['mutant_cellnum'] == 1:
        return "cell_specific"
    elif row['mutant_cellnum'] == 0:
        return "absent"
    elif ((row['prod_somatic_posterior_per_site'] >= pass_tree_cutoff and 
            row['prod_somatic_posterior_per_site_onecell'] < unpass_tree_cutoff) and 
        (row['BB_somatic_posterior_per_site'] >= pass_tree_cutoff and 
                row['BB_somatic_posterior_per_site_onecell'] < unpass_tree_cutoff)):
        return "successful_pass"
    else:
        return "failed_pass"

# Function to determine phylogeny_label based on the criteria
def determine_phylogeny_label_by_one_likelihood(row, pass_tree_cutoff, unpass_tree_cutoff):
    if row['mutant_cellnum'] == 1:
        return "cell_specific"
    elif row['mutant_cellnum'] == 0:
        return "absent"
    elif ((row['somatic_posterior_per_site'] >= pass_tree_cutoff and 
            row['somatic_posterior_per_site_onecell'] < unpass_tree_cutoff)):
        return "successful_pass"
    else:
        return "failed_pass"

def compare_elements_vectorized(val1, val2):
    """
    Vectorized comparison of two NumPy arrays and return the count for each flip type.
    """
    # Create a dictionary to store the flip count
    flip_counts = {'flipping_0_to_1': 0, 'flipping_1_to_0': 0,
                   'flipping_NA_to_0': 0, 'flipping_NA_to_1': 0}
    # Count the flip type
    flip_counts['flipping_NA_to_0'] = np.sum((val1 == 3) & (val2 == 0))
    flip_counts['flipping_NA_to_1'] = np.sum((val1 == 3) & (val2 == 1))
    flip_counts['flipping_1_to_0'] = np.sum((val1 == 1) & (val2 == 0))
    flip_counts['flipping_0_to_1'] = np.sum((val1 == 0) & (val2 == 1))
    return flip_counts

def calculate_flip_counts_per_site(df1, df2):
    # Make sure both data boxes have the same shape
    if df1.shape != df2.shape:
        raise ValueError("The shapes of the two dataframes do not match.")
    # Get column name
    column_names = df1.columns
    # Define a function to calculate the flip count of a column
    def count_flips(column):
        val1 = df1[column].values
        val2 = df2[column].values
        # Calculate the flip count using vectorization operations
        return compare_elements_vectorized(val1, val2)
    # Apply the count flips function to each column
    flip_counts_results = df1.columns.to_series().apply(count_flips)
    # Convert the result to a data frame
    flip_counts_df = pd.DataFrame(list(flip_counts_results))
    flip_counts_df.columns = ['flipping_False_Negative', 'flipping_False_Positive', 'flipping_NA_to_0', 'flipping_NA_to_1']
    flip_counts_df.index = df1.columns
    return flip_counts_df

def find_flipping_spots(series_in_bin, series_phylogeny, condition_in_bin, condition_phylogeny):
    """
    Find a list of eligible line names (Spots)
    """
    return series_in_bin[(series_in_bin == condition_in_bin) & (series_phylogeny == condition_phylogeny)].index.tolist()

def WriteTfile(out_prefix, matrix, rownames, colnames, judge): # writes matrix output as an integer matrix, uses the format given in the documentary. 
    matrix_output = matrix.astype(int)
    df_output = pd.DataFrame(matrix_output)
    df_output.index = rownames
    df_output.columns = colnames
    df_output.index.name = "cellIDxmutID"
    is_cf = scp.ul.is_conflict_free_gusfield(df_output)
    if judge == "yes":
        if is_cf:
            print("Current tree is conflict-free and output binary matrix and plot phylogenetic tree !")
            df_output.to_csv(out_prefix+".CFMatrix", sep="\t")
            tree = scp.ul.to_tree(df_output)
            scp.pl.clonal_tree(tree, output_file=out_prefix+".tree_scphylo.pdf")
        else:
            print("Current tree is not conflict-free !")
    else:
        print("Only output binary matrix !")
        df_output.to_csv(out_prefix+".CFMatrix", sep="\t")


############################## Main ##############################
print(" ===== The input data is already succussfully loaded! ===== ")
# Log information for each dataframe
log_dataframe_info(in_posterior, "Posterior Matrix")
log_dataframe_info(in_llmut, "Likelihood Mut Matrix")
log_dataframe_info(in_llunmut, "Likelihood Unmut Matrix")
log_dataframe_info(in_reads, "Allele Count")
log_dataframe_info(in_baseq, "Base Quality Dictionary")
log_dataframe_info(in_features, "Preprocessed Features")
log_dataframe_info(in_alpha_beta, "Alpha Beta")


# Select project as "phylo" or "valid" to decide following prcessing
##### Usage1: Construct phylogenetic tree #####
if args.phylo_or_valid=="phylo":
    ### Step1: Pre-scan ###
    cleaning_result = data_cleaning(in_posterior, in_reads, in_features, cutoff, args.mut_filter)
    if cleaning_result is not None:
        [pending_sites, pending_posterior, pending_reads, pending_features, pending_bin, no_zeros_ones_posterior, no_zeros_ones_reads, no_zeros_ones_features, no_mut_cellname] = cleaning_result
        if 'in_reads' in globals():
            ## get binary matrix from posterior matrix
            in_bin_withNA3 = posterior2ter_NAto3_bothPosteriorMutallele(pending_posterior, pending_reads, cutoff)
            in_bin_withoutNA = posterior2bin_bothPosteriorMutallele(pending_posterior, pending_reads, cutoff)
            in_bin_withNA3.to_csv(outputpath+"/df_binary_withNA3_for_circosPlot.txt", sep="\t")
            in_bin_withoutNA.to_csv(outputpath+"/df_binary_withoutNA.txt", sep="\t")
        else:
            ## get binary matrix from posterior matrix
            in_bin_withNA3 = posterior2ter_NAto3_onlyPosterior(pending_posterior, cutoff)
            in_bin_withoutNA = posterior2bin_onlyPosterior(pending_posterior, cutoff)
            in_bin_withNA3.to_csv(outputpath+"/df_binary_withNA3_for_circosPlot.txt", sep="\t")
            in_bin_withoutNA.to_csv(outputpath+"/df_binary_withoutNA.txt", sep="\t")
        ## generate conflict matrix 
        conflict_results = find_all_conflicting_mutations(in_bin_withoutNA)
        conflict_df = convert_conflict_results_to_dataframe(conflict_results)
        conflict_df.to_csv(outputpath+"/df_conflict_situation.txt", sep="\t")
        ### Step2: Pivot tree ###
        ## Parameter theta: It is possible to determine whether the mutation is a key parameter of the germ line, set to the outlier value of the cell number.
        mutant_cellnum = pending_features.T["mutant_cellnum"]
        mean_mutant_cellnum = np.mean(mutant_cellnum)
        std_mutant_cellnum = np.std(mutant_cellnum)
        theta = mean_mutant_cellnum + 2 * std_mutant_cellnum
        ### Manually flipping false positive
        correct_bin_withoutNA = in_bin_withoutNA.copy()
        correct_bin_withNA3 = in_bin_withNA3.copy()
        correct_posterior = pending_posterior.copy()
        if args.correct_fp_manually != "no":
            df_correct_unit = pd.read_csv(args.correct_fp_manually, sep='\t', index_col=False)
            for idx,unit in df_correct_unit.iterrows():
                correct_bin_withoutNA.loc[unit['cellID'], unit['mutID']] = 0
                correct_bin_withNA3.loc[unit['cellID'], unit['mutID']] = 0
                correct_posterior.loc[unit['cellID'], unit['mutID']] = 0
        correct_bin_withoutNA.to_csv(outputpath+"/correct_bin_withoutNA.txt", sep="\t")
        correct_bin_withNA3.to_csv(outputpath+"/correct_bin_withNA3_for_circosPlot.txt", sep="\t")
        correct_posterior.to_csv(outputpath+"/correct_posterior.txt", sep="\t")
        # Get shared mutations
        shared_muts_list = pending_features.loc['mutant_cellnum'][pending_features.loc['mutant_cellnum'] >= 2].index.tolist()
        shared_bin_withoutNA = correct_bin_withoutNA.loc[:,shared_muts_list]
        shared_bin_withNA3 = correct_bin_withNA3.loc[:,shared_muts_list]
        shared_posterior = correct_posterior.loc[:,shared_muts_list]
        shared_features = pending_features.loc[:,shared_muts_list]
        shared_bin_withoutNA.to_csv(outputpath+"/shared_bin_withoutNA.txt", sep="\t")
        shared_bin_withNA3.to_csv(outputpath+"/shared_bin_withNA3_for_circosPlot.txt", sep="\t")
        shared_posterior.to_csv(outputpath+"/shared_posterior.txt", sep="\t")
        shared_features.to_csv(outputpath+"/shared_features.txt", sep="\t")
        print("The number of candidate sites for building pivot tree is: "+str(len(shared_muts_list)))
        # Ultimate mutations for tree building: pivot mutations + shared mutations - germline mutaions
        [M_lowR, df_subclone, flipping_penalty_total, pivot_sites, pivot_mutid, df_traversal_early_recoder, filtered_useful_cells_indices] = minFlippingPtree(correct_bin_withoutNA.values, correct_posterior, theta, pending_features, cutoff, fnfp_ratio_subclone, fnfp_ratio_late, general_weight_na)
        # [M_lowR, df_subclone, flipping_penalty_total, pivot_sites, pivot_mutid, df_traversal_early_recoder, filtered_useful_cells_indices] = minFlippingPtree(shared_bin_withoutNA.values, shared_posterior, theta, shared_features, cutoff, fnfp_ratio_subclone, fnfp_ratio_late, general_weight_na)
        print("The number of mutations in the pivot tree is: "+str(M_lowR.shape[1]))
        df_subclone.to_csv(outputpath+"/df_subclone.txt", sep="\t")
        df_traversal_early_recoder.to_csv(outputpath+"/df_traversal_early_recoder.txt", sep="\t")
        out_posterior = correct_posterior.iloc[filtered_useful_cells_indices]
        WriteTfile(outputpath+"/M_lowR_basedPivots.filtered_sites_inferred", M_lowR, out_posterior.index.tolist(), pivot_mutid, judge="yes")
        ##### clean M_lowR: remove all zeros columns(muts) or rows(cells)
        # clean columns
        zero_columns = np.all(M_lowR == 0, axis=0)
        zero_columns_indices = np.where(zero_columns)[0]
        keep_columns = np.ones(M_lowR.shape[1], dtype=bool)
        keep_columns[zero_columns_indices] = False
        # clean rows
        zero_rows = np.all(M_lowR == 0, axis=1)
        zero_rows_indices = np.where(zero_rows)[0]
        keep_rows = np.ones(M_lowR.shape[0], dtype=bool)
        keep_rows[zero_rows_indices] = False
        # out matrix
        final_cleaned_M_lowR = M_lowR[keep_rows, :][:, keep_columns]
        pivot_mutid_clean = [col for idx, col in enumerate(pivot_mutid) if idx not in zero_columns_indices]
        pivot_cellid_clean = [row for idx, row in enumerate(out_posterior.index.tolist()) if idx not in zero_rows_indices]
        print("===================================================================")
        print("The final cleaned tree structrue is: "+str(len(pivot_cellid_clean))+" cells x "+str(len(pivot_mutid_clean))+" mutations.")
        print("===================================================================")
        # The number of mutations in the pivot tree is: 59
        WriteTfile(outputpath+"/final_cleaned_M_lowR_basedPivots.filtered_sites_inferred", final_cleaned_M_lowR, pivot_cellid_clean, pivot_mutid_clean, judge="yes")
        ### refine finale_cleaned binary matrix
        final_cleaned_bin_withoutNA = correct_bin_withoutNA.loc[pivot_cellid_clean, pivot_mutid_clean]
        final_cleaned_bin_withNA3 = correct_bin_withNA3.loc[pivot_cellid_clean, pivot_mutid_clean]
        final_cleaned_posterior = out_posterior.loc[pivot_cellid_clean, pivot_mutid_clean]
        final_cleaned_bin_withoutNA.to_csv(outputpath+"/final_cleaned_bin_withoutNA.txt", sep="\t")
        final_cleaned_bin_withNA3.to_csv(outputpath+"/final_cleaned_bin_withNA3_for_circosPlot.txt", sep="\t")
        final_cleaned_posterior.to_csv(outputpath+"/final_cleaned_posterior.txt", sep="\t")
        # ##### Get subtree + nodeid + spots
        # refer: /storage/douyanmeiLab/yangqing/tools/PhyloMosaicGenie/pmg/spatial/treecode_49muts_0904/spatial_phylo_49muts.2featuresExtract.debug.py
        ### Regenerate M_lowR_withNone by adding some spots/cells with greater than or equal to 3 total alleles but no mutant allele.
        df_phylogeny_pivot = pd.DataFrame(final_cleaned_M_lowR)
        df_phylogeny_pivot.index = final_cleaned_bin_withNA3.index.tolist()
        df_phylogeny_pivot.columns = pivot_mutid
        if args.extra_NoneCluster_cutoff == "no":
            df_phylogeny = df_phylogeny_pivot.copy()
            df_bin_withNA3_for_flipping = final_cleaned_bin_withNA3.copy()
        else:
            extra_NoneCluster_cutoff = int(args.extra_NoneCluster_cutoff)
            # Get >=3 dp but no mutAllele spots
            df_spots_high_cov_but_no_mut = get_valid_spot_ids(in_reads.iloc[1:,], extra_NoneCluster_cutoff)
            all_NoneCluster_spots = list(set().union(*df_spots_high_cov_but_no_mut))
            unique_spots_only_NoneCluster = [s for s in all_NoneCluster_spots if s not in final_cleaned_posterior.index.tolist()]
            # extract
            df_NoneCluster_posterior = in_posterior.loc[unique_spots_only_NoneCluster, :]
            df_NoneCluster_reads = in_reads.loc[['bulk'] + unique_spots_only_NoneCluster, :]
            # get binary matrix from NoneCluster posterior matrix
            df_NoneCluster_bin_withNA3 = posterior2ter_NAto3_bothPosteriorMutallele(df_NoneCluster_posterior, df_NoneCluster_reads, cutoff)
            df_NoneCluster_bin_withoutNA = posterior2bin_bothPosteriorMutallele(df_NoneCluster_posterior, df_NoneCluster_reads, cutoff)
            ## Combine final_cleaned_M_lowR and NoneCluster
            final_cleaned_bin_withNA3_addCol = final_cleaned_bin_withNA3.copy()
            final_cleaned_bin_withNA3_addCol['unmutant_cluster'] = np.zeros(len(final_cleaned_bin_withNA3_addCol), dtype=int)
            final_cleaned_bin_withoutNA_addCol = final_cleaned_bin_withoutNA.copy()
            final_cleaned_bin_withoutNA_addCol['unmutant_cluster'] = np.zeros(len(final_cleaned_bin_withNA3_addCol), dtype=int)
            df_NoneCluster_bin_withNA3_addCol = df_NoneCluster_bin_withNA3.copy()
            df_NoneCluster_bin_withNA3_addCol['unmutant_cluster'] = np.ones(len(df_NoneCluster_bin_withNA3_addCol), dtype=int)
            df_NoneCluster_bin_withoutNA_addCol = df_NoneCluster_bin_withoutNA.copy()
            df_NoneCluster_bin_withoutNA_addCol['unmutant_cluster'] = np.ones(len(df_NoneCluster_bin_withoutNA_addCol), dtype=int)
            # combine
            df_binary_withNA3_correct_with_NoneCluster = pd.concat([final_cleaned_bin_withNA3, df_NoneCluster_bin_withNA3], axis=0)
            df_binary_withoutNA_correct_with_NoneCluster = pd.concat([final_cleaned_bin_withoutNA, df_NoneCluster_bin_withoutNA], axis=0)
            df_binary_withNA3_correct_with_NoneCluster.to_csv(outputpath+"/correct_with_NoneCluster_binary_withNA3_for_circosPlot.txt", sep="\t")
            df_binary_withoutNA_correct_with_NoneCluster.to_csv(outputpath+"/correct_with_NoneCluster_binary_withoutNA.txt", sep="\t")
            # combine (addCol)
            df_binary_withNA3_correct_with_NoneCluster_addCol = pd.concat([final_cleaned_bin_withNA3_addCol, df_NoneCluster_bin_withNA3_addCol], axis=0)
            df_binary_withoutNA_correct_with_NoneCluster_addCol = pd.concat([final_cleaned_bin_withoutNA_addCol, df_NoneCluster_bin_withoutNA_addCol], axis=0)
            df_binary_withNA3_correct_with_NoneCluster_addCol.to_csv(outputpath+"/correct_with_NoneCluster_addCol_binary_withNA3_for_circosPlot.txt", sep="\t")
            df_binary_withoutNA_correct_with_NoneCluster_addCol.to_csv(outputpath+"/correct_with_NoneCluster_addCol_binary_withoutNA.txt", sep="\t")
            ## Combined tree
            df_phylogeny_pivot_addCol = df_phylogeny_pivot.copy()
            df_phylogeny_pivot_addCol['unmutant_cluster'] = np.zeros(len(df_phylogeny_pivot_addCol), dtype=int)
            df_phylogeny = pd.concat([df_phylogeny_pivot_addCol, df_NoneCluster_bin_withoutNA_addCol], axis=0)
            df_phylogeny.index.name = "cellIDxmutID"
            # df_phylogeny.to_csv(outputpath+"/final_cleaned_M_lowR_correct_with_NoneCluster_basedPivots.filtered_sites_inferred.CFMatrix", sep="\t")
            WriteTfile(outputpath+"/final_cleaned_M_lowR_correct_with_NoneCluster_basedPivots.filtered_sites_inferred", df_phylogeny.values, df_binary_withNA3_correct_with_NoneCluster.index.tolist(), pivot_mutid+['unmutant_cluster'], judge="yes")
            df_bin_withNA3_for_flipping = df_binary_withNA3_correct_with_NoneCluster_addCol.copy()
        ### Step3: Calculate somatic posterior based on tree ###
        if args.is_update_somatic_posterior == "yes":
            # [df_subclone_lessCol, M_pivot, df_newSP_allNodes_prod, df_newSP_prod_out, df_newSP_allNodes_BB, df_newSP_BB_out, withoutTree_posterior_prod, withoutTree_posterior_BB] = all_newSomaticPosterior(df_likelihoods_mut_prod, df_likelihoods_mut_BB, df_likelihoods_unmut, df_grnd.values.astype(int))
            # df_subclone_lessCol.to_csv(outputpath+"/df_subclone_lessCol.txt", sep="\t")
            df_newSP_allNodes.to_csv(outputpath+"/df_newSP_allNodes.txt", sep="\t")
            WriteTfile(outputpath+"/M_pivot_basedPivots.filtered_sites_inferred", M_pivot, no_zeros_ones_posterior.index.tolist(), [str(i+1) for i in range(M_pivot.shape[1])], judge="yes")
        ### Step4: Features generation
        out_features = pending_features.T.drop(['somatic_posterior_persite'], axis=1)
        out_features['somatic_posterior_per_site_old'] = no_zeros_ones_features.T['somatic_posterior_persite']
        out_features.to_csv(outputpath+"/out_features.somatic_posterior_basedTree.txt", sep="\t")
        ### Step4: Add flipping counts (flipping:0->1; flipping:1->0; flipping:NA->0; flipping:NA->1)
        ## compare to in_bin_withNA3
        df_flip_counts_tree = calculate_flip_counts_per_site(df_bin_withNA3_for_flipping, df_phylogeny)
        df_flip_counts_tree = df_flip_counts_tree.reindex(index=list(out_features.index))
        df_flip_counts_tree.columns = [f'tree_{flip_type}' for flip_type in df_flip_counts_tree.columns.tolist()]
        assert df_flip_counts_tree.index.equals(out_features.index), "Indexes of df_flip_counts_tree and out_features do not match."
        # 合并 df_flip_counts_tree 和 df_flip_counts_prod
        df_combined = pd.concat([out_features, df_flip_counts_tree], axis=1)
        df_combined.to_csv(outputpath+"/out_features.somatic_posterior_basedTree.txt", sep="\t")
        ### Step5: Output fliiping cells or spots
        # get false_negative_flipping spots
        false_negative_flipping_spots = df_bin_withNA3_for_flipping.apply(
            lambda col: find_flipping_spots(col, df_phylogeny[col.name], condition_in_bin=0, condition_phylogeny=1)
        )
        # get NAto1 spots
        NAto1_flipping_spots = df_bin_withNA3_for_flipping.apply(
            lambda col: find_flipping_spots(col, df_phylogeny[col.name], condition_in_bin=3, condition_phylogeny=1)
        )
        # get false_positive_flipping spots
        false_positive_flipping_spots = df_bin_withNA3_for_flipping.apply(
            lambda col: find_flipping_spots(col, df_phylogeny[col.name], condition_in_bin=1, condition_phylogeny=0)
        )
        # get NAto0 spots
        NAto0_flipping_spots = df_bin_withNA3_for_flipping.apply(
            lambda col: find_flipping_spots(col, df_phylogeny[col.name], condition_in_bin=3, condition_phylogeny=0)
        )
        # process na list
        if false_negative_flipping_spots.empty:
            false_negative_flipping_spots = {col: [] for col in df_bin_withNA3_for_flipping.columns}
        if false_positive_flipping_spots.empty:
            false_positive_flipping_spots = {col: [] for col in df_bin_withNA3_for_flipping.columns}
        if NAto1_flipping_spots.empty:
            NAto1_flipping_spots = {col: [] for col in df_bin_withNA3_for_flipping.columns}
        if NAto0_flipping_spots.empty:
            NAto0_flipping_spots = {col: [] for col in df_bin_withNA3_for_flipping.columns}
        # get flipping_spots dataframe
        df_flipping_spots = pd.DataFrame({
            'Mutation': df_bin_withNA3_for_flipping.columns,
            'false_negative_flipping_spots': [', '.join(false_negative_flipping_spots.get(col, [])) for col in df_bin_withNA3_for_flipping.columns],
            'false_positive_flipping_spots': [', '.join(false_positive_flipping_spots.get(col, [])) for col in df_bin_withNA3_for_flipping.columns],
            'NAto1_flipping_spots': [', '.join(NAto1_flipping_spots.get(col, [])) for col in df_bin_withNA3_for_flipping.columns],
            'NAto0_flipping_spots': [', '.join(NAto0_flipping_spots.get(col, [])) for col in df_bin_withNA3_for_flipping.columns]
        })
        df_flipping_spots.to_csv(outputpath+"/df_flipping_spots.txt", sep="\t", index=False)
        ### Step5: total flipping count
        total_FN_flipping = ((df_bin_withNA3_for_flipping == 0) & (df_phylogeny == 1)).sum().sum()
        total_FP_flipping = ((df_bin_withNA3_for_flipping == 1) & (df_phylogeny == 0)).sum().sum()
        total_NAto0 = ((df_bin_withNA3_for_flipping == 3) & (df_phylogeny == 0)).sum().sum()
        total_NAto1 = ((df_bin_withNA3_for_flipping == 3) & (df_phylogeny == 1)).sum().sum()
        print("Total False Negative flipping count: " + str(total_FN_flipping))
        print("Total False Positive flipping count: " + str(total_FP_flipping))
        print("Total NA to 0 flipping count: " + str(total_NAto0))
        print("Total NA to 1 flipping count: " + str(total_NAto1))
        df_total_flipping_count = pd.DataFrame({
            'total_flipping_False_Negative': [total_FN_flipping],
            'total_flipping_False_Positive': [total_FP_flipping],
            'total_flipping_NA_to_0': [total_NAto0],
            'total_flipping_NA_to_1': [total_NAto1]
            })
        df_total_flipping_count.to_csv(outputpath+"/df_total_flipping_count.txt", sep="\t", index=False)
##### Usage2: Validate candidate somatic sites #####
elif args.phylo_or_valid=="valid":
    df_grnd = pd.read_csv(args.grnd, sep='\t', index_col=0)
    orderd_index = df_grnd.index
    in_posterior = in_posterior.reindex(index=orderd_index)
    in_llmut = in_llmut.reindex(index=orderd_index)
    in_llunmut = in_llunmut.reindex(index=orderd_index)
    in_reads = in_reads.reindex(index=['bulk'] + list(orderd_index))
    in_reads_noBulk = in_reads.reindex(index=orderd_index)
    in_baseq = in_baseq.reindex(index=orderd_index)
    ### normalize input likelihoods
    combined_df = pd.concat([in_llmut, in_llunmut], axis=1)
    normalized_result = apply_normalization(combined_df, is_log=args.is_log_value_for_likelihoods)
    normalized_llmut_BB = normalized_result.filter(like='norm_llmut')
    normalized_llmut_BB.columns = in_llmut.columns
    normalized_llunmut_BB = normalized_result.filter(like='norm_llunmut')
    normalized_llunmut_BB.columns = in_llunmut.columns
    ### Step0: calculate likelihoods by product and beta-binomial methods
    # per site
    df_mutid = pd.DataFrame(in_features.columns, columns=["mutid"])
    # Caculate types of likelihoods
    df_alpha_beta = in_alpha_beta[['alpha', 'beta']]
    df_alpha_beta.index = in_alpha_beta['scid']
    df_alpha_beta = df_alpha_beta.reindex(index=orderd_index)
    results = df_mutid['mutid'].apply(lambda mutid: calculate_likelihood_per_site(mutid, in_posterior, orderd_index, in_reads_noBulk, in_features, in_baseq, df_alpha_beta))
    ## Collect likelihood dataframe
    # mut by product
    df_likelihoods_mut_prod = pd.DataFrame([result[0] for result in results]).T
    df_likelihoods_mut_prod.columns = df_mutid['mutid']
    df_likelihoods_mut_prod.index = in_posterior.index
    df_likelihoods_mut_prod.to_csv(outputpath+"/df_likelihoods_mut_prod.txt", sep="\t")
    # unmut by product
    df_likelihoods_unmut_prod = pd.DataFrame([result[1] for result in results]).T
    df_likelihoods_unmut_prod.columns = df_mutid['mutid']
    df_likelihoods_unmut_prod.index = in_posterior.index
    df_likelihoods_unmut_prod.to_csv(outputpath+"/df_likelihoods_unmut_prod.txt", sep="\t")
    # mut by beta-binomial
    df_likelihoods_mut_BB = pd.DataFrame([result[2] for result in results]).T
    df_likelihoods_mut_BB.columns = df_mutid['mutid']
    df_likelihoods_mut_BB.index = in_posterior.index
    df_likelihoods_mut_BB.to_csv(outputpath+"/df_likelihoods_mut_BB.txt", sep="\t")
    # unmut by beta-binomial
    df_likelihoods_unmut_BB = pd.DataFrame([result[3] for result in results]).T
    df_likelihoods_unmut_BB.columns = df_mutid['mutid']
    df_likelihoods_unmut_BB.index = in_posterior.index
    df_likelihoods_unmut_BB.to_csv(outputpath+"/df_likelihoods_unmut_BB.txt", sep="\t")
    ### Step1: Calculate somatic posterior based on tree ###
    # [df_subclone, M_pivot, df_newSP_allNodes_prod, df_newSP_prod_out, df_newSP_allNodes_BB, df_newSP_BB_out, withoutTree_posterior_prod, withoutTree_posterior_BB]
    [df_subclone_lessCol, M_pivot, df_newSP_allNodes_prod, df_newSP_prod_out, df_newSP_allNodes_BB, df_newSP_BB_out, withoutTree_posterior_prod, withoutTree_posterior_BB] = all_newSomaticPosterior(df_likelihoods_mut_prod, df_likelihoods_unmut_prod, normalized_llmut_BB, normalized_llunmut_BB, df_grnd.values.astype(int))
    df_subclone_lessCol.to_csv(outputpath+"/df_subclone_lessCol.txt", sep="\t")
    df_newSP_allNodes_prod.to_csv(outputpath+"/df_newSP_allNodes_prod.txt", sep="\t")
    df_newSP_allNodes_BB.to_csv(outputpath+"/df_newSP_allNodes_BB.txt", sep="\t")
    WriteTfile(outputpath+"/M_pivot_basedPivots.filtered_sites_inferred", M_pivot, in_posterior.index.tolist(), [str(i+1) for i in range(M_pivot.shape[1])], judge="yes")
    ### Step2: Features generation
    # Rename df_newSP in order to combine dataframes
    df_newSP_prod_out.rename(columns=lambda x: 'prod_' + x, inplace=True)
    df_newSP_BB_out.rename(columns=lambda x: 'BB_' + x, inplace=True)
    # Process and combine features dataframe
    out_features = in_features.T.drop(['somatic_posterior_persite'], axis=1)
    out_features['somatic_posterior_per_site_old'] = in_features.T['somatic_posterior_persite']
    out_features['withoutTree_posterior_prod'] = withoutTree_posterior_prod
    out_features['withoutTree_posterior_BB'] = withoutTree_posterior_BB
    out_features = pd.concat([out_features, df_newSP_prod_out, df_newSP_BB_out], axis=1)
    ### Step3: Make phylogeny labels
    # Ensure that all relevant columns are numeric
    out_features['prod_somatic_posterior_per_site'] = out_features['prod_somatic_posterior_per_site'].astype(float)
    out_features['BB_somatic_posterior_per_site'] = out_features['BB_somatic_posterior_per_site'].astype(float)
    out_features['prod_somatic_posterior_per_site_onecell'] = out_features['prod_somatic_posterior_per_site_onecell'].astype(float)
    out_features['BB_somatic_posterior_per_site_onecell'] = out_features['BB_somatic_posterior_per_site_onecell'].astype(float)
    out_features['mutant_cellnum'] = out_features['mutant_cellnum'].astype(int)
    out_features['phylogeny_label'] = out_features.apply(determine_phylogeny_label, axis=1, args=(args.pass_tree_cutoff,args.unpass_tree_cutoff,))
    ### Step4: Add flipping counts (flipping:0->1; flipping:1->0; flipping:NA->0; flipping:NA->1)
    df_binary_withNA = posterior2ter_NAto3(in_posterior, in_reads, args.threshold)[3]
    df_binary_withNA.to_csv(outputpath+"/df_binary_withNA3_for_circosPlot.txt", sep="\t")
    # prod
    df_prod_cluster = pd.DataFrame(out_features['prod_cluster'].tolist(), index=out_features.index).T
    df_prod_cluster.index = in_posterior.index
    df_flip_counts_prod = calculate_flip_counts_per_site(df_binary_withNA, df_prod_cluster)
    df_flip_counts_prod.columns = [f'prod_{flip_type}' for flip_type in df_flip_counts_prod.columns.tolist()]
    # BB
    df_BB_cluster = pd.DataFrame(out_features['BB_cluster'].tolist(), index=out_features.index).T
    df_BB_cluster.index = in_posterior.index
    df_flip_counts_BB = calculate_flip_counts_per_site(df_binary_withNA, df_BB_cluster)
    df_flip_counts_BB.columns = [f'BB_{flip_type}' for flip_type in df_flip_counts_BB.columns.tolist()]
    # 确保 df_flip_counts_BB 和 df_flip_counts_prod 的索引与 out_features 的索引一致
    assert df_flip_counts_BB.index.equals(out_features.index), "Indexes of df_flip_counts_BB and out_features do not match."
    assert df_flip_counts_prod.index.equals(out_features.index), "Indexes of df_flip_counts_prod and out_features do not match."
    # 合并 df_flip_counts_BB 和 df_flip_counts_prod
    df_combined = pd.concat([out_features, df_flip_counts_prod, df_flip_counts_BB], axis=1)
    # Save file
    df_combined.to_csv(outputpath+"/out_features.somatic_posterior_basedTree.txt", sep="\t")


if args.labelfile:
    # If given labels, test features between groups.
    # args.labelfile="/storage/douyanmeiLab/yangqing/tools/PhyloMosaicGenie/pmg/input/label_78true_withAllele.txt"
    label_df = pd.read_csv(args.labelfile, sep='\t', index_col=0)
    # out_features = out_features.merge(label_df, left_index=True, right_index=True, how='left')
    tp_in_df_label = label_df[label_df['label'] == 1].index
    out_features['label'] = None
    for mutid in tp_in_df_label:
        out_features.loc[out_features.index.str.contains(mutid), 'label'] = label_df.loc[mutid, 'label']
    out_features['label'].fillna(0, inplace=True)
    out_features.to_csv(outputpath+"/out_features.somatic_posterior_basedTree.txt", sep="\t")
    # Feature test
    tp_mutid = list(out_features[out_features['label'] == 1].index)
    fp_mutid = [x for x in in_posterior.columns if x not in tp_mutid]
    # out_features.columns
    # Index(['indid', 'chr', 'pos', 'ref', 'mut', 'genotype', 'bulk_cov', 'bulk_VAF',
    #        'avg_cov_percell', 'avg_mutAF_percell', 'max_mutAF_percell',
    #        'cov_in_maxmutAFcell', 'max_mutAllele_num', 'total_mutAllele_num',
    #        'mutant_cell_fraction', 'mutant_cellnum', 'mutAllele_cellnum',
    #        'moreCutoff_cellnum', 'somatic_posterior_per_site_old',
    #        'withoutTree_posterior_prod', 'withoutTree_posterior_BB',
    #        'prod_node_position', 'prod_cluster', 'prod_branch_state',
    #        'prod_somatic_posterior_per_site_max',
    #        'prod_somatic_posterior_per_site_onecell',
    #        'prod_somatic_posterior_per_site', 'prod_unmutant_posterior_per_site',
    #        'prod_wrong_posterior_per_site', 'BB_node_position', 'BB_cluster',
    #        'BB_branch_state', 'BB_somatic_posterior_per_site_max',
    #        'BB_somatic_posterior_per_site_onecell', 'BB_somatic_posterior_per_site',
    #        'BB_unmutant_posterior_per_site', 'BB_wrong_posterior_per_site',
    #        'phylogeny_label', 'label'],
    #       dtype='object')
    non_test_terms = ['indid', 'chr', 'pos', 'ref', 'mut', 'genotype', 'prod_node_position', 'prod_cluster', 'prod_branch_state', 'BB_node_position', 'BB_cluster', 'BB_branch_state', 'phylogeny_label', 'label']
    df_feature_test = test_features(out_features.drop(non_test_terms, axis=1), tp_mutid, fp_mutid, kind="violinplot")


##### Time #####
finish_time = time.perf_counter()
print("Program finished in {:.4f} seconds".format(finish_time-start_time))

