# Date: 2025/04/09
# Author: Qing Yang, Mengdie Yao
# Work: Calculate R square between mutant allele number and expression.


##### Time #####
import time
start_time = time.perf_counter()


##### input para
import multiprocessing as mp
import argparse
from argparse import ArgumentParser
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--feature_file", default="/storage/douyanmeiLab/yangqing/tools/PhyloMosaicGenie/Benchmark/data_A549/PhyloSOLID_tree/mosaic_mutations/features_test/test.identifier.feature_final.txt", type=str, help="The feature file.")
parser.add_argument("-r", "--reads_filepath", default="/storage/douyanmeiLab/yangqing/tools/PhyloMosaicGenie/Benchmark/data_A549/PhyloSOLID_tree/mosaic_mutations/features_100k/depth_in_spots/", type=str, help="The reads information file path.")
parser.add_argument("-o", "--output_file", default="/storage/douyanmeiLab/yangqing/tools/PhyloMosaicGenie/Benchmark/data_A549/PhyloSOLID_tree/mosaic_mutations/features_test/test.identifier.feature_mutation_vs_expression.txt", type=str, help="The outputpath.")
parser.add_argument("-t", "--thread", default=mp.cpu_count(), type=int, help="Cpu count.")
args = parser.parse_args()


##### Load libraies
import os
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import wilcoxon
from sklearn.preprocessing import StandardScaler


##### Parameters
feature_file = args.feature_file
reads_filepath = args.reads_filepath
output_file = args.output_file


##### Function
def check_allele_dropout(df, thr_r2=0.8, alpha=0.05):
    """
    Assess whether a site may be caused by allele dropout.
    Parameters:
        df (pd.DataFrame): data with 'total_dp' and 'alt_dp' columns
        thr_r2 (float): R-squared threshold; values below this indicate a poor fit
        alpha (float): significance level for the paired Wilcoxon test
    Returns:
        no_alleledrop (bool): True means no allele dropout, False means dropout is possible
        r_squared (float): standardized R-squared
        wilcoxon_pval (float): Wilcoxon test p-value
    """
    depth_data = pd.DataFrame({
        'total_dp': df['total_dp'],
        'alt_dp': df['alt_dp']
    })
    
    # Empty-data check
    if depth_data.empty:
        return "undefined", np.nan, np.nan
    
    # Standardization
    standard_scaler = StandardScaler()
    depth_std = standard_scaler.fit_transform(depth_data)
    depth_std = pd.DataFrame(depth_std, columns=['total_dp', 'alt_dp'])
    
    # Difference
    diff = depth_std['total_dp'] - depth_std['alt_dp']
    
    # R-squared
    if len(depth_std['total_dp'].unique()) > 1:
        _, _, r_value, _, _ = stats.linregress(depth_std['total_dp'], depth_std['alt_dp'])
        r_squared = r_value**2
    else:
        r_squared = 1.0
    
    # Paired Wilcoxon test
    if len(diff.unique()) > 1:
        try:
            wilcoxon_stat, wilcoxon_pval = wilcoxon(depth_std['total_dp'], depth_std['alt_dp'])
        except ValueError:
            # Fallback when the Wilcoxon test fails
            wilcoxon_stat = 0
            wilcoxon_pval = 1.0
    else:
        wilcoxon_stat = 0
        wilcoxon_pval = 1.0
    
    # Decide allele dropout (False means dropout is suspected)
    no_alleledrop = (r_squared < thr_r2) or (wilcoxon_pval < alpha)
    # # Convert the boolean to a string label
    # dropout_qc_result = "pass" if no_alleledrop else "fail"
    
    return no_alleledrop, r_squared, wilcoxon_pval


##### Load feature file
df_features = pd.read_csv(feature_file, sep="\t")
identifier_list = list(df_features['identifier'])


##### depth info and cal R² 
is_no_alleledrop_list = []
r_squared_list = []
wilcoxon_pval_list = []

# for identifier in identifier_list:
#     reads_file = reads_filepath+"/"+identifier+".mut.spots.txt"
#     df_reads = pd.read_csv(reads_file, sep="\t", header=None, names=["barcode", "total_dp", "alt_dp"])
#     is_no_alleledrop_persite, r_squared_persite, wilcoxon_pval_persite = check_allele_dropout(df_reads, thr_r2=0.8, alpha=0.05)
#     is_no_alleledrop_list.append(is_no_alleledrop_persite)
#     r_squared_list.append(r_squared_persite)
#     wilcoxon_pval_list.append(wilcoxon_pval_persite)

def process_identifier(identifier):
    reads_file = reads_filepath + "/" + identifier + ".mut.spots.txt"
    
    # Check that the file exists
    if not os.path.exists(reads_file):
        print(f"Warning: File {reads_file} does not exist!")
        return None, None, None
    
    df_reads = pd.read_csv(reads_file, sep="\t", header=None, names=["barcode", "total_dp", "alt_dp"])
    is_no_alleledrop_persite, r_squared_persite, wilcoxon_pval_persite = check_allele_dropout(df_reads, thr_r2=0.8, alpha=0.05)
    
    return is_no_alleledrop_persite, r_squared_persite, wilcoxon_pval_persite

# Parallel processing
with mp.Pool(processes=args.thread) as pool:
    results = pool.map(process_identifier, identifier_list)

# Collect results
for result in results:
    if result != (None, None, None):
        is_no_alleledrop_list.append(result[0])
        r_squared_list.append(result[1])
        wilcoxon_pval_list.append(result[2])
    else:
        is_no_alleledrop_list.append(None)
        r_squared_list.append(None)
        wilcoxon_pval_list.append(None)


# Add the newly computed feature columns to df_features
out_features = df_features.copy()

out_features['is_no_alleledrop_based_on_expression'] = is_no_alleledrop_list
out_features['r_squared_mutant_vs_expression'] = r_squared_list
out_features['wilcoxon_pval_mutant_vs_expression'] = wilcoxon_pval_list

out_features.to_csv(output_file, sep="\t", index=False)


