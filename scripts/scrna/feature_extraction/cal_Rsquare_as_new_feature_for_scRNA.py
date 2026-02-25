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
    判断一个位点是否可能由于 allele dropout 引起。
    参数：
        df (pd.DataFrame): 包含 'total_dp' 和 'alt_dp' 两列的数据
        thr_r2 (float): R² 的阈值，小于此值表示拟合不好
        alpha (float): Wilcoxon 配对检验的显著性水平
    返回：
        no_alleledrop (bool): True 表示无 allele dropout，False 表示可能存在
        r_squared (float): 标准化后的 R² 值
        wilcoxon_pval (float): Wilcoxon 检验的 p 值
    """
    depth_data = pd.DataFrame({
        'total_dp': df['total_dp'],
        'alt_dp': df['alt_dp']
    })
    
    # 空数据判断
    if depth_data.empty:
        return "undefined", np.nan, np.nan
    
    # 标准化
    standard_scaler = StandardScaler()
    depth_std = standard_scaler.fit_transform(depth_data)
    depth_std = pd.DataFrame(depth_std, columns=['total_dp', 'alt_dp'])
    
    # 差值计算
    diff = depth_std['total_dp'] - depth_std['alt_dp']
    
    # R² 计算
    if len(depth_std['total_dp'].unique()) > 1:
        _, _, r_value, _, _ = stats.linregress(depth_std['total_dp'], depth_std['alt_dp'])
        r_squared = r_value**2
    else:
        r_squared = 1.0
    
    # Wilcoxon 配对检验
    if len(diff.unique()) > 1:
        try:
            wilcoxon_stat, wilcoxon_pval = wilcoxon(depth_std['total_dp'], depth_std['alt_dp'])
        except ValueError:
            # Wilcoxon 检验异常时 fallback
            wilcoxon_stat = 0
            wilcoxon_pval = 1.0
    else:
        wilcoxon_stat = 0
        wilcoxon_pval = 1.0
    
    # 判断是否为 allele dropout（False 表示疑似有 dropout）
    no_alleledrop = (r_squared < thr_r2) or (wilcoxon_pval < alpha)
    # # 将布尔值转换为字符串标签
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
    
    # 文件路径存在性检查
    if not os.path.exists(reads_file):
        print(f"Warning: File {reads_file} does not exist!")
        return None, None, None
    
    df_reads = pd.read_csv(reads_file, sep="\t", header=None, names=["barcode", "total_dp", "alt_dp"])
    is_no_alleledrop_persite, r_squared_persite, wilcoxon_pval_persite = check_allele_dropout(df_reads, thr_r2=0.8, alpha=0.05)
    
    return is_no_alleledrop_persite, r_squared_persite, wilcoxon_pval_persite

# 并行处理
with mp.Pool(processes=args.thread) as pool:
    results = pool.map(process_identifier, identifier_list)

# 结果处理
for result in results:
    if result != (None, None, None):
        is_no_alleledrop_list.append(result[0])
        r_squared_list.append(result[1])
        wilcoxon_pval_list.append(result[2])
    else:
        is_no_alleledrop_list.append(None)
        r_squared_list.append(None)
        wilcoxon_pval_list.append(None)


# 将新生成的特征列添加到 df_features 中
out_features = df_features.copy()

out_features['is_no_alleledrop_based_on_expression'] = is_no_alleledrop_list
out_features['r_squared_mutant_vs_expression'] = r_squared_list
out_features['wilcoxon_pval_mutant_vs_expression'] = wilcoxon_pval_list

out_features.to_csv(output_file, sep="\t", index=False)


