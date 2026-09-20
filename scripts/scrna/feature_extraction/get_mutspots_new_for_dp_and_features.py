import argparse
from collections import defaultdict
from functools import reduce
from math import ceil, log10
import pysam
import os
import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import ttest_ind  # T-test for differential analysis
import multiprocessing
from UMI_combine import calculate_UMI_combine_phred, get_most_candidate_allele
from myutils import check_dir
from concurrent.futures import ProcessPoolExecutor, as_completed
from scipy.stats import mannwhitneyu
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


def handle_cigar(ciagr_symbol):
    '''
    ## handel cigar
    # [(0, 76), (2, 1), (0, 33), (3, 139241), (0, 11)]
    # '76M1D33M139241N11M'
    # the 1st is symbol; and the 2nd is count
    # 0: Match; 1: Insertion; 2: deletion; 3: N; 4: S; 5: H; 6: P; 7: =; 8: X
    output:
    set_cut may be a turple, contain the softclip information; pos_cut are the indel information: [(before_length, insert number)]
    '''
    seq_length_before = 0
    pos_length_before = 0
    
    seq_cut_start = None; seq_cut_end = None
    pos_cut=[]
    for cigars, i in zip(ciagr_symbol,range(1,len(ciagr_symbol)+1)):
        symbol = cigars[0]
        count = cigars[1]
        if symbol in [5,6,7,8]:
            # an api for handeling mapping issues "HP=X"
            print(ciagr_symbol)  ## LOG
        elif symbol in [0, 1, 4]:
            # measure the seq length 
            seq_length_before += count
            if symbol == 0:
                pos_length_before += count
            elif symbol == 4:
                # whether "S" is in this read
                if i == 1:
                    # whether the "S" is in the head or tail
                    seq_cut_start = seq_length_before
                elif i ==len(ciagr_symbol):
                    seq_cut_end = seq_length_before
                else:
                    print(ciagr_symbol) ## LOG
            elif symbol == 1:
                # whether the "I" is in the cigar
                pos_cut.append((pos_length_before,count))
        else:
            pass
    seq_cut = (seq_cut_start, seq_cut_end)
    return seq_cut, pos_cut


def handle_seq(seq, seq_cut):
    # only support one time for cut
    cut_seq=seq[seq_cut[0]:seq_cut[1]]
    return cut_seq


def handle_pos(pos_matrix,pos_cut):
    if len(pos_cut) == 0:
        cut_pos_matrix = pos_matrix
    elif len(pos_cut) == 1:
        times = pos_cut[0][1]
        pos = pos_cut[0][0]
        cut_pos_matrix = pos_matrix[0:pos] + [""] * times + pos_matrix[pos:]
    else:
        start = 0
        cut_pos_matrix = []
        for item in range(len(pos_cut)):
            pos = pos_cut[item][0]; times = pos_cut[item][1]
            cut_pos_matrix = cut_pos_matrix + pos_matrix[start:pos] + [""] * times
            start=pos
            #print(cut_pos_matrix)
        last_pos = pos_cut[-1][0]
        cut_pos_matrix = cut_pos_matrix + pos_matrix[last_pos:]
    return cut_pos_matrix


def handle_quality_matrix(mutation_in_cutseq_index,seq,cut_seq):
    if len(cut_seq[mutation_in_cutseq_index:]) >= len(cut_seq[:mutation_in_cutseq_index]):
        query_str = cut_seq[mutation_in_cutseq_index:]
        raw_index = seq.index(query_str)
    else:
        query_str = cut_seq[:mutation_in_cutseq_index]
        raw_index = seq.index(query_str) + len(query_str)
    return raw_index

def get_identifier_info_from_bam(identifier,bam_file,barcode_list,CBtag,UBtag,run_type="UMI"):
    # geno_index={"A":0,"T":1,"C":2,"G":3}
    bam_handle=pysam.AlignmentFile(bam_file,"r")
    chrom,pos,ref,alt=identifier.split("_")
    reads=bam_handle.fetch(chrom,int(pos)-1,int(pos))
    pos_index = int(pos)-1
    
    # all_reads_count=0
    # effective_DP=0
    # mapping_reads=0;skip_reads=0
    # barcode_low_quality=0   
    
    if barcode_list!=[]:
        filteration="yes"
    else:
        filteration="no"
    
    barcode_name=[]
    site_barcode_UMI_dict={}
    # counts=defaultdict(int)
    
    for item in reads:
        try:
            # a part of reads didn't have the information of "CB", bacause the "CR" didn't pass QC
            CB=item.get_tag(CBtag).strip()
            UB=item.get_tag(UBtag).strip()
        except:
                continue
        
        if filteration=="no" or (filteration=="yes" and CB in barcode_list):
            
            try:
                item.get_reference_positions().index(pos_index)
            except:
                # skip_reads += 1
                continue
            # mapping_reads += 1
            
            seq_cut, pos_cut = handle_cigar(item.cigar)
            cut_seq=handle_seq(item.seq, seq_cut)
            cut_pos=handle_pos(item.get_reference_positions(), pos_cut)
            
            if pos_index in cut_pos:
                # effective_DP += 1
                geno = cut_seq[cut_pos.index(pos_index)]
                if geno not in "ATCG":
                    continue
                
                raw_index = handle_quality_matrix(cut_pos.index(pos_index),item.seq,cut_seq)
                quality=item.get_forward_qualities()[raw_index]
                
                UMI_name = str(UB)
                barcode_name = str(CB)
                
                if barcode_name not in site_barcode_UMI_dict.keys():
                    site_barcode_UMI_dict[barcode_name]=defaultdict(dict)
                
                if UMI_name not in site_barcode_UMI_dict[barcode_name].keys():
                    site_barcode_UMI_dict[barcode_name][UMI_name]["count"]=defaultdict(int)
                    site_barcode_UMI_dict[barcode_name][UMI_name]["quality"]={"A":defaultdict(int),"T":defaultdict(int),"C":defaultdict(int),"G":defaultdict(int)}
                
                site_barcode_UMI_dict[barcode_name][UMI_name]["count"][geno]+=1
                site_barcode_UMI_dict[barcode_name][UMI_name]["quality"][geno][quality]+=1
    # print(site_barcode_UMI_dict)
    if run_type=="read":
        out_lists={}
        for CB in site_barcode_UMI_dict.keys():
            out_list=[0,0,0,0]
            for UB in site_barcode_UMI_dict[CB].keys():
                for i,geno in zip([0,1,2,3],"ATCG"):
                    out_list[i]+=site_barcode_UMI_dict[CB][UB]["count"][geno]
            out_lists[CB]=out_list
        return out_lists
    
    elif run_type=="check":
        out_lists={}
        for CB in site_barcode_UMI_dict.keys():
            for UB in site_barcode_UMI_dict[CB].keys():
                out_list=[0,0,0,0]
                for i,geno in zip([0,1,2,3],"ATCG"):
                    out_list[i]=site_barcode_UMI_dict[CB][UB]["count"][geno]
                out_lists[CB+":"+UB]=out_list
        return out_lists
    
    elif run_type=="UMI": 
        out_lists={}
        # consensus_read_count=defaultdict(int)
        for CB in site_barcode_UMI_dict.keys():
            out_list=[0,0,0,0]
            for UB in site_barcode_UMI_dict[CB].keys():
                count_dict=site_barcode_UMI_dict[CB][UB]["count"]
                quality_dict=site_barcode_UMI_dict[CB][UB]["quality"]
                phred_dict=calculate_UMI_combine_phred(count_dict,quality_dict,weigh=0.5)
                candidate_allele,phred=get_most_candidate_allele(phred_dict,ref)
                out_list["ATCG".index(candidate_allele)]+=1
            out_lists[CB]=out_list
        
        return out_lists


def handle_barcode(barcode_file,pos=0,in_tissue_choose=0):
    '''
    input: (barcode,in_tissue_or_not,simplified_location_x,simplified_location_y,real_location_x,real_location_y) #in_tissue_or_not:0 means not in tissue, 1 means in tissue
    ACGCCTGACACGCGCT-1,0,0,0,323,308
    TACCGATCCAACACTT-1,0,1,1,334,326
    
    argrument:
    pos means which kind of position you want to use, if you given a int not equal to 0, the simplified location will be used.
    
    output:
    dict: {barcode: (in_tissue, pos1, pos2),...}
    '''
    barcode_dict={}
    f=open(barcode_file,"r")
    for line in f.readlines():
        s = line.strip().split(",")
        barcode=s[0]
        if barcode!="barcode":
            in_tissue=int(s[1])
            if pos==0:
                pos1= int(s[2]); pos2=int(s[3])
            else:
                pos1= int(s[4]); pos2=int(s[5])
            
            if in_tissue_choose==0 and in_tissue==1:
                barcode_dict[barcode]= (in_tissue, pos1, pos2)
            elif in_tissue_choose==1:
                barcode_dict[barcode]= (in_tissue, pos1, pos2)
            else:pass
        else: pass
    return barcode_dict


# Function to calculate whole-genome average depth
def process_chunk(args):
    """Process one chromosome region of a BAM file"""
    bam_file, barcodes, chrom = args
    coverage_sums = defaultdict(int)
    base_counts = defaultdict(int)
    
    bam = pysam.AlignmentFile(bam_file, "rb")
    
    for read in bam.fetch(chrom):
        try:
            barcode = read.get_tag("CB")
            if barcode in barcodes and not read.is_unmapped:
                coverage_sums[barcode] += read.reference_length
                base_counts[barcode] += 1
        except KeyError:
            continue
            
    bam.close()
    return dict(coverage_sums), dict(base_counts)

def calculate_coverage_by_barcode(bam_file, barcode_file, threads=None):
    """Calculate average coverage for each barcode in parallel, chunked by chromosome"""
    if not threads:
        threads = multiprocessing.cpu_count()
        
    # Read barcodes
    with open(barcode_file) as f:
        barcodes = set(line.strip() for line in f)
    
    # Get all chromosome names
    bam = pysam.AlignmentFile(bam_file, "rb")
    chromosomes = bam.references
    bam.close()
    
    # Assign chromosomes to different threads
    args = [(bam_file, barcodes, chrom) for chrom in chromosomes]
    
    # Parallel processing
    coverage_sums = defaultdict(int)
    base_counts = defaultdict(int)
    
    with ProcessPoolExecutor(max_workers=threads) as executor:
        results = executor.map(process_chunk, args)
        
        for chunk_covs, chunk_counts in results:
            for barcode, cov in chunk_covs.items():
                coverage_sums[barcode] += cov
            for barcode, count in chunk_counts.items():
                base_counts[barcode] += count
    
    # Calculate average coverage
    avg_coverage = {}
    for barcode in barcodes:
        if base_counts[barcode] > 0:
            avg_coverage[barcode] = coverage_sums[barcode] / base_counts[barcode]
        else:
            avg_coverage[barcode] = 0
            
    return avg_coverage


# Function to update the features file, overwriting the original file
def update_features_file_inplace(features_file, stats_df):
    """
    Update newly computed features onto the existing features file by adding new columns.
    """
    # Load features file and stats file
    features_df = pd.read_csv(features_file, sep="\t", header=0)  # Load features_file
    stats_df
    
    # Print column names for debugging
    print("Features file columns:", features_df.columns)
    print("Stats file columns:", stats_df.columns)
    
    # Ensure column names are consistent: `#identifier` -> `barcode`
    features_df.rename(columns={"#identifier": "identifier"}, inplace=True)
    stats_df.rename(columns={"identifier": "identifier"}, inplace=True)
    
    # Strip whitespace from column names
    features_df.columns = features_df.columns.str.strip()
    stats_df.columns = stats_df.columns.str.strip()
    
    # Check whether the identifier column is present
    if "identifier" not in features_df.columns or "identifier" not in stats_df.columns:
        raise KeyError("Both files must contain a 'identifier' column for merging.")
    
    # Merge data
    updated_features_df = pd.merge(features_df, stats_df, how="left", on="identifier")
    
    # Rename back for consistency (do not add "#")
    # updated_features_df.rename(columns={"identifier": "#identifier"}, inplace=True)
    
    # Overwrite the original file
    feature_out = features_file.replace('.feature.txt', '.feature_depth.txt')
    updated_features_df.to_csv(feature_out, sep="\t", index=False)
    print(f"Features file updated in place: {feature_out}")


# Function for differential analysis (computation)
def process_identifier(identifier, output_folder, wg_data):
    """
    Process computation for a single identifier
    """
    file_name = identifier + ".mut.spots.txt"
    file_path = os.path.join(output_folder, file_name)
    # generate data frame
    data = pd.read_csv(file_path, sep="\t", header=None, names=["barcode", "total_dp", "alt_dp"])
    
    # Check whether the file is empty or missing
    if data.empty:
        print(f"Warning: {identifier} has no reads. Skipping calculate for this identifier.")
        columns = [
            'identifier', 'ref_dp_unmutant_avg', 'alt_dp_unmutant_avg',
            'total_dp_unmutant_avg', 'vaf_unmutant_avg', 'norm_ref_dp_unmutant_avg',
            'norm_alt_dp_unmutant_avg', 'norm_total_dp_unmutant_avg',
            'ref_dp_mutant_avg', 'alt_dp_mutant_avg', 'total_dp_mutant_avg',
            'vaf_mutant_avg', 'norm_ref_dp_mutant_avg', 'norm_alt_dp_mutant_avg',
            'norm_total_dp_mutant_avg', 'ref_dp_unmutant_max',
            'alt_dp_unmutant_max', 'total_dp_unmutant_max', 'vaf_unmutant_max',
            'norm_ref_dp_unmutant_max', 'norm_alt_dp_unmutant_max',
            'norm_total_dp_unmutant_max', 'ref_dp_mutant_max', 'alt_dp_mutant_max',
            'total_dp_mutant_max', 'vaf_mutant_max', 'norm_ref_dp_mutant_max',
            'norm_alt_dp_mutant_max', 'norm_total_dp_mutant_max',
            'ref_dp_unmutant_min', 'alt_dp_unmutant_min', 'total_dp_unmutant_min',
            'vaf_unmutant_min', 'norm_ref_dp_unmutant_min',
            'norm_alt_dp_unmutant_min', 'norm_total_dp_unmutant_min',
            'ref_dp_mutant_min', 'alt_dp_mutant_min', 'total_dp_mutant_min',
            'vaf_mutant_min', 'norm_ref_dp_mutant_min', 'norm_alt_dp_mutant_min',
            'norm_total_dp_mutant_min', 'ref_dp_mean_diff', 'alt_dp_mean_diff',
            'total_dp_mean_diff', 'vaf_mean_diff', 'norm_ref_dp_mean_diff',
            'norm_alt_dp_mean_diff', 'norm_total_dp_mean_diff', 'ref_dp_t_value',
            'alt_dp_t_value', 'total_dp_t_value', 'vaf_t_value',
            'norm_ref_dp_t_value', 'norm_alt_dp_t_value', 'norm_total_dp_t_value',
            'ref_dp_p_value', 'alt_dp_p_value', 'total_dp_p_value', 'vaf_p_value',
            'norm_ref_dp_p_value', 'norm_alt_dp_p_value', 'norm_total_dp_p_value',
            'mutant_cell_number', 'mutant_cell_fraction',
            'mutant_and_unmutant_cell_number', 'mutant_and_unmutant_cell_fraction',
            'ref_dp_in_pseudobulk', 'alt_dp_in_pseudobulk',
            'total_dp_in_pseudobulk', 'vaf_in_pseudobulk',
            'norm_ref_dp_in_pseudobulk', 'norm_alt_dp_in_pseudobulk',
            'norm_total_dp_in_pseudobulk'
        ]
        # Fill columns grouped by column-name suffix
        zero_fill_suffixes = ('_avg', '_max', '_min', '_number', '_fraction', '_in_pseudobulk')
        nan_fill_suffixes = ('_mean_diff', '_t_value', '_p_value')
        zero_fill_columns = [col for col in columns if col.endswith(zero_fill_suffixes)]
        nan_fill_columns = [col for col in columns if col.endswith(nan_fill_suffixes)]
        # Create and fill a DataFrame for identifiers with no reads
        flattened_results = pd.DataFrame(columns=columns)
        # Add one empty row
        flattened_results.loc[0, :] = np.nan  # Fill with NaN first to avoid dtype changes
        # Assign identifier
        flattened_results.loc[0, 'identifier'] = identifier
        # Fill with 0 and NaN
        flattened_results[zero_fill_columns] = flattened_results[zero_fill_columns].fillna(0)
        flattened_results[nan_fill_columns] = flattened_results[nan_fill_columns].fillna(np.nan)
    else:
        # Convert data types to numeric
        data["total_dp"] = pd.to_numeric(data["total_dp"], errors='coerce')
        data["alt_dp"] = pd.to_numeric(data["alt_dp"], errors='coerce')
        
        # Check for values that cannot be converted
        if data["total_dp"].isnull().any() or data["alt_dp"].isnull().any():
            print(f"Warning: NaN values found in {file_path}. These rows will be ignored.")
            data = data.dropna(subset=["total_dp", "alt_dp"])
        
        # Calculate ref_dp
        data["ref_dp"] = data["total_dp"] - data["alt_dp"]
        # Calculate vaf
        data["vaf"] = data["alt_dp"] / data["total_dp"]
        
        # 1. Merge data and wg_data
        merged_data = pd.merge(data, wg_data, on='barcode', how='left')
        
        # Check for missing average_coverage
        missing_coverage = merged_data['average_coverage'].isnull().sum()
        if missing_coverage > 0:
            print(f"Warning: {missing_coverage} barcodes have no matching average_coverage in wg_data.")
        
        # 2. Normalize dp columns
        dp_columns = ['ref_dp', 'alt_dp', 'total_dp']
        for col in dp_columns:
            merged_data[f'norm_{col}'] = merged_data[col] / merged_data['average_coverage']
        
        # Split data into unmutant and mutant groups by alt_dp
        unmutant = merged_data[merged_data['alt_dp'] == 0]
        mutant = merged_data[merged_data['alt_dp'] > 0]
        
        # 4. Compute statistics for unmutant and mutant groups
        statistics_list = ['ref_dp', 'alt_dp', 'total_dp', 'vaf', 'norm_ref_dp', 'norm_alt_dp', 'norm_total_dp']
        unmutant_avg = unmutant[statistics_list].mean()
        mutant_avg = mutant[statistics_list].mean()
        unmutant_max = unmutant[statistics_list].max()
        mutant_max = mutant[statistics_list].max()
        unmutant_min = unmutant[statistics_list].min()
        mutant_min = mutant[statistics_list].min()
        
        # 5. Initialize results dictionary
        results = {}
        
        # 6. For each metric, compute mean difference, t-value, and p-value
        for column in statistics_list:
            # Get data for both groups
            unmutant_data = unmutant[column]
            mutant_data = mutant[column]
            
            # Perform independent t-test (Welch's t-test) -> Mann-Whitney U test
            # t_stat, p_val = stats.ttest_ind(mutant_data, unmutant_data, equal_var=False)
            # t_stat, p_val = stats.mannwhitneyu(mutant_data, unmutant_data, alternative='greater')
            # Check whether data is empty or all zeros
            if len(unmutant_data) == 0 or len(mutant_data) == 0:
                print(f"Warning: No data available for {column}. Skipping Mann-Whitney U test.")
                t_stat, p_val = np.nan, np.nan  # If no data, set t-value and p-value to NaN
            elif unmutant_data.sum() == 0 or mutant_data.sum() == 0:
                # If one group sums to 0, skip the test
                print(f"Warning: One group has all zero values for {column}. Skipping Mann-Whitney U test.")
                t_stat, p_val = np.nan, np.nan
            else:
                # Perform Mann-Whitney U test
                t_stat, p_val = stats.mannwhitneyu(mutant_data, unmutant_data, alternative='greater')
            
            # Calculate mean difference
            mean_diff = mutant_avg[column] - unmutant_avg[column]
            
            # Store all statistics in the results dictionary
            results[column] = {
                'unmutant_avg': unmutant_avg[column],
                'mutant_avg': mutant_avg[column],
                'unmutant_max': unmutant_max[column],
                'mutant_max': mutant_max[column],
                'unmutant_min': unmutant_min[column],
                'mutant_min': mutant_min[column],
                'mean_diff': mean_diff,
                't_value': t_stat,
                'p_value': p_val
            }
        
        # Convert results dictionary to DataFrame and transpose for display
        results_df = pd.DataFrame(results).T
        results_df['p_value'] = np.log10(results_df['p_value'] + 1e-300)
        
        # 7. Flatten results_df and create new column names 'columnname_rowname'
        flattened_results = results_df.unstack()
        flattened_results.index = [f"{col}_{row}" for row, col in flattened_results.index]
        flattened_results = flattened_results.to_frame().T  # Transpose so identifier is one row
        
        # 8. Add identifier column
        flattened_results['identifier'] = identifier
        
        # 9. Reorder columns so identifier comes first
        cols = ['identifier'] + [col for col in flattened_results.columns if col != 'identifier']
        flattened_results = flattened_results[cols]
        
        # 10. Count cell numbers
        flattened_results['mutant_cell_number'] = len(mutant)
        flattened_results['mutant_cell_fraction'] = len(mutant)/len(merged_data)
        flattened_results['mutant_and_unmutant_cell_number'] = len(data)
        flattened_results['mutant_and_unmutant_cell_fraction'] = len(data)/len(merged_data)
        # 11. pseudobulk
        flattened_results['ref_dp_in_pseudobulk'] = sum(merged_data['ref_dp'])
        flattened_results['alt_dp_in_pseudobulk'] = sum(merged_data['alt_dp'])
        flattened_results['total_dp_in_pseudobulk'] = sum(merged_data['total_dp'])
        flattened_results['vaf_in_pseudobulk'] = sum(merged_data['alt_dp'])/sum(merged_data['total_dp'])
        flattened_results['norm_ref_dp_in_pseudobulk'] = sum(merged_data['ref_dp'])/len(wg_data)
        flattened_results['norm_alt_dp_in_pseudobulk'] = sum(merged_data['alt_dp'])/len(wg_data)
        flattened_results['norm_total_dp_in_pseudobulk'] = sum(merged_data['total_dp'])/len(wg_data)
    
    return flattened_results

# Function for differential analysis (parallel)
def calculate_spot_statistics_parallel(output_folder, stats_file, identifier_list):
    """
    Process statistical analysis for all identifiers in parallel
    """
    # Read whole_genome_average_depth file
    wg_data = pd.read_csv(output_folder+"/whole_genome_average_depth.mut.spots.txt", sep="\t", header=0)
    
    # Process all identifiers in parallel with ProcessPoolExecutor
    with ProcessPoolExecutor() as executor:
        futures = [executor.submit(process_identifier, identifier, output_folder, wg_data) for identifier in identifier_list]
        
        all_results = []
        for future in as_completed(futures):
            all_results.append(future.result())
    
    # Merge all results
    final_results_df = pd.concat(all_results, ignore_index=True)
    final_results_df.to_csv(stats_file, sep="\t", index=False)
    print(f"Spot statistics written to: {stats_file}")
    
    return final_results_df


def main():
    file = open(args.mutation_list, "r")
    identifier_list = [i.strip() for i in file.readlines()]
    
    outpath = args.outdir
    check_dir(outpath)
    
    # Load barcode list
    if args.barcode_file != "" and os.path.exists(args.barcode_file):
        barcode_list = pd.read_csv(args.barcode_file, sep="\t", header=None, names=["barcode"])["barcode"].tolist()
    else:
        barcode_list = []
    
    bam_file = args.bam
    
    # Step 1: Calculate whole-genome average depth
    average_depth_file = os.path.join(outpath, "whole_genome_average_depth.mut.spots.txt")
    coverage = calculate_coverage_by_barcode(bam_file, args.barcode_file, args.threads)
    
    with open(average_depth_file, 'w') as f:
        f.write("barcode\taverage_coverage\n")
        for barcode, avg_cov in coverage.items():
            f.write(f"{barcode}\t{avg_cov:.2f}\n")
    
    # Step 2: Loop over each identifier
    for identifier in identifier_list:
        out_file = outpath + "/" + identifier + ".mut.spots.txt"
        with open(out_file, "w") as f:
            print("\n", "====", args.run_type, identifier, "===")
            out_info = get_identifier_info_from_bam(identifier, bam_file, barcode_list, args.CBtag, args.UBtag, args.run_type)
            chrom, pos, ref, alt = identifier.split("_")
            for onekey in out_info.keys():
                countsum = sum(out_info[onekey])
                f.write(f'{onekey}\t{countsum}\t{out_info[onekey]["ATCG".index(alt)]}\n')
    
    # Step 3: Compute mutant/unmutant metrics for the output folder
    stats_file = outpath + "/spot_statistics_summary.txt"
    stats_df_all_identifiers = calculate_spot_statistics_parallel(outpath, stats_file, identifier_list)
    
    # Step 4: Update the features file (if a features file argument is provided)
    if args.features:
        update_features_file_inplace(args.features, stats_df_all_identifiers)


## Argument parser
parser = argparse.ArgumentParser()
parser.add_argument("--bam", required=True, help="bam_file")
parser.add_argument("--CBtag", required=False, default="CB", help="CB tag")
parser.add_argument("--UBtag", required=False, default="UB", help="UB tag")
parser.add_argument("--barcode_file", required=False, default="", help="target barcodes file")
parser.add_argument("--run_type", required=False, choices=["UMI", "read"], default="UMI", help="which type you want to run")
parser.add_argument("--mutation_list", required=False, help="identifier list: chr_pos_ref_alt")
parser.add_argument("--features", required=False, help="existing features file to update")
parser.add_argument("--outdir", "-o", required=False, default="IGV_Plot", help="output dir")
parser.add_argument('--threads', type=int, default=1, help='Number of threads (default: all CPUs)')

args = parser.parse_args()

if __name__ == '__main__':
    main()
