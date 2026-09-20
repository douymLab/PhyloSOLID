#!/usr/bin/env python3

"""
germline_filter.py

Identify putative heterozygous germline variants
based on Methods Section 2 (MCF, J_r, U_r, S_r^FP, gap detection).

Inputs (from data_loader.load_all):
  - P : posterior matrix (cells x muts)
  - V : mutant allele frequency matrix (cells x muts)
  - C : coverage matrix (cells x muts)
  - A : mutant allele count matrix (cells x muts)

Main function:
  identify_germline_variants(P, V, C, A, mcf_cutoff=0.05)

Returns:
  A Python set of mutation IDs (columns in P/V/C/A) identified as germline.
"""
import numpy as np
import pandas as pd
from typing import Set, Dict, List, Tuple, Optional
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.colors import ListedColormap
import logging
from src.reproducibility import set_seed, deterministic_permutation

logger = logging.getLogger(__name__)

# -------------------------
# Helper functions
# -------------------------

def calculate_mutant_fraction(I_detected, df_reads_detected):
    # Count mutant_cell_number in each column (number of 1s)
    mutant_cell_number = (I_detected == 1).sum(axis=0)
    
    # Count non-NaN entries in each column
    coverage_cell_number = df_reads_detected.notna().sum(axis=0)
    
    # Compute mutant_cell_fraction (mutant_cell_number / coverage_cell_number)
    mutant_cell_fraction = mutant_cell_number / coverage_cell_number
    
    return mutant_cell_number, mutant_cell_fraction

def build_binary_I(P: pd.DataFrame, V: pd.DataFrame, C: pd.DataFrame, p_thresh: float = 0.5) -> pd.DataFrame:
    """Binarize matrix I according to formula (1)"""
    I = pd.DataFrame(np.nan, index=P.index, columns=P.columns)
    covered = (C >= 1)
    mask_mut = (V > 0) & (P > p_thresh) & covered
    I[mask_mut] = 1
    mask_ref = (V == 0) & (P <= p_thresh) & covered
    I[mask_ref] = 0
    # Set remaining NA entries with coverage>=1 to 0
    I[covered & I.isna()] = 0
    return I

def pairwise_counts(I: pd.DataFrame, j1: str, j2: str) -> Dict[str,int]:
    """Compute N11, N10, N01, N00 for two mutations (excluding NA)"""
    a = I[j1]; b = I[j2]
    valid = (~a.isna()) & (~b.isna())
    a = a[valid]; b = b[valid]
    N11 = int(((a==1)&(b==1)).sum())
    N10 = int(((a==1)&(b==0)).sum())
    N01 = int(((a==0)&(b==1)).sum())
    N00 = int(((a==0)&(b==0)).sum())
    return dict(N11=N11, N10=N10, N01=N01, N00=N00)

def pairwise_counts_for_two_columns(col1: pd.Series, col2: pd.Series) -> Dict[str, int]:
    """Compute N11, N10, N01, N00 for two columns (excluding NA)"""
    # Ensure the two columns share the same index
    if not col1.index.equals(col2.index):
        raise ValueError("Columns do not have the same index")
    # Filter out NA values
    valid = (~col1.isna()) & (~col2.isna())
    col1 = col1[valid]
    col2 = col2[valid]
    # Compute N11, N10, N01, N00
    N11 = int(((col1 == 1) & (col2 == 1)).sum())
    N10 = int(((col1 == 1) & (col2 == 0)).sum())
    N01 = int(((col1 == 0) & (col2 == 1)).sum())
    N00 = int(((col1 == 0) & (col2 == 0)).sum())
    return dict(N11=N11, N10=N10, N01=N01, N00=N00)

def jaccard_index(I, j1, j2):
    """Symmetric Jaccard"""
    counts = pairwise_counts(I, j1, j2)
    N11, N10, N01 = counts['N11'], counts['N10'], counts['N01']
    denominator = N11 + N10 + N01
    return N11 / denominator if denominator > 0 else 0.0

def f_fraction(I: pd.DataFrame, j1: str, j2: str) -> float:
    """f(j1,j2) = N11 / |S(j1)|"""
    counts = pairwise_counts(I,j1,j2)
    N11 = counts['N11']
    S1 = int((I[j1]==1).sum())
    return N11/S1 if S1>0 else 0.0

def are_mutations_correlated(I: pd.DataFrame, j1: str, j2: str) -> bool:
    """
    Determine whether two mutations are correlated, according to:
    - N11(j1,j2) ≥ 3 ∧ J(j1,j2) ≥ 0.2
    - or N11(j1,j2) ≥ 3 ∧ 0 < J(j1,j2) < 0.2 ∧ max(f(j1,j2), f(j2,j1)) ≥ 0.9
    """
    N11_threshold = 1
    J_val_threshold = 0.08
    f_fraction_threshold = 0.5
    counts = pairwise_counts(I, j1, j2)
    N11 = counts['N11']
    # First check whether N11 meets the minimum cell-count requirement
    if N11 < N11_threshold:
        return False
    # Compute Jaccard index
    J_val = jaccard_index(I, j1, j2)
    # First condition: Jaccard index ≥ 0.2
    if J_val >= J_val_threshold:
        return True
    # Second condition: 0 < Jaccard index < 0.2 and max(f(j1,j2), f(j2,j1)) ≥ 0.9
    if 0 < J_val < J_val_threshold:
        f_j1j2 = f_fraction(I, j1, j2)
        f_j2j1 = f_fraction(I, j2, j1)
        if max(f_j1j2, f_j2j1) >= f_fraction_threshold:
            return True
    return False

def build_J_r(I: pd.DataFrame, r: str) -> Set[str]:
    """
    Build the set J_r of mutations correlated with reference mutation r
    J_r = { j≠r | j is correlated with r }
    """
    J = set()
    for j in I.columns:
        if j == r:
            continue
        if are_mutations_correlated(I, r, j):
            J.add(j)
    return J

def infer_U_r(I: pd.DataFrame, r: str, J_r: Set[str]) -> Set[str]:
    """U_r = { cells with I[r]=1 or q_i ≥ q_threshold }"""
    if len(J_r) < 3:
        return set(I.index[I[r] == 1]), 'unknown'
    # Count mutants for each cell over J_r
    row_sums = I[list(J_r)].apply(lambda x: x[x == 1].count(), axis=1)
    # If the overall upper bound is too low (q_max <= 2) or low values dominate, do not expand
    value_counts = row_sums.value_counts().sort_index()
    low_counts_ratio = (value_counts.get(1,0) + value_counts.get(2,0)) / len(row_sums)
    if row_sums.max() <= 2 or low_counts_ratio > 0.5:
        # Do not expand; keep only cells with r=1
        return set(I.index[I[r] == 1]), low_counts_ratio
    # Otherwise use the cumulative-percentage threshold method
    value_counts_desc = row_sums.value_counts().sort_index(ascending=False)
    cumulative_percent = value_counts_desc.cumsum() / len(row_sums) * 100
    for i, (value, percent) in enumerate(cumulative_percent.items()):
        if percent > 10:
            if i == 0:
                q_threshold = value + 1
            else:
                q_threshold = list(cumulative_percent.index)[i-1]
            break
    else:
        q_threshold = cumulative_percent.index[-1]
    # Build mask
    mask = (I[r] == 1) | (row_sums >= q_threshold)
    return set(I.index[mask]), low_counts_ratio

def compute_S_r_FP(I: pd.DataFrame, r: str) -> float:
    """Compute S_r^FP = mean of S_j^FP(r) over J_r_plus"""
    J_r = build_J_r(I,r)
    J_plus = set(J_r)|{r}
    if not J_plus:
        return 0.0
    U_r, low_counts_ratio = infer_U_r(I,r,J_r)
    scores=[]
    for j in I.columns:
        if j==r: continue
        S_j = set(I.index[I[j]==1])
        if len(S_j)==0: 
            continue
        N_in = len([c for c in S_j if c in U_r])
        N_out = len(S_j)-N_in
        if j in J_plus:
            N_FP = N_out
        else:
            N_FP = N_in
        scores.append(N_FP/len(S_j))
    mean_score = float(np.mean(scores)) if scores else 0.0
    std_score = float(np.std(scores, ddof=0)) if scores else 0.0
    cv_score = std_score / mean_score if mean_score > 0 else 0.0
    return mean_score, std_score, cv_score, U_r, low_counts_ratio


def plot_heatmap_with_germline_mutations(I_raw, germline_mutations, pdf_file):
    """
    Plot a mutation-matrix heatmap with germline_mutations placed on the left,
    x-axis ticks in red, row/column mutation-count bar plots, and a legend below.

    Parameters
    ----------
    I_raw : pd.DataFrame
        Raw mutation matrix (cell x mutation), entries in {0,1,NA}
    germline_mutations : set
        Set of germline mutations, placed at the left of the heatmap
    pdf_file : str
        Path of the PDF file to save
    """
    
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, BoundaryNorm
    from matplotlib.patches import Patch
    
    # -------------------
    # Step 1: Move germline_mutations to the leftmost columns
    # -------------------
    germline_mutations_list = list(germline_mutations)
    # Keep only germline_mutations that are present as columns in I_raw
    germline_mutations_in_data = [mut for mut in germline_mutations_list if mut in I_raw.columns]
    
    # Place germline_mutations at the left of the dataframe
    I_sorted = I_raw[germline_mutations_in_data + [col for col in I_raw.columns if col not in germline_mutations_in_data]]
    
    # -------------------
    # Step 2: Convert matrix values (NA → 2, for a separate color)
    # -------------------
    I_numeric = I_sorted.fillna(np.nan).apply(pd.to_numeric, errors="coerce")
    I_plot = I_numeric.copy()
    I_plot = I_plot.where(~I_plot.isna(), 2)  # Fill NA with 2
    
    # -------------------
    # Step 3: Compute row/column mutation counts (for bar plots)
    # -------------------
    row_sums = I_numeric.sum(axis=1, skipna=True)  # Number of mutations per cell
    col_sums = I_numeric.sum(axis=0, skipna=True)  # Number of cells supporting each mutation
    
    # -------------------
    # Step 4: GridSpec layout
    #   - Left: row bar plot
    #   - Center: heatmap
    #   - Top: column bar plot
    # -------------------
    fig = plt.figure(figsize=(12, 10))
    gs = fig.add_gridspec(5, 6,
                          width_ratios=[0.3, 0.3, 3, 0.05, 0.05, 0.05],
                          height_ratios=[0.5, 0.05, 3, 0.3, 0.3],
                          wspace=0.05, hspace=0.05)
    
    ax_row_bar = fig.add_subplot(gs[2, 0])   # Left row bar plot
    ax_heatmap = fig.add_subplot(gs[2, 2])   # Center heatmap
    ax_col_bar = fig.add_subplot(gs[0, 2])   # Top column bar plot
    ax_dummy = fig.add_subplot(gs[0, 0]); ax_dummy.axis("off")  # Placeholder
    
    # -------------------
    # Step 5: Draw heatmap
    #   - Color map: 0=light blue, 1=dark red, NA=white
    # -------------------
    cmap = ListedColormap(["#D4E8F0", "#7D2224", "white"])
    bounds = [0, 0.5, 1.5, 2.5]
    norm = BoundaryNorm(bounds, cmap.N)
    
    im = ax_heatmap.imshow(I_plot, aspect="auto", cmap=cmap,
                           interpolation="nearest", norm=norm)
    ax_heatmap.set_xlim(-0.5, I_plot.shape[1]-0.5)
    ax_heatmap.set_ylim(I_plot.shape[0]-0.5, -0.5)
    ax_heatmap.set_yticks([])
    
    # Set x-axis mutation names
    ax_heatmap.set_xticks(range(len(I_plot.columns)))
    ax_heatmap.set_xticklabels(I_plot.columns, rotation=90, fontsize=6, ha='center')
    
    # -------------------
    # Step 6: Color x-axis mutation labels (germline_mutations in red)
    # -------------------
    for label in ax_heatmap.get_xticklabels():
        mut_name = label.get_text()
        if mut_name in germline_mutations_in_data:
            label.set_color('red')
        else:
            label.set_color('black')
    
    # -------------------
    # Step 7: Column bar plot (number of cells supporting each mutation)
    # -------------------
    ax_col_bar.bar(range(len(col_sums)), col_sums.values,
                   color="#7D2224", alpha=0.7, align="center")
    ax_col_bar.set_xlim(ax_heatmap.get_xlim())
    ax_col_bar.set_xticks([])
    ax_col_bar.tick_params(axis="y", labelsize=8)
    ax_col_bar.set_ylabel("Cell Number\nper Mutation", fontsize=10)
    
    # -------------------
    # Step 8: Row bar plot (number of mutations per cell)
    # -------------------
    ax_row_bar.barh(range(len(row_sums)), row_sums.values,
                    color="#7D2224", alpha=0.7, align="center")
    ax_row_bar.set_ylim(ax_heatmap.get_ylim())
    ax_row_bar.set_yticks([])
    ax_row_bar.set_xlabel("Mutation\nBurden\nper Cell", fontsize=10)
    ax_row_bar.invert_xaxis()
    
    # -------------------
    # Step 9: Add legend
    #   - Mutation values (0,1,NA)
    # -------------------
    # Mutation-value legend
    heatmap_handles = [Patch(facecolor=c, label=l) 
                       for c, l in zip(["#D4E8F0", "#7D2224", "white"],
                                       ["0 (No Mutation)", "1 (Mutation)", "NA (Missing)"])]
    fig.legend(handles=heatmap_handles, loc="lower center", ncol=3,
               bbox_to_anchor=(0.5, -0.03), frameon=False, fontsize=9,
               title="Mutation Values", title_fontsize=10)
    
    # -------------------
    # Step 10: Save figure
    # -------------------
    plt.suptitle("Heatmap of Mutations with Germline Mutations Highlighted", fontsize=14, y=0.95)
    plt.tight_layout()
    plt.subplots_adjust(bottom=0)
    plt.savefig(pdf_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {pdf_file}")

def update_germline_status(final_df, I, mcn_cutoff=5):
    # Select rows where 'germline_determined' is 'germline'
    germline_mutations = final_df[final_df['germline_determined'] == 'germline']
    
    for idx, row in germline_mutations.iterrows():
        # Get mutation id
        mutation_id = idx
        
        # Get the corresponding mutation column
        mutation_column = I[mutation_id]
        
        # Count entries equal to 1 in this mutation column
        mutant_cell_count = np.sum(mutation_column == 1)
        
        # If below mcn_cutoff, update to 'non-germline'
        if mutant_cell_count <= mcn_cutoff:
            final_df.loc[idx, 'germline_determined'] = 'non-germline'
    
    return final_df


def calculate_prob_threshold(probs_trimmed):
    mean = probs_trimmed.mean()
    std = probs_trimmed.std()
    
    # Compute the initial threshold
    threshold = mean + 2.6 * std
    
    # If the computed value is greater than 0.9, use mean + 2 * std
    if threshold > 0.9:
        print(f"Threshold is greater than 0.9; using {mean + 2 * std} (mean + 2 * std)")
        threshold = mean + 2 * std
    
    # If the computed value is greater than 0.9, use mean + 1.5 * std
    if threshold > 0.9:
        print(f"Threshold is greater than 0.9; using {mean + 1.5 * std} (mean + 1.5 * std)")
        threshold = mean + 1.5 * std
    
    # If the computed value is greater than 0.9, use mean + std
    if threshold > 0.9:
        print(f"Threshold is greater than 0.9; using {mean + std} (mean + std)")
        threshold = mean + std
    
    # If it is still greater than 0.9, set it to 0.9
    if threshold > 0.9:
        print("Threshold is still greater than 0.9; using 0.9")
        threshold = 0.9
    
    # If it is less than 0.1, set it to 0.1
    if threshold < 0.1:
        print("Threshold is less than 0.1; using 0.1")
        threshold = 0.1
    
    return threshold


def update_features_matrix(I, df_reads, df_features, mcf_cutoff):
    """
    Build features matrix and identify candidate founders
    
    Parameters:
    -----------
    I : pandas.DataFrame
        Binary matrix
    df_reads : pandas.DataFrame
        Reads data
    df_features : pandas.DataFrame
        Existing features matrix
    mcf_cutoff : float
        Minimum mutant cell fraction cutoff for candidates
        
    Returns:
    --------
    df_features_new : pandas.DataFrame
        Updated features matrix with mutant cell fraction data
    """
    
    # Step 1. Filter detected mutations and cells
    rows_with_ones = (I == 1).any(axis=1)
    cols_with_ones = (I == 1).any(axis=0)
    I_detected = I.loc[rows_with_ones, cols_with_ones]
    df_reads_detected = df_reads.loc[I_detected.index, I_detected.columns]
    empty_mutations = I.columns[~cols_with_ones]
    if not empty_mutations.empty:
        logger.warning(f"The following mutations are not detected in any cell: {', '.join(empty_mutations)}")
    
    # Step 2. Calculate mutant fractions
    mutant_cell_number_detected, mutant_cell_fraction_detected = calculate_mutant_fraction(I_detected, df_reads_detected)
    mutant_cell_number_input, mutant_cell_fraction_input = calculate_mutant_fraction(I, df_reads)
    
    # Step 3. Build new features matrix
    df_features_new = pd.concat([df_features, 
                             pd.DataFrame([mutant_cell_fraction_detected], index=['mutant_cell_fraction_detected']),
                             pd.DataFrame([mutant_cell_fraction_input], index=['mutant_cell_fraction_input'])])
    
    return df_features_new, empty_mutations


def add_mutation_proportions_to_features(df_features, df_cells):
    """
    Compute the 0/1/NA proportions of each mutation in the cell dataframe and add them to the feature dataframe
    
    Parameters:
    -----------
    df_features : DataFrame
        Mutation feature dataframe (df_features_new), used to store results
    df_cells : DataFrame
        Cell genotype dataframe (I_attached), used for computation
    """
    zero_props = []
    one_props = []
    na_props = []
    
    for mutation in df_features.columns:
        if mutation in df_cells.columns:
            # Get value counts for this mutation column (including NA)
            value_counts = df_cells[mutation].value_counts(dropna=False)
            total_cells = len(df_cells[mutation])
            
            # Compute proportions of each value
            zero_count = value_counts.get(0, 0)
            one_count = value_counts.get(1, 0)
            na_count = value_counts.get(np.nan, 0) if np.nan in value_counts.index else 0
            
            zero_prop = zero_count / total_cells if total_cells > 0 else 0
            one_prop = one_count / total_cells if total_cells > 0 else 0
            na_prop = na_count / total_cells if total_cells > 0 else 0
        else:
            # If the mutation is not in df_cells, set to 0 or NaN
            zero_prop = 0
            one_prop = 0
            na_prop = 1  # Or set to 1, indicating complete missingness
    
        zero_props.append(zero_prop)
        one_props.append(one_prop)
        na_props.append(na_prop)
    
    # Add results to df_features
    df_features.loc['zero_prop_detected'] = zero_props
    df_features.loc['one_prop_detected'] = one_props
    df_features.loc['na_prop_detected'] = na_props
    
    return df_features


def reorder_columns_by_mutant_stats(df_values, df_features_new, 
                                    min_cell_threshold=30, bin_size=5, 
                                    descending=True, return_stats=True):
    """
    Optimized column-reordering function: group by mutant cell number, then sort within groups by mutant cell fraction
    (fully deterministic sorting version)
    
    Parameters:
    -----------
    df_values : DataFrame
        Raw dataframe containing 0,1,NA (rows: cells, columns: mutations)
    df_features_new : DataFrame
        Dataframe containing mutation statistics
    min_cell_threshold : int
        Minimum cell-count threshold; mutations at or above this value form a high-priority group of their own
    bin_size : int
        Bin width for groups below the threshold
    descending : bool
        True: sort from large to small (high mutant cell number first)  
        False: sort from small to large
    return_stats : bool
        Whether to return sorting statistics
    
    Returns:
    --------
    df_reordered : DataFrame
        Dataframe with reordered columns
    sorting_stats : DataFrame (optional)
        Sorting statistics for columns
    """
    
    # 1. Get the intersection of columns from the two dataframes (sorted alphabetically for determinism)
    common_columns = sorted(list(set(df_values.columns) & set(df_features_new.columns)))
    # print(f"Original df_values column count: {len(df_values.columns)}")
    # print(f"Original df_features_new column count: {len(df_features_new.columns)}")
    # print(f"Number of shared columns: {len(common_columns)}")
    
    if len(common_columns) == 0:
        raise ValueError("The two dataframes have no columns in common!")
    
    # 2. Keep shared columns
    df_values_common = df_values[common_columns]
    
    # 3. Extract key statistics (shared columns only)
    mutant_cell_num = df_features_new[common_columns].loc['mutant_cellnum'].astype(int)
    mutant_cell_frac = df_features_new[common_columns].loc['mutant_cell_fraction'].astype(float)
    
    # 4. Create sorting-statistics DataFrame
    stats_df = pd.DataFrame({
        'column_name': mutant_cell_num.index,
        'mutant_cell_num': mutant_cell_num.values,
        'mutant_cell_frac': mutant_cell_frac.values
    })
    
    # 5. Define grouping logic
    def create_mutant_group(num):
        """Create mutant cell number group labels"""
        if num >= min_cell_threshold:
            return f'≥{min_cell_threshold}'
        else:
            lower = (num // bin_size) * bin_size
            upper = lower + bin_size - 1
            return f'{lower:02d}-{upper:02d}'
    
    stats_df['mutant_group'] = stats_df['mutant_cell_num'].apply(create_mutant_group)
    
    # 6. Define group sort order
    # Groups with high mutant cell number come first
    high_priority_groups = [f'≥{min_cell_threshold}']
    
    # Groups with low mutant cell number, from large to small
    low_priority_groups = []
    for i in range(min_cell_threshold - bin_size, -1, -bin_size):
        lower = i
        upper = i + bin_size - 1
        if lower >= 0:
            low_priority_groups.append(f'{lower:02d}-{upper:02d}')
    
    group_order = high_priority_groups + low_priority_groups
    
    # 7. Convert to an ordered categorical variable
    stats_df['mutant_group'] = pd.Categorical(
        stats_df['mutant_group'], 
        categories=group_order, 
        ordered=True
    )
    
    # 8. Fully deterministic sort: group, then mutant cell fraction, then column name
    if descending:
        # Large to small: high mutant number + high fraction first; column names alphabetical
        stats_df_sorted = stats_df.sort_values(
            ['mutant_group', 'mutant_cell_frac', 'column_name'], 
            ascending=[True, False, True]  # Groups by categorical order, fraction descending, column name ascending
        )
    else:
        # Small to large: low mutant number + low fraction first; column names alphabetical
        stats_df_sorted = stats_df.sort_values(
            ['mutant_group', 'mutant_cell_frac', 'column_name'], 
            ascending=[True, True, True]   # Groups by categorical order, fraction ascending, column name ascending
        )
    
    # 9. Get sorted column names
    sorted_columns = stats_df_sorted['column_name'].tolist()
    
    # 10. Reorder dataframe columns (shared columns only)
    df_reordered = df_values_common[sorted_columns]
    
    # 11. Reset index for inspection
    stats_df_sorted = stats_df_sorted.reset_index(drop=True)
    stats_df_sorted['final_order'] = stats_df_sorted.index + 1
    
    # print(f"Final reordered column count: {len(sorted_columns)}")
    # print(f"Group counts:")
    group_counts = stats_df_sorted['mutant_group'].value_counts().sort_index()
    # for group, count in group_counts.items():
    #     print(f"  {group}: {count} mutations")
    
    if return_stats:
        return df_reordered, stats_df_sorted
    else:
        return df_reordered

# # Usage example
# I_attached, sorting_stats_of_I_attached = reorder_columns_by_mutant_stats(
#     I_attached_split, 
#     df_features_new,
#     min_cell_threshold=30,  # ≥30 as the high-priority group
#     bin_size=5,             # Below 30, one group every 5
#     descending=True         # Sort from large to small
# )




# -------------------------
# Main function
# -------------------------

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from typing import Set, Optional, Tuple
from sklearn.linear_model import LogisticRegression
from matplotlib.backends.backend_pdf import PdfPages

def identify_germline_variants(
    P: pd.DataFrame, V: pd.DataFrame, C: pd.DataFrame, df_reads: pd.DataFrame, df_features_new: pd.DataFrame, 
    p_thresh: float = 0.5, mcf_cutoff: float = 0.05, mcn_cutoff: int = 5, 
    outputpath: Optional[str] = None,
    sampleid: Optional[str] = None,
    df_labeled: Optional[pd.DataFrame] = None
) -> Tuple[pd.DataFrame, Set[str]]:
    """
    Identify germline mutations via logistic regression on mean/std/cv,
    and visualize scatter plots.
    
    df_labeled: optional, labeled training dataset for logistic regression.
    """
    # Step 1. Build binary matrix
    I = build_binary_I(P, V, C, p_thresh)
    
    # Step 2. Candidate founder set    
    candidates = df_features_new.loc['mutant_cell_fraction_detected'][df_features_new.loc['mutant_cell_fraction_detected'] > mcf_cutoff].index.tolist()
    if not candidates:
        return pd.DataFrame(), set()
    
    # Step 3. Compute stats
    S_r_scores, S_r_std, S_r_cv, S_r_fn, S_r_lcr = {}, {}, {}, {}, {}
    for r in I.columns:
    # for r in candidates:
        mean_score, std_score, cv_score, U_r, low_counts_ratio = compute_S_r_FP(I, r)
        S_r_scores[r] = mean_score
        S_r_std[r] = std_score
        S_r_cv[r] = cv_score
        S_r_fn[r] = len(U_r) - len(I[I[r] == 1])
        S_r_lcr[r] = low_counts_ratio
    
    stats_df = pd.DataFrame({
        "FP_mean": pd.Series(S_r_scores),
        "FP_std": pd.Series(S_r_std),
        "FP_cv": pd.Series(S_r_cv),
        "FN_num": pd.Series(S_r_fn),
        "low_counts_ratio": pd.Series(S_r_lcr)
    })
    
    # Step 4. Split the dataframe
    # Filter mutations with FP_mean = FP_std = FP_cv = 0 into the non_germline dataframe
    stats_df_non_germline = stats_df[(stats_df['FP_mean'] == 0) &
                                     (stats_df['FP_std'] == 0) &
                                     (stats_df['FP_cv'] == 0)].copy()
    stats_df_candidates = stats_df.drop(stats_df_non_germline.index)  # Remaining mutations
    
    if len(stats_df_candidates)==0:
        return pd.DataFrame(), set()
    
    # Step 5. Label the non-germline dataframe
    stats_df_non_germline['germline_prob'] = np.nan  # Add empty column
    stats_df_non_germline['germline_pred'] = np.nan  # Add empty column
    stats_df_non_germline['germline_determined'] = 'non-germline'
    
    # Step 6. Apply logistic regression to stats_df_candidates
    pred_germline_mutations = set()
    prob_threshold = None
    
    if df_labeled is None:
        # Get the project root (assuming this file is under src/phylosolid/germline_filter/)
        current_file_dir = os.path.dirname(os.path.abspath(__file__))
        # Walk up to the project root: src/phylosolid/germline_filter -> src/phylosolid -> src -> root
        project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(current_file_dir))))
        
        # Build the path under the resource folder
        criteria_path = os.path.join(
            project_root,
            "resource",
            "stats_data_by_merged_3samples_withLabels.csv"
        )
        
        if os.path.exists(criteria_path):
            df_labeled = pd.read_csv(criteria_path)
            print(f"Loaded criteria file from: {criteria_path}")  # Optional debug info
        else:
            print(f"Warning: Criteria file not found at {criteria_path}")  # Optional warning
    
    if df_labeled is not None and not df_labeled.empty:
        if sampleid is not None and "sampleid" in df_labeled.columns:
            df_labeled = df_labeled[df_labeled["sampleid"] != sampleid]
        
        if not df_labeled.empty:
            # Prepare training data
            X_train = df_labeled[['FP_mean', 'FP_std', 'FP_cv']].values
            y_train = (df_labeled['label'] == 'germline').astype(int).values
            
            from sklearn.linear_model import LogisticRegression
            clf = LogisticRegression(class_weight="balanced", random_state=42)
            clf.fit(X_train, y_train)
            
            # Predict on the current dataset
            X_test = stats_df_candidates[['FP_mean', 'FP_std', 'FP_cv']].values
            probs = clf.predict_proba(X_test)[:, 1]
            stats_df_candidates['germline_prob'] = probs
            
            # Automatically compute cutoff
            probs_sorted = np.sort(probs)
            n_remove = int(len(probs_sorted) * 0.05)
            probs_trimmed = probs_sorted[n_remove:-n_remove] if n_remove > 0 else probs_sorted
            # prob_threshold = calculate_prob_threshold(probs_trimmed)
            prob_threshold = max(probs_trimmed.mean() + 1.5 * probs_trimmed.std(), 0.1)
            
            stats_df_candidates['germline_pred'] = (probs > prob_threshold).astype(int)
            pred_germline_mutations = set(stats_df_candidates[stats_df_candidates['germline_pred'] == 1].index.tolist())
            
            print("Identified prob_threshold: ", str(prob_threshold))
            
            # Visualization PDF
            if outputpath is not None:
                os.makedirs(outputpath, exist_ok=True)
                pdf_file = os.path.join(outputpath, "logreg_scatter_plots.pdf")
                from matplotlib.backends.backend_pdf import PdfPages
                import matplotlib.pyplot as plt
                
                with PdfPages(pdf_file) as pdf:  # pdf is only available inside this block
                    pairs = [('FP_mean', 'FP_std'), ('FP_mean', 'FP_cv'), ('FP_std', 'FP_cv')]
                    for xcol, ycol in pairs:
                        plt.figure(figsize=(7, 6))
                        scatter = plt.scatter(stats_df_candidates[xcol], stats_df_candidates[ycol],
                                              c=stats_df_candidates['germline_prob'],
                                              cmap='coolwarm', s=40, edgecolors='k')
                        plt.colorbar(scatter, label="Predicted germline probability")
                        plt.xlabel(xcol)
                        plt.ylabel(ycol)
                        plt.title(f"Logistic regression: {xcol} vs {ycol}\nCutoff={prob_threshold:.3f}")
                        plt.tight_layout()
                        pdf.savefig()  # Must be inside the with block
                        plt.close()
    
    # Add column 'germline_determined', based on whether the row name is in pred_germline_mutations
    stats_df_candidates['germline_determined'] = stats_df_candidates.index.to_series().apply(
        lambda x: 'germline' if x in pred_germline_mutations and candidates else 'non-germline'
    )
    
    # Step 7. Merge the two dataframes    
    merged_df = pd.concat([stats_df_candidates, stats_df_non_germline], axis=0)
    # Step 7. Update the 'germline_determined' column
    print("Updating germline status...")
    final_df = update_germline_status(merged_df, I, mcn_cutoff)
    print("Germline status updated.")
    
    # Step 8. Save output files
    if outputpath is not None:
        os.makedirs(outputpath, exist_ok=True)
        final_df.to_csv(os.path.join(outputpath, "S_r_FP_stats_df.csv"))
    
    final_germline_mutations = set(final_df[final_df['germline_determined'] == 'germline'].index)
    
    print("Identified germline mutations:\n", final_germline_mutations)
    return final_df, final_germline_mutations




# -------------------------
# 3.2 Coverage-based filtration
# -------------------------

def filter_scaffold_muts_by_na_proportion_germline(filtered_sites, df_reads, df_celltype, na_prop_thresh=0.9):
    """
    Identify high-confidence scaffold mutations (shared variants in relatively ubiquitously expressed genes)
    based on cross-cell-type coverage.
    Parameters
    ----------
    filtered_sites : list
        Mutations that have passed per-cell filters (MAF, coverage).
    df_reads : pd.DataFrame
        Rows = cells (first row can be 'bulk'), columns = mutations,
        values = 'mut_count/total_count' (string) or NaN
    df_celltype : pd.DataFrame
        DataFrame containing 'barcode' and 'cell_type' columns
    na_prop_thresh : float or None
        Uniform threshold for NA proportion across cell types.
        If None, use Q3 quantile from the data.
    Returns
    -------
    scaffold_mutations : list
        High-confidence shared mutations for building scaffold phylogeny
    NA_prop : pd.DataFrame
        NA proportion (mutation × cell type)
    """
    # --- 1. Drop the bulk row ---
    reads = df_reads.drop(index='bulk', errors='ignore')
    # --- 2. Get all cell types ---
    cell_types = df_celltype['cell_type'].unique()
    # --- 3. Build coverage matrix (covered = 1, uncovered or NA = 0) ---
    def has_coverage(val):
        if pd.isna(val):
            return 0
        try:
            _, total = val.split('/')
            return 1 if int(total) > 0 else 0
        except:
            return 0
    coverage_matrix = reads.applymap(has_coverage)
    # --- 4. Compute NA proportion ---
    df_NA_prop = pd.DataFrame(index=filtered_sites, columns=cell_types, dtype=float)
    for mut in filtered_sites:
        for t in cell_types:
            cells_in_type = df_celltype.loc[df_celltype['cell_type'] == t, 'barcode']
            valid_cells = [c for c in cells_in_type if c in coverage_matrix.index]
            if len(valid_cells) == 0:
                df_NA_prop.loc[mut, t] = 1.0
            else:
                cov_values = coverage_matrix.loc[valid_cells, mut]
                df_NA_prop.loc[mut, t] = 1.0 - cov_values.sum() / len(valid_cells)
    # --- 5. Compute cutoff ---
    if na_prop_thresh is not None:
        theta = pd.Series(na_prop_thresh, index=cell_types)
    else:
        theta = df_NA_prop.quantile(0.75, axis=0)
    # --- 6. Determine informative ---
    informative = df_NA_prop.lt(theta, axis=1)
    # --- 7. Select high-confidence scaffold mutations ---
    scaffold_mutations = []
    cell_prop = df_celltype['cell_type'].value_counts(normalize=True)
    dominant_ctypes = cell_prop[cell_prop > 0.9].index.tolist()
    has_dominant = len(dominant_ctypes) > 0
    for mut in informative.index:
        n_informative = informative.loc[mut].sum()
        n_celltypes = informative.shape[1]
        if n_celltypes == 1:
            if n_informative >= 1:
                scaffold_mutations.append(mut)
        elif has_dominant:
            if informative.loc[mut, dominant_ctypes].any():
                scaffold_mutations.append(mut)
        else:
            if n_informative >= 2:
                scaffold_mutations.append(mut)
    print("Step1 (coverage-based) mutations:", len(scaffold_mutations))
    return scaffold_mutations, df_NA_prop


def get_total_reads_withoutNAcells_germline(x):
    if pd.isna(x):
        return np.nan
    try:
        _, total = x.split('/')
        return int(total)
    except:
        return np.nan

def get_total_reads_withNAcells_germline(x):
    if pd.isna(x):
        return 0  # Treat NA as 0 reads
    try:
        _, total = x.split('/')
        return int(total)
    except:
        return 0  # In case of any parsing issues, treat as 0


def coverage_filters_germline(kept_mutations, df_reads, df_celltype, params, outputpath):
    """
    High-confidence scaffold mutation filtering (coverage-based)
    Step1: cross-cell-type coverage (NA proportion)
    Step2: CV filter (with median safety)
    The two steps are applied in parallel to the original input; the final result is their union
    
    Parameters
    ----------
    kept_mutations : list
        Candidate mutations
    df_reads : pd.DataFrame
        Rows = cells (first row is 'bulk'), columns = mutations, values = 'mut/total' or NaN
    df_celltype : pd.DataFrame
        DataFrame containing 'barcode' and 'cell_type' columns
    params : dict
        Contains 'na_prop_thresh_global' and 'cv_thresh'
    outputpath : str, optional
        Path to save the summary csv
    
    Returns
    -------
    final_scaffold_mutations : list
        High-confidence shared scaffold mutations
    summary_df : pd.DataFrame
        median / CV / mean / std / pass_filter for each mutation
    df_NA_prop : pd.DataFrame
        NA proportion from Step1
    """
    import numpy as np
    import pandas as pd
    
    if params is None:
        params = DEFAULT_PARAMS
    
    na_prop_thresh = params["na_prop_thresh_global"]
    cv_thresh = params["cv_thresh"]
    
    logger.info("Applying coverage-based filtration (Section 3.2)")
    
    # --- Step1: coverage-based filter ---
    step1_mutations, df_NA_prop = filter_scaffold_muts_by_na_proportion_germline(
        kept_mutations, df_reads, df_celltype, na_prop_thresh
    )
    logger.info("Section 3.2.1) Selection of ubiquitously expressed regions across cell (types)")
    print("=====> Step1 (coverage-based) mutations:", len(step1_mutations))
    
    # --- Step2: CV filter ---
    # Convert the reads to total read counts, treating NA as 0
    reads_matrix_withoutNAcells = df_reads.drop(index='bulk', errors='ignore').applymap(get_total_reads_withoutNAcells_germline)
    reads_matrix_withNAcells = df_reads.drop(index='bulk', errors='ignore').applymap(get_total_reads_withNAcells_germline)
    reads_matrix_withNAcells = reads_matrix_withNAcells.applymap(lambda v: 1 if (not pd.isna(v) and v >= 1) else (0 if not pd.isna(v) else np.nan))
    step2_mutations = []
    median_dict = {}
    cv_dict = {}
    mean_dict = {}
    std_dict = {}
    for mut in kept_mutations:
        values_for_median = reads_matrix_withoutNAcells[mut].dropna() if mut in reads_matrix_withoutNAcells else []
        values_for_cv = reads_matrix_withNAcells[mut].dropna()
        # median
        if len(values_for_median) == 0:
            median_dict[mut] = np.nan
            continue
        median_val = np.median(values_for_median)
        median_dict[mut] = median_val
        # cv
        mean_val = np.mean(values_for_cv)
        std_val = np.std(values_for_cv)
        cv = std_val / mean_val if mean_val > 0 else np.inf
        cv_dict[mut] = cv
        mean_dict[mut] = mean_val
        std_dict[mut] = std_val
        if cv <= cv_thresh:
            step2_mutations.append(mut)
    
    logger.info("Section 3.2.2) Selection of regions with relatively uniform read coverage")
    print("=====> Step2 (CV filter) mutations:", len(step2_mutations))
    
    # --- Union ---
    # final_scaffold_mutations = list(set(step1_mutations) | set(step2_mutations))
    final_scaffold_mutations = list(set(step2_mutations))
    print("=====> Final scaffold mutations (union):", len(final_scaffold_mutations))
    
    # --- Summary ---
    df_cv_stats = pd.DataFrame({
        "cov_median": pd.Series(median_dict),
        "cov_CV": pd.Series(cv_dict),
        "cov_mean": pd.Series(mean_dict),
        "cov_std": pd.Series(std_dict)
    })
    df_cv_stats["pass_CV"] = df_cv_stats.index.isin(step2_mutations)
    df_cv_stats["pass_NA"] = df_cv_stats.index.isin(step1_mutations)
    df_cv_stats["pass_cov"] = df_cv_stats.index.isin(final_scaffold_mutations)
    
    df_summary = pd.concat([df_cv_stats, df_NA_prop], axis=1)
    
    # --- Save ---
    if outputpath is not None:
        os.makedirs(outputpath, exist_ok=True)
        df_summary.to_csv(os.path.join(outputpath, "Summary_df_in_scaffold_filtration.csv"))
    return final_scaffold_mutations, df_summary



# -------------------------
# Compute pairwise Jaccard index then build a Leiden graph
# -------------------------

def compute_clone_and_pair_weights_germline(muts, corr_cache, n_shuffle=100):
    """
    Parameters
    ----------
    muts : list of str
        All mutation IDs
    corr_cache : dict
        {(mut1, mut2): True/False}  whether two mutations are correlated
    n_shuffle : int
        Number of shuffles per mutation
    
    Returns
    -------
    clone_weights : dict
        {tuple(mut_ids): weight}  global weight of each clone, including singleton clones
    pair_weights : dict
        {tuple(m1,m2): weight}  weight of each mutation pair (only clones of length ≥2)
    """
    clone_weights = defaultdict(float)  # Accumulate global clone weights
    
    for ref in muts:
        other_muts = [m for m in muts if m != ref]
        ref_clone_counter = defaultdict(int)  # Count each clone under the current reference
        
        for _ in range(n_shuffle):
            shuffled = list(deterministic_permutation(other_muts))
            remaining = [ref] + shuffled.copy()
            
            while remaining:
                curr_ref = remaining[0]
                current_clone = [curr_ref]
                next_remaining = []
                
                for m in remaining[1:]:
                    key1 = (curr_ref, m)
                    key2 = (m, curr_ref)
                    # Check corr_cache, avoid KeyError
                    is_corr = corr_cache.get(key1, corr_cache.get(key2, False))
                    if is_corr:
                        current_clone.append(m)
                    else:
                        next_remaining.append(m)
                
                # Clone count for the current shuffle
                ref_clone_counter[tuple(sorted(current_clone))] += 1
                remaining = next_remaining
        
        # Step: normalize by n_shuffle → clone proportion for this reference
        for clone, count in ref_clone_counter.items():
            clone_weights[clone] += count / n_shuffle  # Accumulate globally
    
    # Step: compute pair weights, considering only clones of length ≥2
    pair_weights = defaultdict(float)
    for clone, weight in clone_weights.items():
        if len(clone) > 1:
            for m1, m2 in itertools.combinations(sorted(clone), 2):
                pair_weights[(m1, m2)] += weight
    
    return clone_weights, pair_weights


import matplotlib.pyplot as plt
import networkx as nx
def plot_mutation_graph_germline(G_ig, mutation_group, pdf_file, figsize=(8,8), edge_scale=0.2, seed=42):
    """
    Visualize the mutation graph; node color indicates group, edge width indicates weight.
    
    Parameters
    ----------
    G_ig : igraph Graph
        Constructed igraph graph
    mutation_group : dict
        {mutation_id: group_id} group of each mutation
    figsize : tuple
        Figure size
    edge_scale : float
        Edge-weight scaling factor
    seed : int
        Layout random seed
    """
    # 1. Convert to a NetworkX graph
    G_nx = nx.Graph()
    for v in G_ig.vs:
        G_nx.add_node(v['name'])
    for e in G_ig.es:
        m1 = G_ig.vs[e.source]['name']
        m2 = G_ig.vs[e.target]['name']
        G_nx.add_edge(m1, m2, weight=float(e['weight']))
    
    # 2. Node colors
    groups = [mutation_group[n] for n in G_nx.nodes()]
    unique_groups = list(set(groups))
    color_map = plt.cm.get_cmap('tab20', len(unique_groups))
    node_colors = [color_map(g) for g in groups]
    
    # 3. Edge widths
    edges = G_nx.edges()
    edge_weights = [G_nx[u][v]['weight'] for u,v in edges]
    edge_widths = [w*edge_scale for w in edge_weights]
    
    # 4. Layout
    pos = nx.spring_layout(G_nx, seed=seed, k=0.5)  # k controls node spacing
    
    # 5. Draw
    plt.figure(figsize=figsize)
    nx.draw_networkx_nodes(G_nx, pos, node_color=node_colors, node_size=200)  # Slightly smaller nodes
    nx.draw_networkx_edges(G_nx, pos, width=edge_widths, alpha=0.7)
    nx.draw_networkx_labels(G_nx, pos, font_size=10, font_color='black')
    plt.title("Mutation Graph with Leiden Groups")
    plt.axis('off')
    plt.margins(x=0.2, y=0.2)        # Add margin around the plot
    plt.tight_layout(pad=2.0)        # Extra padding
    plt.savefig(pdf_file, dpi=300)
    plt.close()


from typing import List, Tuple, Dict
from collections import defaultdict
import itertools
import numpy as np
import pandas as pd

def cal_jaccard_index_by_pairs_for_graph_elements(I_S: pd.DataFrame):
    """
    Compute clone_weights and pair_jacidx, filtering out low-support singleton mutations.
    
    Args:
        I_S: binary matrix (cells x mutations), values 0/1/NA
        n_shuffle: number of random permutations per reference mutation
        seed: random seed
        min_frac: minimum mutant cell fraction for singleton mutations
        min_cells: minimum mutant cell number for singleton mutations
    
    Returns:
        clone_weights: dict mapping clone tuples to weight
        pair_jacidx: dict mapping mutation pairs to weight
    """
    muts = list(I_S.columns)
    n_mut = len(muts)
    
    # Step1: precompute pairwise jaccard index cache
    jacidx_cache = {}
    for u, v in itertools.combinations(muts, 2):
        jacidx = jaccard_index(I_S, u, v)
        jacidx_cache[(u, v)] = jacidx
        # jacidx_cache[(v, u)] = jacidx
    
    for m in muts:
        jacidx_cache[(m, m)] = 1.0
    
    # Step2: Convert results to Leiden graph input format
    pair_jacidx = defaultdict(float)
    for (var1, var2), weight in jacidx_cache.items():
        if var1 != var2:  # Skip the diagonal
            sorted_pair = tuple(sorted((var1, var2)))
            pair_jacidx[sorted_pair] = weight
                
    return pair_jacidx


import igraph as ig
import leidenalg
def leiden_mutation_groups_using_jaccard_index(pair_jacidx, pdf_file, resolution=1.0, seed=42):
    """
    Build a weighted co-occurrence graph from clone_weights and pair_jacidx, and partition mutation groups with Leiden.
    
    Parameters
    ----------
    clone_weights : dict
        {tuple(mutations): weight}  global weight of each clone, including singleton clones
    pair_jacidx : dict
        {tuple(m1,m2): weight}  weight of each mutation pair (only clones of length>=2)
    resolution : float
        Leiden algorithm resolution parameter
    seed : int
        Random seed
    
    Returns
    -------
    mutation_group : dict
        {mutation_id: group_id} group of each mutation
    partition : leidenalg VertexPartition
        Partition object returned by Leiden (can be used for visualization, etc.)
    G_ig : igraph Graph
        Constructed igraph graph
    """
    # 1. Collect all mutations (including isolated nodes)
    all_mutations = set()
    for clone in pair_jacidx.keys():
        all_mutations.update(clone)
    
    # 2. Build igraph graph
    G_ig = ig.Graph()
    G_ig.add_vertices(list(all_mutations))  # All mutations as nodes
    
    # Add edges (only clones of length>=2)
    for (m1, m2), w in pair_jacidx.items():
        G_ig.add_edge(m1, m2, weight=float(w))
    
    # 3. Run Leiden algorithm
    partition = leidenalg.find_partition(
        G_ig,
        leidenalg.RBConfigurationVertexPartition,
        weights='weight',
        resolution_parameter=resolution,
        seed=seed
    )
    
    # 4. Output mutation -> group dictionary
    mutation_group = {}
    for idx, community in enumerate(partition):
        for v in community:
            mutation_group[G_ig.vs[v]['name']] = idx
    
    # 5. Plot
    plot_mutation_graph_germline(G_ig, mutation_group, pdf_file)
    
    return mutation_group, partition, G_ig


def get_correlation_graph_elements_germline(I_S: pd.DataFrame, n_shuffle: int = 100, seed: int = 42, cutoff_mcf_for_graph: float = 0.05, cutoff_mcn_for_graph: int = 5) -> Tuple[Dict[Tuple[str], float], Dict[Tuple[str,str], float]]:
    """
    Compute clone_weights and pair_weights, filtering out low-support singleton mutations.
    
    Args:
        I_S: binary matrix (cells x mutations), values 0/1/NA
        n_shuffle: number of random permutations per reference mutation
        seed: random seed
        min_frac: minimum mutant cell fraction for singleton mutations
        min_cells: minimum mutant cell number for singleton mutations
    
    Returns:
        clone_weights: dict mapping clone tuples to weight
        pair_weights: dict mapping mutation pairs to weight
    """
    muts = list(I_S.columns)
    n_mut = len(muts)
    
    # Step 1: precompute pairwise correlation cache
    corr_cache = {}
    for u, v in itertools.combinations(muts, 2):
        corr = are_mutations_correlated(I_S, u, v)
        corr_cache[(u, v)] = corr
        corr_cache[(v, u)] = corr
    
    for m in muts:
        corr_cache[(m, m)] = True
    
    # Step 2: compute clone weights and pair weights
    clone_weights, pair_weights = compute_clone_and_pair_weights_germline(muts, corr_cache, n_shuffle=n_shuffle)
    
    # Step 3: Compute mutant fraction and mutant cell number for each mutation
    mutant_cell_fraction = {mut: I_S[mut].mean(skipna=True) for mut in muts}
    mutant_cell_number = {mut: I_S[mut].sum(skipna=True) for mut in muts}
    
    # Step 4: Drop clones corresponding to low-support singleton mutations
    count = 0
    for mut in muts:
        if clone_weights.get((mut,), 0) == n_mut:  # singleton clone
            frac = mutant_cell_fraction.get(mut, 0)
            num = mutant_cell_number.get(mut, 0)
            if frac <= cutoff_mcf_for_graph or num <= cutoff_mcn_for_graph:  # scDNA may need redefinition; this cutoff may not apply
                count +=1
                print(f"Filter out singleton low-support mutation: {mut}, frac={frac:.3f}, num={num}")
                clone_weights.pop((mut,), None)  # Drop this clone
    
    print(f"The number of filtered singleton, low-support mutations is: {count}")
    return clone_weights, pair_weights


import igraph as ig
import leidenalg
def leiden_mutation_groups_germline(clone_weights, pair_weights, pdf_file, resolution=1.0, seed=42):
    """
    Build a weighted co-occurrence graph from clone_weights and pair_weights, and partition mutation groups with Leiden.
    
    Parameters
    ----------
    clone_weights : dict
        {tuple(mutations): weight}  global weight of each clone, including singleton clones
    pair_weights : dict
        {tuple(m1,m2): weight}  weight of each mutation pair (only clones of length>=2)
    resolution : float
        Leiden algorithm resolution parameter
    seed : int
        Random seed
    
    Returns
    -------
    mutation_group : dict
        {mutation_id: group_id} group of each mutation
    partition : leidenalg VertexPartition
        Partition object returned by Leiden (can be used for visualization, etc.)
    G_ig : igraph Graph
        Constructed igraph graph
    """
    # 1. Collect all mutations (including isolated nodes)
    all_mutations = set()
    for clone in clone_weights.keys():
        all_mutations.update(clone)
    
    # 2. Build igraph graph
    G_ig = ig.Graph()
    G_ig.add_vertices(list(all_mutations))  # All mutations as nodes
    
    # Add edges (only clones of length>=2)
    for (m1, m2), w in pair_weights.items():
        G_ig.add_edge(m1, m2, weight=float(w))
    
    # 3. Run Leiden algorithm
    partition = leidenalg.find_partition(
        G_ig,
        leidenalg.RBConfigurationVertexPartition,
        weights='weight',
        resolution_parameter=resolution,
        seed=seed
    )
    
    # 4. Output mutation -> group dictionary
    mutation_group = {}
    for idx, community in enumerate(partition):
        for v in community:
            mutation_group[G_ig.vs[v]['name']] = idx
    
    # 5. Plot
    plot_mutation_graph_germline(G_ig, mutation_group, pdf_file)
    
    return mutation_group, partition, G_ig




##### Find the hub group among groups partitioned from each graph
def detect_hub_clusters_germline(G_ig, mutation_group):
    """
    Detect hub clusters based on weighted degree centrality
    """
    # 1. Build the cluster-level graph (as in the example)
    cluster_graph = build_cluster_graph_germline(G_ig, mutation_group)
    
    # 2. Compute the weighted degree of each cluster
    cluster_degrees = {}
    for cluster_id in set(mutation_group.values()):
        weighted_degree = 0
        for edge in cluster_graph.es:
            source_cluster = cluster_graph.vs[edge.source]['name']
            target_cluster = cluster_graph.vs[edge.target]['name']
            
            if source_cluster == cluster_id or target_cluster == cluster_id:
                weighted_degree += edge['weight']
        
        cluster_degrees[cluster_id] = weighted_degree
    
    # 3. Identify hub clusters (threshold is adjustable)
    hub_threshold = np.percentile(list(cluster_degrees.values()), 75)
    hub_clusters = [cluster_id for cluster_id, degree in cluster_degrees.items() 
                   if degree > hub_threshold]
    
    return hub_clusters, cluster_degrees


def build_cluster_graph_germline(G_ig, mutation_group):
    """
    Build a cluster-level weighted graph (as illustrated)
    """
    clusters = set(mutation_group.values())
    cluster_graph = ig.Graph()
    cluster_graph.add_vertices(list(clusters))
    
    # Compute inter-cluster connection weights
    inter_cluster_weights = {}
    for edge in G_ig.es:
        source_mut = G_ig.vs[edge.source]['name']
        target_mut = G_ig.vs[edge.target]['name']
        
        source_cluster = mutation_group[source_mut]
        target_cluster = mutation_group[target_mut]
        
        if source_cluster != target_cluster:
            pair = tuple(sorted([source_cluster, target_cluster]))
            inter_cluster_weights[pair] = inter_cluster_weights.get(pair, 0) + edge['weight']
    
    # Add edges
    for (cluster1, cluster2), weight in inter_cluster_weights.items():
        cluster_graph.add_edge(cluster1, cluster2, weight=weight)
    
    return cluster_graph




# def identify_germline_variants(
#     P: pd.DataFrame, V: pd.DataFrame, C: pd.DataFrame, df_reads: pd.DataFrame, df_features_new: pd.DataFrame, 
#     p_thresh: float = 0.5, mcf_cutoff: float = 0.05, mcn_cutoff: int = 5, 
#     outputpath: Optional[str] = None,
#     sampleid: Optional[str] = None,
#     df_labeled: Optional[pd.DataFrame] = None
# ):



# pair_jacidx = cal_jaccard_index_by_pairs_for_graph_elements(I_add_germline)
# mutation_group, partition, G_ig = leiden_mutation_groups_using_jaccard_index(pair_jacidx, outputpath + "/" + sampleid + ".graph_for_mut_grouping.pdf")

# group_mutations = list(mutation_group.keys())

# I_selected_and_sorted, mut_df_sorted, group_to_muts, final_order = sort_I_hierarchical_freeze_ones_fixed(I_add_germline, mutation_group)
# df_celltype_sub = df_celltype[df_celltype['barcode'].isin(I_selected_and_sorted.index)].copy()
# plot_heatmap_with_celltype_by_your_sorting(I_selected_and_sorted, df_celltype_sub, mutation_group, list(mut_df_sorted['mutation']), os.path.join(outputpath, sampleid+".heatmap_with_celltype_right_in_I_selected_and_sorted_after_graph_grouping.pdf"))
# logger.info(f"Total mutation groups: {len(mutation_group)}")

# hub_clusters, cluster_degrees = detect_hub_clusters(G_ig, mutation_group)




# logger.info("Performing mutation grouping using Leiden algorithm ...")
# clone_weights, pair_weights = get_correlation_graph_elements(I_add_germline, 100, 42)
# mutation_group, partition, G_ig = leiden_mutation_groups(clone_weights, pair_weights, outputpath + "/" + sampleid + ".graph_for_mut_grouping.pdf")
# group_mutations = list(mutation_group.keys())
# I_selected_and_sorted, mut_df_sorted, group_to_muts, final_order = sort_I_hierarchical_freeze_ones_fixed(I_add_germline, mutation_group)
# df_celltype_sub = df_celltype[df_celltype['barcode'].isin(I_selected_and_sorted.index)].copy()
# plot_heatmap_with_celltype_by_your_sorting(I_selected_and_sorted, df_celltype_sub, mutation_group, list(mut_df_sorted['mutation']), os.path.join(outputpath, sampleid+".heatmap_with_celltype_right_in_I_selected_and_sorted_after_graph_grouping.pdf"))
# logger.info(f"Total mutation groups: {len(mutation_group)}")

# hub_clusters, cluster_degrees = detect_hub_clusters(G_ig, mutation_group)




# -------------------------
# Demo
# -------------------------
if __name__ == "__main__":
    import pandas as pd
    
    # Construct a simple toy matrix (2 cells × 3 sites)
    P = pd.DataFrame([[0.9, 0.2, 0.8], [0.1, 0.8, 0.3]], 
                     index=["cell1", "cell2"], 
                     columns=["mut1", "mut2", "mut3"])
    M = pd.DataFrame([[1, 0, 1], [0, 1, 0]], index=P.index, columns=P.columns)
    C = pd.DataFrame([[10, 10, 10], [10, 10, 10]], index=P.index, columns=P.columns)
    A = (M * C).astype(int)
    
    # Construct a simple labeled dataset for training logistic regression
    df_labeled = pd.DataFrame({
        'FP_mean': [0.8, 0.1, 0.5],
        'FP_std': [0.05, 0.02, 0.1],
        'FP_cv': [0.0625, 0.2, 0.2],
        'label': ['germline', 'mosaic', 'germline']
    }, index=['mut1','mut2','mut3'])
    
    # Call the function
    stats_df, final_germline_mutations = identify_germline_variants(P, M, C, 
                                                       p_thresh=0.5, 
                                                       mcf_cutoff=0.05, 
                                                       outputpath=None, 
                                                       plot=False,
                                                       df_labeled=df_labeled)
    
    print("Stats dataframe:")
    print(stats_df)
    print("Putative germline variants:")
    print(final_germline_mutations)



