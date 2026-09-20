#!/usr/bin/env python3
"""
PhyloSOLID binary-matrix mode: reconstruct a phylogeny from a pre-called 0/1 matrix.

Classifier, germline filter, and coverage-based CV search are skipped because
the input is already a binary mutation matrix. Scaffold and full-tree intermediates
are kept; there is no CV* subdirectory.

Usage:
    python -m src.run_phylosilid_fullTree_binput -s SAMPLE_ID -i matrix.tsv -o /path/to/output

Output Structure:
    outputpath/SAMPLE_ID/
    ├── 01_scaffold_builder/
    ├── 02_mutation_integrator/
    ├── 03_final_results/
    │   ├── phylo/
    │   └── phylo_unpruned/
    └── logs/

Author: Qing
Date: 2026/03/13
Update: 2026/06/18
Update: 2026/09/19 - Numbered output dirs without CV search or unused filter steps
"""

import time
start_time = time.perf_counter()

import os
os.environ['PYTHONHASHSEED'] = '42'
os.environ['MPLBACKEND'] = 'Agg'
os.environ['QT_QPA_PLATFORM'] = 'offscreen'

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

import logging
import copy
import pandas as pd
import numpy as np
from copy import deepcopy
import json
import sys
import argparse

logger = logging.getLogger(__name__)

from src.reproducibility import set_seed
from src.scaffold_builder import *
from src.mutation_integrator import *
from src.full_tree_builder import build_fully_resolved_tree
from src.germline_filter import update_features_matrix, add_mutation_proportions_to_features
from src.utils import save_celltype_table
from src.phylo_export import export_final_phylo_results


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] [PID:%(process)d] %(message)s"
)
root_logger = logging.getLogger()

PRUNING_CONFIDENCE_THRESHOLD = 10.0

SETTING_PARAMS = {
    "models_path": "phylosolid/models/scdna",
    "p_thresh": 0.5,
    "mcf_cutoff": 0.05,
    "mcn_cutoff": 5,
    "pair_N11_min": 0,
    "jaccard_thresh": 0.2,
    "jaccard_low": 0.1,
    "fraction_parent_child_thresh": 0.9,
    "posterior_threshold": 0.5,
    "maf_max_threshold": 0.3,
    "maf_mean_threshold": 0.1,
    "na_prop_thresh_global": 0.95,
    "consensus_runs": 100,
    "consensus_clone_freq_thresh": 0.1,
    "resolution_of_graph": None,
    "min_resolution": 0.5,
    "max_resolution": 2.0,
    "general_weight_NA": 0.001,
    "fnfp_ratio": 0.1,
    "phi": 1.0,
    "pass_tree_cutoff": 0.9,
    "unpass_tree_cutoff": 0.1,
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
    "flipping_to_1_count_per_cells_cutoff": 2,
}


def setup_main_log_file(outputpath):
    """Set up main log file for the root logger."""
    log_dir = os.path.join(outputpath, "logs")
    os.makedirs(log_dir, exist_ok=True)
    main_log_file = os.path.join(log_dir, "run_results.log")

    for handler in root_logger.handlers[:]:
        if isinstance(handler, logging.FileHandler):
            root_logger.removeHandler(handler)

    file_handler = logging.FileHandler(main_log_file, mode='w')
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(logging.Formatter(
        "%(asctime)s [%(levelname)s] [PID:%(process)d] %(message)s"
    ))
    root_logger.addHandler(file_handler)

    for handler in root_logger.handlers[:]:
        if isinstance(handler, logging.StreamHandler) and getattr(handler, "stream", None) == sys.stdout:
            root_logger.removeHandler(handler)

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(logging.Formatter(
        "%(asctime)s [%(levelname)s] %(message)s"
    ))
    root_logger.addHandler(console_handler)

    root_logger.info(f"Main log file: {main_log_file}")
    return main_log_file


def _bin_to_reads(value):
    if pd.isna(value):
        return np.nan
    return "10/10" if int(value) == 1 else "0/10"


def parse_args():
    parser = argparse.ArgumentParser(
        description="PhyloSOLID binary-matrix mode: tree building from a pre-called 0/1 matrix"
    )
    parser.add_argument("-s", "--sampleid", default="", type=str, help="Sample ID")
    parser.add_argument("-o", "--outputpath", default="./output", type=str,
                        help="Output path; results are written to outputpath/sampleid/")
    parser.add_argument("-i", "--inputfile", default="", type=str,
                        help="Input binary matrix file (rows=cells, columns=mutations)")
    parser.add_argument("-c", "--celltype_file", default=None, type=str,
                        help="Optional cell-type file. If omitted, all cells use default_type.")
    parser.add_argument("--seed", default=42, type=int, help="Random seed for reproducibility")
    return parser.parse_args()


def main():
    args = parse_args()
    set_seed(args.seed)

    sampleid = args.sampleid
    inputfile = args.inputfile
    celltype_file = args.celltype_file
    outputpath = os.path.join(args.outputpath, sampleid)
    params = copy.deepcopy(SETTING_PARAMS)
    remove_artifact_mutations = "yes"
    is_filter_quality = "no"

    outputpath_scaffold = os.path.join(outputpath, "01_scaffold_builder")
    outputpath_full = os.path.join(outputpath, "02_mutation_integrator")
    for path in (outputpath, outputpath_scaffold, outputpath_full):
        os.makedirs(path, exist_ok=True)

    main_log_file = setup_main_log_file(outputpath)
    logger.info(f"sampleid: {sampleid}")
    logger.info(f"inputfile: {inputfile}")
    logger.info(f"outputpath: {outputpath}")
    logger.info(f"celltype_file: {celltype_file}")
    logger.info("is_predict_germ: no (binary-matrix mode)")
    logger.info("is_filter_quality: no (binary-matrix mode)")
    logger.info("")
    logger.info("Directory structure:")
    logger.info(f"  01_scaffold_builder: {outputpath_scaffold}")
    logger.info(f"  02_mutation_integrator: {outputpath_full}")
    logger.info(f"  03_final_results/phylo: {os.path.join(outputpath, '03_final_results', 'phylo')}")
    logger.info(f"  Main log file: {main_log_file}")

    # ------------------------------
    # Step 1: Load binary matrix
    # ------------------------------
    logger.info("===== Step1: Loading data ...")
    I_raw = pd.read_csv(inputfile, sep='\t', encoding='utf-8', index_col=0)
    I_raw.columns = I_raw.columns.str.replace(':', '_', regex=False)
    logger.info(f"Loaded data: {len(I_raw)} cells, {len(I_raw.columns)} mutations")

    I_filtered = I_raw[I_raw.eq(1).any(axis=1)]
    df_features = pd.DataFrame([
        (I_raw == 1).sum().astype(int),
        (I_raw == 1).sum() / len(I_raw)
    ], index=['mutant_cellnum', 'mutant_cell_fraction'], columns=I_raw.columns)

    bulk_row = pd.DataFrame(
        [[f"{int((I_raw[col] == 1).sum()) * 10}/{int(I_raw[col].count()) * 10}" for col in I_raw.columns]],
        index=['bulk'],
        columns=I_raw.columns
    )
    cell_reads = I_raw.apply(lambda col: col.map(_bin_to_reads))
    df_reads_raw = pd.concat([bulk_row, cell_reads])

    I = reorder_columns_by_mutant_stats(I_filtered, df_features)[0]
    all_mutations = list(I.columns)
    P = I.copy()
    V = I.copy()
    A = I.replace({0: 0, 1: 10}).fillna(0).astype(int)
    C = I.replace({0: 10, 1: 10}).fillna(0).astype(int)
    df_reads = df_reads_raw.loc[['bulk'] + list(I.index), I.columns]

    df_features_new, empty_mutations = update_features_matrix(I, df_reads, df_features, params["mcf_cutoff"])
    df_features_new = add_mutation_proportions_to_features(df_features_new, I)

    logger.info("===== Step2: Classifier ...")
    logger.info("Skipping classifier (binary-matrix mode: mutations are already called)")

    logger.info("===== Step3: Predict germline mutations ...")
    logger.info("Skipping germline filter (binary-matrix mode)")

    removed_germline_mutations = []
    removed_artifact_mutations = []
    somatic_mutations = list((reorder_columns_by_mutant_stats(I, df_features_new)[0]).columns)
    P_somatic = P[somatic_mutations].copy()
    V_somatic = V[somatic_mutations].copy()
    A_somatic = A[somatic_mutations].copy()
    C_somatic = C[somatic_mutations].copy()
    I_somatic = I[somatic_mutations].copy()
    df_reads_somatic = df_reads[somatic_mutations].copy()

    # ------------------------------
    # Step 4: Scaffold builder
    # ------------------------------
    logger.info("===== Step4: Construct scaffold tree ...")
    if celltype_file is None or celltype_file in ("None", "none", ""):
        barcodes = df_reads_somatic.index.tolist()
        df_celltype = pd.DataFrame({
            "barcode": barcodes,
            "cell_type": ["default_type"] * len(barcodes)
        })
    else:
        df_celltype = pd.read_csv(celltype_file, sep="\t")

    save_celltype_table(df_celltype, os.path.join(outputpath_scaffold, "df_celltype.txt"))
    logger.info(f"Celltype data loaded: {df_celltype.shape[0]} cells")
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
        is_filter_quality=is_filter_quality,
        outputpath_scaffold=outputpath_scaffold,
        sampleid=sampleid,
        immune_mutations=immune_mutations,
        df_celltype=df_celltype,
        logger_obj=logger
    )

    (T_scaffold, M_scaffold, df_flipping_spots, df_total_flipping_count,
     final_cleaned_I_selected_withNA3, final_cleaned_M_scaffold,
     backbone_mutations, mutation_group, spots_to_split, group_mutations,
     no_group_mutations, remained_mutations, conflict_mutations, root_mutations) = results_of_scaffold

    scaffold_mutations = list(M_scaffold.columns)
    logger.info(f"  Scaffold tree built: {len(scaffold_mutations)} mutations")
    print_tree_logger(T_scaffold, logger_obj=logger)

    # ------------------------------
    # Step 5: prepare attached mutations (no DP in binary-matrix mode)
    # ------------------------------
    logger.info("===== Step5: Dynamic programming pass tree ...")
    logger.info("Skipping DP pass-tree (binary-matrix mode)")
    attached_mutations = [
        mut for mut in all_mutations
        if mut not in scaffold_mutations
        and mut not in removed_germline_mutations
        and mut not in removed_artifact_mutations
    ]
    logger.info(f"  Pass tree mutations: {len(all_mutations)}")
    logger.info(f"  Attached mutations: {len(attached_mutations)}")

    I_attached_selected = I[scaffold_mutations + attached_mutations]
    I_attached_selected_sorted = I_attached_selected[
        I_attached_selected.apply(lambda col: (col == 1).sum(), axis=0).sort_values(ascending=False).index
    ]
    I_attached_sorted_non_empty = I_attached_selected_sorted[
        I_attached_selected_sorted.eq(1).any(axis=1)
    ]
    P_attached_sorted_non_empty = P.loc[
        I_attached_sorted_non_empty.index,
        I_attached_sorted_non_empty.columns
    ]
    I_attached_split, P_attached_split = split_spots_by_immune_mutations(
        spots_to_split,
        [mut for mut in immune_mutations if mut in I_attached_sorted_non_empty.columns],
        I_attached_sorted_non_empty,
        P_attached_sorted_non_empty
    )
    I_attached, sorting_stats_of_I_attached = reorder_columns_by_mutant_stats(
        I_attached_split,
        df_features_new,
        min_cell_threshold=30,
        bin_size=5,
        descending=True
    )
    P_attached = P_attached_split[I_attached.columns]
    all_conflict_mutations = list(conflict_mutations)

    # ------------------------------
    # Step 6-8: Fully resolved tree
    # ------------------------------
    logger.info("===== Step 6-8: Building fully resolved tree ...")
    (T_current, M_current, root_mutations, all_conflict_mutations,
     omega_before_qc, to_be_removed_cells, identified_doublet_cells,
     to_be_removed_mutations_by_fp_mutations_cross_all_cells,
     final_remained_mutations, final_conflict_mutations) = build_fully_resolved_tree(
        T_scaffold=T_scaffold,
        M_scaffold=M_scaffold,
        scaffold_mutations=scaffold_mutations,
        I_attached=I_attached,
        P_attached=P_attached,
        df_features_new=df_features_new,
        params=params,
        outputpath_full=outputpath_full,
        sampleid=sampleid,
        attached_mutations=attached_mutations,
        immune_mutations=immune_mutations,
        spots_to_split=spots_to_split,
        conflict_mutations=conflict_mutations,
        remove_artifact_mutations=remove_artifact_mutations,
        logger_obj=logger,
        cv_value=None,
        export_phylo=False
    )
    logger.info(f"  Fully resolved tree built: {M_current.shape[0]} cells, {M_current.shape[1]} mutations")

    # ------------------------------
    # Step 9: Post-processing & output
    # ------------------------------
    logger.info("=" * 80)
    logger.info("Step 9: Post-processing & output")
    logger.info("=" * 80)

    M_current_filtered = M_current.drop(columns=['ROOT'], errors='ignore')
    for mut_on_root in root_mutations:
        M_current_filtered.insert(0, mut_on_root, 1)
    mutations_on_T_current = (
        M_current_filtered.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()
    )
    T_full = copy.deepcopy(T_current)
    M_full = split_merged_columns(M_current_filtered, mutations_on_T_current)

    logger.info("Final full-resolved tree:")
    print_tree_logger(T_full, logger_obj=logger)
    logger.info(f"  Final tree cells: {M_full.shape[0]}")
    logger.info(f"  Final tree mutations: {M_full.shape[1]}")

    phylo_parent = os.path.join(outputpath, "03_final_results")
    phylo_export = export_final_phylo_results(
        phylo_parent, T_full, M_full, I_attached, params['fnfp_ratio'], logger_obj=logger,
    )
    phylo_dir = phylo_export["phylo_dir"]
    final_cleaned_M_full = phylo_export["M_cleaned"]
    omega_final = phylo_export["omega"]
    total_FP_flipping = phylo_export["total_delta_FP"]
    total_FN_flipping = phylo_export["total_delta_FN"]
    total_NAto0 = phylo_export["total_NA_to_0"]
    total_NAto1 = phylo_export["total_NA_to_1"]

    logger.info("")
    logger.info("  ┌─────────────────────────────────────────────────────────────────────┐")
    logger.info("  │              WEIGHTED DISCORDANCE INDEX (FINAL)                    │")
    logger.info("  ├─────────────────────────────────────────────────────────────────────┤")
    logger.info(f"  │  Weighted Discordance Index (Omega)        : {omega_final:>10.4f}      │")
    logger.info(f"  │    - delta_FP discordance                  : {total_FP_flipping:>10}          │")
    logger.info(f"  │    - delta_FN discordance                  : {total_FN_flipping:>10}          │")
    logger.info(f"  │    - NA->0 imputations                     : {total_NAto0:>10}       │")
    logger.info(f"  │    - NA->1 imputations                     : {total_NAto1:>10}       │")
    logger.info(f"  │    - FN/FP weight (lambda)                 : {params['fnfp_ratio']:>10.1f}      │")
    logger.info("  └─────────────────────────────────────────────────────────────────────┘")

    retention_rate = (
        (final_cleaned_M_full.shape[0] + final_cleaned_M_full.shape[1])
        / (I_somatic.shape[0] + I_somatic.shape[1])
        if (I_somatic.shape[0] + I_somatic.shape[1]) else 0.0
    )
    omega_reduced = omega_before_qc - omega_final
    pruning_ratio = omega_reduced / omega_final if omega_final else float('inf')

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
    logger.info(f"  │  Omega (pre-QC)                 : {omega_before_qc:.4f}      │")
    logger.info(f"  │  Omega (final)                  : {omega_final:.4f}      │")
    logger.info(f"  │  Omega (reduced)                : {omega_reduced:.4f}      │")
    logger.info(f"  │  Pruning ratio                  : {pruning_ratio:.4f}      │")
    logger.info(f"  │  Retention rate                 : {retention_rate:.4f}      │")
    logger.info("  │                                   (pruning ratio < 2.0 = confident)│")
    logger.info("  └─────────────────────────────────────────────────────────────────────┘")
    logger.info("")
    logger.info("=" * 80)
    logger.info("PhyloSOLID completed successfully!")
    logger.info("=" * 80)
    logger.info("Program finished in {:.4f} seconds".format(time.perf_counter() - start_time))


if __name__ == "__main__":
    main()
