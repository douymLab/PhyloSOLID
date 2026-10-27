"""
Export cleaned post-QC and pre-QC phylo result bundles.

Both final directories drop all-zero matrix rows/columns and collapse ghost
mutations from the matching tree. phylo_unpruned means before Step 8 QC,
not before Step 7 refinement or before export cleanup.
"""

import json
import logging
import os
from typing import Dict, Optional

import numpy as np
import pandas as pd

from src.mutation_integrator import (
    assign_clone_labels,
    calculate_flip_counts_per_site,
    get_mutation_clone_and_backbone_mut_as_keys_by_first_level_with_frequency,
    prune_mutations_from_tree,
)
from src.scaffold_builder import WriteTfile, find_flipping_spots, tree_to_dict, split_merged_columns

logger = logging.getLogger(__name__)


def prepare_phylo_matrix(M_current: pd.DataFrame, root_mutations) -> pd.DataFrame:
    """Expand merged mutations and restore ROOT-assigned mutations consistently."""
    matrix = M_current.drop(columns=['ROOT'], errors='ignore').copy()
    names = list(dict.fromkeys(m for column in matrix.columns for m in column.split('|')))
    matrix = split_merged_columns(matrix, names)
    for mutation in reversed(list(dict.fromkeys(root_mutations))):
        if mutation not in matrix.columns:
            matrix.insert(0, mutation, 1)
    return matrix


def _align_I_to_M(I_withNA3: pd.DataFrame, M: pd.DataFrame) -> pd.DataFrame:
    missing_rows = M.index.difference(I_withNA3.index)
    missing_cols = M.columns.difference(I_withNA3.columns)
    if len(missing_rows) or len(missing_cols):
        raise ValueError(f'Missing observations for phylo matrix: cells={list(missing_rows)}, mutations={list(missing_cols)}')
    return I_withNA3.loc[M.index, M.columns]


def _flipping_tables(I_bin: pd.DataFrame, M: pd.DataFrame):
    false_negative_flipping_spots = I_bin.apply(
        lambda col: find_flipping_spots(col, M[col.name], condition_in_bin=0, condition_phylogeny=1)
    )
    false_positive_flipping_spots = I_bin.apply(
        lambda col: find_flipping_spots(col, M[col.name], condition_in_bin=1, condition_phylogeny=0)
    )
    NAto1_flipping_spots = I_bin.apply(
        lambda col: find_flipping_spots(col, M[col.name], condition_in_bin=3, condition_phylogeny=1)
    )
    NAto0_flipping_spots = I_bin.apply(
        lambda col: find_flipping_spots(col, M[col.name], condition_in_bin=3, condition_phylogeny=0)
    )
    if false_negative_flipping_spots.empty:
        false_negative_flipping_spots = {col: [] for col in I_bin.columns}
    if false_positive_flipping_spots.empty:
        false_positive_flipping_spots = {col: [] for col in I_bin.columns}
    if NAto1_flipping_spots.empty:
        NAto1_flipping_spots = {col: [] for col in I_bin.columns}
    if NAto0_flipping_spots.empty:
        NAto0_flipping_spots = {col: [] for col in I_bin.columns}

    df_flipping_spots = pd.DataFrame({
        "Mutation": I_bin.columns,
        "delta_FN_spots": [", ".join(false_negative_flipping_spots.get(col, [])) for col in I_bin.columns],
        "delta_FP_spots": [", ".join(false_positive_flipping_spots.get(col, [])) for col in I_bin.columns],
        "NA_to_1_spots": [", ".join(NAto1_flipping_spots.get(col, [])) for col in I_bin.columns],
        "NA_to_0_spots": [", ".join(NAto0_flipping_spots.get(col, [])) for col in I_bin.columns],
    })
    counts = {
        "total_delta_FP": int(((I_bin == 1) & (M == 0)).sum().sum()),
        "total_delta_FN": int(((I_bin == 0) & (M == 1)).sum().sum()),
        "total_NA_to_0": int(((I_bin == 3) & (M == 0)).sum().sum()),
        "total_NA_to_1": int(((I_bin == 3) & (M == 1)).sum().sum()),
    }
    df_total = pd.DataFrame([counts])
    df_per_mut = (calculate_flip_counts_per_site(I_bin, M) if len(M.columns) else
                  pd.DataFrame(columns=['flipping_False_Negative', 'flipping_False_Positive',
                                        'flipping_NA_to_0', 'flipping_NA_to_1']))
    return df_flipping_spots, df_total, df_per_mut, counts


def _write_tree(T, outdir: str, stem: str) -> None:
    with open(os.path.join(outdir, f"{stem}.json"), "w") as handle:
        json.dump(tree_to_dict(T), handle, indent=4)
    T.save_to_file(os.path.join(outdir, f"{stem}.txt"))


def _write_clones(T, M: pd.DataFrame, I_for_freq: pd.DataFrame, outpath: str, logger_obj=None) -> None:
    if M.empty:
        pd.DataFrame(columns=["label", "color", "backbone_mutation"]).to_csv(outpath, index=False)
        return
    mutation_clones = get_mutation_clone_and_backbone_mut_as_keys_by_first_level_with_frequency(
        T, I_for_freq, logger_obj=logger_obj
    )
    present = set(M.columns)
    filtered = {}
    for key, muts in mutation_clones.items():
        kept = [m for m in muts if m in present]
        if not kept:
            continue
        filtered[key if key in present else kept[0]] = kept
    if not filtered:
        pd.DataFrame(columns=["label", "color", "backbone_mutation"]).to_csv(outpath, sep=",", index=False)
        return
    assign_clone_labels(M, filtered).to_csv(outpath, sep=",", index=False)


def _write_inferred_matrix(prefix, matrix, logger_obj=None):
    """Always save the inferred matrix, including empty/non-conflict-free cases."""
    output = matrix.astype(int)
    output.index.name = 'cellIDxmutID'
    output.to_csv(prefix + '.CFMatrix', sep='\t')
    if not matrix.empty:
        WriteTfile(prefix, matrix, matrix.index.tolist(), matrix.columns.tolist(),
                   judge='yes', logger_obj=logger_obj)


def write_phylo_bundle(
    outdir: str,
    T,
    M: pd.DataFrame,
    I_withNA3: pd.DataFrame,
    fnfp_ratio: float,
    *,
    cleaned: bool,
    suffix: str = "",
    logger_obj=None,
) -> Dict[str, object]:
    """Write one phylo bundle. cleaned=True uses final_cleaned_* names."""
    os.makedirs(outdir, exist_ok=True)
    I_aligned = _align_I_to_M(I_withNA3, M)
    M_aligned = M.loc[I_aligned.index, I_aligned.columns]

    if cleaned:
        _write_inferred_matrix(
            os.path.join(outdir, f"final_cleaned_M_full_basedPivots.filtered_sites_inferred{suffix}"),
            M_aligned, logger_obj=logger_obj,
        )
        I_aligned.to_csv(
            os.path.join(outdir, f"final_cleaned_I_full_withNA3_for_circosPlot{suffix}.txt"),
            sep="\t",
        )
        tree_stem = f"final_cleaned_tree_node{suffix}"
    else:
        _write_inferred_matrix(
            os.path.join(outdir, f"M_full_basedPivots.filtered_sites_inferred{suffix}"),
            M_aligned, logger_obj=logger_obj,
        )
        I_aligned.to_csv(os.path.join(outdir, f"I_full_withNA3{suffix}.txt"), sep="\t")
        tree_stem = f"tree_node{suffix}"

    _write_tree(T, outdir, tree_stem)

    df_spots, df_total, df_per_mut, counts = _flipping_tables(I_aligned, M_aligned)
    omega = counts["total_delta_FP"] + fnfp_ratio * counts["total_delta_FN"]
    df_total["weighted_discordance_index_Omega"] = [omega]
    df_spots.to_csv(os.path.join(outdir, f"df_flipping_spots{suffix}.txt"), sep="\t", index=False)
    df_total.to_csv(os.path.join(outdir, f"df_total_flipping_count{suffix}.txt"), sep="\t", index=False)
    df_per_mut.to_csv(os.path.join(outdir, f"df_flipping_count_for_each_mut{suffix}.txt"), sep="\t", index=True)
    _write_clones(
        T, M_aligned, I_aligned,
        os.path.join(outdir, f"df_barcode_clones_from_phylo_tree{suffix}.csv"),
        logger_obj=logger_obj,
    )
    return {
        "cells": int(M_aligned.shape[0]),
        "muts": int(M_aligned.shape[1]),
        "omega": float(omega),
        **counts,
        "outdir": outdir,
    }


def split_cleaned_and_unpruned(T_full, M_full: pd.DataFrame, I_attached: pd.DataFrame):
    """Drop all-zero mutations/cells and collapse those mutations out of the tree."""
    I_withNA3 = I_attached.replace({np.nan: 3}).astype(int)
    ghost_muts = [col for col in M_full.columns if not (M_full[col] != 0).any()]
    T_cleaned = prune_mutations_from_tree(T_full, ghost_muts)
    M_cleaned = M_full.loc[:, (M_full != 0).any(axis=0)]
    M_cleaned = M_cleaned.loc[(M_cleaned != 0).any(axis=1)].copy()
    return T_cleaned, M_cleaned, I_withNA3, ghost_muts


def export_final_phylo_results(
    outputpath_final: str,
    T_full,
    M_full: pd.DataFrame,
    I_attached: pd.DataFrame,
    fnfp_ratio: float,
    logger_obj=None,
    *,
    pre_qc_snapshot,
    output_stage: str = 'final',
) -> Dict[str, object]:
    """
    Write cleaned Step 8 output to phylo/ and cleaned Step 7 snapshot to phylo_unpruned/.

    Only called for the selected tree under final_results, not for every CV folder.
    """
    log = logger_obj if logger_obj is not None else logger
    if output_stage not in ('final', 'pre_fallback'):
        raise ValueError('Unknown phylo output stage: ' + str(output_stage))
    if pre_qc_snapshot is None:
        raise ValueError('A Step 7 pre-QC snapshot is required for phylo_unpruned')
    T_pre, M_pre, root_pre, _ = pre_qc_snapshot
    M_pre_full = prepare_phylo_matrix(M_pre, root_pre)
    T_pre_cleaned, M_pre_cleaned, I_pre, pre_ghost_muts = split_cleaned_and_unpruned(
        T_pre, M_pre_full, I_attached)
    T_cleaned, M_cleaned, I_withNA3, ghost_muts = split_cleaned_and_unpruned(T_full, M_full, I_attached)

    phylo_dir = os.path.join(outputpath_final, "phylo")
    unpruned_dir = os.path.join(outputpath_final, "phylo_unpruned")

    log.info("STEP 9.2: Output result files")
    log.info("-" * 80)
    if ghost_muts:
        log.info("  Collapsing %d all-zero mutation(s) from the cleaned tree: %s", len(ghost_muts), ghost_muts)
    else:
        log.info("  No all-zero mutations to collapse from the tree")

    cleaned = write_phylo_bundle(
        phylo_dir, T_cleaned, M_cleaned, I_withNA3, fnfp_ratio,
        cleaned=True, logger_obj=log,
    )
    unpruned = write_phylo_bundle(
        unpruned_dir, T_pre_cleaned, M_pre_cleaned, I_pre, fnfp_ratio,
        cleaned=False, logger_obj=log,
    )
    log.info("  Cleaned phylo:   %s (%s cells x %s muts)", phylo_dir, cleaned["cells"], cleaned["muts"])
    log.info("  Pre-QC phylo:    %s (%s cells x %s muts)", unpruned_dir, unpruned["cells"], unpruned["muts"])
    for directory, stage, source, stats, ghosts in [
        (phylo_dir, output_stage, 'step8' if output_stage == 'final' else 'step7', cleaned, ghost_muts),
        (unpruned_dir, 'pre_qc', 'step7', unpruned, pre_ghost_muts),
    ]:
        with open(os.path.join(directory, 'output_stage.json'), 'w') as handle:
            json.dump(dict(output_stage=stage, source_step=source,
                           all_zero_rows_columns_removed=True, ghost_mutations=ghosts,
                           cells=stats['cells'], mutations=stats['muts'], omega=stats['omega']),
                      handle, indent=2)
    with open(os.path.join(outputpath_final, 'README_phylo.txt'), 'w') as handle:
        handle.write('phylo/: Step 8 output; stage=' + output_stage + '\n'
                     'phylo_unpruned/: Step 7 snapshot before Step 8 QC.\n'
                     'Both bundles remove all-zero inferred matrix rows/columns and collapse corresponding tree nodes.\n'
                     'pre_fallback means QC produced no usable tree, so phylo also uses the Step 7 snapshot.\n')

    tree_muts = [name for name in T_cleaned.all_names_no_root() for name in name.split("|")]
    extra_on_tree = sorted(set(tree_muts) - set(M_cleaned.columns))
    missing_on_tree = sorted(set(M_cleaned.columns) - set(tree_muts))
    if extra_on_tree or missing_on_tree:
        log.warning("  Cleaned tree/matrix still differ: extra_on_tree=%s missing_on_tree=%s", extra_on_tree, missing_on_tree)
    else:
        log.info("  Cleaned tree nodes match cleaned matrix columns (%s mutations)", cleaned["muts"])

    return {
        "phylo_dir": phylo_dir,
        "phylo_unpruned_dir": unpruned_dir,
        "ghost_muts": ghost_muts,
        "T_cleaned": T_cleaned,
        "M_cleaned": M_cleaned,
        "I_withNA3": I_withNA3,
        **cleaned,
        "unpruned": unpruned,
        "output_stage": output_stage,
        "omega_pre_qc": unpruned["omega"],
    }


def export_cv_phylo_results(
    phylo_dir: str,
    T_full,
    M_full: pd.DataFrame,
    I_attached: pd.DataFrame,
    fnfp_ratio: float,
    suffix: str,
    logger_obj=None,
) -> Dict[str, object]:
    """Per-CV export: cleaned tree/matrix in the existing phylo_results folder (no extra unpruned dir)."""
    T_cleaned, M_cleaned, I_withNA3, ghost_muts = split_cleaned_and_unpruned(T_full, M_full, I_attached)
    write_phylo_bundle(
        phylo_dir, T_full, M_full, I_withNA3, fnfp_ratio,
        cleaned=False, suffix=suffix, logger_obj=logger_obj,
    )
    cleaned = write_phylo_bundle(
        phylo_dir, T_cleaned, M_cleaned, I_withNA3, fnfp_ratio,
        cleaned=True, suffix=suffix, logger_obj=logger_obj,
    )
    cleaned["ghost_muts"] = ghost_muts
    return cleaned
