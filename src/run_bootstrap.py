#!/usr/bin/env python3
"""
run_bootstrap_jaccard.py - Cell-resampling bootstrap analysis for mutation trees

This implementation reports one support measure:
percentage: Mean best-matching-clade Jaccard recovery of each mutation-defined branch.

Output files:
- bootstrap_support_percentage.csv

Output directory structure:
{outputpath}/
├── bootstrap_results/
│   ├── boot_0000_scaffold/          # scaffold files for each bootstrap
│   ├── boot_0000_full/              # full tree files for each bootstrap
│   │   └── phylo_results/           # complete results
│   ├── boot_0001_scaffold/
│   ├── boot_0001_full/
│   │   └── phylo_results/
│   └── ...
├── bootstrap_support_percentage.csv
├── bootstrap_support_percentage.json
└── logs/
    └── bootstrap_jaccard.log

Usage:
    python run_bootstrap_jaccard.py -s SAMPLE_ID -i /path/to/data -t tree.json -o /path/to/output -n 100
"""

import os
import sys
import json
import copy
import time
import argparse
import logging
import warnings
import numpy as np
import pandas as pd
from tqdm import tqdm
from typing import Dict, List, Set, Tuple, Optional, Any
from collections import defaultdict
from multiprocessing import Pool, cpu_count
from functools import partial
import matplotlib
matplotlib.use('Agg')

warnings.filterwarnings("ignore")

# Add src to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from src.data_loader import load_all
from src.scrna_classifier import real_time_classifier_predict
from src.germline_filter import identify_germline_variants
from src.germline_filter import *
from src.scaffold_builder import build_scaffold_tree
from src.scaffold_builder import *
from src.mutation_integrator import *
from src.full_tree_builder import build_fully_resolved_tree


# ============================================================
# Configuration
# ============================================================

DEFAULT_PARAMS = {
    "p_thresh": 0.5,
    "mcf_cutoff": 0.05,
    "mcn_cutoff": 5,
    "posterior_threshold": 0.5,
    "maf_max_threshold": 0.3,
    "maf_mean_threshold": 0.1,
    "na_prop_thresh_global": 0.95,
    "cv_rank_thresh": 0.5,
    "general_weight_NA": 0.001,
    "fnfp_ratio": 0.1,
    "phi": 1.0,
    "pass_tree_cutoff": 0.9,
    "unpass_tree_cutoff": 0.1,
    "fp_ratio_cutoff_across_tree": 0.2,
    "fn_ratio_cutoff_across_tree": 0.9,
    "fp_ratio_cutoff_within_subclone": 0.1,
    "fp_ratio_persite_cutoff": 0.1,
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


# ============================================================
# Bootstrap Support Calculator
# ============================================================

class BootstrapSupportCalculator:
    """
    Calculate bootstrap support using FULL PhyloSOLID pipeline.
    Support is reported as mean best-matching-clade Jaccard recovery.
    """
    
    def __init__(
        self,
        data: Dict[str, pd.DataFrame],
        df_features_new: pd.DataFrame,
        df_celltype: pd.DataFrame,
        optimal_tree: TreeNode,
        optimal_tree_mutations: List[str],
        sampleid: str,
        outputpath: str,
        params: Dict[str, Any] = None,
        n_bootstrap: int = 100,
        n_jobs: int = 1,
        random_seed: int = 42,
        logger: logging.Logger = None,
        cv_rank_thresh: float = 0.5,
        is_filter_quality: str = "yes",
        remove_artifact_mutations: str = "yes",
        immune_mutations: List[str] = None,
        min_clade_size: int = 2,
    ):
        self.data = data
        self.df_features_new = df_features_new
        self.df_celltype = df_celltype
        self.optimal_tree = optimal_tree
        self.optimal_tree_mutations = optimal_tree_mutations
        self.sampleid = sampleid
        self.outputpath = outputpath
        self.params = params if params else DEFAULT_PARAMS.copy()
        self.params['cv_rank_thresh'] = cv_rank_thresh
        self.n_bootstrap = n_bootstrap
        self.n_jobs = n_jobs
        self.random_seed = random_seed
        self.logger = logger if logger else logging.getLogger(__name__)
        self.is_filter_quality = is_filter_quality
        self.remove_artifact_mutations = remove_artifact_mutations
        self.immune_mutations = immune_mutations if immune_mutations else []
        
        # Mode parameters
        self.min_clade_size = min_clade_size
        
        # Create output directory
        self.bootstrap_dir = os.path.join(outputpath, "bootstrap_results")
        os.makedirs(self.bootstrap_dir, exist_ok=True)
        
        # Set seed
        set_seed(random_seed)
        
        # Build branch structure from optimal tree
        self._build_branch_structure(optimal_tree)
        
        # Store optimal tree mutations as set for efficient lookup
        self.optimal_tree_mutations_set = set(optimal_tree_mutations)
        
        self.logger.info("=" * 60)
        self.logger.info("Bootstrap Support Calculator Initialized")
        self.logger.info(f"  Optimal tree: {len(optimal_tree_mutations)} mutations")
        self.logger.info(f"  Branches: {len(self.branches)}")
        self.logger.info("  Support mode: best-matching-clade Jaccard percentage")
        self.logger.info(f"  Min clade size: {min_clade_size}")
        self.logger.info(f"  Output directory: {outputpath}")
        self.logger.info(f"  Bootstrap replicates: {n_bootstrap}")
        self.logger.info(f"  Parallel jobs: {n_jobs}")
        self.logger.info("=" * 60)
    
    def _build_branch_structure(self, tree: TreeNode):
        """Build branch structure from optimal tree."""
        self.branches = {}
        self.branch_hierarchy = {}
        self.branch_mutation_to_branch = {}
        
        for node in tree.traverse():
            if node.name == "ROOT":
                continue
            
            branch_mutations = set()
            for descendant in node.traverse():
                if descendant.name != "ROOT":
                    branch_mutations.update(descendant.name.split("|"))
            
            branch_id = node.name
            
            self.branches[branch_id] = {
                'node_name': node.name,
                'mutations': sorted(branch_mutations),
                'node': node,
                'level': len(node.path_to_root()) - 1,
                'node_mutations': node.name.split("|")
            }
            
            for mut in branch_mutations:
                self.branch_mutation_to_branch[mut] = branch_id
            
            if node.parent and node.parent.name != "ROOT":
                self.branch_hierarchy[branch_id] = node.parent.name
        
        self.logger.info(f"Built {len(self.branches)} branches")
        sizes = [len(b['mutations']) for b in self.branches.values()]
        self.logger.info(f"  Branch sizes: min={min(sizes)}, max={max(sizes)}, mean={np.mean(sizes):.1f}")
    
    def _get_mutations_from_tree(self, tree: TreeNode) -> Set[str]:
        """Extract all mutation names from a tree."""
        mutations = set()
        for node in tree.traverse():
            if node.name != "ROOT":
                mutations.update(node.name.split("|"))
        return mutations
    
    def _find_node_for_mutation(self, tree: TreeNode, mutation: str) -> Optional[TreeNode]:
        """Find the tree node containing a specific mutation."""
        for node in tree.traverse():
            if node.name != "ROOT" and mutation in node.name.split("|"):
                return node
        return None
    
    # ============================================================
    # Core: Check branch support
    # ============================================================
    
    def _check_branch_percentage(
        self,
        tree: TreeNode,
        branch_mutations: List[str],
        branch_id: str
    ) -> Dict[str, Any]:
        """
        Calculate best-matching-clade Jaccard recovery for one reference branch.
        
        Definitions:
        - n_total: Total number of mutations in this branch in the optimal tree
        - n_missing: Reference-clade mutations absent from the entire bootstrap
          tree. They are excluded from that replicate's branch comparison and
          therefore do not reduce branch support.
        - n_outgoing: Retained reference-clade mutations outside the selected
          bootstrap clade.
        - n_incoming: Mutations from other optimal-tree clades inside the
          selected bootstrap clade.
        - n_mismatch: n_outgoing + n_incoming.

        The percentage score is the Jaccard overlap between the retained
        reference clade and its best-matching bootstrap clade after both sets
        are restricted to mutations shared by the optimal and bootstrap trees:

            |reference ∩ bootstrap_clade| / |reference ∪ bootstrap_clade|

        Bootstrap-only mutations and optimal-tree mutations missing from the
        entire bootstrap tree are excluded before clade matching.
        """
        tree_mutations = self._get_mutations_from_tree(tree)
        branch_set = set(branch_mutations)
        shared_universe = self.optimal_tree_mutations_set & tree_mutations
        reference_shared = branch_set & shared_universe
        missing_mutations = branch_set - tree_mutations

        n_total = len(branch_set)
        n_existing = len(reference_shared)
        n_missing = len(missing_mutations)
        
        self.logger.debug(f"Branch {branch_id}: n_total={n_total}, "
                         f"shared={sorted(reference_shared)}, "
                         f"missing={sorted(missing_mutations)}")
        
        # At least two retained reference mutations are required to evaluate a
        # mutation-defined clade. Non-evaluable replicates are excluded from the
        # branch-specific mean rather than counted as unsupported.
        if n_existing < self.min_clade_size:
            return {
                'percentage_support': 0.0,
                'evaluable': False,
                'match_ratio': 0.0,
                'n_total': n_total,
                'n_existing': n_existing,
                'n_missing': n_missing,
                'n_matched': 0,
                'n_mismatch': 0,
                'n_outgoing': 0,
                'n_incoming': 0,
                'n_extra': 0,
                'n_union': 0,
                'clade_mutations': [],
                'matched_node': None,
                'lca_node': None,
                'global_extra_mutations_ignored': sorted(
                    tree_mutations - self.optimal_tree_mutations_set
                ),
            }

        # Enumerate every rooted clade in the bootstrap tree, including ROOT.
        # Each candidate is restricted to the shared mutation universe before
        # comparison. This permits both outgoing and incoming mismatches.
        candidates = []
        seen_clades = set()
        for traversal_order, node in enumerate(tree.traverse()):
            candidate_all = set()
            for descendant in node.traverse():
                if descendant.name != "ROOT":
                    candidate_all.update(descendant.name.split("|"))

            candidate_shared = candidate_all & shared_universe
            if not candidate_shared:
                continue

            frozen_candidate = frozenset(candidate_shared)
            if frozen_candidate in seen_clades:
                continue
            seen_clades.add(frozen_candidate)

            intersection = reference_shared & candidate_shared
            union = reference_shared | candidate_shared
            jaccard = len(intersection) / len(union)
            symmetric_difference = len(reference_shared ^ candidate_shared)
            size_difference = abs(len(reference_shared) - len(candidate_shared))
            depth = len(node.path_to_root()) - 1

            candidates.append({
                'node': node,
                'mutations': candidate_shared,
                'all_mutations': candidate_all,
                'jaccard': jaccard,
                'symmetric_difference': symmetric_difference,
                'size_difference': size_difference,
                'depth': depth,
                'traversal_order': traversal_order,
            })

        if not candidates:
            return {
                'percentage_support': 0.0,
                'evaluable': False,
                'match_ratio': 0.0,
                'n_total': n_total,
                'n_existing': n_existing,
                'n_missing': n_missing,
                'n_matched': 0,
                'n_mismatch': 0,
                'n_outgoing': 0,
                'n_incoming': 0,
                'n_extra': 0,
                'n_union': 0,
                'clade_mutations': [],
                'matched_node': None,
                'lca_node': None,
                'global_extra_mutations_ignored': sorted(
                    tree_mutations - self.optimal_tree_mutations_set
                ),
            }

        # Deterministic tie-breaking: maximize Jaccard overlap, then minimize
        # symmetric difference and size difference, prefer the deeper clade,
        # and finally retain traversal order.
        best = min(
            candidates,
            key=lambda item: (
                -item['jaccard'],
                item['symmetric_difference'],
                item['size_difference'],
                -item['depth'],
                item['traversal_order'],
            )
        )

        clade_mutations = best['mutations']
        matched_mutations = reference_shared & clade_mutations
        outgoing_mutations = reference_shared - clade_mutations
        incoming_mutations = clade_mutations - reference_shared
        union_mutations = reference_shared | clade_mutations

        n_matched = len(matched_mutations)
        n_outgoing = len(outgoing_mutations)
        n_incoming = len(incoming_mutations)
        n_mismatch = n_outgoing + n_incoming
        union_size = len(union_mutations)
        match_ratio = n_matched / union_size if union_size > 0 else 0.0

        self.logger.debug(
            f"Branch {branch_id}: matched_node={best['node'].name}, "
            f"matched={sorted(matched_mutations)}, "
            f"outgoing={sorted(outgoing_mutations)}, "
            f"incoming={sorted(incoming_mutations)}, "
            f"score={match_ratio:.4f}"
        )
        
        return {
            'percentage_support': match_ratio,
            'evaluable': True,
            'match_ratio': match_ratio,
            'n_total': n_total,
            'n_existing': n_existing,
            'n_missing': n_missing,
            'n_matched': n_matched,
            'n_mismatch': n_mismatch,
            'n_outgoing': n_outgoing,
            'n_incoming': n_incoming,
            'n_extra': n_incoming,
            'n_union': union_size,
            'matched_mutations': sorted(matched_mutations),
            'outgoing_mutations': sorted(outgoing_mutations),
            'incoming_mutations': sorted(incoming_mutations),
            'clade_mutations': sorted(clade_mutations),
            'matched_node': best['node'].name,
            # Retained for backward compatibility with earlier detailed output.
            'lca_node': best['node'].name,
            'global_extra_mutations_ignored': sorted(
                tree_mutations - self.optimal_tree_mutations_set
            ),
        }
    
    # ============================================================
    # Bootstrap Iteration
    # ============================================================
    
    def _run_single_bootstrap(self, bootstrap_idx: int) -> Dict[str, Any]:
        """Run a single bootstrap iteration."""
        try:
            # Set seed for this iteration
            iter_seed = self.random_seed + bootstrap_idx * 10000
            set_seed(iter_seed)
            
            # 1. Bootstrap resample cells
            P_raw = self.data["P"]
            V_raw = self.data["V"]
            C_raw = self.data["C"]
            A_raw = self.data["A"]
            df_reads = self.data["df_reads"]
            
            n_cells = P_raw.shape[0]
            sampled_indices = np.random.choice(
                P_raw.index,
                size=n_cells,
                replace=True
            )
            
            # Create bootstrap matrices
            P_boot = P_raw.loc[sampled_indices]
            V_boot = V_raw.loc[sampled_indices]
            C_boot = C_raw.loc[sampled_indices]
            A_boot = A_raw.loc[sampled_indices]
            df_reads_boot = df_reads.loc[['bulk'] + list(sampled_indices)]
            
            # ===== STEP 1: Rename all indices to be unique =====
            from collections import Counter
            
            all_indices = sampled_indices.tolist()
            index_counts = Counter(all_indices)
            duplicated_indices = [idx for idx, count in index_counts.items() if count > 1]
            
            new_index_mapping = {}
            new_indices = []
            
            for idx in all_indices:
                if idx in duplicated_indices:
                    if idx not in new_index_mapping:
                        new_index_mapping[idx] = []
                    seq_num = len(new_index_mapping[idx]) + 1
                    new_idx = f"{idx}_boot{bootstrap_idx}_{seq_num}"
                    new_index_mapping[idx].append(new_idx)
                    new_indices.append(new_idx)
                else:
                    new_index_mapping[idx] = [idx]
                    new_indices.append(idx)
            
            flat_mapping = {}
            for orig_idx, new_idx_list in new_index_mapping.items():
                for new_idx in new_idx_list:
                    flat_mapping[new_idx] = orig_idx
            
            self.logger.debug(f"Bootstrap {bootstrap_idx}: {len(all_indices)} cells, "
                             f"{len(duplicated_indices)} duplicated indices")
            
            # Apply new unique indices
            P_boot.index = new_indices
            V_boot.index = new_indices
            C_boot.index = new_indices
            A_boot.index = new_indices
            df_reads_boot.index = ['bulk'] + new_indices
            
            # ===== STEP 2: Build binary matrix =====
            from src.germline_filter import build_binary_I
            I_boot = build_binary_I(P_boot, V_boot, C_boot, 0.5)
            I_boot = I_boot[I_boot.eq(1).any(axis=1)]
            I_boot = I_boot.replace({np.nan: 0})
            
            if I_boot.shape[0] < 5 or I_boot.shape[1] < 3:
                self.logger.debug(f"Bootstrap {bootstrap_idx}: too few cells/mutations "
                                f"({I_boot.shape[0]} cells, {I_boot.shape[1]} mutations)")
                return {'success': False, 'error': 'Too few cells or mutations'}
            
            # ===== STEP 3: Filter all dataframes =====
            mutant_indices = I_boot.index.tolist()
            
            P_boot_filtered = P_boot.loc[mutant_indices].copy()
            V_boot_filtered = V_boot.loc[mutant_indices].copy()
            C_boot_filtered = C_boot.loc[mutant_indices].copy()
            A_boot_filtered = A_boot.loc[mutant_indices].copy()
            I_boot_filtered = I_boot.copy()
            df_reads_boot_filtered = df_reads_boot.loc[['bulk'] + mutant_indices].copy()
            
            mutant_flat_mapping = {new_idx: flat_mapping[new_idx] for new_idx in mutant_indices}
            
            # ===== STEP 4: Update features =====
            from src.germline_filter import update_features_matrix as ufm
            df_features_boot, _ = ufm(I_boot_filtered, df_reads_boot_filtered, self.df_features_new, 0.05)
            
            # ===== STEP 5: Create bootstrap df_celltype =====
            if 'barcode' in self.df_celltype.columns:
                celltype_dict = dict(zip(self.df_celltype['barcode'], self.df_celltype['cell_type']))
            else:
                celltype_dict = dict(zip(self.df_celltype.index, self.df_celltype['cell_type']))
            
            df_celltype_boot_list = []
            for new_idx, orig_idx in mutant_flat_mapping.items():
                if orig_idx in celltype_dict:
                    df_celltype_boot_list.append({
                        'barcode': new_idx,
                        'cell_type': celltype_dict[orig_idx]
                    })
                else:
                    df_celltype_boot_list.append({
                        'barcode': new_idx,
                        'cell_type': 'default_type'
                    })
            
            df_celltype_boot = pd.DataFrame(df_celltype_boot_list)
            
            # Store mapping
            self._boot_index_mapping = {
                'original_to_new': {orig: [idx for idx in new_list if idx in mutant_indices] 
                                   for orig, new_list in new_index_mapping.items()},
                'new_to_original': mutant_flat_mapping
            }
            
            # ===== STEP 6: Run full pipeline =====
            result = self._run_full_pipeline_on_bootstrap(
                P_boot_filtered,
                V_boot_filtered,
                A_boot_filtered,
                C_boot_filtered,
                I_boot_filtered,
                df_reads_boot_filtered,
                df_features_boot,
                df_celltype_boot,
                bootstrap_idx
            )
            
            if result is None:
                return {'success': False, 'error': 'Pipeline failed'}
            
            T_boot, M_boot = result
            
            # ===== STEP 7: Calculate best-clade Jaccard recovery for each branch =====
            branch_results = {}
            for branch_id, branch_info in self.branches.items():
                branch_mutations = branch_info['mutations']
                branch_results[branch_id] = self._check_branch_percentage(
                    T_boot, branch_mutations, branch_id
                )
            
            return {
                'success': True,
                'branches': branch_results,
                'n_cells': I_boot_filtered.shape[0],
                'n_mutations': I_boot_filtered.shape[1],
                'bootstrap_idx': bootstrap_idx,
                'mapping': mutant_flat_mapping
            }
            
        except Exception as e:
            self.logger.debug(f"Bootstrap {bootstrap_idx} failed: {e}")
            import traceback
            self.logger.debug(traceback.format_exc())
            return {'success': False, 'error': str(e)}
    
    def _run_full_pipeline_on_bootstrap(
        self,
        P: pd.DataFrame,
        V: pd.DataFrame,
        A: pd.DataFrame,
        C: pd.DataFrame,
        I: pd.DataFrame,
        df_reads: pd.DataFrame,
        df_features_boot: pd.DataFrame,
        df_celltype_boot: pd.DataFrame,
        bootstrap_idx: int
    ) -> Optional[Tuple[TreeNode, pd.DataFrame]]:
        """Run the full PhyloSOLID pipeline on bootstrap data."""
        try:
            bootstrap_label = f"boot_{bootstrap_idx:04d}"
            
            run_outputpath_scaffold = os.path.join(self.bootstrap_dir, f"{bootstrap_label}_scaffold")
            os.makedirs(run_outputpath_scaffold, exist_ok=True)
            
            # Build scaffold tree
            from src.scaffold_builder import build_scaffold_tree
            
            scaffold_result = build_scaffold_tree(
                P_somatic=P,
                V_somatic=V,
                A_somatic=A,
                C_somatic=C,
                I_somatic=I,
                df_reads_somatic=df_reads,
                df_features_new=df_features_boot,
                params=self.params,
                is_filter_quality=self.is_filter_quality,
                outputpath_scaffold=run_outputpath_scaffold,
                sampleid=f"{self.sampleid}_boot",
                immune_mutations=self.immune_mutations,
                df_celltype=df_celltype_boot,
                logger_obj=self.logger
            )
            
            if scaffold_result is None:
                self.logger.debug(f"Bootstrap {bootstrap_idx}: build_scaffold_tree returned None")
                return None
            
            # Unpack scaffold results
            (T_scaffold, M_scaffold, df_flipping_spots, df_total_flipping_count,
             final_cleaned_I_selected_withNA3, final_cleaned_M_scaffold,
             backbone_mutations, mutation_group, spots_to_split,
             group_mutations, no_group_mutations,
             remained_mutations, conflict_mutations, root_mutations) = scaffold_result
            
            scaffold_mutations = list(M_scaffold.columns)
            if 'ROOT' in scaffold_mutations:
                scaffold_mutations.remove('ROOT')
            
            self.logger.debug(f"Bootstrap {bootstrap_idx}: scaffold built with {len(scaffold_mutations)} mutations")
            
            # Prepare data for full-resolved tree building
            from src.mutation_integrator import split_spots_by_immune_mutations
            from src.germline_filter import reorder_columns_by_mutant_stats
            
            attached_mutations = [m for m in list(I.columns) if m not in scaffold_mutations]
            
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
            
            actual_immune_mutations = [i for i in self.immune_mutations if i in I_attached_sorted_non_empty.columns]
            I_attached_split, P_attached_split = split_spots_by_immune_mutations(
                spots_to_split,
                actual_immune_mutations,
                I_attached_sorted_non_empty,
                P_attached_sorted_non_empty
            )
            
            I_attached, sorting_stats_of_I_attached = reorder_columns_by_mutant_stats(
                I_attached_split,
                df_features_boot,
                min_cell_threshold=30,
                bin_size=5,
                descending=True
            )
            P_attached = P_attached_split[I_attached.columns]
            
            all_conflict_mutations = conflict_mutations.copy()
            
            # Build fully resolved tree
            from src.full_tree_builder import build_fully_resolved_tree
            
            run_outputpath_full = os.path.join(self.bootstrap_dir, f"{bootstrap_label}_full")
            os.makedirs(run_outputpath_full, exist_ok=True)
            
            full_result = build_fully_resolved_tree(
                T_scaffold=T_scaffold,
                M_scaffold=M_scaffold,
                scaffold_mutations=scaffold_mutations,
                I_attached=I_attached,
                P_attached=P_attached,
                df_features_new=df_features_boot,
                params=self.params,
                outputpath_full=run_outputpath_full,
                sampleid=f"{self.sampleid}_boot",
                attached_mutations=attached_mutations,
                immune_mutations=self.immune_mutations,
                spots_to_split=spots_to_split,
                conflict_mutations=conflict_mutations,
                remove_artifact_mutations=self.remove_artifact_mutations,
                logger_obj=self.logger,
                cv_value=self.params.get('cv_rank_thresh', 0.5),
                export_phylo=True
            )
            
            if full_result is None:
                self.logger.debug(f"Bootstrap {bootstrap_idx}: build_fully_resolved_tree returned None")
                return None
            
            T_current = full_result[0]
            M_current = full_result[1]
            
            self.logger.debug(f"Bootstrap {bootstrap_idx}: full tree built")
            
            return T_current, M_current
            
        except Exception as e:
            self.logger.debug(f"Bootstrap {bootstrap_idx} pipeline failed: {e}")
            import traceback
            self.logger.debug(traceback.format_exc())
            return None
    
    # ============================================================
    # Aggregate Results
    # ============================================================
    
    def compute_bootstrap_support(self) -> Dict[str, Dict[str, Any]]:
        """
        Compute mean best-matching-clade Jaccard bootstrap support values.
        
        Returns:
            {'percentage': {branch_id: support_percentage, ...}, 'details': {...}}
        """
        self.logger.info("=" * 80)
        self.logger.info("Starting Bootstrap Analysis")
        self.logger.info("=" * 80)
        self.logger.info(f"Replicates: {self.n_bootstrap}")
        self.logger.info(f"Parallel jobs: {self.n_jobs}")
        self.logger.info(f"Branches: {len(self.branches)}")
        self.logger.info("=" * 80)
        
        # Initialize accumulators
        percentage_sums = defaultdict(float)
        evaluable_counts = defaultdict(int)
        
        # Also track detailed statistics
        missing_counts = defaultdict(list)
        mismatch_counts = defaultdict(list)
        outgoing_counts = defaultdict(list)
        incoming_counts = defaultdict(list)
        extra_counts = defaultdict(list)
        match_ratios = defaultdict(list)
        
        # Run bootstrap iterations
        if self.n_jobs > 1:
            with Pool(processes=self.n_jobs) as pool:
                results = list(tqdm(
                    pool.imap(self._run_single_bootstrap, range(self.n_bootstrap)),
                    total=self.n_bootstrap,
                    desc="Bootstrap iterations"
                ))
        else:
            results = []
            for i in tqdm(range(self.n_bootstrap), desc="Bootstrap iterations"):
                results.append(self._run_single_bootstrap(i))
        
        # Aggregate results
        valid_iterations = 0
        for result in results:
            if not result.get('success', False):
                continue
            
            valid_iterations += 1
            branch_results = result.get('branches', {})
            
            for branch_id, data in branch_results.items():
                # Best-clade Jaccard recovery. Replicates with fewer than the
                # required number of retained reference mutations are not
                # informative for this branch and are excluded from its mean.
                if data.get('evaluable', False):
                    percentage_sums[branch_id] += data.get('percentage_support', 0.0)
                    evaluable_counts[branch_id] += 1
                
                # Detailed statistics
                missing_counts[branch_id].append(data.get('n_missing', 0))
                if data.get('evaluable', False):
                    mismatch_counts[branch_id].append(data.get('n_mismatch', 0))
                    outgoing_counts[branch_id].append(data.get('n_outgoing', 0))
                    incoming_counts[branch_id].append(data.get('n_incoming', 0))
                    extra_counts[branch_id].append(data.get('n_extra', 0))
                    match_ratios[branch_id].append(data.get('match_ratio', 0.0))
        
        # Calculate final support values
        if valid_iterations == 0:
            raise RuntimeError(
                "All bootstrap replicates failed; branch support cannot be calculated."
            )
        n_effective = valid_iterations
        self.valid_iterations = valid_iterations
        self.logger.info(f"Valid iterations: {valid_iterations}")
        
        support_values = {
            'percentage': {},
            # Also store detailed stats for each branch
            'details': {}
        }
        
        for branch_id in self.branches.keys():
            # Percentage: mean Jaccard recovery across branch-evaluable replicates.
            branch_n_evaluable = evaluable_counts.get(branch_id, 0)
            if branch_n_evaluable > 0:
                support = (percentage_sums.get(branch_id, 0.0) / branch_n_evaluable) * 100
            else:
                support = np.nan
            support_values['percentage'][branch_id] = support
            
            # Detailed statistics (averages)
            support_values['details'][branch_id] = {
                'avg_missing': np.mean(missing_counts[branch_id]) if missing_counts[branch_id] else 0,
                'avg_mismatch': np.mean(mismatch_counts[branch_id]) if mismatch_counts[branch_id] else 0,
                'avg_outgoing': np.mean(outgoing_counts[branch_id]) if outgoing_counts[branch_id] else 0,
                'avg_incoming': np.mean(incoming_counts[branch_id]) if incoming_counts[branch_id] else 0,
                'avg_extra': np.mean(extra_counts[branch_id]) if extra_counts[branch_id] else 0,
                'avg_match_ratio': np.mean(match_ratios[branch_id]) * 100 if match_ratios[branch_id] else 0,
                'n_evaluable_bootstrap': branch_n_evaluable,
                'evaluation_rate': (branch_n_evaluable / n_effective) * 100,
            }
        
        self.logger.info(f"Bootstrap analysis completed")
        return support_values
    
    # ============================================================
    # Save Results
    # ============================================================
    
    def save_results(self, support_values: Dict[str, Dict[str, Any]]):
        """Save best-matching-clade Jaccard bootstrap support results."""
        
        # Retain the historical 'percentage' filename for downstream compatibility.
        for mode_name in ['percentage']:
            mode_data = support_values.get(mode_name, {})
            rows = []
            for branch_id, branch_info in self.branches.items():
                support = mode_data.get(branch_id, 0.0)
                
                # Get detailed stats for this branch
                details = support_values.get('details', {}).get(branch_id, {})
                
                row = {
                    'branch_id': branch_id,
                    'branch_node': branch_info['node_name'],
                    'support_percentage': support,
                    'mutations_in_branch': '|'.join(branch_info['mutations']),
                    'n_mutations': len(branch_info['mutations']),
                    'level': branch_info['level'],
                    'avg_missing': details.get('avg_missing', 0),
                    'avg_mismatch': details.get('avg_mismatch', 0),
                    'avg_outgoing': details.get('avg_outgoing', 0),
                    'avg_incoming': details.get('avg_incoming', 0),
                    'avg_extra': details.get('avg_extra', 0),
                    'avg_match_ratio': details.get('avg_match_ratio', 0),
                    'n_valid_bootstrap': getattr(self, 'valid_iterations', 0),
                    'n_evaluable_bootstrap': details.get('n_evaluable_bootstrap', 0),
                    'evaluation_rate': details.get('evaluation_rate', 0),
                }
                rows.append(row)
            
            df_support = pd.DataFrame(rows)
            df_support = df_support.sort_values('support_percentage', ascending=False)
            
            # Save CSV
            csv_file = os.path.join(self.outputpath, f"bootstrap_support_{mode_name}.csv")
            df_support.to_csv(csv_file, index=False)
            self.logger.info(f"Saved: {csv_file}")
            
            # Print summary for this mode
            support_list = df_support['support_percentage'].dropna().tolist()
            self.logger.info("-" * 60)
            self.logger.info(f"Bootstrap Support Summary ({mode_name} mode)")
            self.logger.info("-" * 60)
            self.logger.info(f"  Total branches: {len(support_list)}")
            if support_list:
                self.logger.info(f"  Mean support: {np.mean(support_list):.1f}%")
                self.logger.info(f"  Median support: {np.median(support_list):.1f}%")
                self.logger.info(f"  Std support: {np.std(support_list):.1f}%")
                self.logger.info(f"  Min support: {min(support_list):.1f}%")
                self.logger.info(f"  Max support: {max(support_list):.1f}%")
            else:
                self.logger.info("  No branches had enough retained mutations to be evaluated")
            self.logger.info("-" * 60)
        
        # Save JSON with support values and diagnostic summaries
        json_file = os.path.join(self.outputpath, "bootstrap_support_percentage.json")
        json_safe = {}
        for mode_name in ['percentage']:
            mode_data = support_values.get(mode_name, {})
            json_safe[mode_name] = {
                k: (float(v) if pd.notna(v) else None)
                for k, v in mode_data.items()
            }
        
        # Also save details
        json_safe['details'] = {}
        for branch_id, details in support_values.get('details', {}).items():
            json_safe['details'][branch_id] = {k: float(v) for k, v in details.items()}
        
        with open(json_file, 'w') as f:
            json.dump(json_safe, f, indent=2)
        self.logger.info(f"Saved combined JSON: {json_file}")


# ============================================================
# Helper Functions
# ============================================================

def _dict_to_tree(tree_dict: Dict) -> TreeNode:
    """Convert dictionary to TreeNode."""
    def _build_node(node_dict):
        node = TreeNode(node_dict['name'])
        if 'children' in node_dict:
            for child_dict in node_dict['children']:
                child = _build_node(child_dict)
                node.add_child(child)
        return node
    return _build_node(tree_dict)


def update_features_matrix(I, df_reads, df_features, mcf_cutoff):
    """Helper function to update features matrix."""
    from src.germline_filter import update_features_matrix as ufm
    return ufm(I, df_reads, df_features, mcf_cutoff)


# ============================================================
# Main Function
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="Cell-resampling bootstrap with best-matching-clade Jaccard recovery"
    )
    parser.add_argument("-s", "--sampleid", required=True, type=str,
                       help="Sample ID")
    parser.add_argument("-i", "--inputpath", required=True, type=str,
                       help="Path to input data directory")
    parser.add_argument("-t", "--tree_json", required=True, type=str,
                       help="Path to optimal tree JSON file")
    parser.add_argument("-o", "--outputpath", required=True, type=str,
                       help="Output directory path")
    parser.add_argument("-n", "--n_bootstrap", default=100, type=int,
                       help="Number of bootstrap replicates (default: 100)")
    parser.add_argument("-j", "--n_jobs", default=1, type=int,
                       help="Number of parallel jobs (default: 1)")
    parser.add_argument("--seed", default=42, type=int,
                       help="Random seed (default: 42)")
    parser.add_argument("--cv_rank_thresh", default=0.5, type=float,
                       help="CV rank threshold (default: 0.5)")
    parser.add_argument("--celltype_file", default=None, type=str,
                       help="Cell type file path (optional)")
    parser.add_argument("--is_filter_quality", default="yes", 
                       choices=["yes", "no"], type=str)
    parser.add_argument("--remove_artifact", default="yes",
                       choices=["yes", "no"], type=str)
    parser.add_argument("--min_clade_size", default=2, type=int,
                       help="Minimum mutations to form a clade (default: 2)")
    
    args = parser.parse_args()
    
    # Setup logging
    log_dir = os.path.join(args.outputpath, "logs")
    os.makedirs(log_dir, exist_ok=True)
    log_file = os.path.join(log_dir, "bootstrap_jaccard.log")
    
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_file, mode='w'),
            logging.StreamHandler(sys.stdout)
        ]
    )
    logger = logging.getLogger(__name__)
    
    logger.info("=" * 80)
    logger.info("BOOTSTRAP SUPPORT ANALYSIS - BEST-CLADE JACCARD RECOVERY")
    logger.info("=" * 80)
    logger.info(f"Sample: {args.sampleid}")
    logger.info(f"Replicates: {args.n_bootstrap}")
    logger.info(f"CV threshold: {args.cv_rank_thresh}")
    logger.info(f"Min clade size: {args.min_clade_size}")
    logger.info("=" * 80)
    
    try:
        # 1. Load data
        logger.info("Loading data...")
        data = load_all(args.inputpath)
        
        # Build features
        I_raw = build_binary_I(data["P"], data["V"], data["C"], 0.5)
        I_filtered = I_raw[I_raw.eq(1).any(axis=1)]
        I_filtered = I_filtered.replace({np.nan: 0})
        
        df_features = data['features']
        df_features_new, _ = update_features_matrix(I_filtered, data['df_reads'], df_features, 0.05)
        
        # Cell types
        if args.celltype_file and args.celltype_file != "None":
            df_celltype = pd.read_csv(args.celltype_file, sep="\t")
        else:
            barcodes = I_filtered.index.tolist()
            df_celltype = pd.DataFrame({
                "barcode": barcodes,
                "cell_type": ["default_type"] * len(barcodes)
            })
        
        logger.info(f"Loaded: {I_filtered.shape[0]} cells, {I_filtered.shape[1]} mutations")
        
        # 2. Load optimal tree
        logger.info(f"Loading optimal tree: {args.tree_json}")
        with open(args.tree_json, 'r') as f:
            tree_dict = json.load(f)
        optimal_tree = _dict_to_tree(tree_dict)
        
        # Get mutations on optimal tree
        optimal_mutations = set()
        for node in optimal_tree.traverse():
            if node.name != "ROOT":
                optimal_mutations.update(node.name.split("|"))
        optimal_mutations = list(optimal_mutations)
        
        logger.info(f"Optimal tree has {len(optimal_mutations)} mutations")
        
        # 3. Initialize calculator
        calculator = BootstrapSupportCalculator(
            data=data,
            df_features_new=df_features_new,
            df_celltype=df_celltype,
            optimal_tree=optimal_tree,
            optimal_tree_mutations=optimal_mutations,
            sampleid=args.sampleid,
            outputpath=args.outputpath,
            params=DEFAULT_PARAMS,
            n_bootstrap=args.n_bootstrap,
            n_jobs=args.n_jobs,
            random_seed=args.seed,
            logger=logger,
            cv_rank_thresh=args.cv_rank_thresh,
            is_filter_quality=args.is_filter_quality,
            remove_artifact_mutations=args.remove_artifact,
            min_clade_size=args.min_clade_size
        )
        
        # 4. Compute best-matching-clade Jaccard bootstrap support
        support_values = calculator.compute_bootstrap_support()
        
        # 5. Save results
        calculator.save_results(support_values)
        
        logger.info("=" * 80)
        logger.info("Bootstrap analysis completed successfully!")
        logger.info("Output files:")
        logger.info(f"  - {args.outputpath}/bootstrap_support_percentage.csv")
        logger.info(f"  - {args.outputpath}/bootstrap_support_percentage.json")
        logger.info("=" * 80)
        
    except Exception as e:
        logger.error(f"Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
