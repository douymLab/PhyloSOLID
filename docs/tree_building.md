# Tree-building controls and workflow

[Back to README](../README.md) · [Usage](usage.md) · [Output stages](phylo_output_stages.md)

## Direct entry points versus pipeline CLI

The pipeline CLI (`phylosolid scrna`, `scdna`, `spacetracer`) orchestrates preprocessing
and tree building. It currently does **not** expose `--cv_rank_thresh`,
`--max_workers`, or resolution-backoff flags. Pass those to the underlying Python
entry points when using prepared inputs. Run these commands from the repository root.

```bash
python -m src.run_phylosilid_fullTree_scRNA \
  --sampleid demo_scrna \
  --inputpath demo/scrna/output/demo_scrna/02_treeinput/data \
  --outputpath results/demo_scrna_tree \
  --features_file demo/scrna/output/demo_scrna/01_features/demo_scrna.benchmark_patched.feature.txt \
  --celltype_file demo/scrna/input/celltype.txt \
  --is_predict_germ yes \
  --cv_rank_thresh auto --max_workers 4
```

Other entry modules are `src.run_phylosilid_fullTree_scDNA`,
`src.run_phylosilid_fullTree_spacetracer`, and `src.run_phylosilid_fullTree_binput`.
Use `python -m MODULE --help` for the mode's arguments. For binary input:

```bash
python -m src.run_phylosilid_fullTree_binput \
  --sampleid demo_binput --inputfile demo/binary-matrix/input/demo_input.tsv \
  --outputpath results --disable_resolution_backoff
```

## Coverage CV threshold search

Here CV refers to the coverage-based filtration parameter, not cross-validation.
It is available in scRNA, scDNA and SpaceTracer; binary-matrix mode has no CV search.

| Direct-entry argument | Behavior |
|---|---|
| `--cv_rank_thresh 0.3` | One run |
| `--cv_rank_thresh 0.3,0.5,0.7` | Search the listed thresholds |
| `--cv_rank_thresh 0.3-0.7:0.1` | Search the range |
| `--cv_rank_thresh auto` | Default: eight values, `0.3` through `1.0` in increments of `0.1` |
| `--max_workers 4` | Limit parallel CV workers |

Each candidate repeats scaffold construction, full-tree refinement and QC.
For successful candidates, the current selection policy is:

1. Prefer usable Step 8 (`final`) candidates over `pre_fallback` candidates.
2. Compute `pruning_ratio = (omega_pre_qc - omega_final) / max(omega_final, 1e-8)`.
3. Among final candidates with ratio `< 10.0`, minimize `omega_final`, then
   `omega_pre_qc`, then CV value.
4. If no final candidate passes that guard, minimize `omega_pre_qc`, then
   `omega_final`, then CV value.
5. If every successful candidate is a pre-QC fallback, minimize `omega_pre_qc`,
   then CV value; fallback outputs are not labeled confident.

`omega_pre_qc = omega_final = 0` therefore has ratio zero. A positive pre-QC Omega
with zero final Omega receives a large ratio. Stable CV tie-breaking makes the
choice independent of candidate completion order.

The code retains names such as `is_confident` and `confident` for compatibility.
These indicate an operational pruning guard, **not a confidence probability** or
a guarantee of biological accuracy. Omega measures discordance with observations:
`delta_FP + lambda * delta_FN`. Default lambda is `1` for scDNA and `0.1` for the
other modes. Pre-QC and final scores describe different retained cell/mutation sets.
ROOT-assigned mutations are included consistently in current scoring and exports.

Search output is `cv_search_summary.csv` plus `cv_selection_report.txt`.
The search summary uses `cv`, `confident`, `omega_pre_qc`, `omega_final`,
`omega_reduced`, `pruning_ratio`, `output_stage`, and `zero_final_omega`, along
with dimensions, retention and flip counts. Single-threshold runs instead write
`cv_threshold_search_results.csv`, using `cv_value` and `is_confident`.

## Resolution backoff

scRNA, SpaceTracer and binary-matrix default to automatic resolution with a
backoff conflict budget of `0.10`. Direct entry points accept:

```text
--resolution_backoff_threshold 0.10
--disable_resolution_backoff
```

The options are mutually exclusive. Fixed `resolution_of_graph` takes precedence;
scDNA retains its fixed `0.1` default, so accepting a backoff flag does not enable
auto-resolution for scDNA. Changing that policy currently requires editing its
parameter dictionary. See [resolution_backoff.md](resolution_backoff.md) for the
fixed-graph procedure, audit files and empirical limitations.

## Internal tree-building steps

These numbers refer to the main tree-building workflow, not the CLI stage names
or the scaffold builder's internal numbering. Front-end filtering differs by mode.

| Step | Purpose |
|---|---|
| 1 | Load observations, likelihoods, coverage and mutation features |
| 2 | Identify candidate mosaic mutations when classifier features are available |
| 3 | Germline filtering when enabled |
| 4 | Build scaffold: filter sites, group with Leiden, select backbone and resolve subgroups |
| 5 | Optional DP pass-tree filtering and preparation of accessory mutations |
| 6 | Place accessory mutations on the scaffold in successive integration passes |
| 7 | Discordance-guided prune-and-regraft refinement; calculate pre-QC Omega |
| 8 | Cell/artifact QC, related reattachment and conflict resolution |
| 9 | Prepare and export the selected tree and its Step 7 snapshot |

Step 7 checks clone-level FP, global FP/FN, locus-level FP and ancestral support;
flagged mutations are pruned and reconsidered for placement. Step 8 can remove
orphan/doublet/chimeric cells and, when enabled, artifact mutations, as well as
change mutation placements. A Step 7 snapshot is already refined; it is not raw
input or a tree before all pruning. See [full_tree_builder.py](../src/full_tree_builder.py).

For the complete CLI, mode-dictionary and configuration catalog, see [parameters.md](parameters.md).
