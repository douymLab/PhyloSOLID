# Parameter reference

[Back to README](../README.md) · [Usage](usage.md) · [Workflow and scoring](tree_building.md)

This reference describes the current repository code (2026-10-27). It covers every
explicit argument in `cli/main.py`, the four `run_phylosilid_fullTree_*` entry points,
every field in their `SETTING_PARAMS`, configuration-template paths, and relevant
runtime/demo environment controls. Internal helper function arguments are described
where they implement the key mechanisms; this is not an API reference for every
private function or a parameter reference for external genotypers/PhyloSOLIDvis.

[Key mechanisms](#key-mechanisms) · [Mode dictionaries](#complete-mode-dictionary-catalog) · [Direct-entry CLI](#direct-entry-cli-catalog) · [Pipeline CLI](#pipeline-cli-catalog) · [Configuration](#configuration-paths) · [Environment](#environment-and-demo-controls) · [Fixed controls](#fixed-implementation-controls-and-output-policy)

Start with the key mechanisms, then use the complete catalogs below.
An accepted argument or a dictionary entry does **not** by itself establish that
its value reaches an active computation. Inactive and partially forwarded settings
are identified explicitly. No new command-line knobs are introduced by this document.

## Where to set a value

| Layer | How to set it | Scope and precedence |
|---|---|---|
| Pipeline CLI | `phylosolid [global options] MODE [mode options]` | Preprocessing and orchestration; limited argument forwarding |
| Direct tree entry | `python -m src.run_phylosilid_fullTree_scRNA ...` (or other mode) | Prepared inputs, CV workers/search and backoff switches |
| Mode dictionary | Edit `SETTING_PARAMS` in the relevant entry file | Algorithm weights/thresholds; not automatically exposed through YAML or CLI |
| Configuration | `config/paths.yaml` / selected `--config` | External tools and resources; not a generic algorithm-parameter override |
| Fixed constants/call arguments | Edit the documented consumer function/call site | Includes active graph thresholds, pruning guard and backoff step limits |

Dictionary tables give semantic domains, not promises that the code validates every
range. Use finite numerical values, valid fractions and positive worker counts.
Run from the repository root; choose a new output directory when comparing settings.

## Key mechanisms

### CV candidate generation and selection

`cv_rank_thresh` is a **retained rank fraction**, not a raw coefficient-of-variation
cutoff and not cross-validation. The current coverage filter computes CV over a
binary covered/not-covered indicator, ranks finite CVs ascending, and for N>10
retains `min(N, max(10, ceil(N*r)))`; N<=10 keeps all ranked sites. It also applies
cell-type coverage filtering. Larger r generally admits more scaffold candidates;
N is the number of finite ranked candidates, not the original mutation count.

Direct-entry `--cv_rank_thresh` accepts a number, comma-separated list, range such
as `0.3-0.7:0.1`, or `auto` (default eight values 0.3–1.0). Each candidate repeats
Steps 4–8. Binary-matrix skips coverage filtration/CV search. With
`is_filter_quality=no`, varying this fraction does not activate the skipped filter.

For nonnegative finite Omega values:

```text
Omega = delta_FP + lambda * delta_FN
pruning_ratio = (Omega_pre - Omega_final) / max(Omega_final, 1e-8)
```

`PRUNING_CONFIDENCE_THRESHOLD=10.0` is a fixed entry-module constant, not a CLI flag.
Prefer usable `final` candidates over `pre_fallback`. Among finals satisfying ratio
<10, sort by `(Omega_final, Omega_pre, CV)`; if none satisfy it, use
`(Omega_pre, Omega_final, CV)`. If all successful outputs are pre-QC fallbacks, use
`(Omega_pre, CV)` and do not label them confident. Zero/zero has ratio 0;
positive/zero gives a large ratio. Omega can increase across stages, yielding a
negative ratio; the implementation does not clip it. The guard is not a calibrated
confidence probability or proof of topology accuracy. Compare independent topology
metrics, founder recovery and retained populations in an evaluation.

### Jaccard and directional support

The active consumer is `_build_pairwise_correlation_cache_fast` in
[scaffold_builder.py](../src/scaffold_builder.py). Its current rule is fixed in code:

```text
J = N11 / (N11 + N10 + N01)             # 0 when denominator is 0
F_forward = N11 / positive_calls(u)    # 0 when denominator is 0
F_reverse = N11 / positive_calls(v)
correlated = N11 >= 1 AND (
    J >= 0.08 OR (0 < J < 0.08 AND max(F_forward, F_reverse) >= 0.5)
)
```

Pair counts use observed 0/1 calls. Directional denominators count each mutation's
positive calls across its full vector, including positives whose partner is missing;
they need not equal pairwise-complete `N11+N10`/`N11+N01`. The support rule rescues
asymmetric/co-nested patterns with low Jaccard; it does not independently prove a
parent-child relationship. Diagonal pairs are set true.

**The declared dictionary fields `pair_N11_min=0`, `jaccard_thresh=0.2`,
`jaccard_low=0.1`, and `fraction_parent_child_thresh=0.9` are not wired into this
active rule.** Changing those fields alone does not change graph edges. To alter
the active thresholds, edit the consumer; this documentation does not silently wire
them up or replace the existing algorithm. Increasing the active Jaccard/support
requirements tends to reject more pair relationships, with downstream graph effects.

### FP/FN weights and placement

Here FP means `I=1, M=0`; FN means `I=0, M=1`. They are observed/inferred
disagreements, not independently measured sequencing error probabilities.
`fnfp_ratio=lambda` defaults to 0.1 in scRNA/SpaceTracer/binary and 1 in scDNA.
For 10 FP and 20 FN, Omega is 12 at lambda=0.1 and 30 at lambda=1.
The same weight also changes mutation placement, so a sensitivity experiment must
rerun the pipeline, not merely reweight an old final table. Empirical error summaries
and external topology metrics answer different questions and should be reported separately.

`general_weight_NA=0.001` is the base missing-data cost, dynamically adjusted in
`compute_dynamic_penalty`; NA-to-1 proportions >0.05, >0.1 and >0.2 multiply the
base by 1.5, 2 and 4. Above 0.4, 0.6 and 0.8, weights use 0.5, 0.7 and 0.9 times
lambda; a cap of 0.95*lambda applies. Posterior-dependent NA weights enter the
placement calculation, while exported Omega counts only FP and FN.

`phi=1` controls a node-complexity component. The `consider_ROOT` placement path
uses the original `phi*log(N_nodes)`, not the dynamically adjusted phi, and adds
other complexity, root, merge and intersection penalties. It is not a single
universal topology regularizer with a guaranteed monotonic effect.

### Reattachment and QC

Step 7 prunes flagged mutations and attempts reintegration; it does not equate a
high conflict ratio with permanent artifact removal. Its clone/global/locus FP,
global FN and ancestral-support cutoffs are listed below. Scaffold/backbone
protection and group handling mean threshold crossing is only part of the decision.
Step 8 applies distinct orphan-cell and artifact/chimeric-cell rules and can also
reattach mutations. Turning `remove_artifact_mutations` off does not disable cell QC.

There is currently **no dedicated public ablation flag** for each Step 7 submodule.
A controlled module ablation requires a deliberate code change, fixed inputs/seeds
and independently scored outputs. Setting an arbitrary extreme threshold is not
necessarily equivalent to disabling a module, especially for conjunctions and
`intersection_count == cutoff` tests.

For the per-mutation artifact test, `fp_cells_ratio_per_mutation` is FP cells
divided by observed positive cells for that mutation (zero positives gives zero).
For the per-cell test, `fp_muts_ratio_per_cell` is FP mutations divided by observed
positive mutations for that cell. These denominators differ from all-covered-cell
or all-covered-site fractions. Step 7 uses its own computed metrics rather than
reusing these Step 8 fields.

The orphan-cell route tests `intersection_count == cutoff` AND
`flipping_count_fn >= cutoff` AND `flipping_to_1_count >= cutoff`; it is joined
with separate outgroup-derived cell candidates. The artifact route combines its
ratio and count conditions with AND and protects scaffold/backbone sites.

### Resolution and conflict-budget backoff

`resolution_of_graph=None` selects initial resolution using DQI. With bounds 0.5
and 2.0, DQI levels >=0.7, >=0.5, >=0.3 and below 0.3 map to 20%, 50%, 75% and
100% of the interval. DQI combines NA score (0.4), rare-site score (0.35), and
coverage score (0.25); rare sites have <=3 positives, effective sites >=5 observed
calls. Rare fractions >0.3/>0.2 add 0.2/0.1, then the result is clipped and rounded
to one decimal. These DQI weights/cutoffs are internal constants.

The backoff budget `tau=0.10` is the fraction of informative cells allowed to have
`q_i >= fp_ratio_per_cell_cross_all_muts_cutoff` (normally 0.5), with
`q_i=1-max_g(n_ig)/sum_g(n_ig)`. It is separate from lambda and the pruning guard.
Backoff decreases gamma by 0.1 on the same graph, accepts the first passing
partition, and retains the initial partition if no candidate passes within the
floor/step limit. Lower tau is stricter but may induce merging; one group has zero
conflict by construction. It is not an estimated doublet rate or an optimality proof.
A numeric fixed gamma bypasses this mechanism; scDNA defaults to fixed gamma=0.1.
See [algorithm, audit and evaluation details](resolution_backoff.md).

## Complete mode dictionary catalog

The following defaults are extracted from each entry's `SETTING_PARAMS`.
`all` means identical in all four modes. Otherwise values are in the order
**scRNA / scDNA / SpaceTracer / binary-matrix**. `—` means not declared.
Unless marked inactive, the field is read on the applicable branch; skipped
classifier/germline/quality/DP stages do not become active merely by editing a threshold.

### Genotype and germline inputs

| Field | Declared defaults | Type / semantic domain | Behavior and conditions |
|---|---|---|---|
| `p_thresh` | `0.5 / 0.9 / 0.5 / 0.5` | float, [0,1] | build_binary_I calls a covered site positive only when V>0 and P>p_thresh; other covered calls become 0, uncovered calls stay missing. Increasing it reduces positive calls. Binary input is already called and bypasses this binarization. Distinct from scDNA preprocessing --threshold. |
| `mcf_cutoff` | `all: 0.05` | float, [0,1] | identify_germline_variants considers sites with detected mutant-cell fraction > this cutoff; update_features_matrix also receives it. Raising it reduces germline-test candidates, not a universal deletion of rare sites. Germline effect requires is_predict_germ=yes. |
| `mcn_cutoff` | `all: 5` | integer, >=0 | update_germline_status resets sites with mutant-cell count <= this threshold to non-germline. Raising it makes more low-count sites eligible for rescue from germline status; relevant when germline filtering runs. |

### Declared but inactive fields

| Field | Declared defaults | Type / semantic domain | Behavior and conditions |
|---|---|---|---|
| `models_path` | `all: 'phylosolid/models/scdna'` | path string | Read into a local variable but not passed to the real-time classifier. Changing this field does not select a saved model in the current four entry paths. Classifier routines manage their own fitting/loading behavior. |
| `pair_N11_min` | `all: 0` | integer, >=0 | Not read by the active graph builder. The current fast pair cache uses N11>=1, not the declared value 0. See the Jaccard section. |
| `jaccard_thresh` | `all: 0.2` | float, [0,1] | Not read by the active graph builder. Its high-Jaccard condition currently uses J>=0.08, not the declared 0.2. |
| `jaccard_low` | `all: 0.1` | float, [0,1] | Not read by the active graph builder. Current rescue condition is 0<J<0.08 with directional support >=0.5; the declared 0.1 is not used as its lower bound. |
| `fraction_parent_child_thresh` | `all: 0.9` | float, [0,1] | Not read by the active graph builder. Current directional-support threshold is fixed at 0.5, not the declared 0.9. |
| `consensus_runs` | `all: 100` | integer, >0 | Not read from params on the active graph path. Call sites pass n_shuffle=100; compute_clone_and_pair_weights uses one deterministic permutation per reference instead of repeating random shuffles. Changing this dictionary field does not change repetitions. |
| `consensus_clone_freq_thresh` | `all: 0.1` | float, [0,1] | Declared in all mode dictionaries, with no active consumer on the four entry paths. Changing it does not filter consensus clone frequencies. |
| `fp_count_persite_cutoff` | `all: 0` | integer, >=0 | Declared, but Step 7.3 selects candidates by fp_ratio_persite alone. This field does not impose an additional FP-count minimum. |
| `cv_thresh` | `— / 6.0 / — / —` | float, >=0 | Only declared in scDNA SETTING_PARAMS. The current scaffold path uses coverage_filters and cv_rank_thresh; it does not consume this absolute CV threshold. Older helper functions accept similarly named arguments, which does not make this entry setting active. |

### Initial scaffold filtration

| Field | Declared defaults | Type / semantic domain | Behavior and conditions |
|---|---|---|---|
| `posterior_threshold` | `all: 0.5` | float, [0,1] | initial_filter retains a cell if it has at least one V>0 call and at least one P>threshold call (the two any conditions need not refer to the same site). Higher values make posterior support harder to obtain. Used only with is_filter_quality=yes. |
| `maf_max_threshold` | `all: 0.3` | float, [0,1] | A scaffold-candidate site must have max V across retained cells >= this value, together with the mean-V condition. Raising it makes site inclusion stricter; used only with is_filter_quality=yes. |
| `maf_mean_threshold` | `all: 0.1` | float, [0,1] | The mean V calculated over V>0 cells must be >= this threshold, together with max V. Raising it makes inclusion stricter. It is not the mean across all cells. |
| `na_prop_thresh_global` | `all: 0.95` | float [0,1], or None | Per-cell-type no-coverage proportion must be strictly below the threshold. None derives per-type Q3 thresholds. Higher values tolerate more missing coverage; downstream requirements depend on the number and dominance of cell types. Used with quality filtration. |

### Leiden resolution

| Field | Declared defaults | Type / semantic domain | Behavior and conditions |
|---|---|---|---|
| `resolution_of_graph` | `None / 0.1 / None / None` | positive float, or None | None invokes DQI-based auto resolution, then optional backoff. A numeric fixed resolution bypasses both auto calculation and backoff and is inherited by downstream subgroup calls. Higher gamma tends to split more groups, without guaranteeing monotone partitions or clone counts. |
| `min_resolution` | `all: 0.5` | positive float | Lower bound for auto resolution and backoff; ignored for fixed resolution. Lower values allow coarser candidates. Keep <= max_resolution; that relationship is a semantic requirement, not a comprehensive CLI validator. |
| `max_resolution` | `all: 2.0` | positive float | Upper bound for initial DQI auto resolution. Backoff starts at the selected initial value and only moves down. Raising this permits a finer initial partition; ignored for fixed resolution. |

### Placement penalties

| Field | Declared defaults | Type / semantic domain | Behavior and conditions |
|---|---|---|---|
| `general_weight_NA` | `all: 0.001` | float, >=0 | Base missing-call penalty passed into placement. Actual weights can be dynamically adjusted using NA-to-1 flips and fnfp_ratio; see compute_dynamic_penalty. Increasing the base tends to penalize missing-data imputation, but the dynamic branches mean the final score is not a simple constant times NA count. |
| `fnfp_ratio` | `0.1 / 1 / 0.1 / 0.1` | float, >0 | FN cost relative to FP baseline 1. Used by placement and Omega; scDNA defaults to 1, others to 0.1. Raising it makes an observed 0 inferred as 1 more costly, and also changes placement rankings and CV scores. It is not a measured error probability. |
| `phi` | `all: 1.0` | float, >=0 | Base coefficient for node-complexity penalty. In the consider_ROOT placement path, BIC_penalty=phi*log(N_nodes); the total also has log(N_nodes), root, merge and intersection-related penalties. That path uses original phi even though a dynamically adjusted value is computed. Increasing phi strengthens this complexity component, not necessarily monotonically reducing final nodes. |

### Optional DP filtering

| Field | Declared defaults | Type / semantic domain | Behavior and conditions |
|---|---|---|---|
| `pass_tree_cutoff` | `all: 0.9` | float, [0,1] | For sites with >=2 mutant cells, successful_pass requires somatic_posterior_per_site >= this value AND somatic_posterior_per_site_onecell < unpass_tree_cutoff. Raising this cutoff makes passing stricter; active only with DP enabled. Sites with 0/1 mutant cells are labeled absent/cell_specific separately. |
| `unpass_tree_cutoff` | `all: 0.1` | float, [0,1] | Upper limit on somatic_posterior_per_site_onecell for successful_pass: the comparison is strictly <. Raising this value relaxes that part of pass acceptance. It concerns a different posterior from pass_tree_cutoff, not the lower end of one interval; active only with DP enabled. |

### Step 7 reattachment

| Field | Declared defaults | Type / semantic domain | Behavior and conditions |
|---|---|---|---|
| `fp_ratio_cutoff_within_subclone` | `all: 0.1` | float, [0,1] | Step 7.1 flags fp_ratio_within_subclone_for_in_subclone_muts >= threshold. Lower values nominate more sites; backbone/scaffold protection and downstream group handling still apply. Nomination is for reattachment, not automatic permanent removal. |
| `fp_ratio_cutoff_across_tree` | `all: 0.2` | float, [0,1] | Step 7.2 flags global fp_ratio >= threshold. Lower values nominate more FP-discordant sites and related groups for pruning/reintegration; protected backbone/scaffold sites are treated separately. |
| `fn_ratio_cutoff_across_tree` | `all: 0.9` | float, [0,1] | Step 7.2 flags global fn_ratio >= threshold. Lower values nominate more FN-discordant sites for reconsideration. This threshold does not change fnfp_ratio or score weights. |
| `fp_ratio_persite_cutoff` | `all: 0.1` | float, [0,1] | Step 7.3 flags fp_ratio_persite >= threshold. Lower values make locus-level reattachment more aggressive; the declared fp_count_persite_cutoff is not combined with this test. |

### Step 7 ancestral support

| Field | Declared defaults | Type / semantic domain | Behavior and conditions |
|---|---|---|---|
| `intersection_vs_fn_flipping_ratio_cutoff` | `all: 0.2` | float, [0,1] | Step 7.4 requires parent_retention_ratio <= threshold, AND either the low-count or low-ratio test below. Raising this threshold nominates more sites lacking ancestral support. It uses the returned parent_retention_ratio field, not an arbitrary new FP/FN quotient. |
| `intersection_cell_count_on_mutation_cutoff` | `all: 5` | integer, >=0 | Count branch of Step 7.4: intersection_cell_count_on_mutation <= threshold. Raising it relaxes this low-support test; combined by OR with the ratio branch, then AND with parent retention. |
| `intersection_cell_ratio_on_mutation_cutoff` | `all: 0.2` | float, [0,1] | Ratio branch of Step 7.4: intersection_cell_ratio_on_mutation <= threshold. Raising it nominates more low-support sites, subject to the parent-retention condition. |

### Step 8 artifact and cell QC

| Field | Declared defaults | Type / semantic domain | Behavior and conditions |
|---|---|---|---|
| `fp_ratio_per_mutation_cross_all_cells_cutoff` | `0.2 / 0.4 / 0.2 / 0.2` | float, [0,1] | Artifact candidates need fp_cells_ratio_per_mutation >= this threshold AND fp_cells_count >= the next threshold, excluding protected scaffold/backbone mutations. Lower values nominate more sites. Removal requires remove_artifact_mutations=yes. |
| `fp_count_per_mutation_cross_all_cells_cutoff` | `all: 5` | integer, >=0 | Absolute count branch of artifact nomination: fp_cells_count >= threshold. Higher values require more discordant cells, making artifact nomination stricter. |
| `fp_ratio_per_cell_cross_all_muts_cutoff` | `all: 0.5` | float, (0,1] for backoff | Step 8 removes cells with fp_muts_ratio_per_cell >= threshold. Lower values flag more cells. The same value also defines high-conflict q_i>=cutoff in resolution backoff; changing it changes both mechanisms. These two ratios use different contexts/populations. |

### Step 8 orphan QC

| Field | Declared defaults | Type / semantic domain | Behavior and conditions |
|---|---|---|---|
| `intersection_count_per_cells_cutoff` | `all: 1` | integer, >=0 | Orphan-cell condition uses intersection_count == this value, NOT <=. Raising it changes which exact count is targeted, so the result need not become monotonically more aggressive. Combined with both flip-count tests below; outgroup-derived doublet candidates are added separately. |
| `flipping_count_fn_per_cells_cutoff` | `all: 1` | integer, >=0 | Orphan-cell condition additionally requires flipping_count_fn >= threshold. Raising it makes this route to cell removal stricter. |
| `flipping_to_1_count_per_cells_cutoff` | `all: 2` | integer, >=0 | Orphan-cell condition additionally requires flipping_to_1_count >= threshold, referring to NA-to-1 imputations. Raising it makes this route stricter. |

### Runtime dictionary additions and internal state

| Name | Default / source | Meaning |
|---|---|---|
| `cv_rank_thresh` | Candidate value; helper fallback 0.5 | Overwritten per CV worker; use direct-entry `--cv_rank_thresh` |
| `resolution_backoff_threshold` | Direct CLI default 0.10; None disables | Validated finite [0,1]; used only for auto gamma |
| `random_seed` | Scaffold helper fallback 42 | Optional programmatic dictionary key; not declared in the four SETTING_PARAMS; scaffold calls set_seed with it |
| `_output_stage` | Builder generated: final/pre_fallback | Internal output status, not a user switch |
| `_pre_qc_snapshot` | Builder generated tuple | Saved Step 7 tree/matrix/ROOT/conflict state, not a tuning parameter |
| `omega_before_qc` | Added to per-CV export params | Metric metadata, not a user override |

## Direct-entry CLI catalog

Applies to `python -m src.run_phylosilid_fullTree_MODE`; underscores in these flag
names differ from pipeline CLI hyphens. Placeholder sample/input/output defaults
are not runnable example datasets. `-h`/`--help` is automatically available in every parser.

### scRNA

| Flags | Type / choices | Parser default | Description / implementation note |
|---|---|---|---|
| `-s` / `--sampleid` | `str` | `''` | The sampleid you can set and check. |
| `-i` / `--inputpath` | `str` | `'yourpath/data'` | Directory of prepared posterior/likelihood/allele-count and feature matrices, not a BAM. |
| `-o` / `--outputpath` | `str` | `'yourpath/results'` | scRNA/scDNA/SpaceTracer write numbered steps directly here; binary appends sampleid under this parent. |
| `-c` / `--celltype_file` | `str` | `None` | Precomputed TSV; None uses default_type annotation. Exported annotation is aligned separately. |
| `--features_file` | `str` | `None` | scRNA and scDNA use different feature schemas; None/no bypasses classifier. SpaceTracer pipeline passes None. |
| `--is_predict_germ` | `str; ['yes', 'no']` | `'no'` | Direct default no; scRNA BAM wrapper explicitly enables yes. Binary mode has no flag and skips germline filtering. |
| `--is_detect_passtree_by_dp` | `str; ['yes', 'no']` | `'no'` | Default no; yes activates Step 5 DP cutoffs. |
| `--is_filter_quality` | `str; ['yes', 'no']` | `'yes'` | Default yes; no skips initial/coverage scaffold filtration. Binary mode fixes no. |
| `--cv_rank_thresh` | `str` | `'auto'` | Default auto; scalar/list/range formats and selection rules are described above. |
| `--remove_artifact_mutations` | `str; ['yes', 'no']` | `'yes'` | Controls Step 8 artifact removal only; no still runs cell QC. Binary mode fixes yes and has no switch. |
| `--seed` | `int` | `42` | Entry calls set_seed; scaffold later resets from params.get(random_seed,42), and CV workers use fixed base 42 plus CV-derived seed. This is not an end-to-end guarantee that every phase follows a supplied seed. |
| `--max_workers` | `int` | `None` | CV multiprocessing only; None resolves to max(1,cpu_count()-1), capped by number of candidates. Give a positive integer; distinct from preprocessing --threads. |
| `--resolution_backoff_threshold` | `float` | `0.1` | Finite [0,1], mutually exclusive with disable; fixed gamma still takes precedence. |
| `--disable_resolution_backoff` | `flag` | `'argparse.SUPPRESS'` | Sets budget to None; suppresses parser default override. Does not disable DQI auto-resolution or alter fixed gamma. |

### scDNA

| Flags | Type / choices | Parser default | Description / implementation note |
|---|---|---|---|
| `-s` / `--sampleid` | `str` | `''` | The sampleid you can set and check. |
| `-i` / `--inputpath` | `str` | `'yourpath/data'` | Directory of prepared posterior/likelihood/allele-count and feature matrices, not a BAM. |
| `-o` / `--outputpath` | `str` | `'yourpath/results'` | scRNA/scDNA/SpaceTracer write numbered steps directly here; binary appends sampleid under this parent. |
| `-c` / `--celltype_file` | `str` | `None` | Precomputed TSV; None uses default_type annotation. Exported annotation is aligned separately. |
| `--features_file` | `str` | `None` | scRNA and scDNA use different feature schemas; None/no bypasses classifier. SpaceTracer pipeline passes None. |
| `--is_predict_germ` | `str; ['yes', 'no']` | `'no'` | Direct default no; scRNA BAM wrapper explicitly enables yes. Binary mode has no flag and skips germline filtering. |
| `--is_detect_passtree_by_dp` | `str; ['yes', 'no']` | `'no'` | Default no; yes activates Step 5 DP cutoffs. |
| `--is_filter_quality` | `str; ['yes', 'no']` | `'yes'` | Default yes; no skips initial/coverage scaffold filtration. Binary mode fixes no. |
| `--cv_rank_thresh` | `str` | `'auto'` | Default auto; scalar/list/range formats and selection rules are described above. |
| `--remove_artifact_mutations` | `str; ['yes', 'no']` | `'yes'` | Controls Step 8 artifact removal only; no still runs cell QC. Binary mode fixes yes and has no switch. |
| `--seed` | `int` | `42` | Entry calls set_seed; scaffold later resets from params.get(random_seed,42), and CV workers use fixed base 42 plus CV-derived seed. This is not an end-to-end guarantee that every phase follows a supplied seed. |
| `--max_workers` | `int` | `None` | CV multiprocessing only; None resolves to max(1,cpu_count()-1), capped by number of candidates. Give a positive integer; distinct from preprocessing --threads. |
| `--resolution_backoff_threshold` | `float` | `0.1` | Finite [0,1], mutually exclusive with disable; fixed gamma still takes precedence. |
| `--disable_resolution_backoff` | `flag` | `'argparse.SUPPRESS'` | Sets budget to None; suppresses parser default override. Does not disable DQI auto-resolution or alter fixed gamma. |

### spacetracer

| Flags | Type / choices | Parser default | Description / implementation note |
|---|---|---|---|
| `-s` / `--sampleid` | `str` | `''` | The sampleid you can set and check. |
| `-i` / `--inputpath` | `str` | `'yourpath/data'` | Directory of prepared posterior/likelihood/allele-count and feature matrices, not a BAM. |
| `-o` / `--outputpath` | `str` | `'yourpath/results'` | scRNA/scDNA/SpaceTracer write numbered steps directly here; binary appends sampleid under this parent. |
| `-c` / `--celltype_file` | `str` | `None` | Precomputed TSV; None uses default_type annotation. Exported annotation is aligned separately. |
| `--is_predict_germ` | `str; ['yes', 'no']` | `'no'` | Direct default no; scRNA BAM wrapper explicitly enables yes. Binary mode has no flag and skips germline filtering. |
| `--is_detect_passtree_by_dp` | `str; ['yes', 'no']` | `'no'` | Default no; yes activates Step 5 DP cutoffs. |
| `--is_filter_quality` | `str; ['yes', 'no']` | `'yes'` | Default yes; no skips initial/coverage scaffold filtration. Binary mode fixes no. |
| `--cv_rank_thresh` | `str` | `'auto'` | Default auto; scalar/list/range formats and selection rules are described above. |
| `--remove_artifact_mutations` | `str; ['yes', 'no']` | `'yes'` | Controls Step 8 artifact removal only; no still runs cell QC. Binary mode fixes yes and has no switch. |
| `--seed` | `int` | `42` | Entry calls set_seed; scaffold later resets from params.get(random_seed,42), and CV workers use fixed base 42 plus CV-derived seed. This is not an end-to-end guarantee that every phase follows a supplied seed. |
| `--max_workers` | `int` | `None` | CV multiprocessing only; None resolves to max(1,cpu_count()-1), capped by number of candidates. Give a positive integer; distinct from preprocessing --threads. |
| `--resolution_backoff_threshold` | `float` | `0.1` | Finite [0,1], mutually exclusive with disable; fixed gamma still takes precedence. |
| `--disable_resolution_backoff` | `flag` | `'argparse.SUPPRESS'` | Sets budget to None; suppresses parser default override. Does not disable DQI auto-resolution or alter fixed gamma. |

### binput

| Flags | Type / choices | Parser default | Description / implementation note |
|---|---|---|---|
| `-s` / `--sampleid` | `str` | `''` | Sample ID |
| `-o` / `--outputpath` | `str` | `'./output'` | scRNA/scDNA/SpaceTracer write numbered steps directly here; binary appends sampleid under this parent. |
| `-i` / `--inputfile` | `str` | `''` | Input binary matrix file (rows=cells, columns=mutations) |
| `-c` / `--celltype_file` | `str` | `None` | Precomputed TSV; None uses default_type annotation. Exported annotation is aligned separately. |
| `--seed` | `int` | `42` | Entry calls set_seed; scaffold later resets from params.get(random_seed,42), and CV workers use fixed base 42 plus CV-derived seed. This is not an end-to-end guarantee that every phase follows a supplied seed. |
| `--resolution_backoff_threshold` | `float` | `0.1` | Finite [0,1], mutually exclusive with disable; fixed gamma still takes precedence. |
| `--disable_resolution_backoff` | `flag` | `'argparse.SUPPRESS'` | Sets budget to None; suppresses parser default override. Does not disable DQI auto-resolution or alter fixed gamma. |

## Pipeline CLI catalog

Global options precede the subcommand, e.g. `phylosolid --workdir results scrna ...`.
The `check-annovar` subcommand also has its own required `--config` after its name.
`-h`/`--help` is implicit. The tables cover all explicit parser options, including
compatibility flags that are currently unused. A direct-entry argument is not
implicitly forwarded just because the spelling looks similar.

### Global

| Flags | Type / choices | Default / requirement | Description / implementation note |
|---|---|---|---|
| `--version` | `flag` | `None` |  |
| `--verbose` / `-v` | `flag` | `False` | Enables CLI logging verbosity; does not automatically enable all internal tree snapshots. |
| `--log-file` | `Path` | `None` | Path to log file |
| `--script-dir` | `Path` | `None` | Explicit external scripts root; overrides discovery using PHYLOSOLID_SCRIPTS / package locations. |
| `--workdir` / `-w` | `Path` | `'./phylosolid_work'` | Parent working directory; pipelines append sample and numbered stages. Binary subcommand instead uses its own --outputpath. |
| `--config` / `-c` | `Path` | `None` | Configuration YAML; not an algorithm SETTING_PARAMS override. check-annovar requires its own instance. |
| `--seed` | `int` | `42` | Declared globally with default 42 but not read/forwarded in cli.main; use direct entry and note its seed limitations above. |

### check-annovar

| Flags | Type / choices | Default / requirement | Description / implementation note |
|---|---|---|---|
| `--config` / `-c` | `Path` | `required` | Configuration YAML; not an algorithm SETTING_PARAMS override. check-annovar requires its own instance. |
| `--threads` / `-t` | `int` | `4` | Preprocessing concurrency; tree-input scRNA/SpaceTracer caps at 2. Does not set direct-entry --max_workers. |

### binary-matrix

| Flags | Type / choices | Default / requirement | Description / implementation note |
|---|---|---|---|
| `-s` / `--sampleid` | `str` | `required` | Sample ID |
| `-i` / `--inputfile` | `str` | `required` | Binary matrix file |
| `-o` / `--outputpath` | `str` | `required` | Output directory |

### spacetracer

| Flags | Type / choices | Default / requirement | Description / implementation note |
|---|---|---|---|
| `--sample` / `-s` | `str` | `required` | Sample ID |
| `--mutation-list` / `-m` | `Path` | `required` | Mutation list file |
| `--bam` / `-b` | `Path` | `required` | BAM file |
| `--barcode` / `-bc` | `Path` | `required` | Barcode file |
| `--celltype-file` | `Path` | `None` | Optional explicit annotation takes precedence over metadata helper; absent annotation falls back to default_type in tree entry. |
| `--metadata` | `Path` | `None` | Invokes the mode annotation helper when no celltype-file is supplied; helper setup may be needed. Prefer explicit celltype-file for reproducible scripted runs. |
| `--read-len` | `int` | `91` | scRNA feature extraction receives read length. SpaceTracer currently accepts it but skips feature extraction and does not pass it in active stage kwargs. |
| `--cellnum` | `int` | `155` | Number of cells (default: 155) |
| `--threads` / `-t` | `int` | `4` | Preprocessing concurrency; tree-input scRNA/SpaceTracer caps at 2. Does not set direct-entry --max_workers. |
| `--running-type` | `str; ['benchmark', 'production']` | `'benchmark'` | scRNA feature-script mode and filename prefix; SpaceTracer skips feature extraction and currently does not consume this value. |
| `--ase-filepath` | `str` | `None` | Passed to scRNA feature extraction (None becomes no); currently unused in SpaceTracer. |
| `--steps` | `str; ['tree_input', 'tree_building']` | `None` | Choose pipeline stages, not internal Step 1–9. Later stages require earlier outputs in the same sample directory. |
| `--parallel` | `flag` | `False` | scRNA runs feature extraction/tree_input concurrently; SpaceTracer accepts it but currently executes stages sequentially. |

### scrna

| Flags | Type / choices | Default / requirement | Description / implementation note |
|---|---|---|---|
| `--sample` / `-s` | `str` | `required` | Sample ID |
| `--mutation-list` / `-m` | `Path` | `required` | Mutation list file |
| `--bam` / `-b` | `Path` | `required` | BAM file |
| `--barcode` / `-bc` | `Path` | `required` | Barcode file |
| `--celltype-file` | `Path` | `None` | Optional explicit annotation takes precedence over metadata helper; absent annotation falls back to default_type in tree entry. |
| `--metadata` | `Path` | `None` | Invokes the mode annotation helper when no celltype-file is supplied; helper setup may be needed. Prefer explicit celltype-file for reproducible scripted runs. |
| `--read-len` | `int` | `91` | scRNA feature extraction receives read length. SpaceTracer currently accepts it but skips feature extraction and does not pass it in active stage kwargs. |
| `--cellnum` | `int` | `155` | Number of cells (default: 155) |
| `--threads` / `-t` | `int` | `4` | Preprocessing concurrency; tree-input scRNA/SpaceTracer caps at 2. Does not set direct-entry --max_workers. |
| `--running-type` | `str; ['benchmark', 'production']` | `'benchmark'` | scRNA feature-script mode and filename prefix; SpaceTracer skips feature extraction and currently does not consume this value. |
| `--ase-filepath` | `str` | `None` | Passed to scRNA feature extraction (None becomes no); currently unused in SpaceTracer. |
| `--steps` | `str; ['feature_extraction', 'tree_input', 'tree_building']` | `None` | Choose pipeline stages, not internal Step 1–9. Later stages require earlier outputs in the same sample directory. |
| `--parallel` | `flag` | `False` | scRNA runs feature extraction/tree_input concurrently; SpaceTracer accepts it but currently executes stages sequentially. |

### scdna

| Flags | Type / choices | Default / requirement | Description / implementation note |
|---|---|---|---|
| `--sample` / `-s` | `str` | `required` | Sample ID |
| `--mutation-list` / `-m` | `Path` | `required` | Mutation list (chr_pos_ref_alt, one per line) |
| `--bam-dir` | `Path` | `one of --bam-dir / --bam-list required` | Directory of single-cell BAM files named {cell}.bam |
| `--bam-list` | `Path` | `one of --bam-dir / --bam-list required` | Text file with one single-cell BAM path per line |
| `--sample-list` | `Path` | `None` | Optional cell IDs (one per line). If omitted, BAM file stems are used |
| `--bulk-bam` | `Path` | `None` | Optional bulk BAM. Not required; omitted by default |
| `--keep-bulk` | `flag` | `False` | Retains a BAM named bulk.bam in single-cell BAM discovery; distinct from optional --bulk-bam pseudobulk input. |
| `--reference` | `Path` | `None` | scDNA optional classifier-feature FASTA; overrides scdna.genome_fasta. |
| `--mappability` | `Path` | `None` | scDNA optional classifier-feature track; overrides scdna.mappability_file. |
| `--celltype-file` | `Path` | `None` | Optional explicit annotation takes precedence over metadata helper; absent annotation falls back to default_type in tree entry. |
| `--threshold` | `float` | `0.9` | scDNA somatic posterior threshold for prepared tree-input generation; default 0.9. Separate from tree-entry p_thresh=0.9. |
| `--threads` / `-t` | `int` | `4` | Preprocessing concurrency; tree-input scRNA/SpaceTracer caps at 2. Does not set direct-entry --max_workers. |
| `--steps` | `str; ['feature_extraction', 'tree_input', 'tree_building']` | `None` | Choose pipeline stages, not internal Step 1–9. Later stages require earlier outputs in the same sample directory. |

## Configuration paths

Copy [paths.yaml.template](../config/paths.yaml.template) to `config/paths.yaml`.
The values below are machine-specific paths, not portable defaults. Existing
configuration loaders/steps route these paths to external tools; adding an arbitrary
weight such as fnfp_ratio to YAML does not forward it to SETTING_PARAMS.

| Key | Meaning / conditions |
|---|---|
| `conda.python` | Python interpreter used by configured shell preprocessing/tree wrappers |
| `conda.env_path` | Tree/preprocessing conda environment path exported to tool wrappers |
| `annovar.script_dir` | Directory containing annotate_variation.pl; required for annotation workflows |
| `annovar.humandb` | ANNOVAR database directory |
| `annovar.build` | Database genome build; template hg38; align all resources to the input build |
| `reference.genome_fasta` | Reference FASTA for preprocessing |
| `reference.gff3_file` | Exon/annotation resource; template points to the supplied BED despite the key name |
| `reference.mappability_file` | Preprocessing mappability resource; template uses k24 bedGraph |
| `reference.gnomad_file` | Population frequency annotation resource |
| `reference.rna_editing_file` | RNA editing annotation resource for scRNA filtering |
| `scdna.genotyper_dir` | Directory containing run_unphased_prod_whether_bulk.py; PHYLOSOLID_GENOTYPER_DIR wins |
| `scdna.genome_fasta` | Optional scDNA classifier-feature reference; --reference wins |
| `scdna.mappability_file` | Optional scDNA classifier-feature track; --mappability wins |
| `genotyper_dir` | Legacy top-level fallback if neither environment nor scdna.genotyper_dir is set |

Binary-matrix does not require the BAM preprocessing resource bundle. scDNA bulk
BAM is optional; configuring resource paths alone does not create missing tools.

## Environment and demo controls

| Name | Default / priority | Effect |
|---|---|---|
| `PHYLOSOLID_SCRIPTS` | Used when explicit --script-dir is absent | External scripts discovery |
| `PHYLOSOLID_GENOTYPER_DIR` | Before nested/top-level genotyper config | scDNA genotyper location |
| `PHYLOSOLID_VERBOSE_TREES` | Off; 1/true/yes/on enables | Additional internal tree snapshots in logs; independent of CLI verbosity |
| `PHYLOSOLID_VIS_RSCRIPT` | First visualization choice when executable | Explicit Rscript for circos; tree environments are refused by helper |
| `PHYLOSOLID_VIS_ENV` | Before phylosolid_vis and circos_env | Visualization conda environment name |
| `PHYLOSOLIDVIS_SRC` | Optional | Local visualization source used by install_vis.sh |
| `MPLCONFIGDIR` | Demo fallback /tmp/mplconfig | Writable matplotlib configuration/cache directory |
| `THREADS` | 4 in scRNA/scDNA demo scripts | Demo preprocessing threads, not CV workers |
| `READ_LEN` | 120 in scRNA demo | Demo read length |
| `CELL_NUM` | 3319 in scRNA demo | Demo cell count |
| `PYTHONPATH` | Demo/wrappers add repository modules | Python import resolution |
| `PYTHONHASHSEED` | Entrypoints set 42; workers set CV-derived values | Implementation state; externally setting it is not a substitute for full seed control |
| `MPLBACKEND` / `QT_QPA_PLATFORM` | Entrypoints set Agg / offscreen | Noninteractive plotting backend |
| `OMP_NUM_THREADS` / `OPENBLAS_NUM_THREADS` / `MKL_NUM_THREADS` | Runtime-library defaults | Optional external numerical-library thread limits; distinct from worker count |

## Fixed implementation controls and output policy

These are documented to distinguish a fixed source-level choice from a public flag.

| Control | Current value / policy | Consumer |
|---|---|---|
| `PRUNING_CONFIDENCE_THRESHOLD` | 10.0; strict ratio < threshold | Each entry's select_best_cv |
| `PRUNING_EPSILON` | 1e-8 | pruning_policy.pruning_ratio denominator |
| Backoff `step` / `max_steps` | 0.1 / 20 at scaffold call site | resolution_backoff.refine_resolution; function accepts programmatic overrides |
| Backoff `threshold` | Required finite [0,1] | Function API; direct CLI supplies default 0.10 |
| Backoff `fp_cutoff` | Default 0.5; scaffold passes per-cell QC cutoff | partition_conflicts; (0,1] |
| Backoff `min_resolution` / `initial_resolution` | Positive; initial from DQI | Function API; no fixed-gamma backoff |
| Backoff `seed` | Function default 42; scaffold passes get_seed() | Partition optimization; reproducibility caveats above apply |
| Backoff `partitioner` / `logger_obj` | None (real Leiden/default logging) | API injection/diagnostics, not CLI controls |
| Graph pair thresholds | N11>=1, J>=0.08 or low-J directional support>=0.5 | _build_pairwise_correlation_cache_fast |
| Initial graph singleton filters | MCF 0.05 / MCN 5 with quality filtering; 0/0 otherwise | build_scaffold_tree call; not dictionary consensus settings |
| Main graph n_shuffle / seed arguments | 100 / 42 at call site | Deterministic clone weighting; seed argument does not itself guarantee all phases are seeded |
| Builder `export_phylo` / `cv_value` | API false / None; CV entries pass true / candidate | Controls per-CV export; not user CLI flags |
| Final all-zero cleanup | Always applied to both final bundles | phylo_export; no public disable switch |
| `output_stage` / pre-QC snapshot | Determined by builder; pre_fallback on unusable final | Not a user-selected QC bypass |

`phylo/` is cleaned Step 8 output (or marked pre_fallback); `phylo_unpruned/` is the
cleaned Step 7 snapshot of the same selected run. See [output-stage semantics](phylo_output_stages.md).

## Source map and verification scope

- [Pipeline CLI](../cli/main.py) and [configuration template](../config/paths.yaml.template)
- [scRNA entry](../src/run_phylosilid_fullTree_scRNA.py), [scDNA entry](../src/run_phylosilid_fullTree_scDNA.py), [SpaceTracer entry](../src/run_phylosilid_fullTree_spacetracer.py), [binary entry](../src/run_phylosilid_fullTree_binput.py)
- [Scaffold grouping/filtering](../src/scaffold_builder.py), [germline/genotype helpers](../src/germline_filter.py)
- [Placement scoring](../src/mutation_integrator.py), [Step 7/8](../src/full_tree_builder.py)
- [Backoff](../src/resolution_backoff.py), [pruning guard](../src/pruning_policy.py), [exports](../src/phylo_export.py)

The dictionary/CLI catalog is checked against static source inventories. Descriptions
of activity are based on current consumers and call sites, not parameter names.
Defaults are implementation choices, not claims of universal biological optima.
