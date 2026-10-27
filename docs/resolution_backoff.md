# Resolution backoff

[Back to README](../README.md) · [Complete parameter reference](parameters.md) · [Workflow](tree_building.md)

## Algorithm and parameters

scRNA, SpaceTracer and binary-matrix enable automatic resolution backoff with
budget tau=0.10. Direct entries accept `--resolution_backoff_threshold` or
`--disable_resolution_backoff`. A fixed `resolution_of_graph` takes precedence;
scDNA retains fixed gamma=0.1 by default, so its backoff flags are inactive until
its parameter dictionary uses None for automatic resolution.

Start with the existing DQI/sparsity-based auto resolution. On the initial graph,
let n_ig be observed positive calls for cell i in group g, and define:

```text
q_i = 1 - max_g(n_ig) / sum_g(n_ig)
H = fraction of informative cells with q_i >= fp_cutoff
```

Informative cells have at least one positive among graph mutations; missing calls
do not count as positive or negative evidence. The scaffold passes
`fp_ratio_per_cell_cross_all_muts_cutoff` as fp_cutoff (default 0.5). This is a
grouping-conflict proxy, not an estimated doublet fraction or tree accuracy.

If H exceeds tau, lower resolution by 0.1 on the same graph and accept the first
candidate with H<=tau. The configured floor is normally 0.5, with at most 20
decreases. Graph vertices, positive evidence and evaluation denominator remain
fixed. H need not be monotone. If no candidate passes, retain the original auto
partition and report the unmet budget; do not automatically select the lowest gamma.
With no informative cells, stop without treating the budget as satisfied.
Only initial main grouping is retried; later subgroup operations inherit the
selected gamma. This module does not remove cells.

Distinguish tau (grouping-conflict budget), lambda=`fnfp_ratio` (FN/FP cost), and
kappa=`PRUNING_CONFIDENCE_THRESHOLD` (CV pruning guard). They apply at different
stages and have different numerical meanings. The complete defaults, source-level
controls and inactive dictionary fields are in [parameters.md](parameters.md).

## Output and audit

Each enabled scaffold run writes the sample-prefixed files:

- `.resolution_backoff.json`: initial/selected gamma, budget, threshold_met, stop reason and trace.
- `.resolution_backoff_trajectory.tsv`: candidate gamma, group counts and conflict statistics.
- `.resolution_backoff_initial_cells.tsv` and `.resolution_backoff_selected_cells.tsv`: per-cell support/proxy values.

For CV search these live in the corresponding candidate scaffold directory;
selected scaffold artifacts are copied into final scaffold_summary by the pipeline.
Inspect the acceptance reason and conflict numerator/denominator, not only the
selected gamma. Fixed-resolution/disabled runs do not produce a backoff trajectory.

## Evaluation and limitations

Ten percent is an empirical development tolerance, not a theoretical optimum or
a presumed biological doublet rate. Lower H alone does not establish correct
partitioning: one group has H=0 by construction. Main group count is not the final
clone count. Mode interface support is distinct from validation of biological
performance in every mode.

Evaluate with fixed code/inputs/seeds and controlled with/without-backoff runs.
Report independent topology and founder recovery metrics, retained cells/mutations,
full gamma/H/group-count trajectories, budget sensitivity and unsuccessful cases.
Without ground truth, increased retention alone cannot establish a better topology.
Integration regression checks establish implementation compatibility, not biological
accuracy or equivalence to a controlled ablation.

Earlier development notes recorded that a stricter tau=0.03 helped an older NN
simulation but degraded founder recovery or merged branches in some real samples
(including P9). These are development observations, not an independent validation
cohort or proof that 0.10 is universally optimal. Claims of broad stability need a
fixed protocol and independent evaluation; include regressions as well as benefits.

Backoff is a bounded search over downward candidates, not statistical regression,
a calibrated probability model or a guarantee of the globally best resolution.
