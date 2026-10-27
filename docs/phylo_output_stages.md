# Final phylo output stages

Each mode writes two bundles under its existing final-results directory:

- `phylo/`: Step 8 output, or the Step 7 snapshot when marked `pre_fallback`.
- `phylo_unpruned/`: Step 7 snapshot before Step 8 QC, from the same selected CV.

Both expand compound mutation columns, restore ROOT-assigned mutations, remove
all-zero inferred matrix columns and rows, and collapse corresponding tree nodes.
Observed matrices are aligned to each bundle separately; missing observations
raise an error rather than silently dropping snapshot cells or mutations.
The existing filenames inside the two bundles are retained. Each bundle adds
`output_stage.json` with its source step, dimensions, removed ghost mutations and
Omega. `README_phylo.txt` explains the distinction in the final-results directory.

Step 7 already performs discordance-guided prune-and-regraft refinement.
"Unpruned" does not mean raw input or a tree before all refinement. Step 8 may
remove cells/artifact mutations, reattach mutations and resolve conflicts, so
the bundles may differ in topology as well as their cell and mutation sets.
Export cleanup uses inferred matrix M, not observed matrix I.

ROOT-assigned mutations now use the same preparation in pre-QC/CV scoring and
final exports. Scores and selected CVs can therefore differ from older runs
that omitted these mutations during scoring. The weight and selection rule
are unchanged. Existing result directories are not migrated automatically;
regenerate results to obtain both snapshots under the new semantics.

Final-results parent directories remain mode-specific: scRNA/scDNA use
`05_final_results`, spacetracer uses `04_final_results`, and binary-input uses
`03_final_results`. Both phylo bundles are emitted only for the selected run.
