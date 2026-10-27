# Expected output layout (binary-matrix)

This snapshot uses the **same directory names** as a real run of `../run_demo.sh` (`output/demo_binput/`). Intermediate files are omitted; key finals are kept.

```
demo_binput/
├── 01_scaffold_builder/
│   └── df_celltype.txt          # shipped (no row index)
├── 02_mutation_integrator/      # created by the run; files omitted
└── 03_final_results/
    ├── phylo/                   # shipped: cleaned tree + matrices
    ├── phylo_unpruned/          # shipped: cleaned Step 7 pre-QC snapshot
    └── circos/circle_tree_output_as_point.pdf
```

Shipped post-QC result: **8 cells × 19 mutations**. Pre-QC snapshot: **13 cells × 21 mutations**. Full intermediates: `bash ../run_demo.sh`.

Refreshed on 2026-10-27 using current source. Both phylo bundles remove all-zero inferred matrix rows/columns and include `output_stage.json`. Upstream feature/tree-input files were reused where applicable; tree building and circos were rerun. See `demo_result_refresh.json` for commands, source/input hashes and dimensions.
