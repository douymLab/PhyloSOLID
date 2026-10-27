# Expected output layout (scRNA)

This snapshot uses the **same directory names** as a real run of `../run_demo.sh` (`output/demo_scrna/`). Intermediate files are omitted; key finals are kept.

```
demo_scrna/
├── pipeline_summary.yaml
├── 01_features/                 # feature tables; files omitted
├── 02_treeinput/                # tree-input matrices; files omitted
└── 03_tree_building/
    ├── cv_search_summary.csv    # shipped
    ├── cv_selection_report.txt  # shipped
    ├── 01_classifier_filter/    # files omitted
    ├── 02_germline_filter/      # files omitted
    ├── 03_scaffold_builder/
    │   └── df_celltype.txt      # shipped (no row index)
    ├── 04_mutation_integrator/  # CV folders omitted
    └── 05_final_results/
        ├── phylo/               # shipped: cleaned tree + matrices
        ├── phylo_unpruned/      # shipped: cleaned Step 7 pre-QC snapshot
        └── circos/circle_tree_output_as_point.pdf
```

Shipped post-QC result: **352 cells × 17 mutations**. Pre-QC snapshot: **433 cells × 17 mutations**. Full intermediates: `bash ../run_demo.sh` (scRNA BAM from Figshare).

Refreshed on 2026-10-27 using current source. Both phylo bundles remove all-zero inferred matrix rows/columns and include `output_stage.json`. Upstream feature/tree-input files were reused where applicable; tree building and circos were rerun. See `demo_result_refresh.json` for commands, source/input hashes and dimensions.
