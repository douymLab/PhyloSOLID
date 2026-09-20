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
        ├── phylo_unpruned/      # files omitted
        └── circos/circle_tree_output_as_point.pdf
```

Shipped result: **352 cells × 17 mutations**. Full intermediates: `bash ../run_demo.sh` (scRNA BAM from Figshare).
