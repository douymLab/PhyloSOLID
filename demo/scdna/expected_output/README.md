# Expected output layout (scDNA)

This snapshot uses the **same directory names** as a real run of `../run_demo.sh` (`output/demo_scdna/`). Intermediate files are omitted; key finals are kept.

```
demo_scdna/
├── pipeline_summary.yaml
├── 01_features/                 # BAM pileup / genotyping; files omitted
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

Shipped result: **22 cells × 78 mutations**. Bulk BAM was not used. Full intermediates: `bash ../run_demo.sh`.
