# Expected output layout (binary-matrix)

This snapshot uses the **same directory names** as a real run of `../run_demo.sh` (`output/demo_binput/`). Intermediate files are omitted; key finals are kept.

```
demo_binput/
├── 01_scaffold_builder/
│   └── df_celltype.txt          # shipped (no row index)
├── 02_mutation_integrator/      # created by the run; files omitted
└── 03_final_results/
    ├── phylo/                   # shipped: cleaned tree + matrices
    ├── phylo_unpruned/          # created by the run; files omitted
    └── circos/circle_tree_output_as_point.pdf
```

Shipped result: **7 cells × 9 mutations** (12 ghost sites pruned). Full intermediates: `bash ../run_demo.sh`.
