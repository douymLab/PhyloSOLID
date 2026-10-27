# Visualization and interactive inspection

[Back to README](../README.md) · [Output stages](phylo_output_stages.md)

Commands below run from the repository root.

PhyloSOLID provides a dedicated R package **PhyloSOLIDvis** for generating publication-ready circular (circos-style) phylogenetic tree visualizations.

<div align="center">
  <img src="https://raw.githubusercontent.com/TsingYang1112/PhyloSOLIDvis/master/inst/demo_circos.png" alt="Example circos plot from PhyloSOLIDvis" width="800">
  <br>
  <em>Figure 2: Example circos plot generated from PhyloSOLID output</em>
</div>

### Installation

PhyloSOLIDvis lives in a **separate** conda environment from tree building. `bash install.sh` does **not** install it into `pmg` / `phylosolid_env`.

```bash
# Recommended: dedicated vis env
bash install_vis.sh
conda activate phylosolid_vis
```

If you already have a working visualization env (for example `circos_env`):

```bash
export PHYLOSOLID_VIS_ENV=circos_env
# or
export PHYLOSOLID_VIS_RSCRIPT=/path/to/envs/circos_env/bin/Rscript
# optional: install from a local PhyloSOLIDvis clone
# export PHYLOSOLIDVIS_SRC=/path/to/PhyloSOLIDvis
```

Demos call `scripts/visualization/run_circos.sh` after tree building. They look for `PHYLOSOLID_VIS_RSCRIPT`, then `PHYLOSOLID_VIS_ENV`, then `phylosolid_vis`, then `circos_env`. The tree-building envs `pmg` and `phylosolid_env` are never used for plotting. If no vis R is found, circos is skipped and the demo still finishes.

Load packages inside the visualization environment:

```r
library(ggplot2)
library(PhyloSOLIDvis)
```

### Demo and manual run

Each `demo/*/run_demo.sh` builds the tree, then draws circos into `0*_final_results/circos/` using the vis env.

Use the exported cell annotation at `0*_scaffold_builder/df_celltype.txt` (tab-separated, **no row index**, barcodes already match the tree). If you did not provide `--celltype-file`, that table still exists with `cell_type=default_type`.

To plot an existing `phylo/` directory without re-running the demo:

```bash
bash scripts/visualization/run_circos.sh \
  path/to/05_final_results/phylo \
  path/to/05_final_results/circos \
  path/to/03_scaffold_builder/df_celltype.txt
```

> **⚠️ Important:** `ggplot2` must be loaded **before** `PhyloSOLIDvis` to avoid namespace conflicts.

```r
library(ggplot2)
library(PhyloSOLIDvis)

# Run the complete pipeline
results <- run_all(
  inputpath = "path/to/workdir/sampleid/03_tree_building/05_final_results/phylo/",
  outputpath = "path/to/workdir/sampleid/03_tree_building/05_final_results/circos/",
  annotation_file = "path/to/workdir/sampleid/03_tree_building/03_scaffold_builder/df_celltype.txt"
)
```

### Locating the correct input directory

After running the PhyloSOLID pipeline with `--workdir ./results` and `--sample SAMPLE_ID`, the required files are located in:

```
workdir/
└── SAMPLE_ID/
    └── 03_tree_building/
        ├── 03_scaffold_builder/
        │   └── df_celltype.txt
        └── 05_final_results/
            ├── phylo/
            │   ├── final_cleaned_I_full_withNA3_for_circosPlot.txt
            │   ├── final_cleaned_M_full_basedPivots.filtered_sites_inferred.CFMatrix
            │   ├── final_cleaned_tree_node.txt
            │   ├── df_flipping_count_for_each_mut.txt
            │   └── df_total_flipping_count.txt
            ├── phylo_unpruned/   # Step 7 pre-QC tree; all-zero rows/columns removed
            └── circos/           # written by the demo vis step
```

**Set `inputpath` to this `phylo/` directory.** Binary-matrix mode uses `03_final_results/phylo/` and `01_scaffold_builder/df_celltype.txt`.

### Required Input Files

| File | Description |
|------|-------------|
| `final_cleaned_I_full_withNA3_for_circosPlot.txt` | Input genotype matrix |
| `final_cleaned_M_full_basedPivots.filtered_sites_inferred.CFMatrix` | Conflict-free matrix |
| `df_flipping_count_for_each_mut.txt` | Per-mutation flipping statistics |
| `df_total_flipping_count.txt` | Total flipping statistics |

### Annotation File Format (Optional)

The annotation file (TSV format) should contain a `barcode` column matching cell IDs in the PhyloSOLID output. Optional columns include:

| Column | Description |
|--------|-------------|
| `cluster_info` | Cell cluster assignments |
| `cell_type` | Cell type annotations |
| `sample` | Sample identifiers |
| `tumor_score` | Numeric tumor scores |
| `B_cell_prop` | B cell proportions |

**Note:** `annotation_file` is optional. If not provided, the plot will be generated without annotation layers (tree + heatmap only).

### Output Files

After running `run_all()` or `plot_circos()`, files are written to `0*_final_results/circos/`:

| File | Description |
|------|-------------|
| `No_target.circle_tree_output_as_point.*.svg/pdf` | ★ Main circular plot |
| `heatmap_and_histograms_for_our_tree.pdf` | ★ Mutation heatmap with histograms |
| `legend_components.circos_annotation.svg/pdf` | Circos annotation legend |
| `legend_components.total_flipping_count.svg/pdf` | Flipping count legend |
| `sorted_cf_matrix.txt` | Sorted mutation matrix |
| `ordered_metadata_for_heatmap.txt` | Heatmap ordering data |
| `CNVtree_data_clone_order_tree.rds` | Tree/clone ordering data |

### Usage Examples

**Basic usage with annotation:**
```r
library(ggplot2)
library(PhyloSOLIDvis)

result <- plot_circos(
  inputpath = "results/demo_scrna/03_tree_building/05_final_results/phylo/",
  outputpath = "results/demo_scrna/03_tree_building/05_final_results/circos/",
  annotation_file = "results/demo_scrna/03_tree_building/03_scaffold_builder/df_celltype.txt"
)
```

**With custom parameters:**
```r
library(ggplot2)
library(PhyloSOLIDvis)

result <- plot_circos(
  inputpath = "results/demo_scrna/03_tree_building/05_final_results/phylo/",
  outputpath = "results/demo_scrna/03_tree_building/05_final_results/circos/",
  annotation_file = "results/demo_scrna/03_tree_building/03_scaffold_builder/df_celltype.txt",
  target_mut = "chr11_65426524_T_C",
  tip_label_offset = 8,
  tip_label_size = 3,
  heatmap_width = 0.3,
  flipping_point_size = 1.3,
  plot_height = 12,
  plot_width = 18,
  verbose = TRUE
)
```

**Complete pipeline (matrix ordering + circos plot + heatmap):**
```r
library(ggplot2)
library(PhyloSOLIDvis)

results <- run_all(
  inputpath = "results/demo_scrna/03_tree_building/05_final_results/phylo/",
  outputpath = "results/demo_scrna/03_tree_building/05_final_results/circos/",
  annotation_file = "results/demo_scrna/03_tree_building/03_scaffold_builder/df_celltype.txt",
  target_mut = "chr11_65426524_T_C",
  run_adjusted = TRUE,
  adjusted_tip_label_offset = 6,
  adjusted_flipping_point_size = 1.3,
  verbose = TRUE
)
```

### For Complete Documentation

> **📘 For detailed parameter descriptions, advanced usage, and troubleshooting, please visit the PhyloSOLIDvis GitHub repository:**
>
> **[https://github.com/TsingYang1112/PhyloSOLIDvis](https://github.com/TsingYang1112/PhyloSOLIDvis)**

### Quick Troubleshooting for Visualization

| Issue | Solution |
|-------|----------|
| `could not find function "plot_circos"` | Run `library(ggplot2)` then `library(PhyloSOLIDvis)` |
| `cannot open file '.../phylo/final_cleaned_I_...'` | Verify `inputpath` points to the exact `phylo/` subdirectory |
| Missing `annotation_file` column | Ensure file has a `barcode` column matching cell IDs |
| Plot renders with overlapping labels | Adjust `tip_label_offset` (increase for more space) |
| Network issues during installation | Set CRAN/BioC mirrors and retry |

---

## Interactive Online Inspection Platform

For interactive post-analysis and quality control, you can upload your PhyloSOLIDvis outputs to our dedicated web server:

**Web Server**: [https://phylosolid.westlake.edu.cn](https://phylosolid.westlake.edu.cn)

For access requirements, follow the instructions on the inspection platform.

**Features:**
- Interactive Circos exploration with zoom and cell-level query
- Tree topology evaluation and quality control
- Genotype flip statistics for assessing tree reliability
- Optional upload of IGV snapshots and ANNOVAR annotation results

**Required upload files:**

| File | Location |
|:-----|:---------|
| `ordered_metadata_for_heatmap.txt` | `05_final_results/circos/` |
| `df_muts_corresponding_to_ordered_tiplabels_by_anticlockwise.txt` | `05_final_results/circos/` |
| `No_target.circle_tree_output_as_point.*.svg` | `05_final_results/circos/` |
| `legend_components.circos_annotation.svg` | `05_final_results/circos/` |
| `df_total_flipping_count.txt` | `05_final_results/phylo/` |
| `df_flipping_count_for_each_mut.txt` | `05_final_results/phylo/` |
