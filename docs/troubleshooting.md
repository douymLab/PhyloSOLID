# Troubleshooting

[Back to README](../README.md) · [Installation](installation.md) · [Tree-building controls](tree_building.md)

## Environment or executable not found

Activate the tree-building environment and install the package:

```bash
conda activate phylosolid_env
pip install -e .
phylosolid --help
```

Run from the repository root. `python -m cli.main` is an alternative to the installed
command. If conda itself is missing, install a conda distribution first.

## ANNOVAR or references not found

Follow [installation and resource setup](installation.md). Check that paths in
`config/paths.yaml` point to extracted, readable files. For ANNOVAR diagnostics:

```bash
phylosolid check-annovar --config config/paths.yaml
```

Binary-matrix tree building does not require ANNOVAR. For scDNA, also check the
external genotyper path in the configuration template.

## Read length or input mismatch

Set `--read-len` to the read length used for your data and `--cellnum` to the cell
count. Match mutation chromosome names and cell IDs to the BAM/barcode inputs.
Later `--steps` stages require earlier outputs under the same sample/work directory.

## CV search is slow or the expected CSV is missing

The default `auto` search evaluates eight thresholds. To reduce work, use a single
threshold or a shorter list with the [direct tree-building entry point](tree_building.md).
`--max_workers` controls parallel CV workers; the pipeline's `--threads` flag is
not a substitute for it.

Search runs write `cv_search_summary.csv` and `cv_selection_report.txt`.
Single-threshold runs write `cv_threshold_search_results.csv`.
Check the main and per-CV logs if all candidates fail. A `pre_fallback` result means
QC produced no usable final matrix; inspect the pre-QC bundle and QC reports.

## Circos is skipped or fails

Use a separate visualization environment, created with `bash install_vis.sh`.
Set `PHYLOSOLID_VIS_ENV` or `PHYLOSOLID_VIS_RSCRIPT` if automatic detection fails.
The helper skips circos if the visualization environment is unavailable.

Point the plot input to the exact `phylo/` directory and use the exported
`df_celltype.txt` for annotations. See [visualization troubleshooting](visualization.md#quick-troubleshooting-for-visualization).

## Large BAMs

When suitable for your analysis, extract reads overlapping target intervals before
preprocessing. Keep chromosome naming, coordinate conventions, indexes and required
barcode tags consistent with the original inputs.

## Reporting a problem

Include the mode, command, environment, main/per-CV error logs and the relevant
`output_stage.json`. See [Contact](../README.md#contact).
