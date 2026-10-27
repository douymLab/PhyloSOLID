# PhyloSOLID demos

Three self-contained demos live here. A live run writes to `output/` (gitignored). GitHub ships a **directory skeleton plus key files** under `expected_output/`, using the same folder names as a real run.

| Mode | Script | Sample ID | Input in repo | Shipped result |
|------|--------|-----------|---------------|----------------|
| scRNA | `scrna/run_demo.sh` | `demo_scrna` | mutations, barcodes, cell types; **BAM from Figshare** | 352 cells × 17 mutations |
| scDNA | `scdna/run_demo.sh` | `demo_scdna` | target-site mini BAMs + mutation list | 22 cells × 78 mutations |
| binary-matrix | `binary-matrix/run_demo.sh` | `demo_binput` | `input/demo_input.tsv` | 8 cells × 19 mutations |

```bash
bash demo/scrna/run_demo.sh
bash demo/scdna/run_demo.sh
bash demo/binary-matrix/run_demo.sh
```

## What is uploaded vs what you regenerate

`expected_output/` keeps the output **framework** so you can see where files land:

- numbered step folders (`01_features`, `02_treeinput`, `03_tree_building`, …)
- `df_celltype.txt` (no row index; for circos)
- cleaned post-QC `phylo/` and Step 7 pre-QC `phylo_unpruned/` trees and matrices
- `output_stage.json` for each bundle and `demo_result_refresh.json` provenance
- one circos PDF
- `cv_search_summary.csv` / `cv_selection_report.txt` (scRNA and scDNA)

Omitted (placeholder `README.md` only): feature tables, tree-input matrices, per-CV folders, logs, extra circos SVG/PDF. Re-run `run_demo.sh` to fill those in under `output/`.

Tree building uses `phylosolid_env` / `pmg`. Circos uses a **separate** env (`bash install_vis.sh`, or `PHYLOSOLID_VIS_ENV=circos_env`). Demos skip circos if that env is missing.

## scDNA notes

scDNA is BAM-to-tree: pileup at mutation sites, unphased genotyping, tree-input matrices, then tree building. Bulk BAM is **not** used in the demo.

Set `scdna.genotyper_dir` in `config/paths.yaml` (copy from `config/paths.yaml.template`) or export `PHYLOSOLID_GENOTYPER_DIR`. Mini BAM indexes are created by `run_demo.sh` if missing.

## Latest result refresh (2026-10-27)

Tree building and circos were rerun with the current code, reusing existing upstream feature/tree-input files for scRNA and scDNA. `output/` contains full regenerated tree-building results; `expected_output/` includes both cleaned bundles and refreshed key files.

| Mode | Step 7 pre-QC cells × mutations | Step 8 post-QC cells × mutations |
|---|---|---|
| scrna | 433 × 17 | 352 × 17 |
| scdna | 22 × 78 | 22 × 78 |
| binary-matrix | 13 × 21 | 8 × 19 |
