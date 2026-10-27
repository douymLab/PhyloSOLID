# Usage and input formats

[Back to README](../README.md) · [Installation](installation.md) · [CV and resolution](tree_building.md)

Run commands from the repository root after activating the tree-building environment.
Global options such as `--workdir` and `--config` go **before** the mode name.
Use `python -m cli.main` in place of `phylosolid` if the executable is not on PATH.

## Binary matrix

```bash
phylosolid binary-matrix --sampleid SAMPLE_ID --inputfile matrix.tsv --outputpath results
```

The tab-separated matrix has cell IDs in its first column and mutation IDs in its header.
Values are `0` (absent), `1` (present), or `NA` (missing); do not use the export sentinel `3` as an input genotype.

```text
cell_id	chr1_1000_A_G	chr2_2000_C_T	chr3_3000_G_A
Cell_1	1	0	1
Cell_2	0	1	NA
Cell_3	1	1	0
```

This mode skips classifier, germline filtering and coverage-based CV search.
For seed, cell-type and resolution controls, use the [direct tree-building entry point](tree_building.md).

## scRNA

```bash
phylosolid --config config/paths.yaml --workdir results scrna \
  --sample SAMPLE_ID \
  --mutation-list mutations.txt \
  --bam sample.bam \
  --barcode barcodes.txt \
  --celltype-file celltype.tsv \
  --read-len 120 --cellnum 3319 --threads 4
```

Set read length and cell count for your own data; the values above match the demo.
Cell-type annotation is optional. For scripted runs, provide a ready-made
`--celltype-file`; the `--metadata` route depends on the annotation helper setup.
`--running-type` defaults to `benchmark`.

Run only selected pipeline stages with `--steps feature_extraction`,
`--steps tree_input`, or `--steps tree_building`. Later stages require the earlier
outputs to exist under the same sample/work directory. `--parallel` allows feature
extraction and tree-input generation to run concurrently.
The CLI stages are coarser than the internal Step 1–9 tree-building steps.

For a complete worked example, see [demo/scrna/run_demo.sh](../demo/scrna/run_demo.sh).

## scDNA

```bash
phylosolid --workdir results scdna \
  --sample SAMPLE_ID \
  --mutation-list mutations.txt \
  --bam-dir /path/to/sc_bams \
  --sample-list sample_list.txt --threads 4
```

BAMs are named `{cell}.bam`. Alternatively, replace `--bam-dir` with `--bam-list bam_paths.txt`,
a file containing one BAM path per line. The sample list contains one cell ID per line;
if omitted, BAM file stems supply the IDs. Bulk BAM is optional (`--bulk-bam bulk.bam`).

Configure `scdna.genotyper_dir` or `PHYLOSOLID_GENOTYPER_DIR` to point to the directory
containing `run_unphased_prod_whether_bulk.py`.
Use `--reference` and `--mappability` for optional classifier BAM features;
see [installation](installation.md) and the [configuration template](../config/paths.yaml.template).
The same `--steps` names as scRNA are available. See [demo/scdna/run_demo.sh](../demo/scdna/run_demo.sh).

## SpaceTracer

```bash
phylosolid --workdir results spacetracer \
  --sample SAMPLE_ID \
  --mutation-list mutations.txt \
  --bam sample.bam \
  --barcode barcodes.txt \
  --read-len 120 --cellnum 3319 --threads 4
```

This route skips feature extraction and the mutation classifier, and generates tree-input
matrices before tree building. Use an optional `--celltype-file` for annotation; `--metadata` invokes the annotation helper.
Its available stages are `tree_input` and `tree_building`.

## Input identifiers and annotation

Mutation lists contain one `chromosome_position_reference_alt` identifier per line:

```text
chr1_1000_A_G
chr1_2000_C_T
chr2_3000_G_A
```

Match chromosome naming and coordinate conventions to your BAM/reference inputs.
Barcode files contain one cell barcode per line, matching the BAM cell identifiers.
scRNA/SpaceTracer use one BAM for an individual; merge replicates beforehand if needed.
scDNA uses one BAM per cell.

For a precomputed cell-type table, use tab-separated `barcode` and `cell_type` columns.
See the [scRNA demo annotation](../demo/scrna/input/celltype.txt).
Optional plotting fields include `cluster_info`, `sample`, `tumor_score`, and `B_cell_prop`.
For plotting, prefer the exported `df_celltype.txt`, whose barcodes are aligned to the
pipeline output and which has no extra pandas row-index column.

Check available options with:

```bash
phylosolid scrna --help
phylosolid scdna --help
phylosolid spacetracer --help
phylosolid binary-matrix --help
```

For the complete CLI, mode-dictionary and configuration catalog, see [parameters.md](parameters.md).
