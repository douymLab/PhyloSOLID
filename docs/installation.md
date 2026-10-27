# Installation and resources

[Back to README](../README.md) · [Usage](usage.md)

## Prerequisites

PhyloSOLID is developed and tested primarily on Linux. Windows users can use WSL. Start with the binary-matrix demo to check tree building before configuring BAM preprocessing.

- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Mambaforge](https://github.com/conda-forge/miniforge#mambaforge) (required)
- samtools, bedtools
- ANNOVAR (for variant annotation)

## Step 1: Clone the repository

```bash
git clone https://github.com/douymLab/PhyloSOLID.git
cd PhyloSOLID
```

## Step 2: Install ANNOVAR

PhyloSOLID uses ANNOVAR for variant annotation. You need to install it separately:

1. Register and download ANNOVAR from the ANNOVAR Download Page (free for academic use).

2. Install ANNOVAR:
```bash
tar -xzvf annovar.latest.tar.gz -C /path/to/software/
cd /path/to/software/annovar
```

3. Download required databases (hg38 build):
```bash
# Create humandb directory
mkdir -p /path/to/software/annovar/humandb

# Download refGene database (required)
./annotate_variation.pl -downdb -buildver hg38 -webfrom annovar refGene humandb/

# Download additional databases (recommended)
./annotate_variation.pl -downdb -buildver hg38 -webfrom annovar cytoBand humandb/
./annotate_variation.pl -downdb -buildver hg38 -webfrom annovar snp138 humandb/
```

## Step 3: Download and setup resource files

The hg38 resources used by the scRNA preprocessing workflow are distributed separately on Figshare. Binary-matrix tree building does not require this resource bundle or ANNOVAR. scDNA BAM-to-tree requires the external unphased genotyper; reference and mappability files for optional classifier BAM features are configured separately under `scdna` in the template.

**Figshare URL**: [https://figshare.com/s/2aee3e878688722e9c6f](https://figshare.com/s/2aee3e878688722e9c6f)

The Figshare repository contains the following files:

**Reference and resource files (required for pipeline execution):**

| File | Description | Size |
|:-----|:------------|------:|
| `genome.fa.gz` | hg38 reference genome | 840.4 MB |
| `genome.fa.fai.gz` | FASTA index file for genome.fa | 2.13 KB |
| `wgEncodeGencodeExonSupportV44.sort.bed.gz` | Exon coordinates for variant annotation | 51.07 MB |
| `k24.umap.bedgraph.gz` | Mappability track for filtering (k-mer 24) | 2.57 GB |
| `k100.umap.bedgraph.gz` | Mappability track for filtering (k-mer 100) | 718.29 MB |
| `hg38_gnomad312_genome_only_af_all.txt.gz` | Population allele frequencies from gnomAD v3.1.2 | 4.65 GB |
| `COMBINED_RADAR_REDIprotal_DARNED_hg38_all_sites.bed.gz` | Curated RNA editing sites | 58.26 MB |

**Demo data files (optional, for testing the pipeline):**

| File | Description | Size |
|:-----|:------------|------:|
| `demo_scrna.bam` | Demo BAM file for testing the pipeline | 52.06 MB |
| `demo_scrna.bam.bai` | BAM index file | 314.56 KB |

**Download and Setup Instructions:**

1. Download the files needed for your mode from Figshare using the URL above
2. Place the downloaded files in a dedicated directory on your system (e.g., `/path/to/resource/`)
3. Extract the reference and resource files:

```bash
# Navigate to your resource directory
cd /path/to/resource/

# Extract reference and resource files
gunzip genome.fa.gz
gunzip genome.fa.fai.gz
gunzip wgEncodeGencodeExonSupportV44.sort.bed.gz
gunzip k24.umap.bedgraph.gz
gunzip k100.umap.bedgraph.gz
gunzip hg38_gnomad312_genome_only_af_all.txt.gz
gunzip COMBINED_RADAR_REDIprotal_DARNED_hg38_all_sites.bed.gz
```

**Important:**
- The BAM file (`demo_scrna.bam`) and its index (`demo_scrna.bam.bai`) are **demo data** for testing the pipeline. They are not required for running PhyloSOLID on your own data.
- If you want to test the pipeline with the provided demo data, place the BAM files in the `demo/scrna/input/` directory.

## Step 4: Configure paths.yaml

Copy the template configuration file and update it with your paths:

```bash
cp config/paths.yaml.template config/paths.yaml
```

Edit `config/paths.yaml` to set the correct paths for:
- Conda environment (Python interpreter and environment path)
- ANNOVAR installation directory
- Reference files (genome FASTA, GFF3, etc.)
- Paths to the extracted resource files
- **scDNA:** `scdna.genotyper_dir` (unphased genotyper). Bulk BAM is not required.

Example `config/paths.yaml` configuration:

```yaml
# Conda environment (required for Python scripts)
conda:
  python: "/path/to/phylosolid_env/bin/python"
  env_path: "/path/to/phylosolid_env"

# ANNOVAR configuration
annovar:
  script_dir: "/path/to/software/annovar"
  humandb: "/path/to/software/annovar/humandb"
  build: "hg38"

# Reference files (paths to extracted files)
reference:
  genome_fasta: "/path/to/resource/genome.fa"
  genome_fasta_index: "/path/to/resource/genome.fa.fai"
  gff3_file: "/path/to/resource/wgEncodeGencodeExonSupportV44.sort.bed"
  mappability_file: "/path/to/resource/k24.umap.bedgraph"
  gnomad_file: "/path/to/resource/hg38_gnomad312_genome_only_af_all.txt"
  rna_editing_file: "/path/to/resource/COMBINED_RADAR_REDIprotal_DARNED_hg38_all_sites.bed"

# scDNA end-to-end (bulk BAM optional)
scdna:
  genotyper_dir: "/path/to/scMosaicCaller/scMosaicCaller/genotyper"
```

## Step 5: Run the environment installation script

```bash
bash install.sh
```

This script will:
- Create a conda environment with all Python/R dependencies for **tree building**
- Install the converTree R package from GitHub (needed for Newick conversion)

Circos plots use a **second** conda environment. After `install.sh`:

```bash
bash install_vis.sh
# or, if you already have one: export PHYLOSOLID_VIS_ENV=circos_env
```

## Step 6: Install PhyloSOLID

```bash
conda activate phylosolid_env
pip install -e .
```

## Step 7: Verify installation

```bash
# Check if ANNOVAR is properly configured
phylosolid check-annovar --config config/paths.yaml

# Run a mode-specific demo (tree building, then circos if the vis env is available)
bash demo/binary-matrix/run_demo.sh
# bash demo/scdna/run_demo.sh
# bash demo/scrna/run_demo.sh
```

Paths in commands above are relative to the repository root. See the [configuration template](../config/paths.yaml.template), [environment definition](../environment.yml) and [demo input requirements](../demo/README.md).
