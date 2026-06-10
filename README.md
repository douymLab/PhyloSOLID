# PhyloSOLID

Tree building from single-cell sequencing data (scRNA-seq and scDNA-seq)


> **Important Note for Users:**
> 
> PhyloSOLID is under active development, especially during its preprint stage. The software and its associated resources are continuously being updated and improved.
> 
> Our current priority is optimizing the runtime efficiency of the core tree-building pipeline, particularly for large datasets with high mutation counts. We are actively addressing performance bottlenecks while maintaining full model accuracy and all constraints.
> 
> If you encounter any issues, have feature requests, or need guidance, **please do not hesitate to contact us**. We are committed to providing timely support and would greatly appreciate your feedback to make PhyloSOLID better for the entire community. You can reach us at the emails provided in the Contact section.


## Overview

PhyloSOLID is a comprehensive pipeline for building phylogenetic trees from single-cell sequencing data. It supports both scRNA-seq and scDNA-seq modes, with features for mutation filtering, tree construction, and visualization.


> **System Requirements:**
> PhyloSOLID is primarily developed and tested on Linux. We recommend running it on Linux or macOS systems. Windows users can use WSL (Windows Subsystem for Linux).


## Release Notes

*PhyloSOLID is actively maintained. Check the Release Notes below and GitHub Releases for the latest updates.*

- 2026/02/27: Version 1.0.0 (**changelog**)  
  This is the first public release of PhyloSOLID, rebuilt from the earlier PhyloMosaicGenie prototype.  
  It provides a more streamlined genotyping and phylogenetic analysis workflow, with improved code structure and reproducibility.


## Installation

### Prerequisites

- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Mambaforge](https://github.com/conda-forge/miniforge#mambaforge) (required)
- samtools, bedtools
- ANNOVAR (for variant annotation)

### Step 1: Clone the repository

```bash
git clone https://github.com/douymLab/PhyloSOLID.git
cd PhyloSOLID
```

### Step 2: Install ANNOVAR

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

### Step 3: Download and setup resource files

All required resource files for running PhyloSOLID with the hg38 genome build are available on Figshare. **Due to their large total size (~8.85 GB), these files are not included in the GitHub repository and must be downloaded separately.**

**Figshare URL**: [https://figshare.com/s/2aee3e878688722e9c6f](https://figshare.com/s/2aee3e878688722e9c6f)

The Figshare repository contains the following six compressed (`.gz`) files:

| File | Description | Size |
|:-----|:------------|------:|
| `genome.fa.gz` | hg38 reference genome | 840.4 MB |
| `wgEncodeGencodeExonSupportV44.sort.bed.gz` | Exon coordinates for variant annotation | 51.07 MB |
| `k24.umap.bedgraph.gz` | Mappability track for filtering (k-mer 24) | 2.57 GB |
| `k100.umap.bedgraph.gz` | Mappability track for filtering (k-mer 100) | 718.29 MB |
| `hg38_gnomad312_genome_only_af_all.txt.gz` | Population allele frequencies from gnomAD v3.1.2 | 4.65 GB |
| `COMBINED_RADAR_REDIprotal_DARNED_hg38_all_sites.bed.gz` | Curated RNA editing sites | 58.26 MB |

**Download and Setup Instructions:**

1. Download all six `.gz` files from Figshare using the URL above
2. Place the downloaded files in a dedicated directory on your system (e.g., `/path/to/resource/`)
3. Extract each file:

```bash
# Navigate to your resource directory
cd /path/to/resource/

# Extract all .gz files
for file in *.gz; do
    gunzip "$file"
done
```

### Step 4: Configure paths.yaml

Copy the template configuration file and update it with your paths:

```bash
cp config/paths.yaml.template config/paths.yaml
```

Edit `config/paths.yaml` to set the correct paths for:
- Conda environment (Python interpreter and environment path)
- ANNOVAR installation directory
- Reference files (genome FASTA, GFF3, etc.)
- Paths to the extracted resource files

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
  gff3_file: "/path/to/resource/wgEncodeGencodeExonSupportV44.sort.bed"
  mappability_file: "/path/to/resource/k24.umap.bedgraph"
  gnomad_file: "/path/to/resource/hg38_gnomad312_genome_only_af_all.txt"
  rna_editing_file: "/path/to/resource/COMBINED_RADAR_REDIprotal_DARNED_hg38_all_sites.bed"
```

### Step 5: Run the environment installation script

```bash
bash install.sh
```

This script will:
- Create a conda environment with all Python/R dependencies
- Install the converTree R package from GitHub

### Step 6: Install PhyloSOLID

```bash
conda activate phylosolid_env
pip install -e .
```

### Step 7: Verify installation

```bash
# Check if ANNOVAR is properly configured
phylosolid check-annovar --config config/paths.yaml

# Run the test pipeline
cd demo
./run_demo.sh
```

## Usage

### Binary matrix mode (Direct tree building)

For users who already have a binary mutation matrix, PhyloSOLID provides a direct mode that skips feature extraction and tree input generation:

```bash
# Basic usage
python -m cli.main binary-matrix --sampleid SAMPLE_ID -inputfile matrix.txt -outputpath output_dir

# Or using the installed command
phylosolid binary-matrix --sampleid SAMPLE_ID -inputfile matrix.txt -outputpath output_dir
```

**Input format (tab-separated):**

Rows: cells (first column contains cell barcodes/IDs)
Columns: mutations (column headers are mutation IDs)
Values: 1 (present), 0 (absent), NA (missing data)

Example (matrix.txt):

| cell_id | chr1_1000_A_G | chr2_2000_C_T | chr3_3000_G_A |
|:--------|:--------------|:--------------|--------------:|
| Cell_1 | 1 | 0 | 1 |
| Cell_2 | 0 | 1 | NA |
| Cell_3 | 1 | 1 | 0 |

**Output:**
```
output_dir/
├── phylo/
│   ├── final_cleaned_M_scaffold_basedPivots.filtered_sites_inferred.CFMatrix
│   └── ...
├── processing/
│   ├── df_celltype.txt
│   └── ...
└── ...
```

### scRNA-seq mode

**Basic usage:**
```bash
phylosolid --workdir ./results scrna \
    --sample SAMPLE_ID \
    --mutation-list mutations.txt \
    --bam sample.bam \
    --barcode barcodes.txt
```

**With all options:**
```bash
phylosolid --workdir ./results scrna \
    --sample SAMPLE_ID \
    --mutation-list mutations.txt \
    --bam sample.bam \
    --barcode barcodes.txt \
    --metadata metadata.txt \
    --read-len 100 \
    --cellnum 155 \
    --threads 8 \
    --workdir ./results \
    --config config/paths.yaml
```

### Running specific steps

```bash
# Run only feature extraction
phylosolid --workdir ./results scrna --sample SAMPLE_ID ... --steps feature_extraction

# Run only tree input generation
phylosolid --workdir ./results scrna --sample SAMPLE_ID ... --steps tree_input

# Run only tree building
phylosolid --workdir ./results scrna --sample SAMPLE_ID ... --steps tree_building
```

### Parallel execution

```bash
# Run feature extraction and tree input in parallel
phylosolid --workdir ./results scrna --sample SAMPLE_ID ... --parallel
```

### scDNA-seq mode

```bash
phylosolid --workdir ./results scdna \
    --sample SAMPLE_ID \
    --mutation-list mutations.txt \
    --bam sample.bam \
    --barcode barcodes.txt
```

### SpaceTracer mode (Skip feature extraction and classifier)

SpaceTracer mode is designed for users who want to directly build trees without running feature extraction and the classifier step. This mode:

- Skips the feature extraction step entirely
- Bypasses the mutation classifier
- Uses a dedicated tree building script
- Accepts the same input parameters as scRNA mode

**Basic usage:**
```bash
phylosolid --workdir ./results spacetracer \
    --sample SAMPLE_ID \
    --mutation-list mutations.txt \
    --bam sample.bam \
    --barcode barcodes.txt
```

**With all options:**
```bash
phylosolid --workdir ./results spacetracer \
    --sample SAMPLE_ID \
    --mutation-list mutations.txt \
    --bam sample.bam \
    --barcode barcodes.txt \
    --metadata metadata.txt \
    --read-len 100 \
    --cellnum 155 \
    --threads 8 \
    --config config/paths.yaml
```

## Input File Formats

### Mutation list (mutations.txt)
```
chr1_1000_A_G
chr1_2000_C_T
chr2_3000_G_A
chr3_4000_T_C
```
Format: `chromosome_position_reference_alt_gene`

### Barcode file (barcodes.txt)
```
AAACCTGAGAAACCAT-1
AAACCTGAGAAACCGG-1
AAACCTGAGAAACCTA-1
```

### Metadata file (metadata.txt)
```
cell_barcode    cell_type
AAACCTGAGAAACCAT-1  CD8+T
AAACCTGAGAAACCGG-1  CD4+T
AAACCTGAGAAACCTA-1  Monocyte
```

## Output Structure

After running, results are organized as:

```
workdir/
└── SAMPLE_ID/
    ├── 01_features/                # Feature extraction results
    │   ├── depth_in_spots/
    │   ├── SAMPLE_ID.benchmark_patched.feature.txt
    │   └── ...
    ├── 02_treeinput/               # Tree input files
    │   ├── treeinput/
    │   │   ├── treeinput_spot_c_155.csv
    │   │   ├── treeinput_scid_barcode.txt
    │   │   └── features_file.txt
    │   └── data/                   # Preprocessed data
    ├── 03_tree_building/           # Tree building results
    │   └── results/
    │       └── PhyloSOLID/
    │           ├── cell_by_mut.CFMatrix
    │           └── celltree.newick
    └── pipeline_summary.yaml       # Pipeline execution summary
```

## Demo Example

A complete demo with test data is available in the `demo/` directory:

```bash
cd demo
./run_demo.sh
```

### Run individual steps with demo data

**1. Feature extraction**

```bash
phylosolid --workdir demo/expected_output scrna \
    --sample Org4S15D63 \
    --mutation-list demo/input/Org4S15D63/02_identifier/identifier.txt \
    --bam demo/input/Org4S15D63/01_rawdata/Org4S15D63.bam \
    --barcode demo/input/Org4S15D63/01_rawdata/Org4S15D63_CB.txt \
    --threads 4 \
    --read-len 100 \
    --steps feature_extraction
```

**2. Tree input generation**

```bash
phylosolid --workdir demo/expected_output scrna \
    --sample Org4S15D63 \
    --mutation-list demo/input/Org4S15D63/02_identifier/identifier.txt \
    --bam demo/input/Org4S15D63/01_rawdata/Org4S15D63.bam \
    --barcode demo/input/Org4S15D63/01_rawdata/Org4S15D63_CB.txt \
    --cellnum 155 \
    --steps tree_input
```

**3. Tree building**

```bash
phylosolid --workdir demo/expected_output scrna \
    --sample Org4S15D63 \
    --mutation-list demo/input/Org4S15D63/02_identifier/identifier.txt \
    --bam demo/input/Org4S15D63/01_rawdata/Org4S15D63.bam \
    --barcode demo/input/Org4S15D63/01_rawdata/Org4S15D63_CB.txt \
    --celltype-file demo/input/Org4S15D63/03_celltype/celltype_file_for_Org4S15D63.txt \
    --steps tree_building
```

**Complete pipeline with demo data**

```bash
phylosolid --workdir demo/test_output scrna \
    --sample Org4S15D63 \
    --mutation-list demo/input/Org4S15D63/02_identifier/identifier.txt \
    --bam demo/input/Org4S15D63/01_rawdata/Org4S15D63.bam \
    --barcode demo/input/Org4S15D63/01_rawdata/Org4S15D63_CB.txt \
    --threads 4 \
    --read-len 100 \
    --cellnum 155
```

## Visualization

PhyloSOLID provides a dedicated R package **PhyloSOLIDvis** for generating publication-ready circular (circos-style) phylogenetic tree visualizations based on the pipeline outputs.

### Installation

```r
remotes::install_github("TsingYang1112/PhyloSOLIDvis", dependencies = TRUE)
```

### Usage

After running the PhyloSOLID pipeline, use the following command to generate the circos plot:

```r
library(PhyloSOLIDvis)

plot_circos(
  inputpath = "path/to/PhyloSOLID/output/",
  outputpath = "path/to/figures/",
  annotation_file = "path/to/annotations.txt"
)
```

The function will create:
- Main circular tree plot (SVG/PDF)
- Separate legend files for annotations and flipping counts
- Sorted mutation matrix and heatmap metadata

### Interactive inspection

After generating the circos plot and its associated output files, you can upload the results to the **interactive inspection platform** for manual review and quality control:

👉 **http://10.28.0.22:32300/home**

Upload the following files to the platform:
- Main circos plot (SVG/PDF)
- Sorted CF matrix (`sorted_cf_matrix.txt`)
- Heatmap metadata (`ordered_metadata_for_heatmap.txt`)
- Tree structure data (`CNVtree_data_clone_order_tree.rds`)

For detailed documentation of the visualization package, visit:  
[https://github.com/TsingYang1112/PhyloSOLIDvis](https://github.com/TsingYang1112/PhyloSOLIDvis)

## Troubleshooting

### Conda not found

If you see "conda: command not found":

1. Download and install Miniconda: https://docs.conda.io/en/latest/miniconda.html
2. Or install Mambaforge: https://github.com/conda-forge/miniforge#mambaforge
3. Restart your terminal after installation

### ANNOVAR not found

If you see "ANNOVAR not found" error:

1. Check that ANNOVAR is installed
2. Verify the paths in config/paths.yaml
3. Run `phylosolid check-annovar` to diagnose

### Reference files not found

If you see errors about missing reference files:

1. Download the resource files package from Figshare using the URL above
2. Ensure files are extracted to the correct location
3. Verify paths in config/paths.yaml point to the extracted files
4. Check file permissions (ensure files are readable)

### Read length mismatch

If you see warnings about read length:

1. Check your BAM file's actual read length:
```bash
samtools view your.bam | head -n1 | awk '{print length($10)}'
```
2. Use the `--read-len` parameter to specify the correct value

### Permission denied errors

Ensure scripts are executable:
```bash
chmod +x scripts/scrna/**/*.sh
chmod +x scripts/scrna/**/*.py
chmod +x scripts/scrna/**/*.R
```

## Citation

If you use PhyloSOLID in your research, please cite:

1. Yang, Q. et al. PhyloSOLID: Robust phylogeny reconstruction from single-cell data despite inherent error and sparsity. (2026) doi:10.64898/2026.02.04.703905.

## License

PhyloSOLID is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

### Code Availability

The source code, documentation and examples are available on GitHub at [https://github.com/douymLab/PhyloSOLID](https://github.com/douymLab/PhyloSOLID).

**Note**: The PhyloSOLID source code is currently being prepared for public release. 
It will be made publicly available on GitHub by **March 1, 2026**, or upon the manuscript's 
formal acceptance, whichever comes first. Until then, the code is available for review 
purposes upon request.

### Third-party Dependencies

PhyloSOLID uses several third-party tools and libraries:
- **ANNOVAR**: Users need to install ANNOVAR separately (free for academic use, registration required)
- **R packages**: Various R packages under GPL/MIT licenses
- **Python packages**: See environment.yml for details

Please respect the licenses of these dependencies when using PhyloSOLID.

## Contact

For questions and support, please contact:  
Qing Yang: yangqing@westlake.edu.cn  
Yanmei Dou: yanmeidou@westlake.edu.cn
