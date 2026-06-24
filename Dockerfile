# Use Mambaforge with specific version for reproducibility
FROM condaforge/mambaforge:24.3.0-0

# Set environment variables
ENV LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    PATH=/opt/conda/envs/phylosolid_env/bin:$PATH \
    CONDA_DEFAULT_ENV=phylosolid_env \
    DEBIAN_FRONTEND=noninteractive \
    # R related environment variables
    R_HOME=/opt/conda/envs/phylosolid_env/lib/R \
    R_LIBS_USER=/opt/conda/envs/phylosolid_env/lib/R/library \
    # Increase R download timeout for large packages
    R_DOWNLOAD_TIMEOUT=600

# Set working directory
WORKDIR /opt

# Install system dependencies (required for R packages and bioinformatics tools)
RUN apt-get update && apt-get install -y \
    build-essential \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libgit2-dev \
    libmagick++-dev \
    libudunits2-dev \
    libgdal-dev \
    cmake \
    wget \
    git \
    vim \
    libglpk-dev \
    libnlopt-dev \
    libgmp-dev \
    libmpfr-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Copy environment.yml file
COPY environment.yml /tmp/environment.yml

# Create conda environment and install all dependencies using mamba
# Set conda channels to use Tsinghua mirror for faster downloads
RUN conda config --set show_channel_urls yes && \
    # Configure conda to use Tsinghua mirror
    conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/ && \
    conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/ && \
    conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/ && \
    conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/ && \
    mamba env create -f /tmp/environment.yml && \
    # Clean conda cache to reduce image size
    mamba clean -afy && \
    # Remove temporary files
    rm -rf /tmp/* && \
    # Verify key Python packages are installed successfully
    conda run -n phylosolid_env python -c "import sys, numpy, pandas, scanpy, torch; print(f'Python: {sys.version}\nNumPy: {numpy.__version__}\nPandas: {pandas.__version__}\nScanpy: {scanpy.__version__}\nPyTorch: {torch.__version__}')"

# Clone PhyloSOLID repository - using v1.0.0 tag for stability
RUN git clone --branch v1.0.0 --depth 1 https://github.com/douymLab/PhyloSOLID.git /opt/PhyloSOLID

# Set working directory to PhyloSOLID project
WORKDIR /opt/PhyloSOLID

# Run PhyloSOLID installation script
RUN bash install.sh && \
    # Install PhyloSOLID package in development mode using pip
    conda run -n phylosolid_env pip install -e . && \
    # Verify phylosolid command is available
    conda run -n phylosolid_env phylosolid --help

# ============ Install PhyloSOLIDvis R package ============
# Set R CRAN and Bioconductor mirrors (using Tsinghua mirror for faster access)
RUN conda run -n phylosolid_env R -e "\
    # Set CRAN and Bioconductor mirrors using Tsinghua mirror\
    options(repos = c(CRAN = 'https://mirrors.tuna.tsinghua.edu.cn/CRAN/'));\
    options(BioC_mirror = 'https://mirrors.tuna.tsinghua.edu.cn/bioconductor');\
    # Increase timeout for large package downloads\
    options(timeout = 600);\
    # Install BiocManager if not already installed\
    if (!requireNamespace('BiocManager', quietly = TRUE)) {\
        install.packages('BiocManager', quiet = TRUE);\
    }\
    # Install remotes if not already installed\
    if (!requireNamespace('remotes', quietly = TRUE)) {\
        install.packages('remotes', quiet = TRUE);\
    }\
    # Install converTree (a dependency of PhyloSOLIDvis)\
    remotes::install_github('xiayh17/converTree', quiet = TRUE);\
    # Install PhyloSOLIDvis itself with all dependencies\
    remotes::install_github('TsingYang1112/PhyloSOLIDvis', dependencies = TRUE, quiet = TRUE);\
    " \
    # Verify PhyloSOLIDvis installation
    && conda run -n phylosolid_env R -e "\
    library(ggplot2);\
    library(PhyloSOLIDvis);\
    print('PhyloSOLIDvis installed successfully!');\
    "

# Configure conda environment to activate automatically at container startup
RUN echo "conda activate phylosolid_env" >> ~/.bashrc && \
    echo "cd /workspace" >> ~/.bashrc && \
    # Set R library path
    echo "export R_LIBS_USER=/opt/conda/envs/phylosolid_env/lib/R/library" >> ~/.bashrc

# Create mount point directories
RUN mkdir -p /workspace /data /resources /output

# Set default working directory
WORKDIR /workspace

# Set entrypoint script to ensure conda environment is always activated
COPY --chmod=755 <<-"EOF" /usr/local/bin/entrypoint.sh
#!/bin/bash
source /opt/conda/etc/profile.d/conda.sh
conda activate phylosolid_env
export R_LIBS_USER=/opt/conda/envs/phylosolid_env/lib/R/library
exec "$@"
EOF

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["/bin/bash"]