# Base image with conda
FROM continuumio/miniconda3:latest

# Set working directory
WORKDIR /app

# Set environment variables
ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    CONDA_AUTO_UPDATE_CONDA=false \
    PATH="/opt/conda/bin:${PATH}"

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gcc \
    g++ \
    make \
    cmake \
    wget \
    curl \
    git \
    vim \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libudunits2-dev \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    libsqlite3-dev \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Copy environment configuration
COPY environment.yml /tmp/environment.yml

# Create conda environment
RUN conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda config --add channels r && \
    conda config --set channel_priority flexible && \
    conda env create -f /tmp/environment.yml && \
    conda clean -afy

# Set conda environment path
ENV PATH="/opt/conda/envs/phylosolid_env/bin:${PATH}"

# Set shell to use conda environment
SHELL ["conda", "run", "-n", "phylosolid_env", "/bin/bash", "-c"]

# Install additional R packages
RUN R -e "install.packages(c('Seurat', 'SeuratDisk', 'scATOMIC', 'Polychrome', 'ggnewscale', 'ggtext', 'gsubfn', 'paletteer', 'ggforce'), repos='https://cloud.r-project.org/')" && \
    R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install(c('ggtreeExtra', 'ComplexHeatmap'))"

# Create directory structure
RUN mkdir -p /app/src /app/data /app/output /app/logs /app/module /app/utils

# Copy project files (uncomment and modify as needed)
# COPY ./src /app/src/
# COPY ./module /app/module/
# COPY ./utils /app/utils/
# COPY ./*.py /app/

# Set working directory
WORKDIR /app

# Expose ports
EXPOSE 8888
EXPOSE 8050

# Default command
CMD ["python", "main.py"]

# Alternative: interactive shell
# CMD ["/bin/bash"]

# Build instructions:
# docker build -t phylosolid:latest .
# docker run -it --rm -v $(pwd)/data:/app/data -v $(pwd)/output:/app/output phylosolid:latest

# GPU support (uncomment if needed)
# FROM nvidia/cuda:11.8.0-cudnn8-runtime-ubuntu20.04
# Then add GPU-specific configurations

