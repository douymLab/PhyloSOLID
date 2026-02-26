#!/bin/bash
# install_with_mamba.sh - One-click installation script using mamba

set -e  # 出错时停止

echo "Installing PhyloSOLID dependencies..."

# 检查是否安装了 mamba，如果没有则安装
if ! command -v mamba &> /dev/null; then
    echo "Installing mamba for faster dependency resolution..."
    conda install -n base -c conda-forge mamba -y
fi

# 用 mamba 创建环境（比 conda 快很多）
mamba env create -f environment.yml -y

echo "Installation complete!"
echo "To activate: conda activate phylosolid_env"
