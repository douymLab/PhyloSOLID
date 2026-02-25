#!/usr/bin/env Rscript

library(converTree)
library(ape)
library(treeio)
library(purrr)

# 设置命令行参数
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("用法: Rscript convert_tree.R <输入文件> <输出文件>")
}

tree_file <- args[1]
output_file <- args[2]

# 主要处理函数
convert_tree <- function(input_file, output_file) {
  # 读取数据并转换
  treedata_truth <- cf2treedata(input_file)
  phylo_tree <- treedata_truth %>% ape::as.phylo()
  
  # 索引标签函数
  index_label <- function(node) {
    treedata_truth[treedata_truth$node == node, ]$label
  }
  
  # 更新叶节点标签
  phylo_tree[["tip.label"]] <- map(phylo_tree[["tip.label"]], index_label) %>% unlist()
  
  # 生成Newick格式树
  nwk_tree <- treeio::write.tree(phylo_tree)
  
  # 去掉内部节点数字（新增这行）
  nwk_tree <- gsub(")[0-9]+", ")", nwk_tree)
  
  # 保存到文件
  writeLines(nwk_tree, output_file)
  
  cat("转换完成！结果保存到:", output_file, "\n")
  cat("生成的树:", nwk_tree, "\n")
}

# 执行转换
convert_tree(tree_file, output_file)