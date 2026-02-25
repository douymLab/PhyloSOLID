# Date: 2025/03/18
# Update: 2025/03/19
# Author: Qing
# Work: Cell type annotation for 10x scRNA-seq data on your sample.


# salloc --time=24:00:00 --mem=128G
# module load R/4.2.1
# /home/douyanmeiLab/xiayonghe/.conda/envs/insight/bin/R
# module load cmake/3.28.0



##### Python and packages path
library(reticulate)
use_python("/home/douyanmeiLab/xiayonghe/.conda/envs/insight/bin/python")
print("===== Test python packages =====")
Rmagic::pymagic_is_available()

.libPaths(c("/home/douyanmeiLab/yangqing/R/x86_64-conda-linux-gnu-library/for_cellanno", "/home/douyanmeiLab/yangzhirui/R/x86_64-pc-linux-gnu-library/4.2_for_seurat"))

##### Time #####
library(stringr)
start_time <- Sys.time()
print(str_c("Start time: ", start_time))

##### library
library(optparse)
option_list = list(
    make_option(c("-i", "--inputpath_10xdata"), type="character", default="/storage/douyanmeiLab/yangqing/tools/PhyloMosaicGenie/Benchmark/data_10xCI792/bam/CI792_output/outs/filtered_feature_bc_matrix/", help="The datapath that need in the program."),
    make_option(c("-s", "--sampleid"), type="character", default="CI792", help="The sample id."),
    make_option(c("-o", "--outputpath"), type="character", default="/storage/douyanmeiLab/yangqing/tools/PhyloMosaicGenie/Benchmark/data_10xCI792/cellanno/results_test/", help="The outputpath, generally set to 'results'.")
    )
parseobj = OptionParser(option_list=option_list,
            usage="usage: %prog [options]",
            add_help_option=TRUE,
            description="Valid the converted tree and calculate the accuray of the phylogenetic tree.")
opt = parse_args(parseobj)

# library(SingleR)
# 使用HumanPrimaryCellAtlasData()函数加载参考数据集
#hpca.se <- HumanPrimaryCellAtlasData()
library(Seurat)
# library(anndata)
library(SeuratDisk)
library(devtools)
library(gridExtra)
library(scATOMIC)
library(AnnotationDbi)
library(org.Hs.eg.db)
# library(biomaRt)
library(dplyr)
library(ggplot2)


##### Parameters
sampleid <- as.character(opt$sampleid)
inputpath_10xdata <- as.character(opt$inputpath_10xdata)
outputpath <- as.character(opt$outputpath)
dir.create(outputpath, recursive=TRUE)

print(str_c("Parameter sampleid is: ", sampleid))
print(str_c("Parameter inputpath_10xdata is: ", inputpath_10xdata))
print(str_c("Parameter outputpath is: ", outputpath))


##########################################
##### run

### Load dataset
mat_sample <- Seurat::Read10X(inputpath_10xdata, gene.column = 1)


### Preprocessing dataset
pct_mt <- colSums(mat_sample[grep("^MT-", row.names(mat_sample)),])/colSums(mat_sample) * 100
nFeatureRNA <- colSums(mat_sample > 0)
mat_sample <- mat_sample[, names(which(pct_mt < 25))]
mat_sample <- mat_sample[, intersect(names(which(nFeatureRNA > 500)), colnames(mat_sample))]
# a <- mat_sample

### 查看数据结构
str(mat_sample)
# Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#   ..@ i       : int [1:15041955] 23 36 60 62 70 73 80 83 86 97 ...
#   ..@ p       : int [1:4476] 0 3805 8054 10015 10596 14455 16855 18943 24362 25659 ...
#   ..@ Dim     : int [1:2] 36601 4475
#   ..@ Dimnames:List of 2
#   .. ..$ : chr [1:36601] "ENSG00000243485" "ENSG00000237613" "ENSG00000186092" "ENSG00000238009" ...
#   .. ..$ : chr [1:4475] "AAACCCAGTAGCACAG-1" "AAACCCAGTATTTCTC-1" "AAACCCAGTTGTCATG-1" "AAACGAAAGCAGTAAT-1" ...
#   ..@ x       : num [1:15041955] 3 2 1 1 1 1 2 1 3 1 ...
#   ..@ factors : list()
dim(mat_sample)
# [1] 36601  4475

length(mat_sample@Dimnames[[1]])
# [1] 36601
head(mat_sample@Dimnames[[1]])
# [1] "ENSG00000243485" "ENSG00000237613" "ENSG00000186092" "ENSG00000238009" "ENSG00000239945" "ENSG00000239906"


##### convert Ensembl ID (such as ENSG00000243485) into Gene Symbol (such as FAM138A)
# ### 获取 merged mapping ids
# # load manual data
# mapping_manual <- read.table("/storage/douyanmeiLab/wuxing/data/0_refGenome/gene_id2name_grh38.txt", sep="\t", header=TRUE)
# colnames(mapping_manual) <- c("ENSEMBL", "SYMBOL")
# dim(mapping_manual)
# # [1] 70711     2
# # load database
# ensembl_ids <- rownames(mat_sample)
# mapping_hsdb <- AnnotationDbi::select(org.Hs.eg.db,
#                                       keys = ensembl_ids,
#                                       keytype = "ENSEMBL",
#                                       columns = "SYMBOL")
# dim(mapping_hsdb)
# # [1] 36825     2
# # 合并两个数据框
# merged_mapping <- full_join(mapping_manual, mapping_hsdb, by = "ENSEMBL", suffix = c("_manual", "_hsdb"))
# # 填充 SYMBOL，按照条件选择
# merged_mapping$SYMBOL <- ifelse(
#   # 如果两个 SYMBOL 都是 NA 或空字符串，结果为 NA
#   is.na(merged_mapping$SYMBOL_manual) & is.na(merged_mapping$SYMBOL_hsdb) | 
#     merged_mapping$SYMBOL_manual == "" & merged_mapping$SYMBOL_hsdb == "", 
#   NA, 
  
#   # 如果 SYMBOL_manual 是 NA 或空字符串，保留 SYMBOL_hsdb
#   ifelse(
#     is.na(merged_mapping$SYMBOL_manual) | merged_mapping$SYMBOL_manual == "", 
#     merged_mapping$SYMBOL_hsdb, 
    
#     # 如果 SYMBOL_hsdb 是 NA 或空字符串，保留 SYMBOL_manual
#     ifelse(
#       is.na(merged_mapping$SYMBOL_hsdb) | merged_mapping$SYMBOL_hsdb == "", 
#       merged_mapping$SYMBOL_manual, 
      
#       # 如果两个都不是 NA 或空字符串，优先保留 SYMBOL_hsdb
#       merged_mapping$SYMBOL_hsdb
#     )
#   )
# )
# # 只保留有 ENSEMBL 和 SYMBOL 两列
# merged_mapping <- merged_mapping %>%
#   select(ENSEMBL, SYMBOL)
# dim(merged_mapping)
# # [1] 71159     2
# write.table(merged_mapping, str_c(outputpath, "/merged_mapping_geneNames_for_human.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
merged_mapping <- read.table("/storage/douyanmeiLab/yangqing/tools/PhyloMosaicGenie/Benchmark/data_10xCI792/cellanno/results_interact/merged_mapping_geneNames_for_human.txt", sep="\t", header=TRUE)
# 计算 ENSEMBL 列中 NA 或空字符串的个数
ensembl_na_empty_count <- sum(is.na(merged_mapping$ENSEMBL) | merged_mapping$ENSEMBL == "")
# 计算 SYMBOL 列中 NA 或空字符串的个数
symbol_na_empty_count <- sum(is.na(merged_mapping$SYMBOL) | merged_mapping$SYMBOL == "")
# 输出结果
cat("ENSEMBL 列中 NA 或空字符串的个数: ", ensembl_na_empty_count, "\n")
# ENSEMBL 列中 NA 或空字符串的个数:  0 
cat("SYMBOL 列中 NA 或空字符串的个数: ", symbol_na_empty_count, "\n")
# SYMBOL 列中 NA 或空字符串的个数:  21757 
# 计算 SYMBOL 列中非 NA 且非空字符串的个数
symbol_non_na_empty_count <- sum(!is.na(merged_mapping$SYMBOL) & merged_mapping$SYMBOL != "")
# 输出结果
cat("SYMBOL 列中非 NA 且非空字符串的个数: ", symbol_non_na_empty_count, "\n")
# SYMBOL 列中非 NA 且非空字符串的个数:  49402 

### 把gene Ensembl_ID(ENSG...)转换成gene_symbol, 需要转换的Ensembl_ID 保留前15位
# 1. 创建映射关系
gene_map <- setNames(merged_mapping$SYMBOL, merged_mapping$ENSEMBL)
# 2. 将 'mat_sample' 的行名替换为对应的基因符号
rownames(mat_sample) <- gene_map[rownames(mat_sample)]
# 3. 查看更新后的行名
length(rownames(mat_sample))
# [1] 36601
sum(rownames(mat_sample) != "" & !is.na(rownames(mat_sample)))
# [1] 26455
sum(is.na(rownames(mat_sample)) | rownames(mat_sample) == "")
# [1] 10146

# 过滤掉基因名为 NA 或空字符串的行
mat_sample <- mat_sample[!is.na(rownames(mat_sample)) & rownames(mat_sample) != "", ]
sum(rownames(mat_sample) != "" & !is.na(rownames(mat_sample)))
# [1] 26455
sum(is.na(rownames(mat_sample)) | rownames(mat_sample) == "")
# [1] 0

### 为 gene name 去重
# 读取原始行名
original_rownames <- rownames(mat_sample)
# 创建一个副本用于存储唯一化的行名
unique_rownames <- original_rownames
# 记录已经遇到的行名及其出现次数
name_counts <- list()
# 遍历所有行名
for (i in seq_along(unique_rownames)) {
  name <- unique_rownames[i]
  # 仅对非 NA 且非空字符串的行名进行处理
  if (!is.na(name) && name != "") {
    if (name %in% names(name_counts)) {
      # 计数加 1，并修改名称
      name_counts[[name]] <- name_counts[[name]] + 1
      unique_rownames[i] <- paste0(name, ".", name_counts[[name]])
    } else {
      # 第一次遇到该名称，初始化计数
      name_counts[[name]] <- 0
    }
  }
}
# 重新赋值给矩阵
rownames(mat_sample) <- unique_rownames
sum(duplicated(rownames(mat_sample)))  # 结果应为 0

### 查看数据结构
str(mat_sample)  
# Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#   ..@ i       : int [1:14672689] 11 20 40 42 50 53 57 59 62 71 ...
#   ..@ p       : int [1:4476] 0 3722 7864 9783 10350 14123 16475 18531 23814 25081 ...
#   ..@ Dim     : int [1:2] 26455 4475
#   ..@ Dimnames:List of 2
#   .. ..$ : Named chr [1:26455] "MIR1302-2HG" "FAM138A" "OR4F5" "OR4F29" ...
#   .. .. ..- attr(*, "names")= chr [1:26455] "ENSG00000243485" "ENSG00000237613" "ENSG00000186092" "ENSG00000284733" ...
#   .. ..$ : chr [1:4475] "AAACCCAGTAGCACAG-1" "AAACCCAGTATTTCTC-1" "AAACCCAGTTGTCATG-1" "AAACGAAAGCAGTAAT-1" ...
#   ..@ x       : num [1:14672689] 3 2 1 1 1 1 2 1 3 1 ...
#   ..@ factors : list()
dim(mat_sample)
# [1] 26455  4475




##### running
cell_predictions <- run_scATOMIC(mat_sample)
# [1] "Starting Layer 1"
# /home/douyanmeiLab/yangqing/.cache/R/reticulate/uv/cache/archive-v0/4vFbFQa7nxtCdbUsxmJ8P/lib/python3.11/site-packages/magic/magic.py:425: UserWarning: Input matrix contains unexpressed genes. Please remove them prior to running MAGIC.
#   warnings.warn(
# [1] "Done Layer 1"
# [1] "Starting Layer 2 Non Blood"
# /home/douyanmeiLab/yangqing/.cache/R/reticulate/uv/cache/archive-v0/4vFbFQa7nxtCdbUsxmJ8P/lib/python3.11/site-packages/magic/magic.py:425: UserWarning: Input matrix contains unexpressed genes. Please remove them prior to running MAGIC.
#   warnings.warn(
# [1] "Done Layer 2 Non Blood"
# [1] "Starting Layer 3 Non Stromal"
# /home/douyanmeiLab/yangqing/.cache/R/reticulate/uv/cache/archive-v0/4vFbFQa7nxtCdbUsxmJ8P/lib/python3.11/site-packages/magic/magic.py:425: UserWarning: Input matrix contains unexpressed genes. Please remove them prior to running MAGIC.
#   warnings.warn(
# [1] "Done Layer 3 Non Stromal"
# [1] "Starting Layer 4 Non GI"
# /home/douyanmeiLab/yangqing/.cache/R/reticulate/uv/cache/archive-v0/4vFbFQa7nxtCdbUsxmJ8P/lib/python3.11/site-packages/magic/magic.py:425: UserWarning: Input matrix contains unexpressed genes. Please remove them prior to running MAGIC.
#   warnings.warn(
# [1] "Done Layer 4 Non GI"
# [1] "Starting Layer 5 Breast Lung Prostate"
# /home/douyanmeiLab/yangqing/.cache/R/reticulate/uv/cache/archive-v0/4vFbFQa7nxtCdbUsxmJ8P/lib/python3.11/site-packages/magic/magic.py:425: UserWarning: Input matrix contains unexpressed genes. Please remove them prior to running MAGIC.
#   warnings.warn(
# [1] "nothing to score in this layer"
# [1] "Done Layer 5 Breast Lung Prostate"
# [1] "Starting Layer 2 Blood"
# /home/douyanmeiLab/yangqing/.cache/R/reticulate/uv/cache/archive-v0/4vFbFQa7nxtCdbUsxmJ8P/lib/python3.11/site-packages/magic/magic.py:425: UserWarning: Input matrix contains unexpressed genes. Please remove them prior to running MAGIC.
#   warnings.warn(
# [1] "Done Layer 2 Blood"
# [1] "Starting Layer 3 TNK"
# [1] "Done Layer 3 TNK"
# [1] "Starting Layer 4 CD4 CD8"
# [1] "nothing to score in this layer"
# [1] "Done Layer 4 CD4 CD8"
# [1] "Starting Layer 4 CD8 NK"
# [1] "nothing to score in this layer"
# [1] "Done Layer 4 CD8 NK"
# [1] "Starting Layer 5 CD4"
# [1] "nothing to score in this layer"
# [1] "Done Layer 5 CD4"
# [1] "Starting Layer 5 CD8"
# [1] "nothing to score in this layer"
# [1] "Done Layer 5 CD8"
# [1] "Starting Layer 3 Myeloid"
# /home/douyanmeiLab/yangqing/.cache/R/reticulate/uv/cache/archive-v0/4vFbFQa7nxtCdbUsxmJ8P/lib/python3.11/site-packages/magic/magic.py:425: UserWarning: Input matrix contains unexpressed genes. Please remove them prior to running MAGIC.
#   warnings.warn(
# /home/douyanmeiLab/yangqing/.cache/R/reticulate/uv/cache/archive-v0/4vFbFQa7nxtCdbUsxmJ8P/lib/python3.11/site-packages/graphtools/base.py:165: RuntimeWarning: Cannot perform PCA to 100 dimensions on data with min(n_samples, n_features) = 74
#   warnings.warn(
# [1] "Done Layer 3 Myeloid"
# [1] "Starting Layer 4 Macrophage Monocyte"
# /home/douyanmeiLab/yangqing/.cache/R/reticulate/uv/cache/archive-v0/4vFbFQa7nxtCdbUsxmJ8P/lib/python3.11/site-packages/magic/magic.py:425: UserWarning: Input matrix contains unexpressed genes. Please remove them prior to running MAGIC.
#   warnings.warn(
# /home/douyanmeiLab/yangqing/.cache/R/reticulate/uv/cache/archive-v0/4vFbFQa7nxtCdbUsxmJ8P/lib/python3.11/site-packages/graphtools/base.py:165: RuntimeWarning: Cannot perform PCA to 100 dimensions on data with min(n_samples, n_features) = 73
#   warnings.warn(
# [1] "nothing to score in this layer"
# [1] "Done Layer 4 Macrophage Monocyte"
# [1] "Starting Layer 5 Macrophage"
# [1] "nothing to score in this layer"
# [1] "Done Layer 5 Macrophage"
# [1] "Starting Layer 3 B Cell"
# /home/douyanmeiLab/yangqing/.cache/R/reticulate/uv/cache/archive-v0/4vFbFQa7nxtCdbUsxmJ8P/lib/python3.11/site-packages/magic/magic.py:425: UserWarning: Input matrix contains unexpressed genes. Please remove them prior to running MAGIC.
#   warnings.warn(
# /home/douyanmeiLab/yangqing/.cache/R/reticulate/uv/cache/archive-v0/4vFbFQa7nxtCdbUsxmJ8P/lib/python3.11/site-packages/graphtools/base.py:165: RuntimeWarning: Cannot perform PCA to 100 dimensions on data with min(n_samples, n_features) = 49
#   warnings.warn(
# [1] "Done Layer 3 B Cell"

results_sample <- create_summary_matrix(prediction_list = cell_predictions, raw_counts = mat_sample)

length(table(results_sample$scATOMIC_pred))
table(results_sample$scATOMIC_pred)

##### Adding cell type information into seurat object


#create seurat object
seurat_sample = CreateSeuratObject(counts = mat_sample, project = sampleid, min.cells = 3,  min.features = 200, meta.data = results_sample)
#run seurat pipeline
seurat_sample <- NormalizeData(seurat_sample)
seurat_sample <- FindVariableFeatures(seurat_sample)
seurat_sample <- ScaleData(seurat_sample)
seurat_sample <- RunPCA(seurat_sample, features = VariableFeatures(object = seurat_sample))
seurat_sample <- RunUMAP(seurat_sample, dims = 1:50)
seurat_sample <- FindNeighbors(seurat_sample)
seurat_sample <- FindClusters(seurat_sample)
# seurat_sample <- AddMetaData(seurat_sample, results_lung)
saveRDS(seurat_sample, file = str_c(outputpath, "/", sampleid, ".seurat_data_with_celltype.rds"))
# seurat_sample <- readRDS(file = str_c("/storage/douyanmeiLab/yangqing/tools/PhyloMosaicGenie/Benchmark/data_10xCI792/cellanno/results_interact/seurat_sample.with_celltype.rds"))


# plot results
p_cellanno <- DimPlot(seurat_sample, group.by = "scATOMIC_pred") + ggtitle(str_c("The ", sampleid, " scRNA-seq data")) + labs(fill="Cell types")
ggsave(str_c(outputpath, "/", sampleid, ".UMAP_plot_with_celltype.pdf"), plot = p_cellanno, width = 10, height = 4)


##### generate annotation file for SComatic
df_metadata <- seurat_sample@meta.data
dim(df_metadata)
# [1] 4475   26
## raw cell type
df_metadata <- df_metadata %>%
  mutate(celltype = str_replace_all(scATOMIC_pred, " ", "_"))
df_metadata$celltype <- gsub("/", "_", df_metadata$celltype)  # 替换 "/" 为 "_"
## output for SComatic
df_anno_for_scomatic <- df_metadata[, c("cell_names", "celltype")]
dim(df_anno_for_scomatic)
# [1] 4475    2
colnames(df_anno_for_scomatic) <- c("Index", "Cell_type")
write.table(df_anno_for_scomatic, str_c(outputpath, "/", sampleid, "_annotation.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

seurat_sample@meta.data$celltype <- df_anno_for_scomatic$Cell_type
saveRDS(seurat_sample, file = str_c(outputpath, "/", sampleid, ".seurat_data_with_celltype_processed.rds"))




