##Author: Zhirui Yang, Qing Yang
##Date: 2024-07-26
##Update: 2024/01/02
##Details: This is one script to calculate hFDR.


##### Time #####
library(stringr)
start_time <- Sys.time()
print(str_c("Start time: ", start_time))


##### library
library(optparse)
option_list = list(
    make_option(c("--feature_file"), type = "character", default = "/storage/douyanmeiLab/yangqing/tools/PhyloMosaicGenie/pmg/pre-classifier/scRNA/all_features/P6_normal.merged_identifier.feature_sigFilter.txt", help = ""),
    make_option(c("--ase_editing"), type = "character", default = "/storage/douyanmeiLab/yangqing/tools/PhyloMosaicGenie/pmg/pre-classifier/scRNA/all_features/P6_normal.ASE_RNAediting/identifer_stat.txt", help = ""),
    make_option(c("--outputfile"), type = "character", default = "/storage/douyanmeiLab/yangqing/tools/PhyloMosaicGenie/pmg/pre-classifier/scRNA/all_features/P6_normal.merged_identifier.feature_sigFilter_ASE_RNAediting_UMIconsistence.txt", help = "")
)
parseobj = OptionParser(option_list=option_list)
opt = parse_args(parseobj)


feature_file <- as.character(opt$feature_file)
ase_editing <- as.character(opt$ase_editing)
outputfile <- as.character(opt$outputfile)

# if (!file.exists(outputpath)) {
#   dir.create(outputpath)
#   print(paste("Path", outputpath, "created successfully"))
# } else {
#   print(paste("Path", outputpath, "already exists"))
# } 


#===========functions=================
library(pracma)
library(dplyr)
library(purrr)
library(ggplot2)


# function to get standard mutation signature matrix


#=============handle data=============

# read files
df_features <- read.table(feature_file, header = TRUE, sep = "\t")
df_stat <- read.table(ase_editing, header = TRUE, sep = "\t")

# Create ASE_filter and RNAediting_filter columns
df_stat$ASE_filter <- ifelse(df_stat$ase %in% c("Unknown", "False"), "pass", "fail")
df_stat$RNAediting_filter <- ifelse(df_stat$editing == "False", "pass", "fail")

# ASE/RNAediting: merge the two data frames
df_out <- merge(df_features, df_stat[, c("identifier", "ASE_filter", "RNAediting_filter")], by = "identifier", all.x = TRUE)

# UMI_consistence_prop
df_out$UMI_consistence_filter <- ifelse(df_out$alt_UMI_consistence_prop >= 0.85, "pass", "fail")


# Save results
write.table(df_out, outputfile, sep = "\t", row.names = FALSE, quote = FALSE)



##### Time #####
end_time <- Sys.time()
print(str_c("End time: ", end_time))
print(str_c("Program finished in ", as.character(round((end_time-start_time), 4)), " seconds"))

