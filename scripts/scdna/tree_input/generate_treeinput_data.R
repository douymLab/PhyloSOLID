# Date: 2023/06/01
# Update: 2024/12/09
# Author: Qing Yang
# Work: Input the original rawdata and complete data preprocessing, which main purpose is extracting somatic posterior matrix, reads matrix and features dataframe.


##### Time #####
library(stringr)
start_time <- Sys.time()
print(str_c("Start time: ", start_time))


##### library
library(optparse)
option_list = list(
    make_option(c("-i", "--inputfile"), type = "character", default = NA, help = "the name of inputfile - posterior raw data & reads information"),
    make_option(c("-n", "--cellnum"), type = "integer", default = 22, help = "Cell number."),
    make_option(c("-o", "--outputpath"), type = "character", default = "data/", help = "The outputpath you can set"),
    make_option(c("-c", "--scid_file"), type = "character", default = "no", help = "Please enter the file containing the cellid used to calculate the genotype (col: sampleid) and the cellid consistent with the ground-truth tree (col: scid_basedTree). If you do not need to replace the cell order, please ignore this parameter directly."),
    make_option(c("-r", "--is_remove_cells"), type = "character", default = "no", help = "Choice option (yes/no). If in order to validate whether the mutation can pass through the phylogenetic tree, it is recommended to set 'no' to perform the operation of not deleting; In order to build a phylogenetic tree from scratch, it is recommended to remove irrelevant cells."),
    make_option(c("-t", "--threshold"), type = "character", default = 0.9, help = "The shreshold of somatic posterior filteration."), 
    make_option(c("-s", "--indid"), type = "character", default = "UMB1465", help = "The individual id of your input.")
)
parseobj = OptionParser(option_list=option_list)
opt = parse_args(parseobj)

library(purrr)
library(ggplot2)


##### Parameters
inputfile <- as.character(opt$inputfile)
cellnum <- as.numeric(as.character(opt$cellnum))
outputpath <- as.character(opt$outputpath)
dir.create(outputpath, recursive=TRUE)
scid_file <- as.character(opt$scid_file)
is_remove_cells <- as.character(opt$is_remove_cells)
valid_choices <- c("yes", "no")
if (!is_remove_cells %in% valid_choices) {
  stop("ParaError: --is_remove_cells must be 'yes' or 'no'")
}
cutoff <- as.numeric(as.character(opt$threshold))
indid <- as.character(opt$indid)

# Display the parameters for verification
cat("Parameters:\n")
cat("Input file: ", inputfile, "\n")
cat("Cell number: ", cellnum, "\n")
cat("Output path: ", outputpath, "\n")
cat("SCID file: ", scid_file, "\n")
cat("Remove cells: ", is_remove_cells, "\n")
cat("Threshold: ", cutoff, "\n")
cat("Individual ID: ", indid, "\n")


##### Functions
get.mut_allele <- function(mut_type){
    str <- mut_type
    split <- unlist(strsplit(str, "[ ,()]"))
    gt_pairs <- split[split != ""]
    ref_gt <- gt_pairs[1:2]
    mut_gt <- gt_pairs[3:4]
    if(length(unique(mut_gt))==1){
        mut_allele <- unique(mut_gt)
    } else {
        mut_allele <- mut_gt[!mut_gt %in% ref_gt]
    }
    return(gsub("'", "", mut_allele))
}

get_alleles <- function(mut_type_str) {
  # Remove extra spaces and parentheses, then split the string into a vector
  mut_type <- gsub("[()]", "", mut_type_str)  # Remove parentheses
  mut_type <- gsub("'", "", mut_type)          # Remove single quotes
  mut_type <- strsplit(mut_type, ",")[[1]]     # Split the string
  # Convert elements to character type
  mut_type <- trimws(mut_type)  # Trim whitespace
  if (length(mut_type) != 4) {
    return(c(NA, "Invalid input: must be a tuple with 4 elements."))
  }
  # Perform comparison
  if (mut_type[1] != mut_type[3]) {  # Compare the 1st and 3rd elements
    raw_allele <- mut_type[1]
    mut_allele <- mut_type[3]
  } else if (mut_type[2] != mut_type[4]) {  # Compare the 2nd and 4th elements
    raw_allele <- mut_type[2]
    mut_allele <- mut_type[4]
  } else {
    return(c(NA, "No mutation found."))  # If no difference is found, return a message
  }
  # Return the mutant allele
  return(c(raw_allele, mut_allele))
}

get_list.base_info <- function(info_set){
    set_str <- gsub("\\{|\\}", "", info_set)
    key_value_pairs <- strsplit(set_str, ", ")[[1]]
    my_list <- list()
    for(pair in key_value_pairs){
      parts <- strsplit(pair, ": ")[[1]]
      key <- parts[1]
      value <- as.numeric(parts[2])
      my_list[[key]] <- value
    }
    return(my_list)
}

stat.allele_info <- function(info_set, mut_allele){
    mutant_idx <- grep(mut_allele, names(get_list.base_info(info_set)))
    mut_allele_count <- 0
    for(idx in mutant_idx){
        mut_allele_count <- mut_allele_count + get_list.base_info(info_set)[[names(get_list.base_info(info_set))[idx]]]
    }
    total_allele_count <- sum(unlist(get_list.base_info(info_set)))
    allele_info <- as.character(str_c(mut_allele_count, "/", total_allele_count))
    return(allele_info)
}

mean_non_zero <- function(x){
    average_non_zero <- mean(x[x != 0], na.rm = TRUE)
    return(average_non_zero)
}

likelihood_split <- function(str_posterior){
    list <- as.numeric(strsplit(strsplit(strsplit(str_posterior, "[", fixed=T)[[1]][2], "]", fixed=T)[[1]][1], ", ")[[1]])
    return(list)
}

extract_read_level_features <- function(basicdata, infodata, output_path) {
  rows_per_output=1000
  ### Get the first few data columns from basicdata
  outputdata <- as.data.frame(cbind(basicdata[, c(1:7)], somatic_posterior_persite=as.data.frame(infodata$somatic_posterior_persite))); dim(outputdata)
  colnames(outputdata) <- c('chr', 'start', 'end', 'ref', 'allele_pool', 'popAF', 'genotype', 'somatic_posterior_persite')
  ### The last mut state column
  patied_extract_posterior_data <- infodata[, remained_scid]
  patied_extract_allelestat_data <- infodata[, paste0(remained_scid, "_AlleleStat")]
  # Create an empty data frame to store results
  result_data <- data.frame(matrix(nrow = nrow(patied_extract_posterior_data), ncol = ncol(patied_extract_posterior_data)))
  colnames(result_data) <- remained_scid
  # Iterate over each cell to apply the rules
  for (i in seq_along(remained_scid)) {
    mut_allele_count <- patied_extract_allelestat_data[, i]
    posterior <- patied_extract_posterior_data[, i]
    # Handle NA values by assigning default values
    mut_allele_count[is.na(mut_allele_count)] <- 0  # Treat NA as 0
    posterior[is.na(posterior)] <- 0               # Treat NA as 0
    # Apply the rules for classification
    result_data[, i] <- ifelse(
      mut_allele_count == 0, 0,  # If mut_allele_count is 0, the result is 0
      ifelse(mut_allele_count > 0 & posterior < 0.5, 1, 2) # Otherwise classify based on posterior
    )
  }
  result_data <- cbind(bulk=rep(1, nrow(result_data)), result_data)
  # Generate the mut_percell column
  mut_percell <- apply(result_data, 1, function(row) {
    paste0("[", paste(row, collapse = ", "), "]") # Join each row's elements with commas and wrap in square brackets
  })
  # Add the new column to the basicdata data frame
  outputdata$mut_percell <- mut_percell; dim(outputdata)
  ### Write results to files in batches
  total_rows <- nrow(outputdata)
  batch_start <- 1
  batch_count <- 0
  while (batch_start <= total_rows) {
    batch_count <- batch_count + 1
    batch_end <- min(batch_start + rows_per_output - 1, total_rows)
    batch_data <- outputdata[batch_start:batch_end, ]
    # Write the current batch to a file
    write.table(batch_data, file=str_c(output_path, "/read_level_features_input.batch_", as.character(batch_count), ".txt"), row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
    # Update the starting row index
    batch_start <- batch_end + 1
  }
  # Return the complete outputdata data frame
  return(outputdata)
}


##########################################################
##### Read in datafile
# 10	100887746	100887747	T	"['T', 'C']"	"{'T': 0.99999, 'C': 1e-05, 'G': 0.0001, 'A': 0.0001}"	"('T', 'T', 'C', 'T')"	1   "[0.6792483610799689, 0.00022161301493556356, 2.4725270340166972e-14, 4.444646655469826e-15, 1.14794505679687e-21, 9.04031190777023e-12, 3.4904489563212095e-17, 3.4693765041291226e-17, 1.849651061910095e-08, 1.1782291811109879e-06, 4.96363963159897e-34, 1.7757789972178067e-14, 6.71625537342966e-23, 2.949279418776493e-30, 1.335025741205901e-22, 3.645471284806071e-11, 7.861860271611663e-22, 5.216365272889868e-14, 4.038071832940647e-30, 9.434140700789552e-08, 2.7213350254134875e-10, 1.0]"	"[-36.79021198351468, -13.882303490116074, -98.73375634103046, -30.534186653291318, -61.74913146145162, -22.895864295072418, -35.396561745256285, -35.3845175469301, -15.269265406890558, -11.101025112762384, -74.26152761593, -29.145552797865538, -48.583891468198914, -73.55097171361723, -47.878526675630056, -21.525636231145583, -54.128875400721995, -36.08435371564084, -65.23676158198413, -67.48971275340719, -75.66722205251457, -34.696232101240554]"	"[-40.107011360822995, -8.03442690950825, -69.9692854124227, -0.053590286578639255, -16.09929760037022, -0.03301651100825664, -0.06912518805259608, -0.05102551701276025, -0.030061342525672167, -0.01600800533733655, -0.1422533774400263, -0.05008034924684311, -0.09544474743495707, -8.121459454463091, -0.07708485820409829, -0.05716551656775981, -8.10050626404293, -8.066442920182924, -0.12145575610712889, -53.87984685774567, -56.208973355883586, -300.8540359115743]"	"[{'T_30': 41, 'T_20': 1, 'C_30': 12, 'C_60': 1, 'C_20': 2, 'T_40': 1}, {'T_30': 45, 'T_20': 3, 'C_30': 5}, {'T_30': 18, 'C_30': 1, 'T_20': 1}, {'T_20': 10, 'T_30': 114, 'A_20': 1, 'T_60': 2, 'C_30': 8}, {'T_29': 2, 'T_30': 41, 'T_20': 1}, {'T_29': 2, 'T_30': 84, 'T_60': 1, 'C_30': 2}, {'T_30': 33}, {'T_30': 49, 'T_20': 2}, {'T_30': 51}, {'T_30': 20, 'T_20': 1, 'T_60': 1}, {'T_30': 16}, {'T_30': 102, 'T_20': 4, 'T_60': 1}, {'T_30': 40, 'T_50': 1, 'T_20': 1}, {'T_29': 1, 'T_20': 3, 'T_30': 64, 'T_60': 2}, {'T_30': 101, 'T_20': 1, 'T_24': 1, 'C_30': 1, 'T_60': 2}, {'T_20': 1, 'T_30': 67, 'T_60': 1}, {'T_30': 27, 'T_60': 1, 'T_20': 3}, {'T_30': 74, 'T_20': 2, 'T_60': 1, 'C_30': 1}, {'T_30': 50, 'C_30': 1, 'T_20': 1}, {'T_29': 1, 'T_30': 90, 'T_20': 3}, {'T_29': 2, 'T_20': 6, 'T_30': 75, 'A_20': 1, 'C_30': 6}, {'T_30': 94, 'T_20': 7, 'T_60': 1, 'C_30': 7}, {'T_29': 1, 'T_30': 12, 'C_30': 35, 'C_60': 1, 'C_20': 1}]"
inputdata <- read.table(inputfile, header=FALSE, sep="\t"); dim(inputdata)
# [1] 327  12
inputdata <- inputdata[inputdata[,8]!="[0]",]; dim(inputdata)
inputdata <- inputdata[inputdata[,8]!=0,]; dim(inputdata)
# [1] 327  12
inputdata <- inputdata[sapply(strsplit(gsub("\\[|\\]", "", inputdata[, 8]), ", "), function(x) x[1]) != "nan", ]; dim(inputdata)
# [1] 327  12
inputdata <- unique(inputdata); dim(inputdata)
# [1] 327  12
alleles_list <- sapply(inputdata[,7], get_alleles, simplify = FALSE)
raw_allele_vector <- sapply(alleles_list, `[`, 1)
mut_allele_vector <- sapply(alleles_list, `[`, 2)

# # Use mutate and rowwise to apply the function and add new columns
# inputdata <- inputdata %>%
#   rowwise() %>%
#   mutate(
#     alleles = list(get_alleles(mutation)),  # Call the function
#     raw_allele = alleles[[1]],               # Extract raw_allele
#     mut_allele = alleles[[2]]                # Extract mut_allele
#   ) %>%
#   select(-alleles)  # Remove the temporary column


##### ID reference (sc).
# You can select wether to replace the cell id.
if(scid_file=="no"){
    scid_data <- as.data.frame(cbind(scid_basedTree=paste0("sc" ,1:cellnum), scid=paste0("sc" ,1:cellnum)))
} else {
    scid_data <- as.data.frame(read.table(scid_file, header=TRUE))
    if(dim(scid_data)[2]==1){
        scid_data <- as.data.frame(scid_data[scid_data$scid_basedTree!="none",])
        colnames(scid_data) <- "scid_basedTree"
    } else {
        scid_data <- scid_data[scid_data$scid_basedTree!="none",]
    }
}
dim(scid_data)
# [1] 22  3


##### posterior data
posterior_data.raw <- inputdata[,c(1,3,4,6,7,8,9)]; dim(posterior_data.raw)
# [1] 327    7
posterior_data.temp <- as.data.frame(cbind(mutid=paste0(posterior_data.raw[,1], "_", posterior_data.raw[,2], "_", raw_allele_vector, "_", mut_allele_vector), posterior_data.raw[,1:dim(posterior_data.raw)[2]])); dim(posterior_data.temp)
# [1] 327   8
posterior_list <- apply(posterior_data.temp, 1, function(x) {likelihood_split(x[dim(posterior_data.temp)[2]])})
posterior_data <- as.data.frame(cbind(mutid=posterior_data.temp[, 1], indid=rep(indid, length(posterior_data.temp[, 1])), chr=posterior_data.temp[, 2], pos=posterior_data.temp[, 3], ref=posterior_data.temp[, 4], mut=mut_allele_vector, popAF=posterior_data.temp[, 5], genotype=posterior_data.temp[, 6], somatic_posterior_persite=posterior_data.temp[, 7], t(as.data.frame(apply(posterior_list, 2, function(x) x))))); dim(posterior_data)
# [1] 327 31
rownames(posterior_data) <- 1:dim(posterior_data)[1]
colnames(posterior_data) <- c("mutid", "indid", "chr", "pos", "ref", "mut", "popAF", "genotype", "somatic_posterior_persite", scid_data$scid_basedTree); dim(posterior_data)
# [1] 327  31
posterior_data[1,]
#              mutid   indid chr       pos ref mut
# 1 10_100004638_C_T UMB1465  10 100004638   C   T
#                                                  popAF             genotype
# 1 {'A': 0.0001, 'T': 0.0001, 'G': 0.0001, 'C': 0.9997} ('C', 'C', 'T', 'C')
#   somatic_posterior_persite                  sc8                  sc4
# 1                         1 2.69318174599167e-21 1.01919475184779e-29
#                    sc3                  sc6                  sc1
# 1 3.63278868728683e-10 5.67230124925702e-12 9.25327628671867e-08
#                    sc9               sc18                 sc14
# 1 9.06738430174915e-11 0.0119161576163528 4.53330164544502e-11
#                   sc20                 sc15                 sc13
# 1 1.11457224985533e-14 3.54992229888129e-13 1.81271521024423e-10
#                    sc7               sc10_2                  sc5
# 1 5.79152160601277e-09 1.08942982283971e-17 2.22240470425053e-14
#                   sc12                sc19                 sc11
# 1 8.88368943993317e-14 9.1232399299247e-11 1.42811659789008e-12
#                    sc2               sc16_1               sc16_2
# 1 1.18124904274177e-05 8.88368943993317e-14 1.82372535992193e-10
#                 sc10_1 sc17
# 1 2.75696779105544e-18    1


##### Likelihood_mut data
likelihood_mut_data.raw <- inputdata[,c(1,3,4,10)]; dim(likelihood_mut_data.raw)
# [1] 327    4
likelihood_mut_data.temp <- as.data.frame(cbind(mutid=paste0(likelihood_mut_data.raw[,1], "_", likelihood_mut_data.raw[,2], "_", raw_allele_vector, "_", mut_allele_vector), likelihood_mut_data.raw[,1:dim(likelihood_mut_data.raw)[2]])); dim(likelihood_mut_data.temp)
# [1] 327   5
likelihood_mut_list <- apply(likelihood_mut_data.temp, 1, function(x) {likelihood_split(x[dim(likelihood_mut_data.temp)[2]])})
likelihood_mut_data <- as.data.frame(cbind(mutid=likelihood_mut_data.temp[, 1], t(as.data.frame(apply(likelihood_mut_list, 2, function(x) x))))); dim(likelihood_mut_data)
# [1] 327 23
rownames(likelihood_mut_data) <- 1:dim(likelihood_mut_data)[1]
colnames(likelihood_mut_data) <- c("mutid", paste0(scid_data$scid_basedTree, "_mut")); dim(likelihood_mut_data)
# [1] 327  23
likelihood_mut_data[1,]
#              mutid          sc8_mut           sc4_mut          sc3_mut
# 1 10_100004638_C_T -44.421500526218 -63.8610048941992 -18.737669532308
#             sc6_mut           sc1_mut           sc9_mut         sc18_mut
# 1 -22.8951980727516 -13.1824673214053 -20.1207807638977 -1.3876281390953
#            sc14_mut          sc20_mut          sc15_mut          sc13_mut
# 1 -20.8137558641087 -29.1455467978435 -25.6704543509422 -19.4267939473342
#             sc7_mut        sc10_2_mut           sc5_mut          sc12_mut
# 1 -15.9578963466118 -36.0788472468341 -28.4463768514536 -27.0587487123583
#            sc19_mut          sc11_mut           sc2_mut        sc16_1_mut
# 1 -20.1329950985486 -24.2955366324939 -8.32510261225097 -27.0587487123583
#          sc16_2_mut        sc10_1_mut          sc17_mut
# 1 -19.4388381456604 -37.4893819299046 -15.2699316292114


##### Likelihood_unmut data
likelihood_unmut_data.raw <- inputdata[,c(1,3,4,11)]; dim(likelihood_unmut_data.raw)
# [1] 327    4
likelihood_unmut_data.temp <- as.data.frame(cbind(mutid=paste0(likelihood_unmut_data.raw[,1], "_", likelihood_unmut_data.raw[,2], "_", raw_allele_vector, "_", mut_allele_vector), likelihood_unmut_data.raw[,1:dim(likelihood_unmut_data.raw)[2]])); dim(likelihood_unmut_data.temp)
# [1] 327   5
likelihood_unmut_list <- apply(likelihood_unmut_data.temp, 1, function(x) {likelihood_split(x[dim(likelihood_unmut_data.temp)[2]])})
likelihood_unmut_data <- as.data.frame(cbind(mutid=likelihood_unmut_data.temp[, 1], t(as.data.frame(apply(likelihood_unmut_list, 2, function(x) x))))); dim(likelihood_unmut_data)
# [1] 327 23
rownames(likelihood_unmut_data) <- 1:dim(likelihood_unmut_data)[1]
colnames(likelihood_unmut_data) <- c("mutid", paste0(scid_data$scid_basedTree, "_unmut")); dim(likelihood_unmut_data)
# [1] 327  23
likelihood_unmut_data[1,]
#              mutid           sc8_unmut          sc4_unmut           sc3_unmut
# 1 10_100004638_C_T -0.0901820275760165 -0.137295208289275 -0.0340643438605063
#             sc6_unmut           sc1_unmut           sc9_unmut
# 1 -0.0320170106751731 -0.0190095063380872 -0.0292737278644454
#             sc18_unmut          sc14_unmut          sc20_unmut
# 1 -0.00200100066716707 -0.0290155096744225 -0.0500713491973429
#            sc15_unmut         sc13_unmut           sc7_unmut
# 1 -0.0360190120095073 -0.028014009340339 -0.0232707258629442
#          sc10_2_unmut           sc5_unmut          sc12_unmut
# 1 -0.0528002813813817 -0.0410205136769249 -0.0390195130097578
#            sc19_unmut          sc11_unmut           sc2_unmut
# 1 -0.0476292265582733 -0.0531171827152595 -0.0110065036699189
#          sc16_1_unmut        sc16_2_unmut       sc10_1_unmut        sc17_unmut
# 1 -0.0390195130097578 -0.0461136803801748 -0.089226859760099 -173.837501395311


##### Reads info data with baseq => Allele_info (mutant allele count/total allele count)
reads_data.raw <- inputdata[,c(1,3,4,12)]; dim(reads_data.raw)
# [1] 327   4
reads_data.temp <- as.data.frame(cbind(mutid=paste0(reads_data.raw[,1], "_", reads_data.raw[,2], "_", raw_allele_vector, "_", mut_allele_vector), reads_data.raw[,4])); dim(reads_data.temp)
# [1] 327   2
reads_list <- t(as.data.frame(apply(reads_data.temp, 1, function(x) {str_c("{", strsplit(strsplit(strsplit(x[dim(reads_data.temp)[2]], "[{", fixed=T)[[1]][2], "}]", fixed=T)[[1]][1], "}, {", fixed=T)[[1]], "}")})))
reads_data <- as.data.frame(cbind(mutid=reads_data.temp[, 1], bulk_info=reads_list[, 1], reads_list[, 2:(cellnum+1)])); dim(reads_data)
# [1] 327  24
rownames(reads_data) <- 1:dim(reads_data)[1]
colnames(reads_data) <- c("mutid", paste0(c("bulk", scid_data$scid_basedTree), "_baseq")); dim(reads_data)
# [1] 327  24
reads_data <- merge(posterior_data[, c("mutid", "chr", "pos", "ref", "mut")], reads_data, by="mutid", all=FALSE); dim(reads_data)
# [1] 327  28
reads_data[1,]
#             mutid chr       pos ref mut                         bulk_baseq
# 1 1_100069266_G_A   1 100069266   G   A {'G_29': 2, 'G_30': 37, 'G_20': 4}
#                                       sc8_baseq   sc4_baseq
# 1 {'G_29': 1, 'G_60': 1, 'G_30': 60, 'G_20': 1} {'G_30': 3}
#                 sc3_baseq               sc6_baseq               sc1_baseq
# 1 {'G_30': 39, 'G_20': 3} {'G_20': 1, 'G_30': 23} {'G_30': 21, 'G_60': 1}
#                 sc9_baseq                         sc18_baseq
# 1 {'G_30': 52, 'G_20': 2} {'G_30': 58, 'G_60': 2, 'G_20': 1}
#                sc14_baseq                         sc20_baseq
# 1 {'G_20': 3, 'G_30': 78} {'G_30': 49, 'G_20': 1, 'G_50': 1}
#                                       sc15_baseq
# 1 {'G_29': 2, 'G_20': 7, 'G_30': 175, 'G_60': 2}
#                           sc13_baseq                          sc7_baseq
# 1 {'G_30': 91, 'G_60': 1, 'G_20': 1} {'G_30': 57, 'G_20': 1, 'G_60': 2}
#                                    sc10_2_baseq
# 1 {'G_29': 2, 'G_30': 91, 'A_20': 1, 'G_60': 2}
#                            sc5_baseq                         sc12_baseq
# 1 {'G_30': 96, 'G_60': 1, 'G_20': 1} {'G_30': 63, 'G_20': 1, 'G_60': 1}
#                           sc19_baseq
# 1 {'G_20': 4, 'G_30': 80, 'G_60': 2}
#                                       sc11_baseq
# 1 {'G_30': 29, 'A_30': 32, 'G_20': 1, 'G_60': 1}
#                                                  sc2_baseq
# 1 {'G_20': 2, 'G_30': 52, 'G_60': 4, 'G_40': 1, 'G_23': 1}
#              sc16_1_baseq                       sc16_2_baseq
# 1 {'G_30': 46, 'G_20': 2} {'G_30': 53, 'G_29': 1, 'G_20': 2}
#               sc10_1_baseq                          sc17_baseq
# 1 {'G_30': 111, 'G_20': 6} {'G_30': 165, 'G_20': 7, 'G_60': 2}

### Allele_info (mutant allele count/total allele count)
allele_data <- as.data.frame(reads_data[, "mutid"]); dim(allele_data)
# [1] 327   1
suppressWarnings({
    for(sample_id in colnames(reads_data)[grep("_baseq", colnames(reads_data))]){
        sample.allele_info <- apply(reads_data, 1, function(x) {stat.allele_info(x[sample_id], x["mut"])})
        allele_data <- as.data.frame(cbind(allele_data, sample.allele_info))
    }
})
colnames(allele_data) <- c("mutid", paste0(c("bulk", scid_data$scid_basedTree), "_AlleleStat")); dim(allele_data)
# [1] 327   24
allele_data[1,]
#             mutid bulk_AlleleStat sc8_AlleleStat sc4_AlleleStat sc3_AlleleStat
# 1 1_100069266_G_A            0/43           0/63            0/3           0/42
#   sc6_AlleleStat sc1_AlleleStat sc9_AlleleStat sc18_AlleleStat sc14_AlleleStat
# 1           0/24           0/22           0/54            0/61            0/81
#   sc20_AlleleStat sc15_AlleleStat sc13_AlleleStat sc7_AlleleStat
# 1            0/51           0/186            0/93           0/60
#   sc10_2_AlleleStat sc5_AlleleStat sc12_AlleleStat sc19_AlleleStat
# 1              1/96           0/98            0/65            0/86
#   sc11_AlleleStat sc2_AlleleStat sc16_1_AlleleStat sc16_2_AlleleStat
# 1           32/63           0/60              0/48              0/56
#   sc10_1_AlleleStat sc17_AlleleStat
# 1             0/117           0/174


##### Merge posterior and read into a df
likelihoods_merged_data <- merge(likelihood_mut_data, likelihood_unmut_data, by="mutid"); dim(likelihoods_merged_data)
# [1] 327  45
posterior_likelihoods_merged_data <- merge(posterior_data, likelihoods_merged_data, by="mutid"); dim(posterior_likelihoods_merged_data)
# [1] 327  75
posterior_likelihoods_allele_merged_data <- merge(posterior_likelihoods_merged_data, allele_data, by="mutid"); dim(posterior_likelihoods_allele_merged_data)
# [1] 327  98
all_merged_data <- merge(posterior_likelihoods_allele_merged_data, reads_data[,c(1,6:dim(reads_data)[2])], by="mutid"); dim(all_merged_data)
# [1] 327 121
all_merged_data[all_merged_data == "0/0"] <- NA
all_merged_data[1,]
#             mutid   indid chr       pos ref mut
# 1 1_100069266_G_A UMB1465   1 100069266   G   A
#                                                  popAF             genotype
# 1 {'A': 0.0001, 'T': 0.0001, 'G': 0.9997, 'C': 0.0001} ('G', 'G', 'A', 'G')
#   somatic_posterior_persite                  sc8                 sc4
# 1                         1 5.32002725850405e-21 0.00595959515725274
#                    sc3                  sc6                  sc1
# 1 1.11490590276031e-14 2.88762196665295e-09 1.15040532516014e-08
#                    sc9                 sc18                 sc14
# 1 2.72459505733053e-18 2.12569907206745e-20 2.05456043697789e-26
#                   sc20                 sc15                 sc13
# 1 2.17019584485403e-17 5.30673870050823e-58 5.00406567224878e-30
#                    sc7               sc10_2                 sc5
# 1 4.24998006603477e-20 1.85932148094195e-28 1.5663811572592e-31
#                   sc12                 sc19 sc11                  sc2
# 1 1.33077944243199e-21 6.44642327530266e-28    1 4.26447559342928e-20
#                 sc16_1               sc16_2               sc10_1
# 1 1.74025393729901e-16 6.81662335491857e-19 3.05351280153619e-37
#                   sc17         sc8_mut           sc4_mut           sc3_mut
# 1 2.16458137420951e-54 -43.71581500536 -2.08144220864295 -29.1582572184905
#             sc6_mut           sc1_mut           sc9_mut          sc18_mut
# 1 -16.6575597683067 -15.2632433077275 -37.4780039538992 -42.3273478969281
#            sc14_mut          sc20_mut          sc15_mut          sc13_mut
# 1 -56.2170059308488 -35.3898794237944 -129.090584679394 -64.5300643447736
#             sc7_mut        sc10_2_mut           sc5_mut          sc12_mut
# 1 -41.6335338273804 -66.6111858251274 -67.9991346925119 -45.1032703974395
#            sc19_mut          sc11_mut           sc2_mut        sc16_1_mut
# 1 -59.6907659331086 -43.7156422583442 -41.6403032154745 -33.3151195366133
#          sc16_2_mut        sc10_1_mut          sc17_mut         sc8_unmut
# 1 -38.8658048400104 -81.2123787320535 -120.764470350791 -0.07134107439312
#             sc4_unmut           sc3_unmut           sc6_unmut
# 1 -0.0030015010007506 -0.0691705205702622 -0.0330618435259227
#             sc1_unmut           sc9_unmut          sc18_unmut        sc14_unmut
# 1 -0.0210115070057542 -0.0721266890533467 -0.0680813552023465 -0.10819003358002
#            sc20_unmut         sc15_unmut         sc13_unmut          sc7_unmut
# 1 -0.0590848522490949 -0.247961346400842 -0.101096866210103 -0.067080854868763
#        sc10_2_unmut          sc5_unmut          sc12_unmut         sc19_unmut
# 1 -5.79734944206152 -0.106099367878021 -0.0730828568697641 -0.120243370101689
#          sc11_unmut           sc2_unmut        sc16_1_unmut        sc16_2_unmut
# 1 -256.242828010336 -0.0772551679466653 -0.0661236870518455 -0.0743869079110367
#         sc10_1_unmut         sc17_unmut bulk_AlleleStat sc8_AlleleStat
# 1 -0.171357552148781 -0.235436906016793            0/43           0/63
#   sc4_AlleleStat sc3_AlleleStat sc6_AlleleStat sc1_AlleleStat sc9_AlleleStat
# 1            0/3           0/42           0/24           0/22           0/54
#   sc18_AlleleStat sc14_AlleleStat sc20_AlleleStat sc15_AlleleStat
# 1            0/61            0/81            0/51           0/186
#   sc13_AlleleStat sc7_AlleleStat sc10_2_AlleleStat sc5_AlleleStat
# 1            0/93           0/60              1/96           0/98
#   sc12_AlleleStat sc19_AlleleStat sc11_AlleleStat sc2_AlleleStat
# 1            0/65            0/86           32/63           0/60
#   sc16_1_AlleleStat sc16_2_AlleleStat sc10_1_AlleleStat sc17_AlleleStat
# 1              0/48              0/56             0/117           0/174
#                           bulk_baseq
# 1 {'G_29': 2, 'G_30': 37, 'G_20': 4}
#                                       sc8_baseq   sc4_baseq
# 1 {'G_29': 1, 'G_60': 1, 'G_30': 60, 'G_20': 1} {'G_30': 3}
#                 sc3_baseq               sc6_baseq               sc1_baseq
# 1 {'G_30': 39, 'G_20': 3} {'G_20': 1, 'G_30': 23} {'G_30': 21, 'G_60': 1}
#                 sc9_baseq                         sc18_baseq
# 1 {'G_30': 52, 'G_20': 2} {'G_30': 58, 'G_60': 2, 'G_20': 1}
#                sc14_baseq                         sc20_baseq
# 1 {'G_20': 3, 'G_30': 78} {'G_30': 49, 'G_20': 1, 'G_50': 1}
#                                       sc15_baseq
# 1 {'G_29': 2, 'G_20': 7, 'G_30': 175, 'G_60': 2}
#                           sc13_baseq                          sc7_baseq
# 1 {'G_30': 91, 'G_60': 1, 'G_20': 1} {'G_30': 57, 'G_20': 1, 'G_60': 2}
#                                    sc10_2_baseq
# 1 {'G_29': 2, 'G_30': 91, 'A_20': 1, 'G_60': 2}
#                            sc5_baseq                         sc12_baseq
# 1 {'G_30': 96, 'G_60': 1, 'G_20': 1} {'G_30': 63, 'G_20': 1, 'G_60': 1}
#                           sc19_baseq
# 1 {'G_20': 4, 'G_30': 80, 'G_60': 2}
#                                       sc11_baseq
# 1 {'G_30': 29, 'A_30': 32, 'G_20': 1, 'G_60': 1}
#                                                  sc2_baseq
# 1 {'G_20': 2, 'G_30': 52, 'G_60': 4, 'G_40': 1, 'G_23': 1}
#              sc16_1_baseq                       sc16_2_baseq
# 1 {'G_30': 46, 'G_20': 2} {'G_30': 53, 'G_29': 1, 'G_20': 2}
#               sc10_1_baseq                          sc17_baseq
# 1 {'G_30': 111, 'G_20': 6} {'G_30': 165, 'G_20': 7, 'G_60': 2}


######################################################################
##### Filter cells that do not contain any mutations
if(is_remove_cells == "yes") {
  
  remained_scid <- c()
  removed_scid <- c()
  
  for(sc in scid_data$scid_basedTree) {
    sc_values <- all_merged_data[, sc]
    
    # Remove if all values are below cutoff or all are NA
    if(all(sc_values < cutoff, na.rm = TRUE) || all(is.na(sc_values))) {
      removed_scid <- c(removed_scid, sc)
    } else {
      # Find sites greater than cutoff
      site_pass <- which(!is.na(sc_values) & sc_values > cutoff)
      
      if(length(site_pass) == 0) {
        removed_scid <- c(removed_scid, sc)
      } else {
        # Extract mutant allele count
        reads_list <- all_merged_data[site_pass, paste0(sc, "_AlleleStat")]
        mutant_dp <- sapply(reads_list, function(x) as.numeric(strsplit(x, "/")[[1]][1]))
        names(mutant_dp) <- rownames(all_merged_data)[site_pass]
        
        # Keep if any mutant allele count >= 1
        if(any(mutant_dp >= 1, na.rm = TRUE)) {
          remained_scid <- c(remained_scid, sc)
        } else {
          removed_scid <- c(removed_scid, sc)
        }
      }
    }
  }
  
  # ===========================
  # Update merged_data, keeping only existing columns
  # ===========================
  base_cols <- c("mutid", "indid", "chr", "pos", "ref", "mut", "somatic_posterior_persite")
  scid_cols <- remained_scid
  scid_unmut_cols <- paste0(remained_scid, "_unmut")
  scid_BB_cols <- paste0(remained_scid, "_mut_by_BB")
  scid_prod_cols <- paste0(remained_scid, "_mut_by_prod")
  scid_AlleleStat_cols <- paste0(remained_scid, "_AlleleStat")
  
  all_cols <- c(base_cols, scid_cols, scid_unmut_cols, scid_BB_cols, scid_prod_cols, scid_AlleleStat_cols)
  
  # Keep only columns that actually exist
  cols_to_keep <- intersect(all_cols, colnames(all_merged_data))
  
  merged_data <- all_merged_data[, cols_to_keep]
  
  # Print information
  print(str_c(
    "Raw cell number is: ", length(scid_data$scid_basedTree), 
    "; After filtering cells that do not contain significant mutations, current cell number is: ", 
    length(remained_scid)
  ))
  
} else if(is_remove_cells == "no") {
  merged_data <- all_merged_data
  remained_scid <- scid_data$scid_basedTree
  print(str_c(
    "Raw cell number is: ", length(scid_data$scid_basedTree), 
    "; And no cells are removed!"
  ))
}

# Print the shape of merged_data
print(str_c(
  "The shape of merged_data is : ", 
  dim(merged_data)[1], " rows X ", 
  dim(merged_data)[2], " columns"
))
# [1] "The shape of merged_data is : 142 rows X 13284 columns"


##### Calculate average mutant cells allele frequency
### Step1: bulk features
cov_in_pseudobulk <- apply(merged_data, 1, function(row) {sum(as.numeric(sapply(strsplit(na.omit(unlist(row[paste0(remained_scid, "_AlleleStat")])), "/"), function(x) as.numeric(x[2]))))})
### Step2
avg_cov_percell <- apply(merged_data, 1, function(row) {mean_non_zero(as.numeric(sapply(strsplit(na.omit(unlist(row[paste0(remained_scid, "_AlleleStat")])), "/"), function(x) as.numeric(x[2]))))})
### => step3: Mutant cells are defined as those that contain at least one mutant allele.
# Detect if mutant allele is present and calculate VAF cross all cells
avg_mutAF_percell <- apply(merged_data, 1, function(row) {mean_non_zero(as.numeric(sapply(strsplit(na.omit(unlist(row[paste0(remained_scid, "_AlleleStat")])), "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))))})
max_mutAF_percell <- apply(merged_data, 1, function(row) {max(as.numeric(sapply(strsplit(na.omit(unlist(row[paste0(remained_scid, "_AlleleStat")])), "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))))})
# The coverage of the max mutant AF cell
cov_in_maxmutAFcell <- apply(merged_data, 1, function(row) {
    alleles <- na.omit(unlist(row[paste0(remained_scid, "_AlleleStat")]))
    if (length(alleles) == 0) {
        return(NA)  # Return NA if all AlleleStat values are NA
    }
    max_ratio <- max(sapply(strsplit(alleles, "/"), function(x) as.numeric(x[1])/as.numeric(x[2])))
    max_value_index <- which(sapply(strsplit(alleles, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]) == max_ratio), arr.ind = TRUE)[1]
    max_value_after_slash <- as.numeric(strsplit(alleles[max_value_index], "/")[[1]][2])
    return(max_value_after_slash)
})
# Count the maximum allele number cross all cells
total_mutAllele_num <- apply(merged_data, 1, function(row) {sum(as.numeric(sapply(strsplit(na.omit(unlist(row[paste0(remained_scid, "_AlleleStat")])), "/"), function(x) as.numeric(x[1]))))})
max_mutAllele_num <- apply(merged_data, 1, function(row) {max(as.numeric(sapply(strsplit(na.omit(unlist(row[paste0(remained_scid, "_AlleleStat")])), "/"), function(x) as.numeric(x[1]))))})
mutAllele_cellnum <- apply(merged_data, 1, function(row) {sum(as.numeric(sapply(strsplit(na.omit(unlist(row[paste0(remained_scid, "_AlleleStat")])), "/"), function(x) as.numeric(x[1])))>0)})
### => Step4: identify mutant cell by somatic_posterior_percell posterior more than a set value
# Count shared cell number by cell posterior value
moreCutoff_cellnum <- apply(merged_data, 1, function(row) {sum(as.numeric(sapply(na.omit(unlist(row[remained_scid])), function(x) as.numeric(x)>cutoff)))})
### => step5: Consider the number of cells with cutoff and with or without mutant allele
mutant_cellnum <- apply(merged_data, 1, function(row) {
    consider_scid <- remained_scid[!is.na(unlist(row[paste0(remained_scid, "_AlleleStat")])) & !is.na(unlist(row[remained_scid]))]
    sum(
        sapply(strsplit(na.omit(unlist(row[paste0(consider_scid, "_AlleleStat")])), "/"), function(x) as.numeric(x[1])) > 0 &
        sapply(na.omit(unlist(row[consider_scid])), function(x) as.numeric(x)) >= cutoff
    )
})
mutant_cell_fraction <- mutant_cellnum/length(remained_scid)
# out to dataframe
df_calAF_data <- as.data.frame(cbind(merged_data[,1:8], cov_in_pseudobulk, avg_cov_percell, avg_mutAF_percell, max_mutAF_percell, cov_in_maxmutAFcell, max_mutAllele_num, total_mutAllele_num, mutant_cell_fraction, mutant_cellnum, mutAllele_cellnum, moreCutoff_cellnum, merged_data[9:dim(merged_data)[2]])); dim(df_calAF_data)
# [1]   142 13294
filter_calAF_data <- df_calAF_data[df_calAF_data$somatic_posterior_persite>=cutoff, ]; dim(filter_calAF_data)
# [1]   142 13294
filter_calAF_data <- filter_calAF_data[filter_calAF_data$total_mutAllele_num!=0, ]; dim(filter_calAF_data)
# [1]   142 13294
# First output: Smmary_data (including only present in one cells)
write.table(filter_calAF_data, str_c(outputpath, "/Summary_data_", as.character(dim(filter_calAF_data)[1]), ".filtered_and_calculatedAF_and_mutantAnno.txt"), sep="\t", quote=FALSE, col.names=T, row.names=F)


##### Extract somatic_posterior_persite and posterior file from raw calculated results
final_data <- filter_calAF_data; dim(final_data)
# [1]   142 13294
# Following output: (except only present in one cells)

### 1. posterior
output.posterior_data <- as.data.frame(lapply(final_data[, remained_scid], as.numeric)); dim(output.posterior_data)
# [1]  142 3319
rownames(output.posterior_data) <- final_data$mutid
colnames(output.posterior_data) <- remained_scid
write.table(output.posterior_data, str_c(outputpath, "/data.posterior_matrix.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)

### 2. likelihood_unmut
output.likelihood_unmut_data <- as.data.frame(lapply(final_data[, paste0(remained_scid, "_unmut")], as.numeric)); dim(output.likelihood_unmut_data)
# [1]  142 3319
rownames(output.likelihood_unmut_data) <- final_data$mutid
colnames(output.likelihood_unmut_data) <- remained_scid
write.table(output.likelihood_unmut_data, str_c(outputpath, "/data.likelihood_unmut_matrix.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)

### 3. likelihood_mut_by_BB
output.likelihood_mut_by_BB_data <- as.data.frame(lapply(final_data[, paste0(remained_scid, "_mut")], as.numeric)); dim(output.likelihood_mut_by_BB_data)
# [1]  142 3319
rownames(output.likelihood_mut_by_BB_data) <- final_data$mutid
colnames(output.likelihood_mut_by_BB_data) <- remained_scid
write.table(output.likelihood_mut_by_BB_data, str_c(outputpath, "/data.likelihood_mut_matrix.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)

### 5. allele_info
output.allele_data <- final_data[, paste0(c("bulk", remained_scid), "_AlleleStat")]; dim(output.allele_data)
# [1]  142 3320
rownames(output.allele_data) <- final_data$mutid
colnames(output.allele_data) <- c("bulk", remained_scid)
write.table(output.allele_data, str_c(outputpath, "/data.allele_count.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)

### 6. features: somatic_posterior_persite + average mutant cell AF + avg_cov_percell
output.features <- as.data.frame(lapply(final_data[, c("somatic_posterior_persite", "cov_in_pseudobulk", "avg_cov_percell", "avg_mutAF_percell", "max_mutAF_percell", "cov_in_maxmutAFcell", "max_mutAllele_num", "total_mutAllele_num", "mutant_cell_fraction", "mutant_cellnum", "mutAllele_cellnum", "moreCutoff_cellnum")], as.numeric)); dim(output.features)
# [1] 142  11
output.features <- as.data.frame(cbind(mutid=final_data[, "mutid"], output.features)); dim(output.features)
# [1] 142  12
output.features <- merge(final_data[, c("mutid", "indid", "chr", "pos", "ref", "mut")], output.features, by="mutid", all=FALSE); dim(output.features)
# [1] 142  18
rownames(output.features) <- output.features$mutid
output.features <- output.features[, -1]; dim(output.features)
# [1] 142  17
write.table(output.features, str_c(outputpath, "/features.preprocess_items.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)

### 7.sequencing data
mutant_cellIndex <- apply(final_data, 1, function(row) {which(na.omit(sapply(strsplit(unlist(row[paste0(remained_scid, "_AlleleStat")]), "/"), function(x) as.numeric(x[1]))>0 & sapply(na.omit(unlist(row[remained_scid])), function(x) as.numeric(x))>cutoff))})
fa_df <- data.frame(cellid=remained_scid)
for(pos in 1:length(mutant_cellIndex)){
    fa_vec_pos <- rep(final_data[pos,"ref"], length(remained_scid))
    cellIndex <- mutant_cellIndex[[pos]]
    fa_vec_pos[cellIndex] <- final_data[pos,"mut"]
    fa_df <- as.data.frame(cbind(fa_df, fa_vec_pos))
}
fa_df_only <- fa_df[,(1:length(mutant_cellIndex))+1]; dim(fa_df_only)
# [1] 22 76
rownames(fa_df_only) <- fa_df[,1]
colnames(fa_df_only) <- final_data[,"mutid"]
fa_str_list <- apply(fa_df_only, 1, function(x) {paste(x, collapse="")})
output.fasta <- data.frame(content=str_c(length(remained_scid), " ", length(mutant_cellIndex)))
for(i in 1:length(fa_str_list)){
    content_i <- str_c(names(fa_str_list[i]), " ", fa_str_list[i])
    output.fasta <- as.data.frame(rbind(output.fasta, content_i))
}
write.table(output.fasta, str_c(outputpath, "/data.allcells_fasta.phy"), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)


######################### plot #########################
##### posterior percell distributuion
plot_df <- as.data.frame(unlist(output.posterior_data)); dim(plot_df)
colnames(plot_df) <- "value"
plot_df <- na.omit(plot_df)
plot_df_stat <- as.data.frame(table(plot_df)); dim(plot_df_stat)
# [1] 8629    2
ratio <- sum(plot_df_stat[as.numeric(as.character(plot_df_stat$value))>0.25&as.numeric(as.character(plot_df_stat$value))<0.75, "Freq"])/sum(plot_df_stat$Freq)
# [1] 0.08128342
ratio <- sprintf("%.2f%%", ratio * 100)
p1 <- ggplot(plot_df, aes(x=value)) +
    geom_histogram(bins=50, fill="steelblue", color="white") +
    labs(x="Value", y="Frequency", title="Histplot.posterior_percellpersite") +
    theme_bw() + 
    theme(plot.title=element_text(hjust=0.5, size=15)) +
    geom_text(x=0.5, y=500, label=str_c("0.25~0.75: ", as.character(ratio)), hjust=0, vjust=-1, size=5)
ggsave(str_c(outputpath, "/Histplot.somatic_posterior_percell.pdf"), p1, width=6, height=4)

##### Generate plot_features dataframe.
plot_features <- as.data.frame(cbind(mutid=rownames(output.features), output.features)); dim(plot_features)
plot_features <- na.omit(plot_features)

##### somatic posterior persite distributuion
p2 <- ggplot(plot_features, aes(x=somatic_posterior_persite)) +
    geom_histogram(bins=50, fill="steelblue", color="white") +
    labs(x="Value", y="Frequency", title="Histplot.somatic_posterior_persite") +
    theme_bw() + 
    theme(plot.title=element_text(hjust=0.5, size=15))
ggsave(str_c(outputpath, "/Histplot.somatic_posterior_persite.pdf"), p2, width=6, height=4)

##### average coverage distribution percell
p3 <- ggplot(plot_features, aes(x=avg_cov_percell)) +
    geom_histogram(bins=50, fill="steelblue", color="white") +
    labs(x="Value", y="Frequency", title="Histplot.avg_cov_percell") +
    theme_bw() + 
    theme(plot.title=element_text(hjust=0.5, size=15))
ggsave(str_c(outputpath, "/Histplot.avg_cov_percell.pdf"), p3, width=6, height=4)

##### mutant AF distributuion
p4 <- ggplot(plot_features, aes(x=avg_mutAF_percell)) +
    geom_histogram(bins=50, fill="steelblue", color="white") +
    labs(x="Value", y="Frequency", title="Histplot.avg_mutAF_percell") +
    theme_bw() + 
    theme(plot.title=element_text(hjust=0.5, size=15))
ggsave(str_c(outputpath, "/Histplot.avg_mutAF_percell.pdf"), p4, width=6, height=4)
p5 <- ggplot(plot_features, aes(x=max_mutAF_percell)) +
    geom_histogram(bins=50, fill="steelblue", color="white") +
    labs(x="Value", y="Frequency", title="Histplot.max_mutAF_percell") +
    theme_bw() + 
    theme(plot.title=element_text(hjust=0.5, size=15))
ggsave(str_c(outputpath, "/Histplot.max_mutAF_percell.pdf"), p5, width=6, height=4)

##### mutant cell number distributuion
p6 <- ggplot(plot_features, aes(x=mutAllele_cellnum)) +
    geom_histogram(bins=50, fill="steelblue", color="white") +
    labs(x="Value", y="Frequency", title="Histplot.mutAllele_cellnum") +
    theme_bw() + 
    theme(plot.title=element_text(hjust=0.5, size=15))
ggsave(str_c(outputpath, "/Histplot.mutAllele_cellnum.pdf"), p6, width=6, height=4)

##### mutant cell number barplot
mutAllele_cellnum_count <- table(plot_features$mutAllele_cellnum)/length(plot_features$mutAllele_cellnum)
mutAllele_cellnum_count_order <- mutAllele_cellnum_count[order(as.numeric(names(mutAllele_cellnum_count)), decreasing=T)]
mutAllele_cellnum_rate <- sprintf("%.2f%%", mutAllele_cellnum_count_order * 100)
x_position=c()
for(mut_cn in (unique(plot_features$mutAllele_cellnum)[order(unique(plot_features$mutAllele_cellnum), decreasing=T)])){
    position_list <- plot_features[plot_features$mutAllele_cellnum==mut_cn, "mutid"]
    if(mut_cn==max(unique(plot_features$mutAllele_cellnum)[order(unique(plot_features$mutAllele_cellnum), decreasing=T)])){
        position <- position_list[length(position_list)]
    } else if(mut_cn==min(unique(plot_features$mutAllele_cellnum)[order(unique(plot_features$mutAllele_cellnum), decreasing=T)])){
        position <- position_list[1]
    } else {
        position <- position_list[ceiling(length(position_list)/2)]
    }
    x_position <- c(x_position, position)
}
p7 <- ggplot(plot_features, aes(x=reorder(mutid, -mutAllele_cellnum), y=mutAllele_cellnum)) +
    geom_bar(stat="identity", color="steelblue") +
    labs(x="Mutations", y="The number of mutant cells", title="Histplot.mutAllele_cellnum") +
    theme_bw() + 
    theme(axis.ticks=element_blank(), axis.text.x=element_blank()) +
    theme(axis.text.y = element_text(size = 10)) +
    theme(axis.title = element_text(size = 10)) + 
    theme(plot.title=element_text(hjust=0.5, size=15)) + 
    annotate("text", x=x_position, y=(unique(plot_features$mutAllele_cellnum)[order(unique(plot_features$mutAllele_cellnum), decreasing=T)])+0.1, label=mutAllele_cellnum_rate, color="black", fontface="bold", size=3)
ggsave(str_c(outputpath, "/Barplot.mutAllele_cellnum.pdf"), p7, width=6, height=4)

# ##### bulk converage distribution
# p8 <- ggplot(plot_features, aes(x=bulk_cov)) +
#     geom_histogram(bins=50, fill="steelblue", color="white") +
#     labs(x="Value", y="Frequency", title="Histplot.bulk_cov") +
#     theme_bw() + 
#     theme(plot.title=element_text(hjust=0.5, size=15))
# ggsave(str_c(outputpath, "/Histplot.bulk_cov.pdf"), p8, width=6, height=4)

# ##### bulk VAF distribution
# p9 <- ggplot(plot_features, aes(x=bulk_VAF)) +
#     geom_histogram(bins=50, fill="steelblue", color="white") +
#     labs(x="Value", y="Frequency", title="Histplot.bulk_VAF") +
#     theme_bw() + 
#     theme(plot.title=element_text(hjust=0.5, size=15))
# ggsave(str_c(outputpath, "/Histplot.bulk_VAF.pdf"), p9, width=6, height=4)


##### Time #####
end_time <- Sys.time()
print(str_c("End time: ", end_time))
print(str_c("Program finished in ", as.character(round((end_time-start_time), 4)), " seconds"))


