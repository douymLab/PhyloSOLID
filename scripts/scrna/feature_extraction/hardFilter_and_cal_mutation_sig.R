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
    make_option(c("-i", "--inputfile"), type = "character", default = "/storage/douyanmeiLab/yangqing/tools/PhyloMosaicGenie/Benchmark/data_10xCI792/labeling/identifier.features.CI792_tumor.txt", help = "the name of inputfile - posterior raw data & reads information"),
    make_option(c("-o", "--outputfile"), type = "character", default = "/storage/douyanmeiLab/yangqing/tools/PhyloMosaicGenie/Benchmark/data_10xCI792/labeling/CI792.identifier.feature.sigFilter.txt", help = "The outputfile you can set")
)
parseobj = OptionParser(option_list=option_list)
opt = parse_args(parseobj)


inputfile <- as.character(opt$inputfile)
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
get_default_labels <- function(choice) {
  if (!(choice %in% c("DNA", "RNA", "96", "192"))) {
    stop("Choice must be 'DNA'/'96' or 'RNA'/'192'")
  }
  if (choice == "DNA" | choice == "96") {
    mid_list <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
    value <- 96
  } else {
    mid_list <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G", "A>C", "A>G", "A>T", "G>A", "G>C", "G>T")
    value <- 192
  }
  first <- c("A", "T", "C", "G")
  inner_bracket <- rep(rep(mid_list, each = 16), times = 1)
  outter_bracket <- expand.grid(first, first)
  result <- sapply(1:value, function(f) {
    paste0(outter_bracket[f %% 16 + 1, 1], "[", inner_bracket[f], "]", outter_bracket[f %% 16 + 1, 2])
  })
  return(result)
}


#=============handle data=============

# read features
df_candidate <- read.csv(inputfile,header = T,sep="\t")
if (!("identifier" %in% colnames(df_candidate))) {
  colnames(df_candidate)[colnames(df_candidate) == "X.identifier"] <- "identifier"
}

counts <- table(df_candidate$RNAMutationType)
sig <- as.data.frame(counts)
colnames(sig)<-c("MutationType","Count")

default_list=get_default_labels("RNA")
default_df=data.frame("MutationType" = default_list,"none" = rep(0,length((default_list))))
df_sigProfile<-left_join(default_df,sig,by="MutationType")

df_sigProfile$Count <- ifelse(is.na(df_sigProfile$Count), 0, df_sigProfile$Count)
df_sigProfile["none"]<-NULL
dim(df_sigProfile)
# [1] 192   2


##### test paired mutation type
# Define base pairing rules
complement_pairs <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C")

# Corrected get_pair function
get_pair <- function(mutation) {
  # Extract the complete MutationType
  # Extract prefix, suffix, and middle parts
  prefix <- substr(mutation, 1, 1) # first base
  suffix <- substr(mutation, nchar(mutation), nchar(mutation)) # last base
  middle <- gsub(".*\\[|\\].*", "", mutation) # middle part (e.g. C>A)

  # Decompose prefix, suffix, and middle parts
  prefix_comp <- complement_pairs[prefix] # complementary prefix
  suffix_comp <- complement_pairs[suffix] # complementary suffix
  X <- substr(middle, 1, 1) # first base of the middle part
  Y <- substr(middle, nchar(middle), nchar(middle)) # second base of the middle part
  X_comp <- complement_pairs[X] # complementary first base of the middle
  Y_comp <- complement_pairs[Y] # complementary second base of the middle

  # Check whether all complementary bases exist
  if (!is.na(prefix_comp) && !is.na(suffix_comp) &&
      !is.na(X_comp) && !is.na(Y_comp)) {
    # Generate the paired MutationType
    pair_mutation <- paste0(prefix_comp, "[", X_comp, ">", Y_comp, "]", suffix_comp)
    return(pair_mutation)
  } else {
    return(NA) # If any pairing fails, return NA
  }
}

# Generate pairs for all MutationType values
df_sigProfile$Pair <- sapply(df_sigProfile$MutationType, get_pair)

# Create a unified pair name (order-independent pairing)
get_unique_pair <- function(mutation, pair) {
  # Ensure consistent order of MutationType and Pair; use the lexicographically smaller as the key
  sorted_pair <- sort(c(mutation, pair))
  return(paste(sorted_pair, collapse = " + "))
}

# Compute unique pairs
UniquePair <- mapply(get_unique_pair, df_sigProfile$MutationType, df_sigProfile$Pair)

# Assign a sig_pair name to each unique pair
unique_pairs <- unique(UniquePair) # find all unique pairs
sig_pair_names <- paste0("pair", seq_along(unique_pairs)) # name each unique pair
pair_map <- setNames(sig_pair_names, unique_pairs) # create pair mapping

df_sigProfile$sig_pair <- pair_map[UniquePair]

# Use sapply to compute binomial_test results and return two columns
results <- t(sapply(seq_len(nrow(df_sigProfile)), function(i) {
  # print(i)
  # Count of the current row and the paired row
  a <- df_sigProfile$Count[i]
  pair_idx <- which(df_sigProfile$MutationType == df_sigProfile$Pair[i])
  
  # If no pair data exists, return NA and "no_pair"
  if (length(pair_idx) == 0) {
    return(c(greater_sig = "no_pair", p_value = NA))
  }
  
  b <- df_sigProfile$Count[pair_idx]
  
  # Determine greater_sig and the larger/smaller values
  if (a > b) {
    greater_sig <- df_sigProfile$MutationType[i]
    larger <- a
    smaller <- b
  } else if (a < b) {
    greater_sig <- df_sigProfile$Pair[i]
    larger <- b
    smaller <- a
  } else {
    greater_sig <- "equal"
    larger <- a
    smaller <- b
  }
  
  # Calculate p-value
  if (a == 0 & b == 0) {
    return(c(greater_sig = greater_sig, p_value = 1))
  } else {
    test <- binom.test(smaller, larger, p = 0.5)
    return(c(greater_sig = greater_sig, p_value = test$p.value))
  }
}))

# Save results to df_sigProfile
df_sigProfile$greater_sig <- results[, "greater_sig"]
df_sigProfile$sig_pvalue <- as.numeric(results[, "p_value"])
# write.table(df_sigProfile, str_c(outputpath, "/muts_SigProfile.txt"), row.names=FALSE, col.names=TRUE,quote=FALSE, sep="\t")


##### Add information into site features dataframe
dim(df_candidate)
# [1] 176 182
# 1. Rename MutationType column in df_sigProfile to RNAMutationType for matching
colnames(df_sigProfile)[colnames(df_sigProfile) == "MutationType"] <- "RNAMutationType"
dim(df_sigProfile)
# [1] 192   7
# 2. Add sig_pvalue column to df_features, matching RNAMutationType with df_sigProfile
df_features <- merge(df_candidate, 
                     df_sigProfile[, c("RNAMutationType", "sig_pvalue")], 
                     by = "RNAMutationType", 
                     all.x = TRUE)
dim(df_features)
# [1] 176 183
# 3. Add signature_filter column based on whether RNAMutationType matches greater_sig
df_removed <- df_sigProfile[
  df_sigProfile$RNAMutationType == df_sigProfile$greater_sig & df_sigProfile$sig_pvalue < 0.01, 
]
dim(df_removed)
# [1] 4 6
print(str_c("The number of mutation types would be removed is: ", as.character(nrow(df_removed))))
removed_muttype_by_sig <- df_removed$RNAMutationType
# [1] "G[T>G]A" "G[T>G]G" "T[G>A]A" "G[G>A]A"
df_features$signature_filter <- ifelse(
  df_features$RNAMutationType %in% removed_muttype_by_sig, 
  "fail", 
  "pass"
)
dim(df_features)
# [1] 176 184


##### geneate bed format columns
# Load tidyr package
library(tidyr)
library(dplyr)

# Split identifier column and keep the original column
df_features <- df_features %>%
  mutate(original_identifier = identifier) %>%  # copy identifier column
  separate(
    col = identifier, 
    into = c("chrom", "position", "ref", "alt"), 
    sep = "_", 
    remove = FALSE  # keep original identifier column
  ) %>%
  mutate(
    start = as.numeric(position) - 1,  # second element minus 1 as start
    end = as.numeric(position)         # second element as end
  )

dim(df_features)


##### output and save
df_out <- df_features[, c('chrom', 'start', 'end', colnames(df_candidate), 'sig_pvalue', 'signature_filter')]; dim(df_out)
write.table(df_out, outputfile, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")


##### check results
# # filter
# df_filter <- df_features[df_features$signature_filter=="pass", ]; dim(df_filter)
# # [1]  63 185

# df_true <- df_filter[df_filter$label=='mosaic',]; dim(df_true)
# # [1]  14 185

# dim(df_candidate[df_candidate$label=='mosaic',])
# # [1]  33 182


# # check
# df_check <- df_features[df_features$label == "mosaic" & df_features$signature_filter == "pass",]
# dim(df_check)
# # [1] 14  5
# head(df_check[, c("label", "RNAMutationType", "sig_pvalue", "signature_filter")])
#     label RNAMutationType   sig_pvalue signature_filter
# 5  mosaic         A[C>T]T 5.335212e-04             pass
# 9  mosaic         A[G>A]C 1.000000e+00             pass
# 15 mosaic         C[C>T]T 1.967493e-11             pass
# 20 mosaic         G[A>G]A 1.000000e+00             pass
# 40 mosaic         G[C>T]T 6.250000e-02             pass
# 41 mosaic         G[C>T]T 6.250000e-02             pass


##### Time #####
end_time <- Sys.time()
print(str_c("End time: ", end_time))
print(str_c("Program finished in ", as.character(round((end_time-start_time), 4)), " seconds"))

