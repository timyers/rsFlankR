# Retrieve flanking sequences of rsID, upstream and downstream
# Reference: https://support.bioconductor.org/p/89688/
# Reference: https://chat.openai.com/share/7e2ec69e-a377-4564-a525-0c9e96647499
# NOTE: dataset = "hsapiens_snp" uses genome build GRCh38.p14, do 'listDatasets(snp_mart)

# install biomaRt
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")

# Load libraries
library(biomaRt)
library(dplyr)


# 1.1) Input list of query rsIDs
# Simple
snp_list <- c("rs3", "rs4")

# 1.2) Input list of query rsIDs for this project.
# Read list of SNP rsIDs from text file
snp_list <- readLines("data/input/rcc_snps_39.txt")

# 2) Input number of flanking sequences to retrieve
up_stream <- 10
down_stream <- 10

# 3) Select database
snp_mart<- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")

# 4) Retrieve sequence and other data
snp_sequence <- getBM(attributes = c("refsnp_id", "snp", "allele", "chr_name", "chrom_start", 
                                     "chrom_end", "chrom_strand"), 
                      filters = c("snp_filter", "upstream_flank", "downstream_flank"), 
                      checkFilters = FALSE, 
                      values = list(snp_list, up_stream, down_stream), 
                      mart = snp_mart, 
                      bmHeader = TRUE)

snp_sequence


# 5) Add new column with genome build, GRCh38.p14 (see notes above).

snp_sequence <- snp_sequence %>%
  mutate(`Genome Build` = "GRCh38.p14") %>%
  relocate(`Genome Build`, .before = Strand)


# 6) Write dataframe "snp_sequence" to CSV file

# Get current date and time to use for unique filename
current_time <- format(Sys.time(), "%Y%m%d_%H%M%S")

# Create length of "snp_list" to use in "file_name" below
snp_list_length <- length(snp_list)

# Create a unique file name using the current date and time
file_name <- paste0("data/output/snp_sequence_", snp_list_length, "_", up_stream, "_", down_stream, "_", current_time, ".csv")

# Write the dataframe "snp_sequence" to a CSV file called "file_name"
write.csv(snp_sequence, file = file_name, row.names = FALSE)
