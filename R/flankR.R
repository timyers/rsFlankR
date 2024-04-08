# Retrieve flanking sequences of rsID, upstream and downstream
# Reference: https://support.bioconductor.org/p/89688/
# Reference: https://chat.openai.com/share/7e2ec69e-a377-4564-a525-0c9e96647499
# NOTE: dataset = "hsapiens_snp" uses GRCh38.p14, do 'listDatasets(snp_mart)

# install biomaRt
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")

library(biomaRt)

# 1.1) Input list of query rsIDs
# Simple
snp_list <- c("rs3", "rs4")

# 1.2) Input list of query rsIDs
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
