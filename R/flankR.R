# Retrieve flanking sequences of rsID, upstream and downstream to
# create Duplex Oligos as ordered from IDT for LUCiferase assays which.
# will be used for Luciferase Assays.
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
library(tidyverse)
library(stringi)

# 1) Input SNPs of interest
## 1.1) Input list of query rsIDs
  ## Simple
### snp_list <- c("rs3", "rs4")  ## Used for simple tests
snp_list <- c("rs6466948", "rs4732060")  ## also used for simple tests

## or

## 1.2) Input list of query rsIDs for this project .
  ## Read list of SNP rsIDs from text file
snp_list <- readLines("data/input/rcc_snps_39.txt")

# 2) Input additional information (required)
## 2.1) Enter number of flanking sequences to retrieve
up_stream <- 10
down_stream <- 10

## 2.2) Enter additional overlapping sequences (required).
  ## For the "Winter" RCC project, we will be double digesting the
  ## Luciferase Vector (pGL4.23) with XhoI and BglII and cloning
  ## the duplex oligos (ordered from IDT) into the vector using the
  ## Takara In-Fusion Snap Assembly (Cat. # 638946) which requires
  ## a 15-bp homologous overlap with the vector ends to which it
  ## will be cloned.
overlap_seq_begin <- "GCTCGCTAGCCTCGA"
overlap_seq_last <- "GATCTGGCCTCGGCG"

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


# 6) Separate the 'Variant sequences' column into three new columns
snp_sequence <- snp_sequence %>%
  tidyr::separate(col = "Variant sequences", 
                  into = c("Upstream", "Alleles", "Downstream"), 
                  sep = "%"
                 )

# 7) Concatenate the upstream and downstream sequences with each allele
snp_sequence_expanded <- snp_sequence %>%
                         mutate(Alleles = str_split(Alleles, "/")) %>%
                         unnest(Alleles) %>%
                         mutate(Full_target_sequence = paste0(Upstream, Alleles, Downstream))

snp_sequence_expanded


# 8) These are the `Forward Orientation` sequences. Add column to label.
snp_sequence_expanded <- snp_sequence_expanded %>%
  mutate(`Orientation` = "Forward") 

# 9) Add two columns with 15-bp overlapping sequences (see notes in 2.2)
snp_sequence_expanded <- snp_sequence_expanded %>%
      mutate(`Overlap_beg` = overlap_seq_begin) %>%
      relocate(`Overlap_beg`, .before = Full_target_sequence) %>%
      mutate(`Overlap_last` = overlap_seq_last) %>%
      relocate(`Overlap_last`, .before = Orientation)


# 10) Duplicate rows and modify with `Reverse` Orientation of `Full_target_sequence`
snp_sequence_modified <- snp_sequence_expanded %>%
  # Ensure Orientation is set correctly before starting
  mutate(Orientation = "Forward") %>%
  # Create duplicate rows with modifications
  slice(rep(1:n(), each = 2)) %>%
  # Apply conditions to modify only the duplicated rows
  mutate(
    Full_target_sequence = if_else(row_number() %% 2 == 0, stringi::stri_reverse(Full_target_sequence), Full_target_sequence),
    Orientation = if_else(row_number() %% 2 == 0, "Reverse", "Forward")
  )


# 11) Concatenate 'Overlap_beg', 'Full_target_sequence', and
# 'Overlap_last' to get the 'Full_sequence_to_order' from IDT

snp_sequence_modified <- snp_sequence_modified %>%
    mutate(Full_sequence_to_order = paste0(Overlap_beg, Full_target_sequence, Overlap_last)) %>%
    relocate(`Full_sequence_to_order`, .before = Orientation)


# 12) Write dataframe "snp_sequence" to CSV file

# Get current date and time to use for unique filename
current_time <- format(Sys.time(), "%Y%m%d_%H%M%S")

# Create length of "snp_list" to use in "file_name" below
snp_list_length <- length(snp_list)

# Create a unique file name using the current date and time
file_name <- paste0("data/output/snp_sequence_", snp_list_length, "_", up_stream, "_", down_stream, "_", current_time, ".csv")

# Write the dataframe "snp_sequence" to a CSV file called "file_name"
write.csv(snp_sequence_expanded, file = file_name, row.names = FALSE)
