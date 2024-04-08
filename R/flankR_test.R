# Retrieve flanking sequences of rsID, 20 bp upstream and downstream
# Reference: https://support.bioconductor.org/p/89688/
# Reference: https://chat.openai.com/share/7e2ec69e-a377-4564-a525-0c9e96647499
# NOTE: dataset = "hsapiens_snp" uses GRCh38.p14, do 'listDatasets(snp_mart)

library(biomaRt)

snp_list <- c("rs25", "rs16944")
snp_list <- c("rs3", "rs4")

# Select database
snp_mart<- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")

# Retrieve sequence and other data
snp_sequence <- getBM(attributes = c("refsnp_id", "snp"), 
                      filters = c("snp_filter", "upstream_flank", "downstream_flank"), 
                      checkFilters = FALSE, 
                      values = list(snp_list, 20, 20), 
                      mart = snp_mart, 
                      bmHeader = TRUE)

snp_sequence
