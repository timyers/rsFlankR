library(dplyr)
library(tidyr)

snp_sequence <- snp_sequence %>%
  tidyr::separate(col = "Variant sequences", into = c("Upstream", "Alleles", "Downstream"), sep = "%")
  
tidyr::separate(col = "Alleles", into = c("Allele 1", "Alleles", "Downstream"), sep = "%")


# View the updated data frame
print(snp_sequence)


library(tidyverse)

# Assuming your dataframe is named df
df_expanded <- snp_sequence %>%
  mutate(Alleles = str_split(Alleles, "/")) %>%
  unnest(Alleles) %>%
  mutate(Full_sequence = paste0(Upstream, Alleles, Downstream))

df_expanded



