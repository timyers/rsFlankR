library(Biostrings)

##### Reverse complement of a DNA sequence #####
# Given DNA sequence
dna_sequence <- "ACTGGCCGGTACCTGAATCTTTTAAATATTCATTTCCGCTAGCCTCGAGGA"

# Creating a DNAString object
dna_string <- DNAString(dna_sequence)

# Getting the reverse complement
reverse_complement <- reverseComplement(dna_string)

# Print the reverse complement
print(reverse_complement)



##### Reverse DNA sequence without taking the complement #####
# Given DNA sequence
dna_sequence <- "ACTGGCCGGTACCTGAATCTTTTAAATATTCATTTCCGCTAGCCTCGAGGA"

# Reverse the DNA sequence
reversed_dna_sequence <- sapply(strsplit(dna_sequence, ""), function(x) paste(rev(x), collapse = ""))

# Print the reversed DNA sequence
print(reversed_dna_sequence)
