#####################################################################
## Sript to create oligos for LUC-iferase assay
#####################################################################

########## Information necessary ##################

region = "3q26" ## only for the name
snp_name = "rs60792873"     ## SNP name
alleles = c("A", "G") ## Which are the alleles for the SNP?
snp_10seq_begore = "AATCTTTTAA"  ## 10 bp before the SNP ## DO NOT INCLUDE THE ALLELE
snp_10seq_after = "TATTCATTTC" ## 10 bp after the SNP ## DO NOT INCLUDE THE ALLELE

######## Functions and Values #####################
inser_seq_beg = "ACTGGCCGGTACCTG"
inser_seq_last = "CGCTAGCCTCGAGGA"

rev.comp<-function(x,rev=TRUE) {
x<-toupper(x)
y<-rep("N",nchar(x))
xx<-unlist(strsplit(x,NULL))
for (bbb in 1:nchar(x)) {
		if(xx[bbb]=="A") y[bbb]<-"T"		
		if(xx[bbb]=="C") y[bbb]<-"G"		
		if(xx[bbb]=="G") y[bbb]<-"C"		
		if(xx[bbb]=="T") y[bbb]<-"A" }
if(rev==FALSE) {
	 	for(ccc in (1:nchar(x))) {
		if(ccc==1) yy<-y[ccc] else yy<-paste(yy,y[ccc],sep="") } }
if(rev==T) {
	zz<-rep(NA,nchar(x))
	for(ccc in (1:nchar(x))) {
		zz[ccc]<-y[nchar(x)+1-ccc]
		if(ccc==1) yy<-zz[ccc] else yy<-paste(yy,zz[ccc],sep="") }}
	return(yy)	}

rev.only <- function(x)
        sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")

########## Script: creating the sequences ############


if (length(alleles)==2) {
	
table = data.frame(matrix(vector(), 8, 2))
colnames(table) <- c("Name", "Sequence")

table$Name[1] <- paste0("LUC_", region, "_", snp_name, "_", alleles[1], "_", "F.F")
table$Sequence[1] <- paste0(inser_seq_beg, snp_10seq_begore, alleles[1], snp_10seq_after, inser_seq_last)
table$Name[2] <- paste0("LUC_", region, "_", snp_name, "_", alleles[1], "_", "F.R")
table$Sequence[2] <- rev.comp(paste0(inser_seq_beg, snp_10seq_begore, alleles[1], snp_10seq_after, inser_seq_last))

table$Name[3] <- paste0("LUC_", region, "_", snp_name, "_", alleles[2], "_", "F.F")
table$Sequence[3] <- paste0(inser_seq_beg, snp_10seq_begore, alleles[2], snp_10seq_after, inser_seq_last)
table$Name[4] <- paste0("LUC_", region, "_", snp_name, "_", alleles[2], "_", "F.R")
table$Sequence[4] <- rev.comp(paste0(inser_seq_beg, snp_10seq_begore, alleles[2], snp_10seq_after, inser_seq_last))

table$Name[5] <- paste0("LUC_", region, "_", snp_name, "_", alleles[1], "_", "R.F")
table$Sequence[5] <- paste0(inser_seq_beg, rev.only(paste0(snp_10seq_begore, alleles[1], snp_10seq_after)), inser_seq_last)
table$Name[6] <- paste0("LUC_", region, "_", snp_name, "_", alleles[1], "_", "R.R")
table$Sequence[6] <- rev.comp(paste0(inser_seq_beg, rev.only(paste0(snp_10seq_begore, alleles[1], snp_10seq_after)), inser_seq_last))

table$Name[7] <- paste0("LUC_", region, "_", snp_name, "_", alleles[2], "_", "R.F")
table$Sequence[7] <- paste0(inser_seq_beg, rev.only(paste0(snp_10seq_begore, alleles[2], snp_10seq_after)), inser_seq_last)
table$Name[8] <- paste0("LUC_", region, "_", snp_name, "_", alleles[2], "_", "R.R")
table$Sequence[8] <- rev.comp(paste0(inser_seq_beg, rev.only(paste0(snp_10seq_begore, alleles[2], snp_10seq_after)), inser_seq_last))

} else {

if (length(alleles)==3) {
	
table = data.frame(matrix(vector(), 12, 2))
colnames(table) <- c("Name", "Sequence")

table$Name[1] <- paste0("LUC_", region, "_", snp_name, "_", alleles[1], "_", "F.F")
table$Sequence[1] <- paste0(inser_seq_beg, snp_10seq_begore, alleles[1], snp_10seq_after, inser_seq_last)
table$Name[2] <- paste0("LUC_", region, "_", snp_name, "_", alleles[1], "_", "F.R")
table$Sequence[2] <- rev.comp(paste0(inser_seq_beg, snp_10seq_begore, alleles[1], snp_10seq_after, inser_seq_last))

table$Name[3] <- paste0("LUC_", region, "_", snp_name, "_", alleles[2], "_", "F.F")
table$Sequence[3] <- paste0(inser_seq_beg, snp_10seq_begore, alleles[2], snp_10seq_after, inser_seq_last)
table$Name[4] <- paste0("LUC_", region, "_", snp_name, "_", alleles[2], "_", "F.R")
table$Sequence[4] <- rev.comp(paste0(inser_seq_beg, snp_10seq_begore, alleles[2], snp_10seq_after, inser_seq_last))

table$Name[5] <- paste0("LUC_", region, "_", snp_name, "_", alleles[3], "_", "F.F")
table$Sequence[5] <- paste0(inser_seq_beg, snp_10seq_begore, alleles[3], snp_10seq_after, inser_seq_last)
table$Name[6] <- paste0("LUC_", region, "_", snp_name, "_", alleles[3], "_", "F.R")
table$Sequence[6] <- rev.comp(paste0(inser_seq_beg, snp_10seq_begore, alleles[3], snp_10seq_after, inser_seq_last))

table$Name[7] <- paste0("LUC_", region, "_", snp_name, "_", alleles[1], "_", "R.F")
table$Sequence[7] <- paste0(inser_seq_beg, rev.only(paste0(snp_10seq_begore, alleles[1], snp_10seq_after)), inser_seq_last)
table$Name[8] <- paste0("LUC_", region, "_", snp_name, "_", alleles[1], "_", "R.R")
table$Sequence[8] <- rev.comp(paste0(inser_seq_beg, rev.only(paste0(snp_10seq_begore, alleles[1], snp_10seq_after)), inser_seq_last))

table$Name[9] <- paste0("LUC_", region, "_", snp_name, "_", alleles[2], "_", "R.F")
table$Sequence[9] <- paste0(inser_seq_beg, rev.only(paste0(snp_10seq_begore, alleles[2], snp_10seq_after)), inser_seq_last)
table$Name[10] <- paste0("LUC_", region, "_", snp_name, "_", alleles[2], "_", "R.R")
table$Sequence[10] <- rev.comp(paste0(inser_seq_beg, rev.only(paste0(snp_10seq_begore, alleles[2], snp_10seq_after)), inser_seq_last))

table$Name[11] <- paste0("LUC_", region, "_", snp_name, "_", alleles[3], "_", "R.F")
table$Sequence[11] <- paste0(inser_seq_beg, rev.only(paste0(snp_10seq_begore, alleles[2], snp_10seq_after)), inser_seq_last)
table$Name[12] <- paste0("LUC_", region, "_", snp_name, "_", alleles[3], "_", "R.R")
table$Sequence[12] <- rev.comp(paste0(inser_seq_beg, rev.only(paste0(snp_10seq_begore, alleles[2], snp_10seq_after)), inser_seq_last))

} }

write.table(table, paste0("Oligos_LUC_pGL4.10_for_", snp_name), quote=F, col.names=T, row.names=F)


