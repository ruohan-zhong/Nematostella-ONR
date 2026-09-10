# Bulk RNA-seq TPM processing for Fig. 5c
# Author: Chris W. Seidel

# Packages
library(edgeR)
library(GenomicFeatures)

# Directory containing STAR *_ReadsPerGene.out.tab files generated from the bulk RNA-seq data deposited in GEO under accession GSE316975
path <- "path/to/STAR_counts"
countfiles <- dir(path, pattern="_ReadsPerGene.out.tab")

# Read STAR gene-count files
counts <- c()
for( i in seq_along(countfiles) ){
  x <- read.table(file=paste0(path, "/", countfiles[i]), sep="\t", header=F, as.is=T)
  counts <- cbind(counts, x[,2])
}

rownames(counts) <- x[,1]
colnames(counts) <- sub("_ReadsPerGene.out.tab","", countfiles)

# Remove STAR alignment summary rows
iv <- grep("^N_", rownames(counts))
counts <- counts[-1*iv,]

# Define tissue groups
Group <- factor(gsub("_[0-9]$","",colnames(counts)))

# Filter lowly expressed genes
keep <- filterByExpr(counts, group=Group)
countData <- counts[keep,]

# RefSeq annotation for Nematostella vectensis genome assembly jaNemVect1.1 (GCF_932526225.1)
gtf <- "path/to/GCF_932526225.1_jaNemVect1.1_chr.gtf"

# Calculate gene lengths
txdb <- makeTxDbFromGFF(file=gtf, organism="Nematostella vectensis", dataSource="RefSeq")

allEx <- exonsBy(txdb, by="gene")
gbyx <- GenomicRanges::reduce(allEx)
geneLengths <- lapply(gbyx, function(x) sum(width(x)))
geneLengths <- unlist(geneLengths)

# Calculate TPM
counts2tpm <- function(counts,len) {
  x <- counts/len
  denom <- sum(x)/10^6
  return(x/denom)
}

iv <- match(rownames(countData), names(geneLengths))
geneLengths <- geneLengths[iv]

tpm <- apply(countData, 2, function(x) counts2tpm(x,geneLengths))

# Export TPM matrix
write.table(tpm, file="tpm_RefSeq_jaNemVect1.1.txt", sep="\t", col.names=NA)

# Session information
sessionInfo()
