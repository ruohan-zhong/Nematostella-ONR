# scRNA-seq processing and Seurat object generation
# Author: Chris W. Seidel

# Packages
library(Seurat)

# Cell Ranger filtered feature-barcode matrix generated from the scRNA-seq data deposited in GEO under accession GSE316974
jaNemVect1.10x <- "path/to/filtered_feature_bc_matrix"
jaNemVect1.data <- Read10X(data.dir=jaNemVect1.10x)

# Create Seurat object with initial filtering
nv <- CreateSeuratObject(counts=jaNemVect1.data, project="nvOD", min.cells=3, min.features=200)

# NCBI feature table for Nematostella vectensis genome assembly jaNemVect1.1 (GCF_932526225.1)
anno <- read.delim(file="path/to/GCF_932526225.1_jaNemVect1.1_feature_table.txt", sep="\t")

# Identify mitochondrial genes
mito.id <- unique(anno$symbol[grep("mitoc", anno$name, ignore.case=T)])

# Select mitochondrial features with >500 counts
mt.iv <- which(rownames(jaNemVect1.data) %in% mito.id)
mt.counts <- rowSums(jaNemVect1.data[mt.iv,])
mt.features <- names(mt.counts[mt.counts > 500])

# Calculate mitochondrial transcript percentage
nv[["percent.mt"]] <- PercentageFeatureSet(nv, features=mt.features)

# Apply QC filtering
nv <- subset(nv, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & nCount_RNA < 15000 & percent.mt < 5)

# Normalize, cluster, and generate dimensionality reductions
nv <- SCTransform(nv, vars.to.regress = "percent.mt", verbose=FALSE)
nv <- RunPCA(nv, verbose = FALSE)
nv <- FindNeighbors(nv, dims = 1:15)
nv <- FindClusters(nv, resolution = 0.5)
nv <- RunUMAP(nv, dims = 1:15)
nv <- RunTSNE(object = nv, dims.use = 1:15, do.fast = TRUE)

rm(anno)

# Save Seurat object
saveRDS(nv, file="nv_oral_seurat.RDS")

# Session information
sessionInfo()
