# Fig. 5c: GABAAR-like gene expression across adult Nematostella tissues
# Author: Ruohan Zhong

# Packages
library(pheatmap)

# Load TPM values
tpm <- read.table(
  "tpm_RefSeq_jaNemVect1.1.txt",
  header = TRUE,
  row.names = 1,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Average TPM across biological replicates for each tissue
groups <- sub("_.*", "", colnames(tpm))
tpm.avg <- t(apply(tpm, 1, function(x) tapply(as.numeric(x), groups, mean, na.rm = TRUE)))

# Order tissues: proximal tentacles, oral disc, body wall, pharynx
order_cols <- intersect(c("ttl", "ord", "bdw", "phx"), colnames(tpm.avg))
tpm.avg <- tpm.avg[, order_cols, drop = FALSE]

# Load curated GABAAR-like gene list
gaba_A_df <- read.table(
  "GABAAR_like_gene_list.txt",
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)
colnames(gaba_A_df) <- c("GeneID", "Description")

# Select GABAAR-like genes detected in the bulk RNA-seq dataset
valid_genes <- gaba_A_df$GeneID[gaba_A_df$GeneID %in% rownames(tpm.avg)]
data_subset <- tpm.avg[valid_genes, , drop = FALSE]
data_subset <- data_subset[
  apply(data_subset, 1, function(x) all(is.finite(x))),
  ,
  drop = FALSE
]

# Heatmap color palette
color_palette <- colorRampPalette(c("#d0cdff", "#ffffe8", "#f59b3c"))

# Generate heatmap
fig5c <- pheatmap(
  data_subset,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  scale = "row",
  color = color_palette(100),
  border_color = NA,
  annotation_legend = FALSE,
  fontsize_row = 6
)

fig5c

# Session information
sessionInfo()
