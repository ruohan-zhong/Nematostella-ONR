# Fig. 6a and Fig. S5: Neurotransmitter receptor expression in clusters 5 and 6
# Author: Ruohan Zhong

# Packages
library(Seurat)
library(Matrix)
library(readxl)
library(dplyr)
library(pheatmap)
library(viridisLite)

# Load Seurat object
nv <- readRDS("nv_oral_seurat.RDS")
DefaultAssay(nv) <- "SCT"
Idents(nv) <- "seurat_clusters"

dat <- LayerData(nv[["SCT"]], layer = "data")

# Load shared cluster 5/6 gene table
c56_xlsx <- "figure_inputs/Fig6a_FigS5_cluster56_genes.xlsx"

goi_c56_all <- read_excel(
  c56_xlsx,
  col_names = TRUE
) %>%
  distinct(Gene, .keep_all = TRUE)

# Shared annotation colors
category_colors <- c(
  "fast-acting-ex" = "#fed976",
  "fast-acting-in" = "#FFA500",
  "modulatory" = "#95b0d1",
  "transporter" = "#9577e5"
)

cluster_colors <- c(
  "5" = "#ffac28",
  "6" = "#ff7438"
)


# Fig. S5 - cluster 5 and 6 E/I/M receptor heatmap

goi_s5 <- goi_c56_all %>%
  filter(Category %in% c(
    "fast-acting-ex",
    "fast-acting-in",
    "modulatory"
  ))

c5_ids <- WhichCells(nv, idents = "5")
c6_ids <- WhichCells(nv, idents = "6")

# Extract gene x cell matrices
expr_c5 <- dat[goi_s5$Gene, c5_ids]
expr_c6 <- dat[goi_s5$Gene, c6_ids]

# Sort cells within each cluster by total expression, then number of expressed genes
order_c5 <- order(
  colSums(expr_c5),
  colSums(expr_c5 > 0),
  decreasing = TRUE
)

order_c6 <- order(
  colSums(expr_c6),
  colSums(expr_c6 > 0),
  decreasing = TRUE
)

c5_ids_sorted <- c5_ids[order_c5]
c6_ids_sorted <- c6_ids[order_c6]
c56_ordered_ids <- c(c5_ids_sorted, c6_ids_sorted)

# Column annotations
annotationCol <- data.frame(
  Cluster = rep(
    c("5", "6"),
    times = c(length(c5_ids_sorted), length(c6_ids_sorted))
  )
)
rownames(annotationCol) <- c56_ordered_ids

# Calculate expression and specificity scores
df <- as.data.frame(
  dat[goi_s5$Gene, c56_ordered_ids]
) %>%
  mutate(
    ExpSum5 = rowSums(.[, 1:length(c5_ids_sorted)]),
    ExpSum6 = rowSums(.[, (length(c5_ids_sorted) + 1):ncol(.)]),
    SpecificityScore = ExpSum5 - ExpSum6,
    category = goi_s5$Category,
    gene = rownames(.)
  ) %>%
  group_by(category) %>%
  arrange(
    category,
    if_else(SpecificityScore >= 0, -ExpSum5, ExpSum6),
    desc(SpecificityScore),
    .by_group = TRUE
  ) %>%
  as.data.frame()

# Row annotations
rownames(df) <- df$gene

annotationRow <- data.frame(
  category = df$category
)
rownames(annotationRow) <- rownames(df)

annotation_colors <- list(
  category = category_colors,
  Cluster = cluster_colors
)

# Prepare matrix for heatmap
plot_data <- df %>%
  select(
    -starts_with("ExpSum"),
    -SpecificityScore,
    -category,
    -gene
  )

# Generate heatmap
figs5 <- pheatmap(
  plot_data,
  annotation_row = annotationRow,
  annotation_col = annotationCol,
  annotation_colors = annotation_colors,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  treeheight_col = 0,
  fontsize_row = 6,
  color = viridisLite::magma(100, direction = -1),
  breaks = seq(0, 3, length.out = 100),
  border_color = NA,
  gaps_col = length(c5_ids_sorted)
)

figs5


# Fig. 6a - inhibitory receptor and transporter heatmap

goi_fig6a <- goi_c56_all %>%
  filter(Category %in% c(
    "fast-acting-in",
    "transporter"
  ))

genes <- intersect(goi_fig6a$Gene, rownames(dat))
goi_fig6a <- goi_fig6a[
  match(genes, goi_fig6a$Gene),
  ,
  drop = FALSE
]

c5_ids <- WhichCells(nv, idents = "5")
c6_ids <- WhichCells(nv, idents = "6")

expr_c5 <- dat[genes, c5_ids, drop = FALSE]
expr_c6 <- dat[genes, c6_ids, drop = FALSE]

# Sort cells within each cluster by summed expression and detection rate
order_c5 <- order(
  Matrix::colSums(expr_c5),
  Matrix::colSums(expr_c5 > 0),
  decreasing = TRUE
)

order_c6 <- order(
  Matrix::colSums(expr_c6),
  Matrix::colSums(expr_c6 > 0),
  decreasing = TRUE
)

c5_ids_sorted <- c5_ids[order_c5]
c6_ids_sorted <- c6_ids[order_c6]
c56_ids_sorted <- c(c5_ids_sorted, c6_ids_sorted)

mat <- cbind(
  expr_c5[, order_c5, drop = FALSE],
  expr_c6[, order_c6, drop = FALSE]
)

ExpSum5 <- Matrix::rowSums(
  mat[, seq_along(c5_ids_sorted), drop = FALSE]
)

ExpSum6 <- Matrix::rowSums(
  mat[, (length(c5_ids_sorted) + 1):ncol(mat), drop = FALSE]
)

SpecificityScore <- ExpSum5 - ExpSum6

# Sort genes within each category by cluster-specific abundance and specificity
row_df <- data.frame(
  gene = genes,
  category = factor(
    goi_fig6a$Category,
    levels = c("fast-acting-in", "transporter")
  ),
  ExpSum5 = ExpSum5,
  ExpSum6 = ExpSum6,
  SpecificityScore = SpecificityScore,
  stringsAsFactors = FALSE
) %>%
  arrange(
    category,
    if_else(SpecificityScore >= 0, -ExpSum5, ExpSum6),
    desc(SpecificityScore)
  )

mat <- mat[row_df$gene, , drop = FALSE]

# Column annotations
annotationCol <- data.frame(
  Cluster = rep(
    c("5", "6"),
    times = c(length(c5_ids_sorted), length(c6_ids_sorted))
  )
)
rownames(annotationCol) <- c56_ids_sorted

# Row annotations
annotationRow <- data.frame(
  category = row_df$category
)
rownames(annotationRow) <- row_df$gene

annotation_colors <- list(
  category = category_colors,
  Cluster = cluster_colors
)

magma_colors <- viridisLite::magma(100, direction = -1)

# Generate heatmap
fig6a <- pheatmap(
  mat,
  annotation_row = annotationRow,
  annotation_col = annotationCol,
  annotation_colors = annotation_colors,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  show_rownames = TRUE,
  fontsize_row = 6,
  color = magma_colors,
  border_color = NA,
  treeheight_col = 0,
  gaps_col = length(c5_ids_sorted),
  width = 10,
  height = 2.65
)

fig6a

# Session information
sessionInfo()
