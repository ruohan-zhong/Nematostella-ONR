# Fig. S6: Molecular features of clusters 11, 12, and 14
# Author: Ruohan Zhong

# Packages
library(Seurat)
library(readxl)
library(ggplot2)
library(dplyr)
library(viridis)
library(pheatmap)
library(viridisLite)

# Load Seurat object
nv <- readRDS("nv_oral_seurat.RDS")
DefaultAssay(nv) <- "SCT"
Idents(nv) <- "seurat_clusters"

dat <- LayerData(nv[["SCT"]], layer = "data")

# Define cluster order
custom_order <- c(
  "2", "3", "5", "6", "9", "13", "11", "12",
  "14", "0", "4", "8", "10", "7", "1", "15"
)

Idents(nv) <- factor(Idents(nv), levels = custom_order)


# Fig. S6a - selected markers for clusters 11, 12, and 14

dot_data <- read_excel(
  "figure_inputs/FigS6a_cluster11_12_14_markers.xlsx"
)

process_cluster <- function(cluster_df, cluster_id, seurat_obj) {
  # Filter genes present in Seurat object
  genes_present <- cluster_df$gene_id %in% rownames(seurat_obj)
  cluster_df <- cluster_df[genes_present, ]
  
  # Calculate percentage of cells expressing each gene
  cluster_cells <- WhichCells(
    seurat_obj,
    idents = as.character(cluster_id)
  )
  
  expr_matrix <- FetchData(
    seurat_obj,
    vars = cluster_df$gene_id,
    cells = cluster_cells
  )
  
  cluster_df$pct_expressed <- colMeans(expr_matrix > 0) * 100
  
  # Sort genes by percentage expressed
  cluster_df %>%
    arrange(desc(pct_expressed))
}

c11_dot <- process_cluster(
  filter(dot_data, cluster == 11),
  11,
  nv
)

c12_dot <- process_cluster(
  filter(dot_data, cluster == 12),
  12,
  nv
)

c14_dot <- process_cluster(
  filter(dot_data, cluster == 14),
  14,
  nv
)

combined_dot <- bind_rows(
  c11_dot,
  c12_dot,
  c14_dot
)

combined_dot$gene_order <- factor(
  combined_dot$gene_id,
  levels = rev(
    c(
      c11_dot$gene_id,
      c12_dot$gene_id,
      c14_dot$gene_id
    )
  )
)

figs6a_dot <- DotPlot(
  nv,
  features = combined_dot$gene_id
) +
  coord_flip() +
  scale_x_discrete(
    limits = levels(combined_dot$gene_order),
    labels = setNames(
      combined_dot$description,
      combined_dot$gene_id
    ),
    position = "bottom"
  ) +
  scale_y_discrete(position = "right") +
  labs(x = NULL, y = NULL) +
  ggtitle("Combined Gene Expression") +
  scale_color_viridis_c(
    option = "magma",
    direction = -1,
    limits = c(-1, 3)
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1
    ),
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5
    ),
    axis.text.y = element_text(
      size = 8,
      angle = 0,
      hjust = 1
    ),
    legend.key.size = unit(0.6, "cm")
  )

figs6a_dot


# Fig. S6b-d - receptor and sensory gene expression in clusters 11, 12, and 14

# Annotation colors
annotation_colors <- list(
  category = c(
    "fast-acting-ex" = "#fed976",
    "fast-acting-in" = "#FFA500",
    "modulatory receptor" = "#95b0d1",
    "sensory conduction" = "#9577e5"
  )
)

# Generate heatmap for each cluster
plot_cluster_heatmap <- function(gene_file, cell_ids, cluster_label) {
  
  # Read gene list
  goi <- read_excel(
    gene_file,
    col_names = TRUE
  )
  
  # Extract expression matrix
  mat_raw <- as.matrix(
    dat[goi$Gene, cell_ids, drop = FALSE]
  )
  
  # Calculate expression and abundance per gene
  ExpSum <- rowSums(mat_raw)
  Abundance <- rowSums(mat_raw > 0) / length(cell_ids)
  
  # Build gene metadata
  df_info <- data.frame(
    gene = rownames(mat_raw),
    category = goi$Category,
    ExpSum = ExpSum,
    Abundance = Abundance,
    stringsAsFactors = FALSE
  )
  
  # Sort genes within category by expression and abundance
  df_info_sorted <- df_info %>%
    group_by(category) %>%
    arrange(
      category,
      desc(ExpSum),
      desc(Abundance),
      .by_group = TRUE
    ) %>%
    ungroup()
  
  mat_sorted <- mat_raw[
    df_info_sorted$gene,
    ,
    drop = FALSE
  ]
  
  # Row annotation
  annotationRow <- data.frame(
    category = df_info_sorted$category
  )
  
  rownames(annotationRow) <- df_info_sorted$gene
  
  # Generate heatmap
  heatmap_out <- pheatmap(
    mat_sorted,
    annotation_row = annotationRow,
    annotation_colors = annotation_colors,
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    show_colnames = FALSE,
    treeheight_col = 20,
    fontsize_row = 8,
    color = viridisLite::magma(
      100,
      direction = -1
    ),
    breaks = seq(
      0,
      3,
      length.out = 100
    ),
    border_color = NA,
    main = paste(
      "Cluster",
      cluster_label
    )
  )
  
  return(heatmap_out)
}

# Fig. S6b - cluster 11
go_list_11 <- "figure_inputs/FigS6b_cluster11_heatmap_genes.xlsx"
ids_11 <- WhichCells(nv, idents = "11")

figs6b_heatmap <- plot_cluster_heatmap(
  go_list_11,
  ids_11,
  "11"
)

figs6b_heatmap

# Fig. S6c - cluster 12
go_list_12 <- "figure_inputs/FigS6c_cluster12_heatmap_genes.xlsx"
ids_12 <- WhichCells(nv, idents = "12")

figs6c_heatmap <- plot_cluster_heatmap(
  go_list_12,
  ids_12,
  "12"
)

figs6c_heatmap

# Fig. S6d - cluster 14
go_list_14 <- "figure_inputs/FigS6d_cluster14_heatmap_genes.xlsx"
ids_14 <- WhichCells(nv, idents = "14")

figs6d_heatmap <- plot_cluster_heatmap(
  go_list_14,
  ids_14,
  "14"
)

figs6d_heatmap

# Session information
sessionInfo()
