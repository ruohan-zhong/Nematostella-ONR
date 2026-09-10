# Fig. 3b-g: Single-cell RNA-seq analysis
# Author: Ruohan Zhong

# Packages
library(Seurat)
library(ggplot2)
library(grid)
library(patchwork)

# Load Seurat object
nv <- readRDS("nv_oral_seurat.RDS")
DefaultAssay(nv) <- "SCT"
Idents(nv) <- "seurat_clusters"

custom_order <- c(
  "2", "3", "5", "6", "9", "13", "11", "12",
  "14", "0", "4", "8", "10", "7", "1", "15"
)

Idents(nv) <- factor(Idents(nv), levels = custom_order)

cluster_colors <- c(
  "0" = "#E0F78C",
  "4" = "#8ECC00",
  "8" = "#00A607",
  "10" = "#005700",
  "7" = "#2EC0E6",
  "1" = "#A5DDFF",
  "15" = "#1021FF",
  "3" = "#E5C3FF",
  "2" = "#FFD97D",
  "5" = "#FFAC28",
  "6" = "#FF7438",
  "11" = "#FF00FF",
  "9" = "#B918FF",
  "14" = "#FF6969",
  "12" = "#FF9DD1",
  "13" = "#910054"
)

# Fig. 3b - UMAP
fig3b_umap <- DimPlot(
  nv,
  reduction = "umap",
  cols = cluster_colors
) +
  theme_void()

fig3b_umap

# Fig. 3c - Dot plot
figure3c_genes <- read.delim(
  "figure_inputs/Fig3c_cluster_markers.txt",
  header = FALSE,
  col.names = c("ID", "GeneName", "Note")
)

figure3c_genes <- subset(
  figure3c_genes,
  ID %in% rownames(nv[["SCT"]])
)

gene_labels <- setNames(
  figure3c_genes$GeneName,
  figure3c_genes$ID
)

fig3c_dot <- DotPlot(
  nv,
  features = figure3c_genes$ID,
  dot.min = 0,
  dot.scale = 6
) +
  coord_flip() +
  scale_x_discrete(
    limits = rev(figure3c_genes$ID),
    labels = gene_labels,
    position = "bottom"
  ) +
  scale_y_discrete(position = "right") +
  labs(x = NULL, y = NULL) +
  geom_hline(
    yintercept = 9.5,
    color = "black",
    linewidth = 0.5
  ) +
  scale_color_viridis_c(
    option = "magma",
    direction = -1
  ) +
  theme(
    plot.title = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 1
    ),
    axis.text.x.top = element_text(
      angle = 0,
      hjust = 0.5
    ),
    axis.text.y.right = element_text(
      angle = 0,
      hjust = 0
    ),
    legend.key.size = unit(0.6, "cm"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10)
  )

fig3c_dot

# Fig. 3d-g - Feature plots
genes4 <- c(
  "LOC5502552",
  "LOC5517708",
  "LOC5520912",
  "LOC5515003"
)

fp_list <- FeaturePlot(
  object = nv,
  features = genes4,
  reduction = "umap",
  pt.size = 0.08,
  ncol = 4,
  order = TRUE,
  combine = FALSE
)

fig3b_features_4 <- wrap_plots(
  lapply(
    fp_list,
    \(p) p +
      scale_color_gradient(
        low = "lightgrey",
        high = "#FF8C00"
      ) +
      theme_classic()
  ),
  nrow = 1
)

fig3b_features_4

# Session information
sessionInfo()
