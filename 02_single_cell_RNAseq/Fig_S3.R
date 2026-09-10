# Fig. S3: Neural gene expression across clusters
# Author: Ruohan Zhong

# Packages
library(Seurat)
library(ggplot2)
library(grid)

# Load Seurat object
nv <- readRDS("nv_oral_seurat.RDS")
DefaultAssay(nv) <- "SCT"
Idents(nv) <- "seurat_clusters"

custom_order <- c(
  "2", "3", "5", "6", "9", "13", "11", "12",
  "14", "0", "4", "8", "10", "7", "1", "15"
)

Idents(nv) <- factor(Idents(nv), levels = custom_order)

# Load and prepare gene list
read_gene_list_id_name <- function(file_path) {
  df <- read.delim(file_path, header = FALSE, stringsAsFactors = FALSE)
  stopifnot(ncol(df) >= 3)
  df <- df[, c(1, 3)]
  colnames(df) <- c("ID", "GeneName")
  df <- df[df$ID %in% rownames(nv[["SCT"]]), , drop = FALSE]
  df <- df[!duplicated(df$ID), , drop = FALSE]
  df
}

figs3a_genes <- read_gene_list_id_name("FigS3_neural_genes.txt")

figs3a_labels <- setNames(
  paste(figs3a_genes$GeneName, figs3a_genes$ID),
  figs3a_genes$ID
)

# Fig. S3 - Dot plot
figs3a_dot <- DotPlot(
  nv,
  features = figs3a_genes$ID,
  dot.min = 0,
  dot.scale = 6
) +
  coord_flip() +
  scale_x_discrete(
    limits = rev(figs3a_genes$ID),
    labels = figs3a_labels,
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

figs3a_dot

# Session information
sessionInfo()
