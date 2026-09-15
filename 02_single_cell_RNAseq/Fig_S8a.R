# Fig. S8a: GABA signaling component co-expression in clusters 5 and 6
# Author: Ruohan Zhong

# Packages
library(Seurat)
library(Matrix)
library(readxl)
library(dplyr)
library(ggplot2)

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


# Fig. S8a - GABAAR-like receptor and transporter co-expression

gabaar_genes <- goi_c56_all %>%
  filter(Category == "fast-acting-in") %>%
  pull(Gene)

gat_gene <- "LOC5521615"
vgat_gene <- "LOC5506413"

# Keep genes present in the SCT expression matrix
gabaar_genes <- intersect(gabaar_genes, rownames(dat))

# Extract expression status for each cell
get_marker_status <- function(cluster_id) {
  
  cluster_cells <- WhichCells(
    nv,
    idents = as.character(cluster_id)
  )
  
  expr_sub <- dat[
    ,
    cluster_cells,
    drop = FALSE
  ]
  
  gabaar_status <- Matrix::colSums(
    expr_sub[gabaar_genes, , drop = FALSE] > 0
  ) > 0
  
  gat_status <- as.vector(
    expr_sub[gat_gene, ] > 0
  )
  
  vgat_status <- as.vector(
    expr_sub[vgat_gene, ] > 0
  )
  
  data.frame(
    Cell = cluster_cells,
    GABAARs = as.integer(gabaar_status),
    GAT = as.integer(gat_status),
    VGAT = as.integer(vgat_status),
    Cluster = as.character(cluster_id),
    stringsAsFactors = FALSE
  )
}

marker_status_5 <- get_marker_status("5")
marker_status_6 <- get_marker_status("6")

marker_combined <- bind_rows(
  marker_status_5,
  marker_status_6
) %>%
  mutate(
    combo = case_when(
      GABAARs == 1 & GAT == 1 & VGAT == 1 ~ "GABAARs + VGAT + GAT",
      GABAARs == 1 & GAT == 0 & VGAT == 1 ~ "GABAARs + VGAT",
      GABAARs == 1 & GAT == 1 & VGAT == 0 ~ "GABAARs + GAT",
      GABAARs == 1 & GAT == 0 & VGAT == 0 ~ "GABAARs only",
      GABAARs == 0 & GAT == 1 & VGAT == 1 ~ "GAT + VGAT",
      GABAARs == 0 & GAT == 0 & VGAT == 1 ~ "VGAT only",
      GABAARs == 0 & GAT == 1 & VGAT == 0 ~ "GAT only",
      TRUE ~ "None"
    ),
    combo = factor(
      combo,
      levels = c(
        "None",
        "GAT only",
        "VGAT only",
        "GAT + VGAT",
        "GABAARs only",
        "GABAARs + GAT",
        "GABAARs + VGAT",
        "GABAARs + VGAT + GAT"
      )
    )
  )

# Calculate percentage of cells in each expression category
plot_data <- marker_combined %>%
  count(
    Cluster,
    combo,
    .drop = FALSE
  ) %>%
  group_by(Cluster) %>%
  mutate(
    percent = n / sum(n) * 100
  ) %>%
  ungroup()

# Shared annotation colors
combo_colors <- c(
  "None" = "#e5e5e5",
  "GAT only" = "#aecde0",
  "VGAT only" = "#79accc",
  "GAT + VGAT" = "#48677a",
  "GABAARs only" = "#ffa500",
  "GABAARs + GAT" = "#ffd27f",
  "GABAARs + VGAT" = "#a993c7",
  "GABAARs + VGAT + GAT" = "#54278F"
)

# Generate stacked bar plot
figs8a <- ggplot(
  plot_data,
  aes(
    x = Cluster,
    y = percent,
    fill = combo
  )
) +
  geom_col() +
  scale_fill_manual(
    values = combo_colors,
    drop = FALSE
  ) +
  labs(
    x = "Cluster",
    y = "Percentage of cells (%)",
    fill = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

figs8a

# Session information
sessionInfo()
