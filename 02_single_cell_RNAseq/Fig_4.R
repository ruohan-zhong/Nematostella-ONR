# Fig. 4: Neural cluster analysis
# Author: Ruohan Zhong

# Packages
library(Seurat)
library(Matrix)
library(ggplot2)
library(readxl)
library(dplyr)
library(purrr)
library(tibble)
library(ggtern)
library(patchwork)

# Load Seurat object
nv <- readRDS("nv_oral_seurat.RDS")
DefaultAssay(nv) <- "SCT"
Idents(nv) <- "seurat_clusters"

# Define cluster order
custom_order <- c(
  "2", "3", "5", "6", "9", "13", "11", "12",
  "14", "0", "4", "8", "10", "7", "1", "15"
)

Idents(nv) <- factor(Idents(nv), levels = custom_order)

# Neural clusters used for Fig. 4a-c
cluster_ids <- c("2", "3", "5", "6", "9", "13", "11", "12", "14")

# Gene-list input folders
receptor_list_dir <- "Fig4_neurotransmitter_receptor_gene_lists"
sensory_list_dir <- "Fig4_sensory_gene_lists"

# Output folders for gene-list quantification
receptor_output_dir <- file.path(receptor_list_dir, "results")
sensory_output_dir <- file.path(sensory_list_dir, "results")

dir.create(receptor_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sensory_output_dir, recursive = TRUE, showWarnings = FALSE)

# Quantify genes across neuronal clusters
pct_cutoff <- 10
avg_cutoff <- 0

quantify_gene_lists <- function(input_dir, output_dir, seurat_obj, cluster_ids) {
  input_files <- list.files(
    input_dir,
    pattern = "\\.txt$",
    full.names = TRUE
  )
  
  cluster_cols <- setNames(
    rep(list(numeric()), length(cluster_ids)),
    paste0("Cluster_", cluster_ids)
  )
  
  summary_df <- data.frame(
    GeneList = character(),
    cluster_cols,
    Total = numeric(),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  sct_genes <- rownames(
    LayerData(seurat_obj[["SCT"]], layer = "counts")
  )
  
  for (input_file in input_files) {
    gene_table <- tryCatch(
      read.table(
        input_file,
        header = FALSE,
        sep = "\t",
        col.names = c("GeneID", "Description"),
        colClasses = c("character", "character")
      ),
      error = function(e) NULL
    )
    
    if (is.null(gene_table)) next
    
    gene_table$GeneID <- trimws(gene_table$GeneID)
    gene_list_name <- sub("\\.txt$", "", basename(input_file))
    valid_genes <- gene_table$GeneID[gene_table$GeneID %in% sct_genes]
    
    cluster_counts <- setNames(
      rep(0, length(cluster_ids)),
      paste0("Cluster_", cluster_ids)
    )
    
    if (length(valid_genes) > 0) {
      all_results <- list()
      
      for (cl in cluster_ids) {
        cluster_cells <- subset(seurat_obj, idents = cl)
        expr_matrix <- LayerData(cluster_cells[["SCT"]], layer = "counts")
        
        gene_subset <- expr_matrix[
          rownames(expr_matrix) %in% valid_genes,
          ,
          drop = FALSE
        ]
        
        if (nrow(gene_subset) > 0) {
          avg_expression <- Matrix::rowMeans(gene_subset)
          pct_expressed <- (Matrix::rowSums(gene_subset > 0) / ncol(gene_subset)) * 100
          
          pass <- (avg_expression > avg_cutoff) & (pct_expressed > pct_cutoff)
          
          cluster_counts[paste0("Cluster_", cl)] <- sum(pass)
          
          if (any(pass)) {
            results_df <- data.frame(
              Cluster = cl,
              Gene = rownames(gene_subset),
              Average_Expression = as.numeric(avg_expression),
              Percent_Expressed = as.numeric(pct_expressed),
              Pass_Threshold = pass,
              row.names = NULL,
              check.names = FALSE
            )
            
            results_df <- merge(
              results_df,
              gene_table,
              by.x = "Gene",
              by.y = "GeneID",
              all.x = TRUE,
              sort = FALSE
            )
            
            all_results[[cl]] <- subset(results_df, Pass_Threshold)[
              ,
              c("Cluster", "Gene", "Average_Expression", "Percent_Expressed", "Description")
            ]
          }
        }
      }
      
      if (length(all_results) > 0) {
        combined_results <- do.call(rbind, all_results)
        
        write.table(
          combined_results,
          file.path(output_dir, paste0(gene_list_name, "_results.txt")),
          sep = "\t",
          quote = FALSE,
          row.names = FALSE
        )
      }
    }
    
    summary_entry <- data.frame(
      GeneList = gene_list_name,
      t(cluster_counts),
      Total = sum(cluster_counts),
      check.names = FALSE
    )
    
    summary_df <- rbind(summary_df, summary_entry)
  }
  
  write.table(
    summary_df,
    file.path(output_dir, "all_clusters_summary.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  summary_df
}

# Run gene-list quantification
receptor_summary <- quantify_gene_lists(
  receptor_list_dir,
  receptor_output_dir,
  nv,
  cluster_ids
)

sensory_summary <- quantify_gene_lists(
  sensory_list_dir,
  sensory_output_dir,
  nv,
  cluster_ids
)

# Results from the gene-list analyses were manually combined into the following figure input files
receptors_xlsx <- "figure_inputs/Fig4_neurotransmitter_receptor_genes.xlsx"
sensory_xlsx <- "figure_inputs/Fig4_sensory_genes.xlsx"

# Build neurotransmitter receptor table used for Fig. 4a-c
map_category <- c(
  "fast-acting-ex" = "E",
  "fast-acting-in" = "I",
  "modulatory" = "M"
)

receptors_df <- read_excel(receptors_xlsx) |>
  mutate(
    Cluster = as.numeric(Cluster),
    Type = recode(Category, !!!map_category, .default = NA_character_)
  ) |>
  filter(!is.na(Type)) |>
  select(Cluster, Gene, Category, Type)

# Build sensory gene table
sensory_df <- read_excel(sensory_xlsx) |>
  mutate(
    Cluster = as.numeric(Cluster),
    Type = "Sensory"
  ) |>
  select(Cluster, Gene, Type)

# Fig. 4a - sensory and neurotransmitter receptor gene counts per cell

gene_sets <- bind_rows(
  receptors_df |> select(Cluster, Gene, Type),
  sensory_df
)

annotate_counts_one <- function(seurat_obj, cl, gene_sets) {
  obj_cl <- subset(seurat_obj, idents = as.character(cl))
  expr <- LayerData(obj_cl[["SCT"]], layer = "counts")
  types <- unique(gene_sets$Type)
  
  map_dfr(types, function(tp) {
    genes <- gene_sets %>%
      filter(Cluster == cl, Type == tp) %>%
      pull(Gene) %>%
      intersect(rownames(expr))
    
    if (length(genes) == 0) {
      tibble(
        cell_id = colnames(expr),
        cluster = as.character(cl),
        Type = tp,
        Count = 0L
      )
    } else {
      counts <- Matrix::colSums(expr[genes, , drop = FALSE] > 0)
      
      tibble(
        cell_id = colnames(expr),
        cluster = as.character(cl),
        Type = tp,
        Count = as.integer(counts)
      )
    }
  })
}

counts_long <- map_dfr(
  as.numeric(cluster_ids),
  ~ annotate_counts_one(nv, .x, gene_sets)
) %>%
  mutate(
    cluster = factor(cluster, levels = cluster_ids),
    Type = factor(Type, levels = c("Sensory", "E", "I", "M"))
  )

type_colors <- c(
  "Sensory" = "#9577e5",
  "E" = "#fed976",
  "I" = "#FFA500",
  "M" = "#7d92b5"
)

fig4a_violin <- ggplot(
  counts_long,
  aes(x = cluster, y = Count, fill = Type)
) +
  facet_wrap(~Type, ncol = 1, scales = "free_y") +
  geom_violin(
    trim = FALSE,
    scale = "width",
    color = "black",
    alpha = 0.85
  ) +
  geom_boxplot(
    width = 0.2,
    fill = "white",
    color = "black",
    outlier.shape = NA
  ) +
  scale_fill_manual(values = type_colors) +
  labs(
    x = "Cluster",
    y = "Number of expressed genes per cell by cluster"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    panel.grid.minor = element_blank()
  )

fig4a_violin

# Fig. 4b - cell counts and E/I/M receptor co-expression

# Cell count per cluster
counts_tbl <- as.data.frame(
  table(factor(Idents(nv), levels = cluster_ids))
)
colnames(counts_tbl) <- c("cluster", "count")

cluster_colors <- c(
  "2" = "#FFD97D",
  "3" = "#E5C3FF",
  "5" = "#FFAC28",
  "6" = "#FF7438",
  "9" = "#B918FF",
  "11" = "#FF00FF",
  "12" = "#FF9DD1",
  "13" = "#910054",
  "14" = "#FF6969"
)

fig4b_bar <- ggplot(
  counts_tbl,
  aes(x = cluster, y = count, fill = cluster)
) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(
    aes(label = count),
    vjust = -0.3,
    size = 3.5,
    color = "black"
  ) +
  scale_fill_manual(values = cluster_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Cluster", y = "Number of cells") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text = element_text(color = "grey30", size = 10),
    axis.title = element_text(size = 12),
    legend.position = "none"
  )

fig4b_bar

# E/I/M receptor combinations per cell
goi_eim <- receptors_df |>
  transmute(
    Cluster = as.character(Cluster),
    Gene,
    Category = Type
  )

combination_levels <- c(
  "None", "E", "I", "M", "EI", "EM", "IM", "EIM"
)

combination_colors <- c(
  "None" = "#EAEAEA",
  "E" = "#FED976",
  "I" = "#FFA500",
  "M" = "#7D92B5",
  "EI" = "#FF7438",
  "EM" = "#C4B0A3",
  "IM" = "#82786C",
  "EIM" = "#42392E"
)

get_comb_counts <- function(seurat_obj, cl, goi_df) {
  
  expr <- seurat_obj %>%
    subset(idents = as.character(cl)) %>%
    { LayerData(.[["SCT"]], layer = "counts") }
  
  pick <- function(cat) {
    intersect(
      goi_df %>%
        filter(Cluster == cl, Category == cat) %>%
        pull(Gene),
      rownames(expr)
    )
  }
  
  e_genes <- pick("E")
  i_genes <- pick("I")
  m_genes <- pick("M")
  
  e_any <- if (length(e_genes)) {
    Matrix::colSums(expr[e_genes, , drop = FALSE] > 0) > 0
  } else {
    rep(FALSE, ncol(expr))
  }
  
  i_any <- if (length(i_genes)) {
    Matrix::colSums(expr[i_genes, , drop = FALSE] > 0) > 0
  } else {
    rep(FALSE, ncol(expr))
  }
  
  m_any <- if (length(m_genes)) {
    Matrix::colSums(expr[m_genes, , drop = FALSE] > 0) > 0
  } else {
    rep(FALSE, ncol(expr))
  }
  
  comb <- paste0(
    as.integer(e_any),
    as.integer(i_any),
    as.integer(m_any)
  )
  
  comb <- recode(
    comb,
    "000" = "None",
    "100" = "E",
    "010" = "I",
    "001" = "M",
    "110" = "EI",
    "101" = "EM",
    "011" = "IM",
    "111" = "EIM"
  )
  
  as_tibble(
    table(factor(comb, levels = combination_levels)),
    .name_repair = "minimal"
  ) %>%
    rename(combination = 1, n = 2) %>%
    mutate(
      cluster = cl,
      percentage = as.numeric(n) / sum(n) * 100
    )
}

all_counts <- map_dfr(
  cluster_ids,
  ~ get_comb_counts(nv, .x, goi_eim)
) %>%
  mutate(
    cluster = factor(cluster, levels = cluster_ids),
    combination = factor(combination, levels = combination_levels)
  )

fig4b_sbp <- ggplot(
  all_counts,
  aes(x = cluster, y = percentage, fill = combination)
) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = combination_colors) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = "Cluster",
    y = "Percentage of cells",
    fill = "Combination"
  ) +
  theme_minimal(base_size = 12)

fig4b_sbp

# Combine Fig. 4b panels
fig4b_combined <- fig4b_bar / fig4b_sbp +
  plot_layout(heights = c(0.2, 1))

fig4b_combined

# Fig. 4c - ternary plot of E/I/M receptor counts per cell

pick_genes <- function(goi_df, cl, type_code) {
  goi_df %>%
    dplyr::filter(Cluster == as.numeric(cl), Type == type_code) %>%
    dplyr::pull(Gene)
}

annotate_counts <- function(seurat_obj, cluster_id, goi_df) {
  expr <- seurat_obj %>%
    subset(idents = as.character(cluster_id)) %>%
    { LayerData(.[["SCT"]], layer = "counts") }
  
  E_genes <- intersect(pick_genes(goi_df, cluster_id, "E"), rownames(expr))
  I_genes <- intersect(pick_genes(goi_df, cluster_id, "I"), rownames(expr))
  M_genes <- intersect(pick_genes(goi_df, cluster_id, "M"), rownames(expr))
  
  E_count <- if (length(E_genes)) Matrix::colSums(expr[E_genes, , drop = FALSE] > 0) else rep.int(0L, ncol(expr))
  I_count <- if (length(I_genes)) Matrix::colSums(expr[I_genes, , drop = FALSE] > 0) else rep.int(0L, ncol(expr))
  M_count <- if (length(M_genes)) Matrix::colSums(expr[M_genes, , drop = FALSE] > 0) else rep.int(0L, ncol(expr))
  
  tibble::tibble(
    cell_id = colnames(expr),
    E_count = as.integer(E_count),
    I_count = as.integer(I_count),
    M_count = as.integer(M_count),
    cluster = paste0("Cluster_", cluster_id)
  )
}

all_cells <- purrr::map_dfr(cluster_ids, ~ annotate_counts(nv, .x, receptors_df)) %>%
  dplyr::mutate(cluster = factor(cluster, levels = paste0("Cluster_", cluster_ids)))

tern_colors <- c("E"="#f9c758", "I"="#f99900", "M"="#7ea6ce")

fig4c_tern <- ggtern::ggtern(all_cells, aes(x = E_count, y = I_count, z = M_count)) +
  geom_point(alpha = 0.3, color = "black", size = 1.2) +
  facet_wrap(~ cluster, ncol = 3) +
  labs(x = "Ex", y = "In", z = "Md") +
  theme_bw() +
  theme(
    panel.background = element_blank(),
    plot.background  = element_blank(),
    strip.text       = element_text(face = "bold", size = 12),
    strip.background = element_blank(),
    tern.axis.title.L = element_text(color = tern_colors["E"], size = 12, face = "bold"),
    tern.axis.title.T = element_text(color = tern_colors["I"], size = 12, face = "bold"),
    tern.axis.title.R = element_text(color = tern_colors["M"], size = 12, face = "bold"),
    tern.axis.line.L  = element_line(color = tern_colors["E"], linewidth = 0.8),
    tern.axis.line.T  = element_line(color = tern_colors["I"], linewidth = 0.8),
    tern.axis.line.R  = element_line(color = tern_colors["M"], linewidth = 0.8),
    tern.panel.grid.major.L = element_line(color = "gray80", linewidth = 0.4),
    tern.panel.grid.major.T = element_line(color = "gray80", linewidth = 0.4),
    tern.panel.grid.major.R = element_line(color = "gray80", linewidth = 0.4),
    tern.panel.grid.minor = element_blank(),
    tern.axis.ticks.length.major = grid::unit(0.5, "cm"),
    tern.axis.ticks.major.L = element_line(color = "gray80", linewidth = 1),
    tern.axis.ticks.major.T = element_line(color = "gray80", linewidth = 1),
    tern.axis.ticks.major.R = element_line(color = "gray80", linewidth = 1),
    tern.axis.ticks.minor = element_blank(),
    tern.axis.text.L = element_text(color = "black", size = 13),
    tern.axis.text.T = element_text(color = "black", size = 13),
    tern.axis.text.R = element_text(color = "black", size = 13)
  )

fig4c_tern

# Fig. 4d - clusters 5 and/or 6 enriched markers

c56_dot_xlsx <- "figure_inputs/Fig4d_cluster56_markers.xlsx"

c56_dot <- read_excel(c56_dot_xlsx) %>%
  dplyr::filter(gene_id %in% rownames(LayerData(nv[["SCT"]], layer = "data"))) %>%
  dplyr::mutate(gene_id = factor(gene_id, levels = gene_id))

gene_labels <- setNames(c56_dot$description, c56_dot$gene_id)

fig4d_dot <- DotPlot(nv, features = levels(c56_dot$gene_id), assay = "SCT") +
  coord_flip() +
  scale_x_discrete(
    limits = rev(levels(c56_dot$gene_id)),
    labels = gene_labels,
    position = "bottom"
  ) +
  scale_y_discrete(
    limits = custom_order,
    position = "right"
  ) +
  labs(x = NULL, y = NULL) +
  scale_color_viridis_c(option = "magma", direction = -1) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major = element_line(color = "lightgrey", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(hjust = 0.5),
    axis.text.y = element_text(size = 9, hjust = 1),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10)
  )

fig4d_dot

# Session information
sessionInfo()
