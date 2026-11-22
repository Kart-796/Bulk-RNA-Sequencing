# ==================== BULK RNA-SEQ: DESeq2 ANALYSIS ====================
# Input : <bulk_rna_workdir>/gene_counts.csv (alignment/03_count_reads.sh; columns gene_id,
#         gene_name, one column per sample) and <bulk_rna_workdir>/Metadata.csv (Sample, Condition)
# Output: <deseq2_workdir>/  QC plots, PCA, DE tables, volcano plots, Enrichr results
# ---- Configuration: edit config/config.yaml (run from the repository root) ----
source("utils/load_config.R")
counts_file   <- file.path(cfg_path("bulk_rna_workdir"), "gene_counts.csv")
metadata_file <- file.path(cfg_path("bulk_rna_workdir"), "Metadata.csv")
setwd(cfg_path("deseq2_workdir", create = TRUE))

# ==================Load Required Libraries ==============================================

library(tidyverse)
library(magrittr)
library(dplyr)
library(DESeq2)         
library(ggplot2)
library(ggrepel)
library(reshape2)
library(enrichR)
library(clusterProfiler)
library(msigdbr)
library(enrichplot)
library(gridExtra)
library(fgsea)
library(org.Hs.eg.db)
library(pheatmap)
library(preprocessCore)
library(naniar)
library(imputeLCMD)
library(plotly)
library(patchwork)
library(sva)
library(pcaMethods)
library(STRINGdb)
library(readr)
library(purrr)
library(stringr)
library(VennDiagram)
library(UpSetR)
library(ellipse)
library(combinat)
library(readxl)
library(ReactomePA)
library(DOSE)

# =====================Directory Setup===========================================

dir.create("Boxplots",          showWarnings = FALSE)
dir.create("histograms",        showWarnings = FALSE)
dir.create("violin_plots",      showWarnings = FALSE)
dir.create("volcano_plots",     showWarnings = FALSE)
dir.create("missing_variables", showWarnings = FALSE)
dir.create("PCA",               showWarnings = FALSE)
dir.create("Scree_Plots",       showWarnings = FALSE)
dir.create("Loadings_Plots",    showWarnings = FALSE)
dir.create("Data",              showWarnings = FALSE)
dir.create("PC_loadings",       showWarnings = FALSE)
dir.create("DE_analysis",       showWarnings = FALSE)
dir.create("GSEA",              showWarnings = FALSE)

# =================Load Data & Metadata===============================================

counts   <- read_csv(counts_file)
metadata <- read_csv(metadata_file)

# Sample columns: everything after gene_id and gene_name
sample_idx <- 3:ncol(counts)
stopifnot(identical(colnames(counts)[1:2], c("gene_id", "gene_name")))

# Ensure sample order in metadata matches the count columns
metadata <- metadata %>%
  mutate(Sample = as.character(Sample)) %>%
  slice(match(colnames(counts)[sample_idx], Sample))

# Parse two-factor structure from Condition
metadata <- metadata %>%
  mutate(
    Genotype  = str_extract(Condition, "MutB|MutA|Control"),
    Treatment = str_extract(Condition, "Vehicle|Treated"),
    Condition = factor(Condition, levels = c(
      "Control_Vehicle",  "Control_Treated",
      "MutA_Vehicle",  "MutA_Treated",
      "MutB_Vehicle",  "MutB_Treated"
    ))
  )

# ===================Data Cleaning=============================================
# Remove rows with all-zero mean across samples
counts_filtered <- counts[rowMeans(counts[, sample_idx], na.rm = TRUE) != 0, ]

# ======================Log2 Transformation==========================================
counts_filtered_log <- counts_filtered %>%
  mutate(across(all_of(sample_idx), ~ log2(. + 1)))
write.csv(counts_filtered_log, "Data/Log_transformed_data.csv", row.names = TRUE)

# Missing variable plot
p_miss <- gg_miss_var(counts_filtered_log[, sample_idx])
ggsave("missing_variables/gg_miss_var_logtrans.png", plot = p_miss,
       width = 8, height = 6, dpi = 300)

# Long format for QC plots — join with metadata for colouring
counts_filtered_log_long <- counts_filtered_log %>%
  pivot_longer(cols = all_of(sample_idx), names_to = "Sample", values_to = "Expression") %>%
  left_join(metadata, by = "Sample")

# ===============Colour Palette (used throughout)=================================================

condition_colors <- c(
  "Control_Vehicle"   = "#666666",
  "Control_Treated" = "#CCCCCC",
  "MutA_Vehicle"   = "#000066",
  "MutA_Treated" = "#0033CC",
  "MutB_Vehicle"   = "#660066",
  "MutB_Treated" = "#CC6699"
)

# Shared theme for all plots
theme_custom <- theme_minimal() +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line        = element_line(color = "black"),
    axis.ticks       = element_line(color = "black"),
    legend.position  = "right"
  )

# ======================QC Plots: Boxplot, Violin, Histogram==========================================

# Boxplot
box_log <- ggplot(counts_filtered_log_long,
                  aes(x = Sample, y = Expression, fill = Condition)) +
  geom_boxplot() +
  scale_fill_manual(values = condition_colors) +
  labs(title = "Gene Expression per Sample (Log2)",
       x = "Sample", y = "Expression (Log2)") +
  theme_custom
ggsave("Boxplots/boxplots_log_transformed_data.png",
       box_log, width = 14, height = 6, dpi = 300)

# Violin plot
violin_log <- ggplot(counts_filtered_log_long,
                     aes(x = Sample, y = Expression, fill = Condition)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = condition_colors) +
  labs(title = "Gene Expression Distribution per Sample (Log2)",
       x = "Sample", y = "Expression (Log2)") +
  theme_custom
ggsave("violin_plots/violin_log_transformed_data.png",
       violin_log, width = 14, height = 6, dpi = 300)

# Histogram
hist_log <- ggplot(counts_filtered_log_long,
                   aes(x = Expression, fill = Condition)) +
  geom_histogram(color = "#e9ecef", alpha = 0.6,
                 position = "identity", bins = 30) +
  scale_fill_manual(values = condition_colors) +
  labs(title = "Histogram of Expression per Sample (Log2)",
       x = "Expression (Log2)", y = "Frequency") +
  theme_custom +
  theme(axis.text.x = element_text(angle = 0))
ggsave("histograms/histogram_log_transformed_data.png",
       hist_log, width = 12, height = 6, dpi = 300)

# =======================DESeq2: Build DESeqDataSet=========================================

raw_counts_matrix <- as.matrix(counts_filtered[, sample_idx])
rownames(raw_counts_matrix) <- counts_filtered[[1]]   # Gene IDs in column 1
storage.mode(raw_counts_matrix) <- "integer"           # DESeq2 requires integer counts

col_data <- as.data.frame(metadata)
rownames(col_data) <- col_data$Sample

# Full-model design using the 6-level Condition factor (no intercept equivalent
# is achieved by ~ 0 + Condition; DESeq2 uses intercept models, so we use
# ~ Condition with relevel or handle contrasts manually via results())
dds <- DESeqDataSetFromMatrix(
  countData = raw_counts_matrix,
  colData   = col_data,
  design    = ~ Condition
)

# ---- Filter low-expressed genes (mirror edgeR filterByExpr logic) ----
keep <- rowSums(counts(dds) >= 10) >= min(table(metadata$Condition))
dds  <- dds[keep, ]
cat("Genes retained after filtering:", sum(keep), "\n")

# ---- Run DESeq2 ----
dds <- DESeq(dds, parallel = FALSE)

# ---- Extract VST-normalized counts for QC & PCA ----
vst_data <- vst(dds, blind = TRUE)
vst_mat  <- assay(vst_data)
write.csv(as.data.frame(vst_mat),
          "Data/DESeq2_VST_normalized.csv", row.names = TRUE)

# ---- QC Boxplot of VST-normalized counts ----
vst_df <- as.data.frame(vst_mat) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Sample", values_to = "VST") %>%
  left_join(metadata, by = "Sample")

box_vst <- ggplot(vst_df, aes(x = Sample, y = VST, fill = Condition)) +
  geom_boxplot() +
  scale_fill_manual(values = condition_colors) +
  labs(title = "DESeq2 VST-Normalized Distribution",
       x = "Sample", y = "VST Expression") +
  theme_custom
ggsave("Boxplots/boxplots_VST_normalized.png",
       box_vst, width = 14, height = 6, dpi = 300)

# ---- QC Violin Plot of VST-normalized counts ----
violin_vst <- ggplot(vst_df, aes(x = Sample, y = VST, fill = Condition)) +
  geom_violin(trim = FALSE, alpha = 0.85) +
  geom_boxplot(width = 0.08, fill = "white",
               outlier.size = 0.5, outlier.alpha = 0.4,
               color = "gray30") +
  scale_fill_manual(values = condition_colors) +
  labs(title = "DESeq2 VST-Normalized Distribution (Violin)",
       x = "Sample", y = "VST Expression") +
  theme_custom +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("violin_plots/violin_VST_normalized.png",
       violin_vst, width = 14, height = 6, dpi = 300)

# ---- QC Histogram of VST-normalized counts ----
hist_vst <- ggplot(vst_df, aes(x = VST, fill = Condition)) +
  geom_histogram(bins = 60, alpha = 0.7, color = NA,
                 position = "identity") +
  scale_fill_manual(values = condition_colors) +
  facet_wrap(~ Sample, scales = "free_y") +
  labs(title = "DESeq2 VST-Normalized Distribution (Histogram)",
       x = "VST Expression", y = "Gene Count") +
  theme_custom +
  theme(strip.text    = element_text(size = 7),
        axis.text.x   = element_text(size = 6),
        panel.spacing = unit(0.3, "lines"))
ggsave("histograms/histogram_VST_normalized.png",
       hist_vst, width = 18, height = 12, dpi = 300)

# ===========================PCA Analysis (on VST matrix)=====================================

pca_input <- t(vst_mat)   # samples as rows, genes as columns
pca_res   <- prcomp(pca_input, center = TRUE, scale. = TRUE)

# Variance explained
var_explained <- summary(pca_res)$importance["Proportion of Variance", ]
pct_var       <- round(var_explained * 100, 2)

# PCA scores with metadata
pca_scores <- as.data.frame(pca_res$x) %>%
  rownames_to_column("Sample") %>%
  left_join(metadata, by = "Sample")

# Centroids per Condition
centroids_12 <- pca_scores %>%
  group_by(Condition) %>%
  summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop")

centroids_13 <- pca_scores %>%
  group_by(Condition) %>%
  summarise(PC1 = mean(PC1), PC3 = mean(PC3), .groups = "drop")

# ---- Scree Plot ----
n_pcs    <- min(18, length(pct_var))
scree_df <- data.frame(
  PC         = factor(paste0("PC", 1:n_pcs), levels = paste0("PC", 1:n_pcs)),
  Variance   = pct_var[1:n_pcs],
  Cumulative = cumsum(pct_var[1:n_pcs])
)

scree_plot <- ggplot(scree_df, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "#377EB8",
           color = "black", linewidth = 0.3) +
  geom_line(aes(y = Cumulative, group = 1),
            color = "#E41A1C", linewidth = 0.8, linetype = "dashed") +
  geom_point(aes(y = Cumulative), color = "#E41A1C", size = 2) +
  geom_hline(yintercept = 80, linetype = "dotted", color = "grey40") +
  annotate("text", x = 2, y = 82, label = "80% threshold",
           size = 3, color = "grey40") +
  labs(title = "Scree Plot",
       x = "Principal Component", y = "Variance Explained (%)") +
  theme_custom
ggsave("Scree_Plots/Scree_plot.png", scree_plot, width = 10, height = 6, dpi = 300)

# ---- PCA PC1 vs PC2 ----
pca_pc1_pc2 <- ggplot(pca_scores,
                      aes(x = PC1, y = PC2, color = Condition,
                          shape = Condition, label = Sample)) +
  stat_ellipse(aes(fill = Condition), type = "norm",
               geom = "polygon", alpha = 0.10,
               level = 0.95, show.legend = FALSE) +
  geom_segment(
    data = pca_scores %>%
      left_join(centroids_12, by = "Condition",
                suffix = c("", "_centroid")),
    aes(xend = PC1_centroid, yend = PC2_centroid),
    alpha = 0.25, linewidth = 0.4, show.legend = FALSE
  ) +
  geom_point(size = 3.5, alpha = 0.85) +
  geom_point(data = centroids_12,
             aes(x = PC1, y = PC2, color = Condition),
             inherit.aes = FALSE,
             shape = 23, size = 5, fill = "white",
             stroke = 1.5, show.legend = FALSE) +
  geom_text_repel(size = 2.8, max.overlaps = 20,
                  box.padding = 0.4, segment.color = "grey60",
                  show.legend = FALSE) +
  scale_color_manual(values = condition_colors) +
  scale_fill_manual(values  = condition_colors) +
  scale_shape_manual(values = c(16, 17, 15, 18, 8, 7)) +
  labs(
    title = "PCA — PC1 vs PC2",
    x     = paste0("PC1 (", pct_var[1], "% variance)"),
    y     = paste0("PC2 (", pct_var[2], "% variance)"),
    color = "Condition", shape = "Condition"
  ) +
  theme_custom +
  theme(axis.text.x = element_text(angle = 0))
ggsave("PCA/PCA_PC1_vs_PC2.png", pca_pc1_pc2, width = 11, height = 8, dpi = 300)

# ---- PCA PC1 vs PC3 ----
pca_pc1_pc3 <- ggplot(pca_scores,
                      aes(x = PC1, y = PC3, color = Condition,
                          shape = Condition, label = Sample)) +
  stat_ellipse(aes(fill = Condition), type = "norm",
               geom = "polygon", alpha = 0.10,
               level = 0.95, show.legend = FALSE) +
  geom_segment(
    data = pca_scores %>%
      left_join(centroids_13, by = "Condition",
                suffix = c("", "_centroid")),
    aes(xend = PC1_centroid, yend = PC3_centroid),
    alpha = 0.25, linewidth = 0.4, show.legend = FALSE
  ) +
  geom_point(size = 3.5, alpha = 0.85) +
  geom_point(data = centroids_13,
             aes(x = PC1, y = PC3, color = Condition),
             inherit.aes = FALSE,
             shape = 23, size = 5, fill = "white",
             stroke = 1.5, show.legend = FALSE) +
  geom_text_repel(size = 2.8, max.overlaps = 20,
                  box.padding = 0.4, segment.color = "grey60",
                  show.legend = FALSE) +
  scale_color_manual(values = condition_colors) +
  scale_fill_manual(values  = condition_colors) +
  scale_shape_manual(values = c(16, 17, 15, 18, 8, 7)) +
  labs(
    title = "PCA — PC1 vs PC3",
    x     = paste0("PC1 (", pct_var[1], "% variance)"),
    y     = paste0("PC3 (", pct_var[3], "% variance)"),
    color = "Condition", shape = "Condition"
  ) +
  theme_custom +
  theme(axis.text.x = element_text(angle = 0))
ggsave("PCA/PCA_PC1_vs_PC3.png", pca_pc1_pc3, width = 11, height = 8, dpi = 300)

# ---- Gene Loadings: PC1, PC2 & PC3 ----
loadings_df <- as.data.frame(pca_res$rotation) %>%
  rownames_to_column("gene_id") %>%
  select(gene_id, PC1, PC2, PC3) %>%
  left_join(counts_filtered %>% select(gene_id, gene_name),
            by = "gene_id")

write.csv(loadings_df, "PC_loadings/PCA_loadings_all.csv", row.names = FALSE)
write.csv(pca_scores,  "PCA/PCA_scores_with_metadata.csv",  row.names = FALSE)

plot_loadings <- function(loadings_df, pc_col, n = 30, save_path) {
  df <- loadings_df %>%
    arrange(desc(abs(.data[[pc_col]]))) %>%
    slice_head(n = n) %>%
    mutate(
      gene_name = factor(gene_name,
                         levels = gene_name[order(.data[[pc_col]])]),
      Direction = ifelse(.data[[pc_col]] > 0, "Positive", "Negative")
    )
  
  p <- ggplot(df, aes(x = gene_name, y = .data[[pc_col]], fill = Direction)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
    scale_fill_manual(values = c("Positive" = "#E41A1C",
                                 "Negative" = "#377EB8")) +
    coord_flip() +
    labs(title = paste("Top", n, "Gene Loadings —", pc_col),
         x = "Gene", y = "Loading Score") +
    theme_custom +
    theme(axis.text.x = element_text(angle = 0))
  
  ggsave(save_path, p, width = 8, height = 10, dpi = 300)
}

plot_loadings(loadings_df, "PC1", 30, "PC_loadings/Loadings_PC1.png")
plot_loadings(loadings_df, "PC2", 30, "PC_loadings/Loadings_PC2.png")
plot_loadings(loadings_df, "PC3", 30, "PC_loadings/Loadings_PC3.png")

# ======================DESeq2: Define Contrasts & Extract DEGs==========================================

# Helper: get numeric contrast vector from a named list of +1/-1 weights
make_numeric_contrast <- function(dds, weights) {
  rn  <- resultsNames(dds)
  vec <- setNames(rep(0, length(rn)), rn)
  for (nm in names(weights)) {
    if (!nm %in% rn) stop("Name '", nm, "' not found in resultsNames(dds).\nAvailable: ", paste(rn, collapse=", "))
    vec[nm] <- weights[[nm]]
  }
  vec
}

# ---- Simple contrasts: c("Condition", numerator, denominator) ----
# (This list was referenced but not defined in the original script. These are the
#  pairwise comparisons used in the PyDESeq2 notebook; add or remove as needed.)
simple_contrasts <- list(
  Control_Treated_vs_Vehicle = c("Condition", "Control_Treated", "Control_Vehicle"),
  MutA_Treated_vs_Vehicle    = c("Condition", "MutA_Treated",    "MutA_Vehicle"),
  MutB_Treated_vs_Vehicle    = c("Condition", "MutB_Treated",    "MutB_Vehicle"),
  MutA_vs_Control_Vehicle    = c("Condition", "MutA_Vehicle",    "Control_Vehicle"),
  MutB_vs_Control_Vehicle    = c("Condition", "MutB_Vehicle",    "Control_Vehicle"),
  MutA_vs_MutB_Vehicle       = c("Condition", "MutA_Vehicle",    "MutB_Vehicle")
)

# ---- Interaction contrasts: numeric vectors on resultsNames ----
# resultsNames(dds) are:
# [1] "Intercept"
# [2] "Condition_Control_Treated_vs_Control_Vehicle"
# [3] "Condition_MutA_Vehicle_vs_Control_Vehicle"
# [4] "Condition_MutA_Treated_vs_Control_Vehicle"
# [5] "Condition_MutB_Vehicle_vs_Control_Vehicle"
# [6] "Condition_MutB_Treated_vs_Control_Vehicle"

interaction_contrasts <- list(
  
  # (MutB_Treated - MutB_Vehicle) - (Control_Treated - Control_Vehicle)
  # = Condition_MutB_Treated_vs_Control_Vehicle
  #   - Condition_MutB_Vehicle_vs_Control_Vehicle
  #   - Condition_Control_Treated_vs_Control_Vehicle
  Interaction_MutB_Treated = make_numeric_contrast(dds, list(
    "Condition_MutB_Treated_vs_Control_Vehicle"  =  1,
    "Condition_MutB_Vehicle_vs_Control_Vehicle"    = -1,
    "Condition_Control_Treated_vs_Control_Vehicle"  = -1
  )),
  
  # (MutA_Treated - MutA_Vehicle) - (Control_Treated - Control_Vehicle)
  # = Condition_MutA_Treated_vs_Control_Vehicle
  #   - Condition_MutA_Vehicle_vs_Control_Vehicle
  #   - Condition_Control_Treated_vs_Control_Vehicle
  Interaction_MutA_Treated = make_numeric_contrast(dds, list(
    "Condition_MutA_Treated_vs_Control_Vehicle"  =  1,
    "Condition_MutA_Vehicle_vs_Control_Vehicle"    = -1,
    "Condition_Control_Treated_vs_Control_Vehicle"  = -1
  )),
  
  # (MutA_Treated - MutA_Vehicle) - (Control_Treated - Control_Vehicle)
  #   - [(MutB_Treated - MutB_Vehicle) - (Control_Treated - Control_Vehicle)]
  # = Condition_MutA_Treated - Condition_MutA_Vehicle
  #   - Condition_MutB_Treated + Condition_MutB_Vehicle
  Interaction_MutA_Treated_Unique = make_numeric_contrast(dds, list(
    "Condition_MutA_Treated_vs_Control_Vehicle"  =  1,
    "Condition_MutA_Vehicle_vs_Control_Vehicle"    = -1,
    "Condition_MutB_Treated_vs_Control_Vehicle"  = -1,
    "Condition_MutB_Vehicle_vs_Control_Vehicle"    =  1
  )),
  
  # (MutB_Treated - MutB_Vehicle) - (Control_Treated - Control_Vehicle)
  #   - [(MutA_Treated - MutA_Vehicle) - (Control_Treated - Control_Vehicle)]
  # = Condition_MutB_Treated - Condition_MutB_Vehicle
  #   - Condition_MutA_Treated + Condition_MutA_Vehicle
  Interaction_MutB_Treated_Unique = make_numeric_contrast(dds, list(
    "Condition_MutB_Treated_vs_Control_Vehicle"  =  1,
    "Condition_MutB_Vehicle_vs_Control_Vehicle"    = -1,
    "Condition_MutA_Treated_vs_Control_Vehicle"  = -1,
    "Condition_MutA_Vehicle_vs_Control_Vehicle"    =  1
  )),
  
  # (Control_Treated - Control_Vehicle)
  #   - (MutA_Treated - MutA_Vehicle)
  #   - [(MutB_Treated - MutB_Vehicle) - (MutA_Treated - MutA_Vehicle)]
  # = Condition_Control_Treated
  #   - Condition_MutA_Treated + Condition_MutA_Vehicle
  #   - Condition_MutB_Treated + Condition_MutB_Vehicle
  Interaction_Control_Treated_Unique = make_numeric_contrast(dds, list(
    "Condition_Control_Treated_vs_Control_Vehicle"  =  1,
    "Condition_MutA_Treated_vs_Control_Vehicle"  = -1,
    "Condition_MutA_Vehicle_vs_Control_Vehicle"    =  1,
    "Condition_MutB_Treated_vs_Control_Vehicle"  = -1,
    "Condition_MutB_Vehicle_vs_Control_Vehicle"    =  1
  ))
)

contrast_names <- c(names(simple_contrasts), names(interaction_contrasts))

# =====================DEG Extraction, Export & Volcano Plots===========================================

extract_and_plot_DEGs_deseq2 <- function(dds, contrast_spec, coef_name,
                                         pval_cut = 0.05) {
  # contrast_spec: either a character(3) for simple contrasts,
  #                or a numeric vector for compound contrasts
  if (is.character(contrast_spec)) {
    res_raw <- results(dds, contrast  = contrast_spec,
                       pAdjustMethod = "BH", independentFiltering = TRUE)
  } else {
    res_raw <- results(dds, contrast  = contrast_spec,
                       pAdjustMethod = "BH", independentFiltering = TRUE)
  }
  
  res <- as.data.frame(res_raw) %>%
    rownames_to_column("gene_id") %>%
    left_join(counts_filtered %>% select(gene_id, gene_name), by = "gene_id") %>%
    mutate(
      gene_name    = coalesce(gene_name, gene_id),
      Significance = case_when(
        padj < pval_cut & log2FoldChange >  0 ~ "Upregulated",
        padj < pval_cut & log2FoldChange <  0 ~ "Downregulated",
        TRUE ~ "Not Significant"
      )
    ) %>%
    arrange(padj)
  
  write.csv(res,
            paste0("DE_analysis/DEG_", coef_name, "_full.csv"),
            row.names = FALSE)
  write.csv(res %>% filter(Significance != "Not Significant"),
            paste0("DE_analysis/DEG_", coef_name, "_significant.csv"),
            row.names = FALSE)
  
  cat(coef_name,
      "| UP:",   sum(res$Significance == "Upregulated",   na.rm = TRUE),
      "| DOWN:", sum(res$Significance == "Downregulated", na.rm = TRUE), "\n")
  
  # Top 20 labelled genes
  top_genes <- res %>%
    filter(Significance != "Not Significant") %>%
    slice_head(n = 20)
  
  vol <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj),
                         color = Significance, label = gene_name)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = c(
      "Upregulated"     = "coral2",
      "Downregulated"   = "#4169E1",
      "Not Significant" = "gray70"
    )) +
    geom_hline(yintercept = -log10(pval_cut),
               linetype = "dashed", color = "gray28", linewidth = 0.5) +
    geom_vline(xintercept = 0,
               linetype = "dashed", color = "gray28", linewidth = 0.5) +
    geom_text_repel(data = top_genes, size = 3, max.overlaps = 20,
                    box.padding = 0.4, segment.color = "grey50",
                    show.legend = FALSE) +
    labs(
      title = paste("Volcano Plot:", gsub("_", " ", coef_name)),
      x     = "Log2 Fold Change",
      y     = "-Log10 Adjusted P-Value",
      color = "Regulation"
    ) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line        = element_line(color = "black", linewidth = 0.5)
    )
  
  ggsave(paste0("volcano_plots/volcano_", coef_name, ".png"),
         vol, width = 9, height = 7, dpi = 300)
  
  return(invisible(res))
}

# Run for all contrasts
all_DEG_results <- c(
  lapply(names(simple_contrasts), function(cn) {
    extract_and_plot_DEGs_deseq2(dds, simple_contrasts[[cn]], cn)
  }),
  lapply(names(interaction_contrasts), function(cn) {
    extract_and_plot_DEGs_deseq2(dds, interaction_contrasts[[cn]], cn)
  })
)
names(all_DEG_results) <- contrast_names

# --------- Enrichment Analysis ---------
library(enrichR)
library(ggplot2)

enrichment_groups <- c("GO_Biological_Process", "GO_Cellular_Component",
                       "GO_Molecular_Function", "KEGG", "SynGO")

# ---- 1. Create output directories for every contrast × group ----
for (contrast in contrast_names) {
  for (grp in enrichment_groups) {
    dir.create(file.path("Enrichment_analysis", contrast, grp, "barplots"),
               recursive = TRUE, showWarnings = FALSE)
  }
}

# ---- 2. Resolve the best-matching Enrichr library for each group ----
avail_libs  <- enrichR::listEnrichrDbs()
wanted_libs <- sapply(enrichment_groups, function(group) {
  libs <- avail_libs$libraryName[grep(paste0("^", group), avail_libs$libraryName)]
  if (length(libs) > 0) libs[1] else NA
})
wanted_libs <- na.omit(wanted_libs)

# ---- 3. Core enrichment + plotting function ----
do_enrichment_and_plot <- function(res_df, contrast, direction, fill_color) {
  genes <- res_df$gene_name
  genes <- unique(genes[!is.na(genes) & genes != ""])
  if (length(genes) == 0) {
    message("No genes for ", contrast, " (", direction, ") – skipping.")
    return(invisible(NULL))
  }
  
  enrichr_res <- enrichR::enrichr(genes, wanted_libs)
  
  for (lib in names(enrichr_res)) {
    grp <- enrichment_groups[sapply(enrichment_groups,
                                    function(g) startsWith(lib, g))]
    if (length(grp) == 0) next
    
    enrich_csv <- file.path("Enrichment_analysis", contrast, grp,
                            paste0(contrast, "_", direction, "_", lib, "_enrichr.csv"))
    write.csv(enrichr_res[[lib]], enrich_csv, row.names = FALSE)
    
    enrich_df <- enrichr_res[[lib]]
    if (nrow(enrich_df) > 0 && "Term" %in% colnames(enrich_df)) {
      plot_df      <- head(enrich_df[order(enrich_df$Adjusted.P.value), ], 15)
      plot_df$Term <- factor(plot_df$Term, levels = rev(plot_df$Term))
      
      p <- ggplot(plot_df, aes(x = Term, y = -log10(Adjusted.P.value))) +
        geom_col(fill = fill_color) +
        coord_flip() +
        labs(
          title = paste0("Top 15 ", direction, ": ",
                         gsub("_", " ", contrast), " (", lib, ")"),
          x = "Term",
          y = expression(-log[10](Adjusted ~ P-value))
        ) +
        theme_minimal() +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line        = element_line(color = "black", size = 0.5),
          axis.ticks       = element_line(color = "black"),
          plot.title       = element_text(size = 11)
        )
      
      plot_file <- file.path("Enrichment_analysis", contrast, grp, "barplots",
                             paste0(contrast, "_", direction, "_", lib,
                                    "_enrichr_barplot.png"))
      ggsave(plot_file, p, width = 8, height = 6)
    }
  }
}

# ---- 4. Iterate over all contrasts × 2 directions ----
for (cn in contrast_names) {
  res <- all_DEG_results[[cn]]
  
  up_df   <- res[!is.na(res$Significance) & res$Significance == "Upregulated",   ]
  down_df <- res[!is.na(res$Significance) & res$Significance == "Downregulated", ]
  
  do_enrichment_and_plot(up_df,   cn, "upregulated",   "coral2")
  do_enrichment_and_plot(down_df, cn, "downregulated", "#4169E1")
}
