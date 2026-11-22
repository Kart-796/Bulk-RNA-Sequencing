# ==================== BULK RNA-SEQ: limma-voom ANALYSIS ====================
# Input : <bulk_rna_workdir>/counts.txt (alignment/03_count_reads.sh) and Metadata.csv
#         (Sample, Genotype, Treatment, Batch)
# ---- Configuration: edit config/config.yaml (run from the repository root) ----
source("utils/load_config.R")
setwd(cfg_path("bulk_rna_workdir", create = TRUE))
set.seed(1234)

# ======================== LIBRARIES =============================================
library(tidyverse)
library(magrittr)
library(dplyr)
library(limma)
library(ggplot2)
library(ggrepel)
library(hrbrthemes)
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
library(edgeR)
library(kohonen)
library(RColorBrewer)

# =================== Directory Setup ===================
dir.create("Boxplots", showWarnings = FALSE)
dir.create("Violinplots", showWarnings = FALSE)
dir.create("Histograms", showWarnings = FALSE)
dir.create("PCA", showWarnings = FALSE)
dir.create("Scree_Plots", showWarnings = FALSE)
dir.create("Normalization_QC", showWarnings = FALSE)

# ---- Shared ggplot theme (used by every plot in this script) ----
common_theme <- theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )

# ===========  DATA LOADING & CLEANING ====================
counts <- read.table("counts.txt", header = TRUE,
                     sep = "\t", stringsAsFactors = FALSE)

metadata <- read_csv("Metadata.csv")

# ---- Replace NA with 0 ----
counts[is.na(counts)] <- 0

# ---- Map Ensembl IDs to gene symbols ----
counts$gene_name <- mapIds(
  org.Hs.eg.db,
  keys = counts$Geneid,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Fallback: keep original Geneid if no symbol found
counts$gene_name[is.na(counts$gene_name)] <-
  counts$Geneid[is.na(counts$gene_name)]

# ---- Identify sample columns from metadata ----
sample_cols <- metadata$Sample
sample_cols  # sanity check — should match the count columns in counts

# ---- Inspect duplicate gene names before collapsing ----
dup_genes <- counts$gene_name[duplicated(counts$gene_name)]
length(dup_genes)                                                # number of duplicate entries
sort(table(counts$gene_name), decreasing = TRUE)[1:10]          # top repeated gene names
unique(dup_genes)                                                # unique list of duplicated names

# ---- Collapse counts by gene_name (sum duplicate rows) ----
counts_filtered <- counts %>%
  select(gene_name, all_of(sample_cols)) %>%
  group_by(gene_name) %>%
  summarise(across(all_of(sample_cols), sum), .groups = "drop") %>%
  column_to_rownames(var = "gene_name")

# ---- Sanity checks ----
dim(counts_filtered)
anyDuplicated(rownames(counts_filtered))   # should be 0

# ---- Remove genes with zero counts across all samples ----
sum(rowSums(counts_filtered) == 0)          # rows to be dropped
counts_filtered <- counts_filtered[rowSums(counts_filtered) > 0, ]
dim(counts_filtered)

# ==================  Genotype / Treatment factors =====================
# ---- Verify assignment: one row per sample, every Genotype x Treatment group present, zero NAs ----
metadata %>% select(Sample, Genotype, Treatment, Batch)
stopifnot(sum(is.na(metadata$Genotype)) == 0)
stopifnot(sum(is.na(metadata$Treatment)) == 0)
stopifnot(sum(is.na(metadata$Batch)) == 0)
table(metadata$Genotype, metadata$Treatment, metadata$Batch)

# ==================  Low-Expression Filtering (edgeR::filterByExpr) =====================
dge_raw <- DGEList(counts = counts_filtered)

metadata_ord <- metadata[match(colnames(counts_filtered), metadata$Sample), ]
stopifnot(identical(metadata_ord$Sample, colnames(counts_filtered)))

design_filt <- model.matrix(~ Genotype + Treatment + Batch, data = metadata_ord)

# Identify genes with sufficient expression across samples/groups
keep_genes <- filterByExpr(dge_raw, design = design_filt)
sum(keep_genes)                       # genes retained
sum(!keep_genes)                      # genes dropped

# Apply the filter
counts_filtered <- counts_filtered[keep_genes, ]
dim(counts_filtered)                  # sanity check — should drop substantially

# ==================  Log2 transformation (for QC/diagnostic plots only) =====================
counts_log2 <- log2(counts_filtered + 1)

# Sanity checks
dim(counts_log2)
head(counts_log2)
summary(counts_log2)

# ---- Remove genes with zero counts across all samples ----
sum(rowSums(counts_log2) == 0)          # rows to be dropped
counts_log2 <- counts_log2[rowSums(counts_log2) > 0, ]
dim(counts_log2)

# ---- Convert to long format for plotting ----
gene_log_long <- counts_log2 %>%
  rownames_to_column(var = "gene_name") %>%
  pivot_longer(cols = -gene_name, names_to = "Sample", values_to = "log2_counts")

# ---- Boxplot (Log-transformed Data) ----
box1 <- ggplot(gene_log_long, aes(x = Sample, y = log2_counts, fill = Sample)) +
  geom_boxplot() +
  labs(title = "Gene Count Distribution per Sample (Log2 Data)",
       x = "Sample", y = "counts (Log2)") +
  common_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_discrete(name = "Sample")
ggsave("Boxplots/boxplots_log_transformed_data.png", box1, width = 12, height = 6)

# ---- Histogram (Log-transformed Data) ----
hist1 <- ggplot(gene_log_long, aes(x = log2_counts, fill = Sample)) +
  geom_histogram(color = "#e9ecef", alpha = 0.6, position = "identity", bins = 30) +
  labs(fill = "Sample", x = "counts (Log2)", y = "Frequency",
       title = "Histogram of Gene counts per Sample (Log2 Data)") +
  common_theme
ggsave("Histograms/histogram_log_transformed_data.png", hist1, width = 12, height = 6)

# ---- Violin Plot (Log-transformed Data) ----
violin1 <- ggplot(gene_log_long, aes(x = Sample, y = log2_counts, fill = Sample)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  labs(title = "Gene Count Distribution per Sample (Log2 Data) - ",
       x = "Sample", y = "counts (Log2)") +
  common_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_discrete(name = "Sample")
ggsave("Violinplots/violin_log_transformed_data.png", violin1, width = 12, height = 6)


# =========== NORMALIZATION: edgeR (TMM) + limma/voom (logCPM) ===========
dge <- DGEList(counts = counts_filtered)
dge <- normLibSizes(dge, method = "TMM")   # (renamed from calcNormFactors)

metadata_ord <- metadata[match(colnames(counts_filtered), metadata$Sample), ]
stopifnot(identical(metadata_ord$Sample, colnames(counts_filtered)))

design <- model.matrix(~ Genotype + Treatment + Batch, data = metadata_ord)

png("Normalization_QC/voom_mean_variance.png", width = 1600, height = 1200, res = 200)
voom_fit <- voom(dge, design = design, plot = TRUE)
dev.off()

counts_norm <- voom_fit$E   # log2-CPM, TMM+voom normalized, common genes x bulk samples

dim(counts_norm)
head(counts_norm)

# ---- Convert normalized data to long format for plotting ----
gene_norm_long <- counts_norm %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene_name") %>%
  pivot_longer(cols = -gene_name, names_to = "Sample", values_to = "norm_counts")

# ---- Boxplot (TMM + voom Normalized Data) ----
box2 <- ggplot(gene_norm_long, aes(x = Sample, y = norm_counts, fill = Sample)) +
  geom_boxplot() +
  labs(title = "Gene Count Distribution per Sample (TMM + voom Normalized)",
       x = "Sample", y = "counts (log2-CPM)") +
  common_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_discrete(name = "Sample")
ggsave("Boxplots/boxplots_normalized_data.png", box2, width = 12, height = 6)

# ---- Histogram (TMM + voom Normalized Data) ----
hist2 <- ggplot(gene_norm_long, aes(x = norm_counts, fill = Sample)) +
  geom_histogram(color = "#e9ecef", alpha = 0.6, position = "identity", bins = 30) +
  labs(fill = "Sample", x = "counts (log2-CPM)", y = "Frequency",
       title = "Histogram of Gene counts per Sample (TMM + voom Normalized)") +
  common_theme
ggsave("Histograms/histogram_normalized_data.png", hist2, width = 12, height = 6)

# ---- Violin Plot (TMM + voom Normalized Data) ----
violin2 <- ggplot(gene_norm_long, aes(x = Sample, y = norm_counts, fill = Sample)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  labs(title = "Gene Count Distribution per Sample (TMM + voom Normalized)",
       x = "Sample", y = "counts (log2-CPM)") +
  common_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_discrete(name = "Sample")
ggsave("Violinplots/violin_normalized_data.png", violin2, width = 12, height = 6)

# ============ PCA (on TMM+voom normalized data) =============
pca <- prcomp(t(counts_norm), scale. = TRUE, center = TRUE)

# % variance explained by each PC
pca_var <- pca$sdev^2
pca_percentVar <- round(100 * pca_var / sum(pca_var), 2)

# Build a data frame for plotting (PC1/PC2 + sample metadata)
pca_df <- data.frame(
  Sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
) %>%
  left_join(metadata_ord, by = "Sample")

# ---- PCA scatter plot ----
# Colored by Genotype, shaped by Treatment — reflects the real 3x2 design
# rather than falling back to per-sample coloring.
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, label = Sample)) +
  geom_point(size = 3, aes(color = Genotype, shape = Treatment)) +
  geom_text_repel(size = 3) +
  labs(
    title = " PCA (TMM + voom normalized)",
    x = paste0("PC1 (", pca_percentVar[1], "% variance)"),
    y = paste0("PC2 (", pca_percentVar[2], "% variance)")
  ) +
  common_theme
ggsave("PCA/PCA_PC1_PC2.png", pca_plot, width = 8, height = 6)

# ---- Scree plot (variance explained per PC) ----
scree_df <- data.frame(
  PC = factor(seq_along(pca_percentVar), levels = seq_along(pca_percentVar)),
  PercentVar = pca_percentVar
)
scree_plot <- ggplot(scree_df[1:min(20, nrow(scree_df)), ],
                     aes(x = PC, y = PercentVar)) +
  geom_col(fill = "#4E79A7") +
  labs(title = " PCA Scree Plot", x = "Principal Component", y = "% Variance Explained") +
  common_theme
ggsave("Scree_Plots/PCA_scree_plot.png", scree_plot, width = 6, height = 5)


# ==================== DIFFERENTIAL EXPRESSION (limma-voom) ====================
# Same approach as the proteomics analysis (limma + volcano + Venn + Enrichr),
# but genes are called significant on the ADJUSTED p-value (BH, adj.P.Val < 0.05).
#
# Design: additive model ~ Genotype + Treatment + Batch (same as the voom model
# above). Per-genotype treatment effects and interactions are tested with DESeq2
# in analysis/02_deseq2_analysis.R.
# The factor levels are set explicitly so the coefficients read as
# "vs Control" and "Treated vs Vehicle". This design spans the same model space as
# the one used for voom above, so the voom precision weights stay valid.

de_pval <- 0.05                                   # threshold on adj.P.Val
de_dir  <- "DE"
for (d in c("", "volcano_plots/wo_labels", "Venn_Diagrams", "Enrichment_analysis")) {
  dir.create(file.path(de_dir, d), recursive = TRUE, showWarnings = FALSE)
}

de_meta <- metadata_ord %>%
  mutate(Genotype  = factor(Genotype,  levels = c("Control", "MutA", "MutB")),
         Treatment = factor(Treatment, levels = c("Vehicle", "Treated")),
         Batch     = factor(Batch))
stopifnot(!anyNA(de_meta$Genotype), !anyNA(de_meta$Treatment))

de_design <- model.matrix(~ Genotype + Treatment + Batch, data = de_meta)
colnames(de_design) <- make.names(colnames(de_design))
if (!is.null(limma::nonEstimable(de_design))) {
  stop("Design is not estimable (check Batch vs Genotype/Treatment): ",
       paste(limma::nonEstimable(de_design), collapse = ", "))
}
de_resid_df <- nrow(de_design) - ncol(de_design)
cat("DE design:", ncol(de_design), "coefficients,", de_resid_df, "residual df\n")
if (de_resid_df < 2) warning("Only ", de_resid_df, " residual df: DE results have very little power.")

de_contrasts <- makeContrasts(
  MutA_vs_Control    = GenotypeMutA,
  MutB_vs_Control    = GenotypeMutB,
  MutB_vs_MutA       = GenotypeMutB - GenotypeMutA,
  Treated_vs_Vehicle = TreatmentTreated,
  levels = de_design
)

de_fit  <- lmFit(voom_fit, de_design)
de_fit2 <- eBayes(contrasts.fit(de_fit, de_contrasts), robust = TRUE)

png(file.path(de_dir, "SA_plot.png"), width = 7, height = 6, units = "in", res = 150)
plotSA(de_fit2, main = "Residual standard deviation vs average log-expression")
dev.off()
print(summary(decideTests(de_fit2, p.value = de_pval)))   # BH-adjusted

# ---- DEG extraction (adjusted p-value) ----
de_get_deg <- function(fit2, coef, out_dir = de_dir, pval = de_pval) {
  tbl_all <- topTable(fit2, coef = coef, number = Inf, adjust.method = "BH")
  tbl_sig      <- tbl_all[tbl_all$adj.P.Val < pval, ]
  tbl_sig_up   <- tbl_sig[tbl_sig$logFC > 0, ]
  tbl_sig_down <- tbl_sig[tbl_sig$logFC < 0, ]

  write.csv(tbl_all,      file.path(out_dir, paste0(coef, "_all.csv")),           row.names = TRUE)
  write.csv(tbl_sig,      file.path(out_dir, paste0(coef, "_DEG.csv")),           row.names = TRUE)
  write.csv(tbl_sig_up,   file.path(out_dir, paste0(coef, "_upregulated.csv")),   row.names = TRUE)
  write.csv(tbl_sig_down, file.path(out_dir, paste0(coef, "_downregulated.csv")), row.names = TRUE)

  cat(coef, "- adj.P.Val <", pval, ":", nrow(tbl_sig),
      "(", nrow(tbl_sig_up), "up /", nrow(tbl_sig_down), "down )\n")
  list(all = tbl_all, sig = tbl_sig, up = tbl_sig_up, down = tbl_sig_down)
}

de_results <- lapply(setNames(colnames(de_contrasts), colnames(de_contrasts)),
                     function(cf) de_get_deg(de_fit2, cf))

# ---- Volcano plot (adj.P.Val; labels on/off) ----
de_make_volcano <- function(df, name, out_dir = de_dir, pthresh = de_pval, labels = TRUE) {
  df <- df %>%
    rownames_to_column("Gene") %>%
    mutate(Significance = case_when(
      adj.P.Val < pthresh & logFC > 0 ~ "Upregulated",
      adj.P.Val < pthresh & logFC < 0 ~ "Downregulated",
      TRUE                            ~ "Not Significant"),
      Significance = factor(Significance, levels = c("Upregulated", "Downregulated", "Not Significant")))

  n_up   <- sum(df$Significance == "Upregulated")
  n_down <- sum(df$Significance == "Downregulated")

  p <- ggplot(df, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = c("Upregulated" = "coral2", "Downregulated" = "#4169E1",
                                  "Not Significant" = "gray70"), drop = FALSE) +
    geom_hline(yintercept = -log10(pthresh), linetype = "dashed", color = "gray28", linewidth = 0.5) +
    geom_vline(xintercept = 0,               linetype = "dashed", color = "gray28", linewidth = 0.5) +
    theme_minimal() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black", linewidth = 0.5)) +
    labs(title    = paste("Volcano Plot:", name),
         subtitle = paste0(n_up, " up / ", n_down, " down (adj.P.Val < ", pthresh, ")"),
         x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value")

  if (labels) {
    p <- p + geom_text_repel(data = df %>% filter(Significance != "Not Significant") %>%
                               arrange(adj.P.Val) %>% head(50),        # top 50 to keep it readable
                             aes(label = Gene), size = 3, color = "black",
                             max.overlaps = 20, show.legend = FALSE)
    out_file <- file.path(out_dir, "volcano_plots", paste0(name, ".png"))
  } else {
    out_file <- file.path(out_dir, "volcano_plots/wo_labels", paste0(name, ".png"))
  }
  ggsave(out_file, p, width = 9, height = 7)
}

# ---- Venn diagram: overlap of DE genes between two contrasts ----
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

de_make_venn <- function(set_a, set_b, name_a, name_b, title, out_file,
                         fill = c("#E69F00", "#56B4E9")) {
  if (length(set_a) == 0 && length(set_b) == 0) {
    cat("Venn skipped (both sets empty):", title, "\n"); return(invisible(NULL))
  }
  v <- venn.diagram(
    x = setNames(list(set_a, set_b), c(name_a, name_b)), filename = NULL,
    fill = fill, alpha = 0.5, col = "black", lwd = 1, cex = 1.6, fontface = "bold",
    cat.cex = 1.2, cat.fontface = "bold", cat.pos = c(-20, 20), cat.dist = 0.05,
    main = title, main.cex = 1.4, main.fontface = "bold", margin = 0.1)
  png(out_file, width = 7, height = 7, units = "in", res = 300)
  grid::grid.newpage(); grid::grid.draw(v)
  dev.off()
}

# All / up / down Venns for two contrasts + common and unique gene lists
de_run_venn <- function(deg_a, deg_b, name_a, name_b, out_dir = de_dir) {
  venn_dir <- file.path(out_dir, "Venn_Diagrams")
  pair_tag <- paste0(name_a, "_vs_", name_b)
  for (direction in c("sig", "up", "down")) {
    dir_label <- c(sig = "All DE", up = "Upregulated", down = "Downregulated")[[direction]]
    genes_a <- rownames(deg_a[[direction]])
    genes_b <- rownames(deg_b[[direction]])

    de_make_venn(genes_a, genes_b, name_a, name_b,
                 title    = paste0(dir_label, " genes: ", name_a, " vs ", name_b),
                 out_file = file.path(venn_dir, paste0("Venn_", pair_tag, "_", direction, ".png")))

    common <- intersect(genes_a, genes_b)
    common_df <- data.frame(
      Gene        = common,
      logFC_A     = deg_a$all[common, "logFC"],
      P.Value_A   = deg_a$all[common, "P.Value"],
      adj.P.Val_A = deg_a$all[common, "adj.P.Val"],
      logFC_B     = deg_b$all[common, "logFC"],
      P.Value_B   = deg_b$all[common, "P.Value"],
      adj.P.Val_B = deg_b$all[common, "adj.P.Val"])
    colnames(common_df) <- sub("_A$", paste0("_", name_a), colnames(common_df))
    colnames(common_df) <- sub("_B$", paste0("_", name_b), colnames(common_df))
    write.csv(common_df, file.path(venn_dir, paste0("Common_", pair_tag, "_", direction, ".csv")),
              row.names = FALSE)
    write.csv(data.frame(Gene = setdiff(genes_a, genes_b)),
              file.path(venn_dir, paste0("Unique_", name_a, "_in_", pair_tag, "_", direction, ".csv")),
              row.names = FALSE)
    write.csv(data.frame(Gene = setdiff(genes_b, genes_a)),
              file.path(venn_dir, paste0("Unique_", name_b, "_in_", pair_tag, "_", direction, ".csv")),
              row.names = FALSE)
    cat(dir_label, ":", name_a, "=", length(genes_a), "|", name_b, "=", length(genes_b),
        "| common =", length(common), "\n")
  }
}

# Both mutant genotypes against Control (the bulk equivalent of the two proteomics tissues)
de_run_venn(de_results$MutA_vs_Control, de_results$MutB_vs_Control,
            name_a = "MutA_vs_Control", name_b = "MutB_vs_Control")

# ---- Enrichment (Enrichr; newest version of each library) ----
setEnrichrSite("Enrichr")
de_enrichment_groups <- c("GO_Biological_Process", "GO_Cellular_Component",
                          "GO_Molecular_Function", "KEGG", "SynGO")
de_avail_libs <- enrichR::listEnrichrDbs()$libraryName

de_pick_latest_lib <- function(base) {
  libs <- de_avail_libs[grepl(paste0("^", base, "_\\d{4}"), de_avail_libs)]
  if (base == "KEGG") libs <- libs[grepl("Human", libs)]
  if (length(libs) == 0) return(NA_character_)
  libs[which.max(as.numeric(str_extract(libs, "\\d{4}")))]
}
de_libs <- sapply(de_enrichment_groups, de_pick_latest_lib)
if (any(is.na(de_libs))) warning("No Enrichr library found for: ",
                                 paste(names(de_libs)[is.na(de_libs)], collapse = ", "))
de_libs <- de_libs[!is.na(de_libs)]
print(de_libs)

de_do_enrichment <- function(deg_tbl, contrast, direction, fill_color, out_dir = de_dir) {
  genes <- unique(rownames(deg_tbl))
  genes <- genes[!grepl("^ENSG", genes)]          # Enrichr needs symbols; drop unmapped Ensembl IDs
  if (length(genes) == 0) {
    cat(contrast, direction, ": no genes, skipped\n"); return(invisible(NULL))
  }
  cat(contrast, direction, ":", length(genes), "genes\n")
  enrichr_res <- enrichr(genes, unname(de_libs))
  Sys.sleep(1)                                     # be polite to the Enrichr server

  for (grp in names(de_libs)) {
    lib <- de_libs[[grp]]
    enrich_df <- enrichr_res[[lib]]
    if (is.null(enrich_df) || nrow(enrich_df) == 0) next

    grp_dir <- file.path(out_dir, "Enrichment_analysis", contrast, grp)
    dir.create(file.path(grp_dir, "barplots"), recursive = TRUE, showWarnings = FALSE)
    write.csv(enrich_df, file.path(grp_dir, paste0(contrast, "_", direction, "_", lib, "_enrichr.csv")),
              row.names = FALSE)

    plot_df <- enrich_df %>% arrange(Adjusted.P.value) %>% head(15) %>%
      mutate(Term = factor(Term, levels = rev(unique(Term))))
    p <- ggplot(plot_df, aes(x = Term, y = -log10(Adjusted.P.value))) +
      geom_col(fill = fill_color) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray28", linewidth = 0.5) +
      coord_flip() +
      labs(title = paste0("Top 15 ", direction, ": ", contrast, " (", lib, ")"),
           x = "Term", y = expression(-log[10](Adjusted~P-value))) +
      theme_minimal() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.line = element_line(color = "black", linewidth = 0.5),
            axis.ticks = element_line(color = "black"), plot.title = element_text(size = 11))
    ggsave(file.path(grp_dir, "barplots", paste0(contrast, "_", direction, "_", lib, "_enrichr_barplot.png")),
           p, width = 8, height = 6)
  }
}

# ---- Volcano plots + enrichment for every contrast ----
for (contrast in names(de_results)) {
  deg <- de_results[[contrast]]
  de_make_volcano(deg$all, contrast, labels = TRUE)
  de_make_volcano(deg$all, contrast, labels = FALSE)
  de_do_enrichment(deg$up,   contrast, "upregulated",   "coral2")
  de_do_enrichment(deg$down, contrast, "downregulated", "#4169E1")
}


# ==================== SELF-ORGANIZING MAP (SOM) GENE CLUSTERING ====================
dir.create("SOM", showWarnings = FALSE)

# ---- 1. Select variable genes to cluster (top 25% by variance) ----
gene_var <- apply(counts_norm, 1, var)
n_top <- ceiling(0.25 * length(gene_var))
top_var_genes <- names(sort(gene_var, decreasing = TRUE))[1:n_top]

som_input <- counts_norm[top_var_genes, ]
dim(som_input)

# ---- 2. Mean-center and variance-scale each gene across samples ----
som_scaled <- t(scale(t(som_input)))

# Drop any genes that became NA (e.g. zero-variance rows slipping through)
som_scaled <- som_scaled[complete.cases(som_scaled), ]
dim(som_scaled)

# ---- 3. PCA on the scaled, variable-gene matrix (for visual context) ----
pca_som <- prcomp(som_scaled, scale. = TRUE)
pca_som_scores <- data.frame(pca_som$x) %>%
  rownames_to_column("gene") %>%
  select(gene, PC1, PC2)

# ---- 4. Run the SOM ----
set.seed(2)
som_matrix <- as.matrix(som_scaled)

som_grid <- somgrid(6, 6, "hexagonal")
som_model <- som(X = som_matrix, grid = som_grid)

summary(som_model)

# ---- 5. Diagnostic SOM plots ----
png("SOM/SOM_training_changes.png", width = 1600, height = 1200, res = 200)
plot(som_model, type = "changes")
dev.off()

png("SOM/SOM_codes.png", width = 1600, height = 1200, res = 200)
plot(som_model, type = "codes")
dev.off()

png("SOM/SOM_counts.png", width = 1600, height = 1200, res = 200)
plot(som_model, type = "counts")
dev.off()

png("SOM/SOM_dist_neighbours.png", width = 1600, height = 1200, res = 200)
plot(som_model, type = "dist.neighbours")
dev.off()

# Per-sample property heatmaps across the SOM grid
n_samples_som <- ncol(som_model$codes[[1]])
png("SOM/SOM_property_heatmaps.png", width = 2400, height = 1600, res = 200)
par(mfrow = c(ceiling(n_samples_som / 3), min(3, n_samples_som)))
for (i in seq_len(n_samples_som)) {
  plot(som_model, type = "property", property = som_model$codes[[1]][, i],
       main = colnames(som_model$codes[[1]])[i])
}
dev.off()

# ---- 6. Further cluster the SOM nodes via hierarchical clustering ----
n_clusters <- 5   # adjust based on how distinct the SOM regions look visually
som_cluster <- cutree(hclust(dist(som_model$codes[[1]])), n_clusters)

png("SOM/SOM_clusters_mapped.png", width = 1600, height = 1200, res = 200)
plot(som_model, type = "mapping", bgcol = som_cluster, main = "SOM Clusters")
add.cluster.boundaries(som_model, som_cluster)
dev.off()

# ---- 7. Attach gene-level SOM results to a results table ----
som_results <- data.frame(
  gene         = rownames(som_scaled),
  unit.classif = som_model$unit.classif,   # which of the 36 nodes each gene landed in
  distance     = som_model$distances       # distance to that node's codebook vector
)

node_cluster_key <- data.frame(
  unit.classif = seq_len(nrow(som_model$codes[[1]])),
  som_cluster  = som_cluster
)

som_results <- som_results %>%
  left_join(node_cluster_key, by = "unit.classif")

# Merge back with scaled expression values + PCA scores for full context
som_results_full <- som_results %>%
  left_join(
    som_scaled %>% as.data.frame() %>% rownames_to_column("gene"),
    by = "gene"
  ) %>%
  left_join(pca_som_scores, by = "gene")

write.csv(som_results_full, "SOM/SOM_gene_cluster_assignments.csv", row.names = FALSE)

# ---- 8. Visualize genes on PCA, colored by SOM cluster ----
som_pca_plot <- ggplot(som_results_full, aes(x = PC1, y = PC2, color = factor(som_cluster))) +
  geom_point(alpha = 0.7, size = 1.5) +
  labs(title = "Genes in PCA space, colored by SOM cluster",
       x = "PC1", y = "PC2", color = "SOM Cluster") +
  common_theme
ggsave("SOM/SOM_PCA_by_cluster.png", som_pca_plot, width = 8, height = 6)

# ==================== MINIMUM SPANNING TREE OF SAMPLES ====================
library(igraph)

dir.create("SOM/MST", showWarnings = FALSE)

# ---- 1. Distance matrix between the 6 samples ----
# Uses the full TMM+voom normalized expression matrix (counts_norm) --
# the same data used for the sample-level PCA earlier in this script.
sample_dist_mat <- as.matrix(dist(t(counts_norm)))

# ---- 2. Build a fully connected weighted graph, then compute the MST ----
g_full_samples <- graph_from_adjacency_matrix(
  sample_dist_mat, mode = "undirected", weighted = TRUE, diag = FALSE
)
g_mst_samples <- mst(g_full_samples)

# ---- 3. Attach Genotype/Treatment metadata for coloring/shaping ----
sample_meta <- metadata_ord[match(V(g_mst_samples)$name, metadata_ord$Sample), ]
V(g_mst_samples)$Genotype  <- as.character(sample_meta$Genotype)
V(g_mst_samples)$Treatment <- as.character(sample_meta$Treatment)

genotype_colors <- setNames(
  RColorBrewer::brewer.pal(length(unique(sample_meta$Genotype)), "Set2"),
  unique(sample_meta$Genotype)
)
treatment_shapes <- setNames(c("circle", "square"), unique(sample_meta$Treatment))

# ---- 4. Plot ----
png("SOM/MST/Sample_MST.png", width = 1600, height = 1600, res = 200)
set.seed(42)
plot(
  g_mst_samples,
  vertex.label       = V(g_mst_samples)$name,
  vertex.size        = 25,
  vertex.color       = genotype_colors[V(g_mst_samples)$Genotype],
  vertex.shape       = treatment_shapes[V(g_mst_samples)$Treatment],
  vertex.label.cex   = 0.8,
  vertex.label.color = "black",
  vertex.label.dist  = 2.5,
  edge.width         = 2,
  edge.label         = round(E(g_mst_samples)$weight, 1),
  edge.label.cex     = 0.7,
  edge.color         = "gray50",
  layout             = layout_with_fr(g_mst_samples),
  main               = "Minimum Spanning Tree of Samples (TMM+voom normalized)"
)
legend("topright", legend = names(genotype_colors), fill = genotype_colors,
       border = NA, bty = "n", cex = 0.8, title = "Genotype")
legend("bottomright", legend = names(treatment_shapes),
       pch = ifelse(treatment_shapes == "circle", 21, 22), pt.bg = "gray70",
       bty = "n", cex = 0.8, title = "Treatment")
dev.off()

# ==================== MINIMUM SPANNING TREE OF SAMPLES (FROM SOM NODE PROFILES) ====================
library(igraph)

dir.create("SOM/MST", showWarnings = FALSE)

# ---- 1. For each SOM node, compute each sample's mean (scaled) expression ----
# across the genes assigned to that node -- i.e. summarize each sample through
# the SOM's own 36-node structure, using the same som_scaled matrix that was
# actually fed into the SOM (not the raw counts_norm).
node_profiles <- som_scaled %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  left_join(som_results %>% select(gene, unit.classif), by = "gene") %>%
  group_by(unit.classif) %>%
  summarise(across(all_of(colnames(counts_norm)), mean), .groups = "drop")

# node_profiles is nodes x samples; transpose to get samples x nodes
node_profiles_mat <- node_profiles %>%
  column_to_rownames("unit.classif") %>%
  as.matrix()

sample_node_mat <- t(node_profiles_mat)   # rows = samples, cols = SOM nodes

dim(sample_node_mat)   # sanity check: should be 6 samples x (up to 36) nodes

# ---- 2. Distance between samples, based on their SOM node profiles ----
sample_dist_mat <- as.matrix(dist(sample_node_mat))

# ---- 3. Build a fully connected weighted graph, then compute the MST ----
g_full_samples <- graph_from_adjacency_matrix(
  sample_dist_mat, mode = "undirected", weighted = TRUE, diag = FALSE
)
g_mst_samples <- mst(g_full_samples)

# ---- 4. Attach Genotype/Treatment metadata for coloring/shaping ----
sample_meta <- metadata_ord[match(V(g_mst_samples)$name, metadata_ord$Sample), ]
V(g_mst_samples)$Genotype  <- as.character(sample_meta$Genotype)
V(g_mst_samples)$Treatment <- as.character(sample_meta$Treatment)

genotype_colors <- setNames(
  RColorBrewer::brewer.pal(length(unique(sample_meta$Genotype)), "Set2"),
  unique(sample_meta$Genotype)
)
treatment_shapes <- setNames(c("circle", "square"), unique(sample_meta$Treatment))

# ---- 5. Plot ----
png("SOM/MST/Sample_MST_from_SOM.png", width = 1600, height = 1600, res = 200)
set.seed(42)
plot(
  g_mst_samples,
  vertex.label        = V(g_mst_samples)$name,
  vertex.size         = 25,
  vertex.color        = genotype_colors[V(g_mst_samples)$Genotype],
  vertex.shape        = treatment_shapes[V(g_mst_samples)$Treatment],
  vertex.label.cex    = 0.8,
  vertex.label.color  = "black",
  vertex.label.dist   = 2.5,
  edge.width          = 2,
  edge.label          = round(E(g_mst_samples)$weight, 1),
  edge.label.cex      = 0.7,
  edge.color          = "gray50",
  layout              = layout_with_fr(g_mst_samples),
  main                = "Minimum Spanning Tree of Samples (from SOM node profiles)"
)
legend("topright", legend = names(genotype_colors), fill = genotype_colors,
       border = NA, bty = "n", cex = 0.8, title = "Genotype")
legend("bottomright", legend = names(treatment_shapes),
       pch = ifelse(treatment_shapes == "circle", 21, 22), pt.bg = "gray70",
       bty = "n", cex = 0.8, title = "Treatment")
dev.off()

# ==================== EXTRACT METAGENES & PATHWAY ENRICHMENT (PER NODE) ====================
org_db   <- org.Hs.eg.db   # already loaded at the top of this script
organism <- "human"

# ---- 1. Genes per SOM node (BMU) ----
# Built from som_results (already computed above) rather than re-deriving
# from som_model$data, since som_results already has gene + node mapping.
gene_nodes <- som_results %>%
  select(gene, BMU = unit.classif)

# ---- 2. Reshape into "metagenes": one column per node, genes listed as rows ----
metagenes <- gene_nodes %>%
  group_by(BMU) %>%
  mutate(row = row_number()) %>%
  ungroup() %>%
  pivot_wider(
    names_from  = BMU,
    values_from = gene
  )
metagenes <- metagenes[, -1]   # drop the helper "row" column

# ---- 3. Global gene -> Entrez ID lookup (SYMBOL, with ENSEMBL fallback) ----
# IMPORTANT: a large share of "gene_name" values in this pipeline are raw
# Ensembl IDs (e.g. "ENSG00000273906"), not real gene symbols -- these are
# genes where the mapIds(ENSEMBL -> SYMBOL) step at the TOP of this script
# returned NA, and the fallback line kept the Ensembl ID in place of a symbol:
#   counts$gene_name[is.na(counts$gene_name)] <- counts$Geneid[is.na(counts$gene_name)]
# bitr(..., fromType = "SYMBOL") can never resolve these (they were never
# symbols to begin with) -- this is almost certainly the main driver of the
# "30.38% of input gene IDs are fail to map" warning you saw on the universe.
#
# Fix: build ONE lookup table for every gene in counts_filtered, trying
# SYMBOL -> ENTREZID first, and falling back to ENSEMBL -> ENTREZID for
# anything that still looks like a raw Ensembl ID. Doing this once, globally,
# also avoids calling bitr() 36 times (once per node) and re-triggering the
# same failed-mapping warnings for genes that repeat across nodes.
all_genes     <- rownames(counts_filtered)
is_ensembl_id <- grepl("^ENSG", all_genes)

symbol_genes  <- all_genes[!is_ensembl_id]
ensembl_genes <- all_genes[is_ensembl_id]

sym_to_entrez <- bitr(symbol_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org_db)
names(sym_to_entrez)[names(sym_to_entrez) == "SYMBOL"] <- "gene"

# Only attempt bitr() on Ensembl IDs that actually exist as keys in org_db.
# NOTE: these are the exact same Geneid values that already failed the
# org.Hs.eg.db ENSEMBL->SYMBOL lookup at the TOP of this script (that's WHY
# they fell back to a raw Ensembl ID in the first place) -- so it's expected
# that many of them will also fail here. This isn't a new bug: it means
# org.Hs.eg.db's annotation doesn't include these specific gene IDs at all,
# most likely because counts.txt was generated against a newer/different
# Ensembl gene annotation than the one built into your installed org.Hs.eg.db
# (Bioconductor annotation packages often lag behind Ensembl's own releases).
# Calling bitr() directly on IDs org_db has never heard of throws a hard
# error ("None of the keys entered are valid keys") rather than a warning --
# so we filter to valid keys first to avoid crashing the whole script.
valid_ensembl_keys  <- keys(org_db, keytype = "ENSEMBL")
ensembl_genes_valid <- intersect(ensembl_genes, valid_ensembl_keys)

message(
  length(ensembl_genes), " genes are Ensembl-ID fallbacks (unmapped symbols); ",
  length(ensembl_genes_valid), " of those exist in org_db and will be mapped; ",
  length(ensembl_genes) - length(ensembl_genes_valid),
  " are not in org_db at all and will remain unmapped (excluded from enrichment)."
)

if (length(ensembl_genes_valid) > 0) {
  ens_to_entrez <- bitr(ensembl_genes_valid, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org_db)
  names(ens_to_entrez)[names(ens_to_entrez) == "ENSEMBL"] <- "gene"
} else {
  ens_to_entrez <- data.frame(gene = character(0), ENTREZID = character(0))
}

gene_to_entrez <- bind_rows(sym_to_entrez, ens_to_entrez) %>%
  distinct(gene, .keep_all = TRUE)   # keep first ENTREZID if a gene mapped >1

# Sanity check: how much more do we recover vs. the SYMBOL-only approach?
length(all_genes)                          # total genes
nrow(gene_to_entrez)                       # genes successfully mapped
round(mean(all_genes %in% gene_to_entrez$gene) * 100, 1)  # % mapped

# ---- 4. Map each node's genes to Entrez IDs via the lookup table ----
metagenes_entrezIDs <- lapply(colnames(metagenes), function(node_nr) {
  message("Node ", node_nr)
  
  gene_list <- metagenes[[node_nr]]
  gene_list <- gene_list[!is.na(gene_list)]
  
  gene_df <- gene_to_entrez %>%
    filter(gene %in% gene_list) %>%
    select(SYMBOL = gene, ENTREZID)
  
  if (nrow(gene_df) == 0) {
    message("Node ", node_nr, ": no mappable genes")
    return(NULL)
  }
  return(gene_df)
})
names(metagenes_entrezIDs) <- colnames(metagenes)

# ---- 5. Background/universe for enrichment ----
# Uses the full filtered gene set from this pipeline (post-filterByExpr),
# i.e. every gene that could plausibly have been detected -- not just the
# top-25%-variable subset that went into the SOM. (Previously this pointed
# at an undefined variable, region_old_opc$X, left over from a different
# script -- replaced with the actual gene set from this analysis.)
universe_entrez <- gene_to_entrez$ENTREZID
length(universe_entrez)

# ---- 6. Run Reactome ORA per node ----
# NOTE: nodes with very few genes (e.g. 1-2) are unlikely to yield any
# significant enrichment -- ORA needs a reasonably sized gene set to have
# power. Expect many nodes to come back NULL or empty; that's not a bug.
metagenes_enrichment <- lapply(names(metagenes_entrezIDs), function(node) {
  message(node)
  gene_entrez <- metagenes_entrezIDs[[node]]$ENTREZID
  
  if (is.null(gene_entrez) || length(gene_entrez) == 0) {
    message("Node ", node, ": no mappable genes, skipping")
    return(NULL)
  }
  
  enrich_node <- tryCatch(
    enrichPathway(
      gene          = gene_entrez,
      universe      = universe_entrez,
      organism      = organism,
      pAdjustMethod = "BH",
      pvalueCutoff  = 1,
      qvalueCutoff  = 1,
      readable      = TRUE
    ),
    error = function(e) {
      message("Node ", node, ": enrichPathway() failed (", conditionMessage(e), ")")
      return(NULL)
    }
  )
  return(enrich_node)
})
names(metagenes_enrichment) <- names(metagenes_entrezIDs)

saveRDS(metagenes_enrichment, "SOM/MetageneEnrichment.rds")

# ---- 7. Keep only significantly enriched pathways per node ----
pathways <- readRDS("SOM/MetageneEnrichment.rds")

# Drop nodes with NULL results (no mappable genes, or enrichPathway
# failed/returned nothing) -- determined dynamically rather than
# hardcoding a specific node number, since which nodes end up empty
# depends on your data and grid size, not a fixed node index.
pathways <- pathways[!sapply(pathways, is.null)]

significantMetagenePathways <- lapply(names(pathways), function(node) {
  message(node)
  pathways_node <- pathways[[node]]@result
  pathways_node <- pathways_node %>% filter(p.adjust < 0.05)
  return(pathways_node)
})
names(significantMetagenePathways) <- names(pathways)

# ---- 8. Keep only nodes that had at least one significant pathway ----
significantMetagenePathways.filtered <- lapply(names(significantMetagenePathways), function(node) {
  message(node)
  pw <- significantMetagenePathways[[node]]
  if (nrow(pw) > 0) {
    return(pw)
  } else {
    return(NULL)
  }
})
names(significantMetagenePathways.filtered) <- names(significantMetagenePathways)
significantMetagenePathways.filtered <- significantMetagenePathways.filtered[!sapply(significantMetagenePathways.filtered, is.null)]

saveRDS(significantMetagenePathways.filtered, "SOM/SignificantMetagenePathways.rds")

# ==================== PER-SAMPLE ORA (TOP-EXPRESSED GENES) ====================
dir.create("SOM/PerSample_ORA", showWarnings = FALSE)

n_top_per_sample <- 500   # adjust as desired; 500 is a common starting point

# ---- 1. Top-N expressed genes per sample ----
sample_top_genes <- lapply(colnames(counts_norm), function(s) {
  vals <- counts_norm[, s]
  names(sort(vals, decreasing = TRUE))[1:n_top_per_sample]
})
names(sample_top_genes) <- colnames(counts_norm)

# ---- 2. Map each sample's gene list to Entrez IDs via the existing lookup ----
# (gene_to_entrez was already built above in the metagene section -- reused
# here rather than re-running bitr())
sample_entrezIDs <- lapply(names(sample_top_genes), function(s) {
  gene_list <- sample_top_genes[[s]]
  gene_to_entrez %>%
    filter(gene %in% gene_list) %>%
    select(SYMBOL = gene, ENTREZID)
})
names(sample_entrezIDs) <- names(sample_top_genes)

# Quick check: how many of each sample's top genes actually mapped?
sapply(sample_entrezIDs, nrow)

# ---- 3. Run Reactome ORA per sample ----
sample_enrichment <- lapply(names(sample_entrezIDs), function(s) {
  message(s)
  gene_entrez <- sample_entrezIDs[[s]]$ENTREZID
  
  if (is.null(gene_entrez) || length(gene_entrez) == 0) {
    message(s, ": no mappable genes, skipping")
    return(NULL)
  }
  
  enrich_s <- tryCatch(
    enrichPathway(
      gene          = gene_entrez,
      universe      = universe_entrez,
      organism      = organism,
      pAdjustMethod = "BH",
      pvalueCutoff  = 1,
      qvalueCutoff  = 1,
      readable      = TRUE
    ),
    error = function(e) {
      message(s, ": enrichPathway() failed (", conditionMessage(e), ")")
      return(NULL)
    }
  )
  return(enrich_s)
})
names(sample_enrichment) <- names(sample_entrezIDs)

saveRDS(sample_enrichment, "SOM/PerSample_ORA/SampleEnrichment.rds")

# ---- 4. Keep only significantly enriched pathways per sample ----
sample_enrichment_nonnull <- sample_enrichment[!sapply(sample_enrichment, is.null)]

significantSamplePathways <- lapply(names(sample_enrichment_nonnull), function(s) {
  message(s)
  res <- sample_enrichment_nonnull[[s]]@result
  res %>% filter(p.adjust < 0.05)
})
names(significantSamplePathways) <- names(sample_enrichment_nonnull)

# ---- 5. Keep only samples that had at least one significant pathway ----
significantSamplePathways.filtered <- lapply(names(significantSamplePathways), function(s) {
  pw <- significantSamplePathways[[s]]
  if (nrow(pw) > 0) pw else NULL
})
names(significantSamplePathways.filtered) <- names(significantSamplePathways)
significantSamplePathways.filtered <- significantSamplePathways.filtered[!sapply(significantSamplePathways.filtered, is.null)]

saveRDS(significantSamplePathways.filtered, "SOM/PerSample_ORA/SignificantSamplePathways.rds")

# ---- 6. Combine into one long table for easy comparison across samples ----
significantSamplePathways.combined <- bind_rows(
  lapply(names(significantSamplePathways.filtered), function(s) {
    significantSamplePathways.filtered[[s]] %>% mutate(Sample = s)
  })
)
write.csv(significantSamplePathways.combined, "SOM/PerSample_ORA/SignificantSamplePathways_combined.csv", row.names = FALSE)

# ==================== BAR PLOTS: TOP ENRICHED TERMS PER SAMPLE ====================
library(enrichplot)

dir.create("SOM/PerSample_ORA/BarPlots", showWarnings = FALSE)

n_terms_to_show <- 10   # adjust as desired

for (s in names(sample_enrichment)) {
  enrich_s <- sample_enrichment[[s]]
  
  # skip samples with no enrichment result, or zero significant terms
  if (is.null(enrich_s) || nrow(enrich_s@result) == 0) {
    message(s, ": no enrichment results, skipping barplot")
    next
  }
  
  bar_plot <- barplot(enrich_s, showCategory = n_terms_to_show) +
    ggtitle(paste("Top Enriched Pathways -", s)) +
    common_theme
  
  ggsave(
    filename = paste0("SOM/PerSample_ORA/BarPlots/", s, "_top_pathways.png"),
    plot = bar_plot,
    width = 10, height = 6
  )
}

# ==================== IDENTIFY & PLOT UNIQUE ENRICHED TERMS PER SAMPLE ====================
dir.create("SOM/PerSample_ORA/BarPlots_Unique", showWarnings = FALSE)

# ---- 1. Collect all significant terms across samples into one table ----
all_terms_df <- bind_rows(lapply(names(sample_enrichment), function(s) {
  enrich_s <- sample_enrichment[[s]]
  if (is.null(enrich_s) || nrow(enrich_s@result) == 0) return(NULL)
  enrich_s@result %>%
    filter(p.adjust < 0.05) %>%
    mutate(Sample = s)
}))

# ---- 2. Count how many samples each pathway term appears in ----
term_sample_counts <- all_terms_df %>%
  distinct(Description, Sample) %>%
  count(Description, name = "n_samples")

# ---- 3. Terms significant in exactly one sample = "unique" to that sample ----
unique_terms_df <- all_terms_df %>%
  left_join(term_sample_counts, by = "Description") %>%
  filter(n_samples == 1)

write.csv(unique_terms_df, "SOM/PerSample_ORA/BarPlots_Unique/unique_terms_per_sample.csv", row.names = FALSE)

# ---- 4. Bar plot of unique terms, one per sample ----
for (s in unique(unique_terms_df$Sample)) {
  sample_unique <- unique_terms_df %>%
    filter(Sample == s) %>%
    arrange(p.adjust) %>%
    slice_head(n = n_terms_to_show)
  
  if (nrow(sample_unique) == 0) next
  
  # order bars so the most significant term appears at the top
  sample_unique <- sample_unique %>%
    mutate(Description = factor(Description, levels = rev(Description)))
  
  unique_plot <- ggplot(sample_unique, aes(x = Description, y = -log10(p.adjust))) +
    geom_col(fill = "#4E79A7") +
    coord_flip() +
    labs(title = paste("Unique Enriched Pathways -", s),
         x = NULL, y = expression(-log[10](adjusted~p))) +
    common_theme
  
  ggsave(
    filename = paste0("SOM/PerSample_ORA/BarPlots_Unique/", s, "_unique_pathways.png"),
    plot = unique_plot,
    width = 10, height = 6
  )
}

# ==================== CIRCOS/CHORD PLOT: SHARED TERMS BETWEEN SAMPLES ====================
library(circlize)

dir.create("SOM/PerSample_ORA/Circos", showWarnings = FALSE)

terms_by_sample <- split(all_terms_df$Description, all_terms_df$Sample)
samples <- names(sample_enrichment)

n <- length(samples)
shared_matrix <- matrix(0, nrow = n, ncol = n, dimnames = list(samples, samples))

for (i in seq_len(n)) {
  for (j in seq_len(n)) {
    if (i != j) {
      shared_matrix[i, j] <- length(intersect(
        terms_by_sample[[samples[i]]],
        terms_by_sample[[samples[j]]]
      ))
    }
  }
}

shared_matrix

sample_colors <- setNames(RColorBrewer::brewer.pal(n, "Set2"), samples)

png("SOM/PerSample_ORA/Circos/shared_terms_circos.png", width = 2200, height = 1800, res = 200)

# Extra margin on the right reserves room for the legend; expanded canvas
# limits give the outer sample-name labels room so they don't get clipped.
par(mar = c(1, 1, 2, 8))
circos.clear()
circos.par(canvas.xlim = c(-1.3, 1.3), canvas.ylim = c(-1.3, 1.3))

chordDiagram(
  shared_matrix,
  symmetric         = TRUE,
  grid.col          = sample_colors,
  transparency      = 0.4,
  annotationTrack   = c("grid"),
  preAllocateTracks = list(track.height = 0.15)
)

circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(
    CELL_META$xcenter, CELL_META$ylim[1] + 1,
    CELL_META$sector.index,
    facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.8
  )
}, bg.border = NA)

title("Shared Enriched Pathway Terms Between Samples")

# ---- Legend: which color = which sample ----
legend(
  x = "right",
  inset = -0.25,        # pushes legend outside the circle, into the margin reserved above
  legend = samples,
  fill   = sample_colors,
  border = NA,
  bty    = "n",
  cex    = 0.8,
  title  = "Sample",
  xpd    = TRUE          # allows drawing outside the normal plot region
)

circos.clear()
dev.off()

# ==================== PER-SAMPLE ORA: UP- AND DOWN-REGULATED GENES ====================
# "Up"/"down" here = relative to the mean of the OTHER 5 samples -- with no
# control/baseline group and no replicates, this is the natural stand-in for
# "what's distinctive about this sample" rather than "highest absolute
# expression" (which mostly just surfaces housekeeping/ribosomal genes in
# every sample, as the earlier top-expressed version tended to do).

dir.create("SOM/PerSample_ORA/Up", showWarnings = FALSE)
dir.create("SOM/PerSample_ORA/Down", showWarnings = FALSE)

n_top_per_direction <- 250   # per direction, per sample; adjust as desired
samples <- colnames(counts_norm)

# ---- 1. log2FC of each sample vs. the mean of all OTHER samples ----
# (counts_norm is already log2-CPM, so a simple subtraction IS the log2FC)
sample_fc_list <- lapply(samples, function(s) {
  other_mean <- rowMeans(counts_norm[, setdiff(samples, s), drop = FALSE])
  data.frame(gene = rownames(counts_norm), log2FC = counts_norm[, s] - other_mean) %>%
    arrange(desc(log2FC))
})
names(sample_fc_list) <- samples

# ---- 2. Top up- and down-regulated gene lists per sample ----
sample_up_genes   <- lapply(sample_fc_list, function(df) head(df, n_top_per_direction)$gene)
sample_down_genes <- lapply(sample_fc_list, function(df) tail(df, n_top_per_direction)$gene)

# ---- 3. Map to Entrez IDs via the existing lookup (gene_to_entrez) ----
map_genes_to_entrez <- function(gene_list) {
  gene_to_entrez %>% filter(gene %in% gene_list) %>% select(SYMBOL = gene, ENTREZID)
}
sample_up_entrez   <- lapply(sample_up_genes,   map_genes_to_entrez)
sample_down_entrez <- lapply(sample_down_genes, map_genes_to_entrez)

sapply(sample_up_entrez, nrow)    # sanity check: how many mapped, up
sapply(sample_down_entrez, nrow)  # sanity check: how many mapped, down

# ---- 4. Run Reactome ORA per sample, per direction ----
run_ora_list <- function(entrez_list, label) {
  result <- lapply(names(entrez_list), function(s) {
    message(label, " - ", s)
    gene_entrez <- entrez_list[[s]]$ENTREZID
    if (is.null(gene_entrez) || length(gene_entrez) == 0) {
      message(s, ": no mappable genes, skipping")
      return(NULL)
    }
    tryCatch(
      enrichPathway(
        gene = gene_entrez, universe = universe_entrez, organism = organism,
        pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE
      ),
      error = function(e) { message(s, ": enrichPathway() failed (", conditionMessage(e), ")"); NULL }
    )
  })
  names(result) <- names(entrez_list)
  result
}

sample_enrichment_up   <- run_ora_list(sample_up_entrez,   "UP")
sample_enrichment_down <- run_ora_list(sample_down_entrez, "DOWN")

saveRDS(sample_enrichment_up,   "SOM/PerSample_ORA/Up/SampleEnrichment_Up.rds")
saveRDS(sample_enrichment_down, "SOM/PerSample_ORA/Down/SampleEnrichment_Down.rds")

# ---- 5. Significant pathways per sample, per direction ----
get_significant <- function(enrichment_list) {
  nonnull <- enrichment_list[!sapply(enrichment_list, is.null)]
  sig <- lapply(names(nonnull), function(s) nonnull[[s]]@result %>% filter(p.adjust < 0.05))
  names(sig) <- names(nonnull)
  sig[sapply(sig, nrow) > 0]
}

significant_up   <- get_significant(sample_enrichment_up)
significant_down <- get_significant(sample_enrichment_down)

saveRDS(significant_up,   "SOM/PerSample_ORA/Up/SignificantSamplePathways_Up.rds")
saveRDS(significant_down, "SOM/PerSample_ORA/Down/SignificantSamplePathways_Down.rds")

combine_direction <- function(sig_list, direction) {
  bind_rows(lapply(names(sig_list), function(s) sig_list[[s]] %>% mutate(Sample = s, Direction = direction)))
}

significant_up_combined     <- combine_direction(significant_up,   "Up")
significant_down_combined   <- combine_direction(significant_down, "Down")
significant_updown_combined <- bind_rows(significant_up_combined, significant_down_combined)

write.csv(significant_updown_combined,
          "SOM/PerSample_ORA/SignificantSamplePathways_UpDown_combined.csv", row.names = FALSE)


# ==================== BAR PLOTS: TOP ENRICHED TERMS PER SAMPLE, BY DIRECTION ====================
library(enrichplot)
n_terms_to_show <- 10

plot_direction_barplots <- function(enrichment_list, direction, outdir) {
  for (s in names(enrichment_list)) {
    enrich_s <- enrichment_list[[s]]
    if (is.null(enrich_s) || nrow(enrich_s@result) == 0) {
      message(s, " (", direction, "): no enrichment results, skipping barplot")
      next
    }
    bar_plot <- barplot(enrich_s, showCategory = n_terms_to_show) +
      ggtitle(paste0("Top Enriched Pathways (", direction, ") - ", s)) +
      common_theme
    ggsave(paste0(outdir, "/", s, "_", tolower(direction), "_pathways.png"),
           bar_plot, width = 10, height = 6)
  }
}

plot_direction_barplots(sample_enrichment_up,   "Up",   "SOM/PerSample_ORA/Up")
plot_direction_barplots(sample_enrichment_down, "Down", "SOM/PerSample_ORA/Down")


# ==================== UNIQUE ENRICHED TERMS PER SAMPLE, BY DIRECTION ====================
dir.create("SOM/PerSample_ORA/Up/Unique", showWarnings = FALSE)
dir.create("SOM/PerSample_ORA/Down/Unique", showWarnings = FALSE)

find_unique_terms <- function(combined_df) {
  term_counts <- combined_df %>% distinct(Description, Sample) %>% count(Description, name = "n_samples")
  combined_df %>% left_join(term_counts, by = "Description") %>% filter(n_samples == 1)
}

unique_up_df   <- find_unique_terms(significant_up_combined)
unique_down_df <- find_unique_terms(significant_down_combined)

write.csv(unique_up_df,   "SOM/PerSample_ORA/Up/Unique/unique_terms_up.csv",   row.names = FALSE)
write.csv(unique_down_df, "SOM/PerSample_ORA/Down/Unique/unique_terms_down.csv", row.names = FALSE)

plot_unique_barplots <- function(unique_df, direction, outdir) {
  bar_color <- ifelse(direction == "Up", "#E15759", "#4E79A7")
  for (s in unique(unique_df$Sample)) {
    sample_unique <- unique_df %>% filter(Sample == s) %>% arrange(p.adjust) %>% slice_head(n = n_terms_to_show)
    if (nrow(sample_unique) == 0) next
    sample_unique <- sample_unique %>% mutate(Description = factor(Description, levels = rev(Description)))
    unique_plot <- ggplot(sample_unique, aes(x = Description, y = -log10(p.adjust))) +
      geom_col(fill = bar_color) +
      coord_flip() +
      labs(title = paste0("Unique Enriched Pathways (", direction, ") - ", s),
           x = NULL, y = expression(-log[10](adjusted~p))) +
      common_theme
    ggsave(paste0(outdir, "/", s, "_unique_", tolower(direction), "_pathways.png"),
           unique_plot, width = 10, height = 6)
  }
}

plot_unique_barplots(unique_up_df,   "Up",   "SOM/PerSample_ORA/Up/Unique")
plot_unique_barplots(unique_down_df, "Down", "SOM/PerSample_ORA/Down/Unique")


# ==================== CIRCOS/CHORD PLOTS: SHARED TERMS, BY DIRECTION ====================
library(circlize)
dir.create("SOM/PerSample_ORA/Circos", showWarnings = FALSE)

make_shared_matrix <- function(combined_df, samples) {
  terms_by_sample <- split(combined_df$Description, combined_df$Sample)
  n <- length(samples)
  m <- matrix(0, nrow = n, ncol = n, dimnames = list(samples, samples))
  for (i in seq_len(n)) for (j in seq_len(n)) if (i != j) {
    a <- terms_by_sample[[samples[i]]]; b <- terms_by_sample[[samples[j]]]
    m[i, j] <- if (is.null(a) || is.null(b)) 0 else length(intersect(a, b))
  }
  m
}

plot_circos <- function(shared_matrix, samples, title_text, outfile) {
  sample_colors <- setNames(RColorBrewer::brewer.pal(length(samples), "Set2"), samples)
  png(outfile, width = 2200, height = 1800, res = 200)
  par(mar = c(1, 1, 2, 8))
  circos.clear()
  circos.par(canvas.xlim = c(-1.3, 1.3), canvas.ylim = c(-1.3, 1.3))
  chordDiagram(shared_matrix, symmetric = TRUE, grid.col = sample_colors, transparency = 0.4,
               annotationTrack = c("grid"), preAllocateTracks = list(track.height = 0.15))
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1] + 1, CELL_META$sector.index,
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.8)
  }, bg.border = NA)
  title(title_text)
  legend(x = "right", inset = -0.25, legend = samples, fill = sample_colors,
         border = NA, bty = "n", cex = 0.8, title = "Sample", xpd = TRUE)
  circos.clear()
  dev.off()
}

shared_matrix_up   <- make_shared_matrix(significant_up_combined,   samples)
shared_matrix_down <- make_shared_matrix(significant_down_combined, samples)
shared_matrix_up     # inspect before plotting
shared_matrix_down

plot_circos(shared_matrix_up,   samples, "Shared Enriched Pathway Terms (Up-regulated)",
            "SOM/PerSample_ORA/Circos/shared_terms_circos_up.png")
plot_circos(shared_matrix_down, samples, "Shared Enriched Pathway Terms (Down-regulated)",
            "SOM/PerSample_ORA/Circos/shared_terms_circos_down.png")

# ==================== GO / KEGG / SynGO ENRICHMENT (PER SAMPLE, UP/DOWN) ====================

dbs <- c(
  "GO_Biological_Process_2026",
  "GO_Molecular_Function_2026",
  "GO_Cellular_Component_2026",
  "KEGG_2026",
  "SynGO_2024"
)

dir.create("SOM/PerSample_ORA/EnrichR", showWarnings = FALSE)
dir.create("SOM/PerSample_ORA/EnrichR/Up", showWarnings = FALSE)
dir.create("SOM/PerSample_ORA/EnrichR/Down", showWarnings = FALSE)

# ---- 1. Download each library's gene sets from Enrichr (once, reused for all samples) ----
get_enrichr_geneset <- function(library_name) {
  url <- paste0("https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=", library_name)
  lines <- readLines(url, warn = FALSE)
  bind_rows(lapply(lines, function(line) {
    parts <- strsplit(line, "\t")[[1]]
    term  <- parts[1]
    genes <- parts[-c(1, 2)]   # parts[2] is a description/blank field in Enrichr's GMT format
    genes <- genes[genes != ""]
    if (length(genes) == 0) return(NULL)
    data.frame(term = term, gene = genes)
  }))
}

message("Downloading Enrichr gene set libraries...")
enrichr_genesets <- lapply(dbs, get_enrichr_geneset)
names(enrichr_genesets) <- dbs
sapply(enrichr_genesets, function(x) length(unique(x$term)))   # sanity check: terms per library

# ---- 2. Universe: post-filterByExpr gene set, restricted per-database to ----
# genes that actually appear in that library's own gene sets (genes absent
# from a library's TERM2GENE would distort enricher()'s background math if
# left in).
universe_symbols_all <- rownames(counts_filtered)
length(universe_symbols_all)

build_universe_for_db <- function(term2gene) {
  intersect(universe_symbols_all, unique(term2gene$gene))
}
universe_per_db <- lapply(enrichr_genesets, build_universe_for_db)
sapply(universe_per_db, length)   # sanity check: usable universe size per library

# ---- 3. Run enricher() per sample, per direction, per database ----
run_enricher_for_gene_list <- function(gene_list, db, label) {
  message(label, " - ", db)
  tryCatch(
    clusterProfiler::enricher(
      gene          = gene_list,
      universe      = universe_per_db[[db]],
      TERM2GENE     = enrichr_genesets[[db]],
      pAdjustMethod = "BH",
      pvalueCutoff  = 1,
      qvalueCutoff  = 1
    ),
    error = function(e) {
      message(label, " - ", db, ": enricher() failed (", conditionMessage(e), ")")
      return(NULL)
    }
  )
}

run_enricher_all_dbs <- function(gene_lists, direction) {
  result <- lapply(names(gene_lists), function(s) {
    per_db <- lapply(dbs, function(db) run_enricher_for_gene_list(gene_lists[[s]], db, paste(direction, "-", s)))
    names(per_db) <- dbs
    per_db
  })
  names(result) <- names(gene_lists)
  result
}

# Reuses sample_up_genes / sample_down_genes from the earlier up/down section
enrichr_up   <- run_enricher_all_dbs(sample_up_genes,   "Up")
enrichr_down <- run_enricher_all_dbs(sample_down_genes, "Down")

saveRDS(enrichr_up,   "SOM/PerSample_ORA/EnrichR/Up/EnrichR_Up.rds")
saveRDS(enrichr_down, "SOM/PerSample_ORA/EnrichR/Down/EnrichR_Down.rds")

# ---- 4. Keep only significantly enriched terms per sample, per database ----
enrichr_up   <- readRDS("SOM/PerSample_ORA/EnrichR/Up/EnrichR_Up.rds")
enrichr_down <- readRDS("SOM/PerSample_ORA/EnrichR/Down/EnrichR_Down.rds")

flatten_enricher <- function(enricher_list, direction) {
  bind_rows(lapply(names(enricher_list), function(s) {
    message(s)
    bind_rows(lapply(names(enricher_list[[s]]), function(db) {
      res <- enricher_list[[s]][[db]]
      if (is.null(res) || nrow(res@result) == 0) return(NULL)
      res@result %>%
        filter(p.adjust < 0.05) %>%
        mutate(Sample = s, Direction = direction, Database = db)
    }))
  }))
}

enrichr_up_combined   <- flatten_enricher(enrichr_up,   "Up")
enrichr_down_combined <- flatten_enricher(enrichr_down, "Down")
enrichr_all_combined  <- bind_rows(enrichr_up_combined, enrichr_down_combined)

write.csv(enrichr_all_combined, "SOM/PerSample_ORA/EnrichR/EnrichR_UpDown_AllDatabases_combined.csv", row.names = FALSE)


# ==================== BAR PLOTS: TOP TERMS PER SAMPLE x DATABASE (ALL / UP / DOWN) ====================
n_terms_to_show <- 10

# "All" bars are colored by Direction (two-color legend); Up/Down bars keep
# their existing single color per plot.
plot_enricher_barplots <- function(combined_df, label, outdir, color_by_direction = FALSE) {
  for (s in unique(combined_df$Sample)) {
    for (db in unique(combined_df$Database)) {
      df_sub <- combined_df %>%
        filter(Sample == s, Database == db) %>%
        arrange(p.adjust) %>%
        slice_head(n = n_terms_to_show)
      if (nrow(df_sub) == 0) next
      df_sub <- df_sub %>% mutate(Description = factor(Description, levels = rev(Description)))
      
      if (color_by_direction) {
        bar_plot <- ggplot(df_sub, aes(x = Description, y = -log10(p.adjust), fill = Direction)) +
          geom_col() +
          scale_fill_manual(values = c(Up = "#E15759", Down = "#4E79A7")) +
          coord_flip() +
          labs(title = paste0(db, " (", label, ") - ", s),
               x = NULL, y = expression(-log[10](adjusted~p))) +
          common_theme
      } else {
        bar_color <- ifelse(label == "Up", "#E15759", "#4E79A7")
        bar_plot <- ggplot(df_sub, aes(x = Description, y = -log10(p.adjust))) +
          geom_col(fill = bar_color) +
          coord_flip() +
          labs(title = paste0(db, " (", label, ") - ", s),
               x = NULL, y = expression(-log[10](adjusted~p))) +
          common_theme
      }
      
      ggsave(paste0(outdir, "/", s, "_", db, "_", tolower(label), ".png"),
             bar_plot, width = 10, height = 6)
    }
  }
}

dir.create("SOM/PerSample_ORA/EnrichR/All", showWarnings = FALSE)

plot_enricher_barplots(enrichr_all_combined,  "All",  "SOM/PerSample_ORA/EnrichR/All",  color_by_direction = TRUE)
plot_enricher_barplots(enrichr_up_combined,   "Up",   "SOM/PerSample_ORA/EnrichR/Up")
plot_enricher_barplots(enrichr_down_combined, "Down", "SOM/PerSample_ORA/EnrichR/Down")


# ==================== UNIQUE TERMS PER SAMPLE, WITHIN EACH DATABASE (ALL / UP / DOWN) ====================
dir.create("SOM/PerSample_ORA/EnrichR/All/Unique",  showWarnings = FALSE)
dir.create("SOM/PerSample_ORA/EnrichR/Up/Unique",   showWarnings = FALSE)
dir.create("SOM/PerSample_ORA/EnrichR/Down/Unique", showWarnings = FALSE)

find_unique_enricher_terms <- function(combined_df) {
  combined_df %>%
    group_by(Database) %>%
    group_modify(~ {
      term_counts <- .x %>% distinct(Description, Sample) %>% count(Description, name = "n_samples")
      .x %>% left_join(term_counts, by = "Description") %>% filter(n_samples == 1)
    }) %>%
    ungroup()
}

unique_enrichr_all  <- find_unique_enricher_terms(enrichr_all_combined)
unique_enrichr_up   <- find_unique_enricher_terms(enrichr_up_combined)
unique_enrichr_down <- find_unique_enricher_terms(enrichr_down_combined)

write.csv(unique_enrichr_all,  "SOM/PerSample_ORA/EnrichR/All/Unique/unique_terms_all.csv",   row.names = FALSE)
write.csv(unique_enrichr_up,   "SOM/PerSample_ORA/EnrichR/Up/Unique/unique_terms_up.csv",     row.names = FALSE)
write.csv(unique_enrichr_down, "SOM/PerSample_ORA/EnrichR/Down/Unique/unique_terms_down.csv", row.names = FALSE)

plot_unique_enricher_barplots <- function(unique_df, label, outdir, color_by_direction = FALSE) {
  for (s in unique(unique_df$Sample)) {
    for (db in unique(unique_df$Database)) {
      df_sub <- unique_df %>%
        filter(Sample == s, Database == db) %>%
        arrange(p.adjust) %>%
        slice_head(n = n_terms_to_show)
      if (nrow(df_sub) == 0) next
      df_sub <- df_sub %>% mutate(Description = factor(Description, levels = rev(Description)))
      
      if (color_by_direction) {
        unique_plot <- ggplot(df_sub, aes(x = Description, y = -log10(p.adjust), fill = Direction)) +
          geom_col() +
          scale_fill_manual(values = c(Up = "#E15759", Down = "#4E79A7")) +
          coord_flip() +
          labs(title = paste0(db, " unique (", label, ") - ", s),
               x = NULL, y = expression(-log[10](adjusted~p))) +
          common_theme
      } else {
        bar_color <- ifelse(label == "Up", "#E15759", "#4E79A7")
        unique_plot <- ggplot(df_sub, aes(x = Description, y = -log10(p.adjust))) +
          geom_col(fill = bar_color) +
          coord_flip() +
          labs(title = paste0(db, " unique (", label, ") - ", s),
               x = NULL, y = expression(-log[10](adjusted~p))) +
          common_theme
      }
      
      ggsave(paste0(outdir, "/", s, "_", db, "_unique_", tolower(label), ".png"),
             unique_plot, width = 10, height = 6)
    }
  }
}

plot_unique_enricher_barplots(unique_enrichr_all,  "All",  "SOM/PerSample_ORA/EnrichR/All/Unique",  color_by_direction = TRUE)
plot_unique_enricher_barplots(unique_enrichr_up,   "Up",   "SOM/PerSample_ORA/EnrichR/Up/Unique")
plot_unique_enricher_barplots(unique_enrichr_down, "Down", "SOM/PerSample_ORA/EnrichR/Down/Unique")


# ==================== CIRCOS PLOTS: SHARED TERMS BETWEEN SAMPLES (ALL / UP / DOWN), PER DATABASE ====================
library(circlize)
dir.create("SOM/PerSample_ORA/EnrichR/Circos", showWarnings = FALSE)

make_shared_matrix_enricher <- function(combined_df, samples) {
  terms_by_sample <- split(combined_df$Description, combined_df$Sample)
  n <- length(samples)
  m <- matrix(0, nrow = n, ncol = n, dimnames = list(samples, samples))
  for (i in seq_len(n)) for (j in seq_len(n)) if (i != j) {
    a <- terms_by_sample[[samples[i]]]; b <- terms_by_sample[[samples[j]]]
    m[i, j] <- if (is.null(a) || is.null(b)) 0 else length(intersect(a, b))
  }
  m
}

plot_circos_enricher <- function(shared_matrix, samples, title_text, outfile) {
  sample_colors <- setNames(RColorBrewer::brewer.pal(length(samples), "Set2"), samples)
  png(outfile, width = 2200, height = 1800, res = 200)
  par(mar = c(1, 1, 2, 8))
  circos.clear()
  circos.par(canvas.xlim = c(-1.3, 1.3), canvas.ylim = c(-1.3, 1.3))
  chordDiagram(shared_matrix, symmetric = TRUE, grid.col = sample_colors, transparency = 0.4,
               annotationTrack = c("grid"), preAllocateTracks = list(track.height = 0.15))
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1] + 1, CELL_META$sector.index,
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.8)
  }, bg.border = NA)
  title(title_text)
  legend(x = "right", inset = -0.25, legend = samples, fill = sample_colors,
         border = NA, bty = "n", cex = 0.8, title = "Sample", xpd = TRUE)
  circos.clear()
  dev.off()
}

samples <- colnames(counts_norm)

# NOTE: "All" here means shared TERMS regardless of direction -- if a term
# is Up in one sample and Down in another, it still counts as shared here.
# That's an intentional choice (same biological term implicated either way),
# not a bug -- Up/Down circos plots below are direction-specific.
circos_inputs <- list(
  All  = enrichr_all_combined,
  Up   = enrichr_up_combined,
  Down = enrichr_down_combined
)

for (label in names(circos_inputs)) {
  combined_df_full <- circos_inputs[[label]]
  for (db in dbs) {
    df_db <- combined_df_full %>% filter(Database == db)
    if (nrow(df_db) == 0) {
      message(db, " (", label, "): no significant terms, skipping circos")
      next
    }
    shared_mat <- make_shared_matrix_enricher(df_db, samples)
    
    # circlize::chordDiagram() errors internally on an all-zero matrix
    # (no pair of samples shares any term) -- skip rather than crash.
    if (sum(shared_mat) == 0) {
      message(db, " (", label, "): terms exist per-sample but none are shared, skipping circos")
      next
    }
    
    plot_circos_enricher(
      shared_mat, samples,
      paste0("Shared ", db, " Terms (", label, ")"),
      paste0("SOM/PerSample_ORA/EnrichR/Circos/", db, "_", tolower(label), "_circos.png")
    )
  }
}