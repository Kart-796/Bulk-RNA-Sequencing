# =============================================================================
# Install the R packages used by analysis/01-03.
#   Rscript envs/install_r_packages.R
# Tip: afterwards run renv::init() / renv::snapshot() and commit renv.lock.
# =============================================================================
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("remotes", quietly = TRUE))     install.packages("remotes")

cran <- c("yaml", "tidyverse", "magrittr", "ggplot2", "ggrepel", "gridExtra", "hrbrthemes",
          "RColorBrewer", "pheatmap", "plotly", "VennDiagram", "UpSetR", "ellipse", "combinat",
          "readxl", "naniar", "kohonen", "enrichR", "msigdbr", "igraph", "reshape2", "patchwork",
          "circlize", "futile.logger", "imputeLCMD", "WGCNA")
bioc <- c("limma", "edgeR", "sva", "fgsea", "clusterProfiler", "enrichplot", "ReactomePA",
          "org.Hs.eg.db", "preprocessCore", "pcaMethods", "STRINGdb", "DESeq2", "DOSE",
          "GO.db", "impute", "preprocessCore")          # GO.db / impute: WGCNA dependencies

BiocManager::install(setdiff(unique(bioc), rownames(installed.packages())), update = FALSE, ask = FALSE)
install.packages(setdiff(cran, rownames(installed.packages())))
if (!requireNamespace("CorLevelPlot", quietly = TRUE)) remotes::install_github("kevinblighe/CorLevelPlot")
