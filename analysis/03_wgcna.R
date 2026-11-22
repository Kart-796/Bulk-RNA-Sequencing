# ==================== BULK RNA-SEQ: WGCNA ====================
# Input : <bulk_rna_workdir>/gene_counts.csv (alignment/03_count_reads.sh) and
#         <bulk_rna_workdir>/Metadata.csv (Sample, Condition, Age, Sex, Genotype, Type)
# Output: <wgcna_workdir>/  module assignments, module-trait heatmap, gene significance tables
# ---- Configuration: edit config/config.yaml (run from the repository root) ----
source("utils/load_config.R")
counts_file   <- file.path(cfg_path("bulk_rna_workdir"), "gene_counts.csv")
metadata_file <- file.path(cfg_path("bulk_rna_workdir"), "Metadata.csv")
setwd(cfg_path("wgcna_workdir", create = TRUE))

library(readr)
library(WGCNA)
library(DESeq2)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)

# 1. Fetch Data ------------------------------------------------
data <- read_csv(counts_file)
data <- data %>% 
  select(-gene_name) %>% 
  column_to_rownames(var = 'gene_id')

# 2. QC - outlier detection ------------------------------------------------
# detect outlier genes

gsg <- goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)

# remove genes that are detectd as outliers
data <- data[gsg$goodGenes == TRUE,]

# detect outlier samples - hierarchical clustering - method 1
htree <- hclust(dist(t(data)), method = "average")
plot(htree)


# 3. Normalization ----------------------------------------------------------------------
# create a deseq2 dataset

coldata <- read_csv(metadata_file)
# Put the metadata rows in the same order as the count columns
coldata <- coldata[match(colnames(data), coldata$Sample), ]
stopifnot(identical(coldata$Sample, colnames(data)))

dds <- DESeqDataSetFromMatrix(countData = data, colData = coldata, design = ~Condition)
dds <- DESeq(dds)
dds_norm <- vst(dds)
norm.counts <- assay(dds_norm) %>% 
  t()

# 4. Network Construction  ---------------------------------------------------
# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)


sft.data <- sft$fitIndices

# visualization to pick power

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()

print(a1)

a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()

print(a2)

grid.arrange(a1, a2, nrow = 2)

# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

soft_power <- 10
temp_cor <- cor
cor <- WGCNA::cor


# memory estimate w.r.t blocksize
bwnet <- blockwiseModules(norm.counts,
                          maxBlockSize = 16000, #based on RAM size of computer 16GB could be  approx 14000
                          TOMType = "signed", #both positive and negative correlations unsigned is only +ve
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)


cor <- temp_cor

# 5. Module Eigengenes ---------------------------------------------------------
# Merge close modules (those with highly correlated eigengenes)
mergeCutHeight = 0.5  # Adjust this threshold as needed

mergedModules <- mergeCloseModules(
  exprData = norm.counts,            # The normalized counts matrix
  colors = bwnet$colors,             # The initial module assignments from blockwiseModules
  cutHeight = mergeCutHeight,        # Height threshold for merging (modules with eigengene correlation > 0.75 will be merged)
  verbose = 3                        # Verbose output
)

# Update the colors and eigengenes to the merged ones
mergedColors <- mergedModules$colors
mergedMEs <- mergedModules$newMEs

# Plot the dendrogram with merged module colors
plotDendroAndColors(
  bwnet$dendrograms[[1]],            # Dendrogram from the first block
  mergedColors[bwnet$blockGenes[[1]]], # Colors for the genes in the first block after merging
  "Merged Module colors",            # Title for the color strip
  dendroLabels = FALSE,              # Do not display labels on the dendrogram
  hang = 0.03,                       # Control for the position of the leaves
  addGuide = TRUE,                   # Add a color guide on the plot
  guideHang = 0.05                   # Control for the position of the guide
)

# Create a new data frame with the merged module colors
merged_module_df <- data.frame(
  gene_id = names(mergedColors),
  colors = labels2colors(mergedColors)
)

# Preview the merged module data frame
head(merged_module_df)

# Save the merged gene-module assignments to a file
write_delim(merged_module_df,
            file = "merged_gene_modules.csv",
            delim = "\t")


module_eigengenes <- bwnet$MEs   # (unmerged) eigengenes, used for module membership below


# Print out a preview
#head(module_eigengenes)


# get number of genes for each module
#table(bwnet$colors)

#table(bwnet$unmergedColors)

# Convert labels to colors for plotting
#mergedColors = labels2colors(bwnet$colors)
# Plot the dendrogram and the module colors underneath
#plotDendroAndColors(
  #bwnet$dendrograms[[1]],
  #mergedColors[bwnet$blockGenes[[1]]],
 #"Module colors",
 # dendroLabels = FALSE,
  #hang = 0.03,
 # addGuide = TRUE,
 # guideHang = 0.05 )

#module_df <- data.frame(
  #gene_id = names(bwnet$colors),
  #colors = labels2colors(bwnet$colors)
#)

#module_df[1:500000,]

#write_delim(module_df,
        #    file = "gene_modules.csv",
       #    delim = "\t")


# 6A. Relate modules to traits --------------------------------------------------
# module trait associations



# create traits file - binarize categorical variables
design <- model.matrix(~0 + as.factor(coldata$Condition))
original_colnames <- colnames(design)
new_colnames <- gsub("^as.factor\\(coldata\\$Condition\\)", "", original_colnames)
colnames(design) <- new_colnames

Age <- model.matrix(~0 + as.factor(coldata$Age))
original_colnames <- colnames(Age)
new_colnames <- gsub("^as.factor\\(coldata\\$Age\\)", "", original_colnames)
colnames(Age) <- new_colnames

Sex <- model.matrix(~0 + as.factor(coldata$Sex))
original_colnames <- colnames(Sex)
new_colnames <- gsub("^as.factor\\(coldata\\$Sex\\)", "", original_colnames)
colnames(Sex) <- new_colnames

Genotype <- model.matrix(~0 + as.factor(coldata$Genotype))
original_colnames <- colnames(Genotype)
new_colnames <- gsub("^as.factor\\(coldata\\$Genotype\\)", "", original_colnames)
colnames(Genotype) <- new_colnames

Type <- model.matrix(~0 + as.factor(coldata$Type))
original_colnames <- colnames(Type)
new_colnames <- gsub("^as.factor\\(coldata\\$Type\\)", "", original_colnames)
colnames(Type) <- new_colnames

traits <- cbind(design, Age, Sex, Type, Genotype)
new_row_names <- rownames(t(data))
rownames(traits) <- new_row_names

# Define numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)

#module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
#module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

# Calculate module-trait correlations with merged modules
module.trait.corr <- cor(mergedMEs, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

# visualize module-trait association as a heatmap

#heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')

# Combine merged module eigengenes with traits for heatmap visualization
# Combine merged module eigengenes with traits for heatmap visualization
heatmap.data <- merge(mergedMEs, traits, by = 'row.names')

# Preview the combined data
head(heatmap.data)

# Set the row names to the sample names
heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')

# Define the columns for traits and merged module eigengenes
trait_columns <- names(heatmap.data)[(ncol(mergedMEs) + 1):ncol(heatmap.data)]
module_columns <- names(heatmap.data)[1:ncol(mergedMEs)]

# Visualize module-trait association as a heatmap
png("merged_plot.png", width = 1600, height = 2600)  # Set width and height in pixels
Corrplot <- CorLevelPlot(heatmap.data,
                         x = trait_columns,    # Traits go on the x-axis
                         y = module_columns,   # Merged module eigengenes on the y-axis
                         col = c("blue1", "skyblue", "white", "pink", "red"))
dev.off()
print(Corrplot)

Control.Vehicle1.module.gene.mapping <- as.data.frame(bwnet$colors)
Control.Vehicle1.module.gene.mapping %>% 
  filter(`bwnet$colors` == 'darkgreen') %>% 
  rownames()

Control.Vehicle1.module.membership.measure <- cor(module_eigengenes, norm.counts, use = 'p')
Control.Vehicle1.module.membership.measure.pvals <- corPvalueStudent(Control.Vehicle1.module.membership.measure, nSamples)

Control.Vehicle1.module.membership.measure.pvals[1:10,1:10]

Control.Vehicle1.gene.signf.corr <- cor(norm.counts, traits[, "Control_Vehicle"], use = 'p')
Control.Vehicle1.gene.signf.corr.pvals <- corPvalueStudent(Control.Vehicle1.gene.signf.corr, nSamples)

Control.Vehicle1.gene.signf.corr.pvals %>% 
  as.data.frame() %>% 
  filter(V1 < 0.05) %>% 
  arrange(V1)
  head(500000)

write.csv(Control.Vehicle1.gene.signf.corr.pvals, file = "Control_Vehicle1_gene_signf_corr_pvals.csv", row.names = TRUE)

Control.Treated1.module.gene.mapping <- as.data.frame(bwnet$colors)
Control.Treated1.module.gene.mapping %>% 
  filter(bwnet$colors == 'violet') %>% 
  rownames()

Control.Treated1.module.membership.measure <- cor(module_eigengenes, norm.counts, use = 'p')
Control.Treated1.module.membership.measure.pvals <- corPvalueStudent(Control.Treated1.module.membership.measure, nSamples)

Control.Treated1.module.membership.measure.pvals[1:10,1:10]

Control.Treated1.gene.signf.corr <- cor(norm.counts, traits[, "Control_Treated"], use = 'p')
Control.Treated1.gene.signf.corr.pvals <- corPvalueStudent(Control.Treated1.gene.signf.corr, nSamples)

Control.Treated1.gene.signf.corr.pvals %>% 
  as.data.frame() %>% 
  filter(V1 < 0.05) %>% 
  arrange(V1)
  head(500000)

write.csv(Control.Treated1.gene.signf.corr.pvals, file = "Control_Treated1_gene_signf_corr_pvals.csv", row.names = TRUE)

MutA.Vehicle.module.gene.mapping <- as.data.frame(bwnet$colors)
MutA.Vehicle.module.gene.mapping %>% 
  filter(`bwnet$colors` == 'yellowgreen') %>% 
  rownames()

MutA.Vehicle.module.membership.measure <- cor(module_eigengenes, norm.counts, use = 'p')
MutA.Vehicle.module.membership.measure.pvals <- corPvalueStudent(MutA.Vehicle.module.membership.measure, nSamples)

MutA.Vehicle.module.membership.measure.pvals[1:10,1:10]

MutA.Vehicle.gene.signf.corr <- cor(norm.counts, traits[, "MutA_Vehicle"], use = 'p')
MutA.Vehicle.gene.signf.corr.pvals <- corPvalueStudent(MutA.Vehicle.gene.signf.corr, nSamples)

MutA.Vehicle.gene.signf.corr.pvals %>% 
  as.data.frame() %>% 
  filter(V1 < 0.05) %>% 
  arrange(V1)
  head(500000)

write.csv(MutA.Vehicle.gene.signf.corr.pvals, file = "MutA_Vehicle_gene_signf_corr_pvals.csv", row.names = TRUE)


MutA.Treated.module.gene.mapping <- as.data.frame(bwnet$colors)
MutA.Treated.module.gene.mapping %>% 
  filter(`bwnet$colors` == 'lightyellow') %>% 
  rownames()

MutA.Treated.module.membership.measure <- cor(module_eigengenes, norm.counts, use = 'p')
MutA.Treated.module.membership.measure.pvals <- corPvalueStudent(MutA.Treated.module.membership.measure, nSamples)


MutA.Treated.module.membership.measure.pvals[1:10,1:10]

MutA.Treated.gene.signf.corr <- cor(norm.counts, traits[, "MutA_Treated"], use = 'p')
MutA.Treated.gene.signf.corr.pvals <- corPvalueStudent(MutA.Treated.gene.signf.corr, nSamples)

MutA.Treated.gene.signf.corr.pvals %>% 
  as.data.frame() %>% 
  filter(V1 < 0.05) %>% 
  arrange(V1)
  head(500000)

write.csv(MutA.Treated.gene.signf.corr.pvals, file = "MutA_Treated_gene_signf_corr_pvals.csv", row.names = TRUE)


MutB.Vehicle.module.gene.mapping <- as.data.frame(bwnet$colors)
MutB.Vehicle.module.gene.mapping %>% 
  filter(`bwnet$colors` == 'darkorange2') %>% 
  rownames()

MutB.Vehicle.module.membership.measure <- cor(module_eigengenes, norm.counts, use = 'p')
MutB.Vehicle.module.membership.measure.pvals <- corPvalueStudent(MutB.Vehicle.module.membership.measure, nSamples)

MutB.Vehicle.module.membership.measure.pvals[1:10,1:10]

MutB.Vehicle.gene.signf.corr <- cor(norm.counts, traits[, "MutB_Vehicle"], use = 'p')
MutB.Vehicle.gene.signf.corr.pvals <- corPvalueStudent(MutB.Vehicle.gene.signf.corr, nSamples)

MutB.Vehicle.gene.signf.corr.pvals %>% 
  as.data.frame() %>% 
  filter(V1 < 0.05) %>% 
  arrange(V1)
  head(500000)

write.csv(MutB.Vehicle.gene.signf.corr.pvals, file = "MutB_Vehicle_gene_signf_corr_pvals.csv", row.names = TRUE)


MutB.Treated.module.gene.mapping <- as.data.frame(bwnet$colors)
MutB.Treated.module.gene.mapping %>% 
  filter(`bwnet$colors` == 'brown4') %>% 
  rownames()

MutB.Treated.module.membership.measure <- cor(module_eigengenes, norm.counts, use = 'p')
MutB.Treated.module.membership.measure.pvals <- corPvalueStudent(MutB.Treated.module.membership.measure, nSamples)

MutB.Treated.module.membership.measure.pvals[1:10,1:10]

MutB.Treated.gene.signf.corr <- cor(norm.counts, traits[, "MutB_Treated"], use = 'p')
MutB.Treated.gene.signf.corr.pvals <- corPvalueStudent(MutB.Treated.gene.signf.corr, nSamples)

MutB_Treated_gene_signf_corr_pvals <- MutB.Treated.gene.signf.corr.pvals %>% 
  as.data.frame() %>% 
  filter(V1 < 0.05) %>% 
  arrange(V1)
  head(500000)

write.csv(MutB_Treated_gene_signf_corr_pvals, file = "MutB_Treated_gene_signf_corr_pvals.csv", row.names = TRUE)



