## Activate conda environment externally if needed (not in R script)
## mamba activate R4.4 (R version 4.3)

# Description:
# WM_line_Lesion_1 : Iron Positive
# WM_line_Lesion_2 : Iron Negative


## mamba activate R4.4

library(edgeR)
library(dplyr)
library(ComplexHeatmap)

f <- "RawCounts_MS177-6-5Only_4slidesOnly_ROIbyGene_CORRECTED_AnnotAdded.txt"

# Read data and preprocess
lines <- readLines(f, n = 3)
g <- strsplit(gsub(" ", "_", lines[1]), "\t")[[1]][-1]
pa <- strsplit(gsub(" ", "_", lines[2]), "\t")[[1]][-1]

tt <- read.table(f, header = TRUE, skip = 2)
colnames(tt) <- gsub(" ", "_", colnames(tt))
rownames(tt) <- tt$gene

# Aggregate duplicate gene entries by mean
x <- aggregate(. ~ gene, data = tt, FUN = mean)
rownames(x) <- x$gene
x <- x[, -1]

# Convert group labels to factor
g <- as.factor(g)

design <- model.matrix(~ 0 + g)  # No intercept, one column per group
colnames(design) <- levels(g)   # Optional: name columns based on group levels

# Prepare the DGEList and filter
y = DGEList(counts=x, group=g)
keep = filterByExpr(y, group=g)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)  # Correct normalization step
y <- estimateDisp(y, design)

m = cpm(y, normalized.lib.sizes=TRUE)


fit <- glmQLFit(y, design)

# Function to perform contrast and extract significant genes using QL test
get_significant_genes <- function(contrast_name) {
  contrast <- makeContrasts(contrasts = contrast_name, levels = design)
  qlf <- glmQLFTest(fit, contrast = contrast)
  topTags(qlf, adjust.method = "fdr", n = Inf)$table %>%
    as.data.frame() %>%
    filter(FDR < 0.05)
}

# Define all pairwise comparisons
comparisons <- list(
  "WM_line_Lesion_1 - WM_line_Lesion_2",
  "NAWM - WM_line_Lesion_1",
  "NAWM - WM_line_Lesion_2",
  "WML_Lesion_1 - WML_Lesion_2",
  "NAWM - WML_Lesion_1",
  "NAWM - WML_Lesion_2",
  "WM_line_Lesion_1 - WML_Lesion_1",
  "WM_line_Lesion_1 - WML_Lesion_2",
  "WM_line_Lesion_2 - WML_Lesion_1",
  "WM_line_Lesion_2 - WML_Lesion_2"
)

# Perform all pairwise comparisons
results <- lapply(comparisons, get_significant_genes)

# Combine results and include comparison labels
combined_results <- do.call(rbind, lapply(seq_along(results), function(i) {
  df <- results[[i]]
  df$comparison <- comparisons[i]
  df$rowname <- row.names(df)
  row.names(df) <- NULL  # Reset row names to default numeric indices
  df
}))

# View the combined results
print(combined_results)

write.csv(combined_results, file = "Combined_DE_Results.csv", row.names = FALSE)

#### All gene's results.

 
get_all_genes <- function(contrast_name) {
  contrast <- makeContrasts(contrasts = contrast_name, levels = design)
  qlf <- glmQLFTest(fit, contrast = contrast)
  topTags(qlf, adjust.method = "fdr", n = Inf)$table %>%
    as.data.frame() %>%
    mutate(gene = rownames(.), comparison = contrast_name)
}
# Get full differential expression results for all comparisons
all_results <- lapply(comparisons, get_all_genes)

# Combine into a single data frame
combined_all_results <- do.call(rbind, all_results)
row.names(combined_all_results) <- NULL
write.csv(combined_all_results, file = "All_DE_Results.csv", row.names = FALSE)

 


#######################################
#### 
####        USE THIS !!!!  Log2-Transformed Group Means Correlation Heatmap
####
#######################################

# Load required library
library(pheatmap)

# Step 1: Compute log2-transformed CPM values
log_cpm <- cpm(y, log = TRUE)  # y should be TMM-normalized using calcNormFactors(y)

# Step 2: Compute log2-transformed group means
group_means_log <- t(sapply(levels(g), function(group) {
  group_samples <- which(g == group)
  rowMeans(log_cpm[, group_samples])
}))
group_means_log <- t(group_means_log)  # Transpose so rows are genes, columns are groups
colnames(group_means_log) <- levels(g)

# Step 3: Compute correlation matrix
group_correlations_log <- cor(group_means_log)

# Step 4: Plot heatmap
pheatmap(group_correlations_log,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "complete",
         color = colorRampPalette(c("white", "yellow", "red"))(100),  # Gradient from low to high
         na_col = "grey",         # Set NA values (if any) to grey
         breaks = seq(0, 1, length.out = 101),  # For correlation values 0–1
         main = "Correlation Heatmap of Group Means")


 
 ########### PCA Plot ###########
 # Define color mapping
 
 group_colors <- c(
  "NAWM" = "#1B9E77", 
  "WM_line_Lesion_1" = "#377EB8", 
  "WM_line_Lesion_2" = "#E41A1C",
  "WML_Lesion_1" = "#984EA3", 
  "WML_Lesion_2" = "#FF7F00"
)
 # Perform PCA on normalized expression values (CPM for each sample)
#y <- calcNormFactors(y)  # Ensure normalization factors are applied
cpm_data <- cpm(y, log = TRUE)  # Use log-transformed CPM values for PCA

# Perform PCA on the normalized CPM values for each sample
pca_result <- prcomp(t(cpm_data), scale. = TRUE)  # Transpose to get samples as rows and genes as columns

# Create a data frame with PCA results for plotting
pca_data <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  Sample = colnames(cpm_data),
  Group = g  # Group information (factor)
)

# Plot PCA using ggplot2
library(ggplot2)
ggplot(pca_data, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(size = 5) +
#  geom_text(vjust = 1.5, size = 3) +  # Optionally label points with sample names
  labs(title = "PCA of Individual Samples", x = "PC1", y = "PC2") +
  theme_minimal() +
  scale_color_manual(values = group_colors)

###############################
#### PCA plot only with NAWM, WM_line_Lesion1 and WM_line_Lesion2

# Filter the samples to include only NAWM, WM_line_Lesion1, and WM_line_Lesion2
selected_groups <- c("NAWM", "WM_line_Lesion_1", "WM_line_Lesion_2")
selected_samples <- which(g %in% selected_groups)

# Subset the DGEList to include only the selected samples
y_selected <- y[, selected_samples]

# Normalize and log-transform the counts
y_selected <- calcNormFactors(y_selected)
cpm_data_selected <- cpm(y_selected, log = TRUE)

# Perform PCA on the normalized CPM values for selected samples
pca_result_selected <- prcomp(t(cpm_data_selected), scale. = TRUE)

# Create a data frame with PCA results for plotting
pca_data_selected <- data.frame(
  PC1 = pca_result_selected$x[, 1],
  PC2 = pca_result_selected$x[, 2],
  Sample = colnames(cpm_data_selected),
  Group = g[selected_samples]  # Corresponding group labels
)

# Plot PCA using ggplot2
library(ggplot2)
ggplot(pca_data_selected, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(size = 5) +
#  geom_text(vjust = 1.5, size = 3) +  # Optionally label points with sample names
  labs(title = "PCA of Selected Groups (NAWM, WM_line_Lesion_1, WM_line_Lesion_2)", x = "PC1", y = "PC2") +
  theme_minimal() +
  scale_color_manual(values = group_colors)


###############################
#### PCA plot only with NAWM, WML_Lesion1 and WML_Lesion2

# Filter the samples to include only NAWM, WM_line_Lesion1, and WM_line_Lesion2
selected_groups <- c("NAWM", "WML_Lesion_1", "WML_Lesion_2")
selected_samples <- which(g %in% selected_groups)

# Subset the DGEList to include only the selected samples
y_selected <- y[, selected_samples]

# Normalize and log-transform the counts
y_selected <- calcNormFactors(y_selected)
cpm_data_selected <- cpm(y_selected, log = TRUE)

# Perform PCA on the normalized CPM values for selected samples
pca_result_selected <- prcomp(t(cpm_data_selected), scale. = TRUE)

# Create a data frame with PCA results for plotting
pca_data_selected <- data.frame(
  PC1 = pca_result_selected$x[, 1],
  PC2 = pca_result_selected$x[, 2],
  Sample = colnames(cpm_data_selected),
  Group = g[selected_samples]  # Corresponding group labels
)

# Plot PCA using ggplot2
library(ggplot2)
ggplot(pca_data_selected, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(size = 5) +
#  geom_text(vjust = 1.5, size = 3) +  # Optionally label points with sample names
  labs(title = "PCA of Selected Groups (NAWM, WML_Lesion1, WML_Lesion2)", x = "PC1", y = "PC2") +
  theme_minimal() +
  scale_color_manual(values = group_colors)

 
####################################################
############
############            Clustering within Gene Categories (common DEGs, mutual exclusive DEGs)
###########
####################################################

library(pheatmap)
library(dplyr)

#y <- calcNormFactors(y)  # Ensure normalization factors are applied
#cpm_data <- cpm(y, log = TRUE)  # Use log-transformed CPM values for PCA

# Extract DEGs and add gene names as a column
degs_nawm_wm1 <- get_significant_genes("NAWM - WM_line_Lesion_1") %>%
  mutate(gene = rownames(.))

degs_nawm_wm2 <- get_significant_genes("NAWM - WM_line_Lesion_2") %>%
  mutate(gene = rownames(.))

# Identify common and unique DEGs
common_degs <- intersect(degs_nawm_wm1$gene, degs_nawm_wm2$gene)
unique_nawm_wm1 <- setdiff(degs_nawm_wm1$gene, degs_nawm_wm2$gene)
unique_nawm_wm2 <- setdiff(degs_nawm_wm2$gene, degs_nawm_wm1$gene)

# Subset expression matrix to include only NAWM, WM_line_Lesion_1, WM_line_Lesion_2
selected_groups <- c("NAWM", "WM_line_Lesion_1", "WM_line_Lesion_2")
selected_samples <- which(g %in% selected_groups)

# Order columns by group for visual clarity
ordered_samples <- order(factor(g[selected_samples], levels = selected_groups))
expr_subset <- cpm_data[, selected_samples[ordered_samples]]

# Extract expression for DEGs and order rows by category
heatmap_data_common <- expr_subset[common_degs, ]
heatmap_data_nawm_wm1 <- expr_subset[unique_nawm_wm1, ]
heatmap_data_nawm_wm2 <- expr_subset[unique_nawm_wm2, ]

# Scale expression values row-wise
heatmap_data_common_scaled <- t(scale(t(heatmap_data_common)))
heatmap_data_nawm_wm1_scaled <- t(scale(t(heatmap_data_nawm_wm1)))
heatmap_data_nawm_wm2_scaled <- t(scale(t(heatmap_data_nawm_wm2)))

# Perform hierarchical clustering within each category
if (nrow(heatmap_data_common_scaled) > 1) {
  common_clust <- hclust(dist(heatmap_data_common_scaled))$order
  heatmap_data_common_scaled <- heatmap_data_common_scaled[common_clust, ]
}

if (nrow(heatmap_data_nawm_wm1_scaled) > 1) {
  nawm_wm1_clust <- hclust(dist(heatmap_data_nawm_wm1_scaled))$order
  heatmap_data_nawm_wm1_scaled <- heatmap_data_nawm_wm1_scaled[nawm_wm1_clust, ]
}

if (nrow(heatmap_data_nawm_wm2_scaled) > 1) {
  nawm_wm2_clust <- hclust(dist(heatmap_data_nawm_wm2_scaled))$order
  heatmap_data_nawm_wm2_scaled <- heatmap_data_nawm_wm2_scaled[nawm_wm2_clust, ]
}

# Combine scaled data after clustering
heatmap_data_scaled <- rbind(heatmap_data_common_scaled, 
                             heatmap_data_nawm_wm1_scaled, 
                             heatmap_data_nawm_wm2_scaled)

# Create row annotation for gene category
gene_category <- c(rep("Common DEGs", nrow(heatmap_data_common_scaled)),
                   rep("Unique NAWM vs WM_line_Lesion_1", nrow(heatmap_data_nawm_wm1_scaled)),
                   rep("Unique NAWM vs WM_line_Lesion_2", nrow(heatmap_data_nawm_wm2_scaled)))

gene_annotation <- data.frame(Category = gene_category)
rownames(gene_annotation) <- rownames(heatmap_data_scaled)

# Create column annotation for sample groups
column_annotation <- data.frame(Group = g[selected_samples][ordered_samples])
rownames(column_annotation) <- colnames(heatmap_data_scaled)

# Generate heatmap with clustered rows within each category
pheatmap(heatmap_data_scaled,
         cluster_rows = FALSE,  # Clustering is done manually
         cluster_cols = FALSE,  # No sample clustering
         annotation_row = gene_annotation, 
         annotation_col = column_annotation,  
         gaps_row = c(nrow(heatmap_data_common_scaled), 
                      nrow(heatmap_data_common_scaled) + nrow(heatmap_data_nawm_wm1_scaled)),  # Row splits
         gaps_col = c(sum(g[selected_samples][ordered_samples] == "NAWM"), 
                      sum(g[selected_samples][ordered_samples] %in% c("NAWM", "WM_line_Lesion_1"))),  # Column splits
         scale = "none",
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100))

################ End of Heatmap of Categori clsutering ####################

# Extract FDR and logFC for each gene category
degs_nawm_wm1 <- get_significant_genes("NAWM - WM_line_Lesion_1") %>%
  mutate(gene = rownames(.), Category = "Unique NAWM vs WM_line_Lesion_1")

degs_nawm_wm2 <- get_significant_genes("NAWM - WM_line_Lesion_2") %>%
  mutate(gene = rownames(.), Category = "Unique NAWM vs WM_line_Lesion_2")

common_degs_data <- degs_nawm_wm1 %>%
  filter(gene %in% common_degs) %>%
  mutate(Category = "Common DEGs")

unique_nawm_wm1_data <- degs_nawm_wm1 %>%
  filter(gene %in% unique_nawm_wm1)

unique_nawm_wm2_data <- degs_nawm_wm2 %>%
  filter(gene %in% unique_nawm_wm2)


####  still working on
common1 <- degs_nawm_wm1 %>% filter(gene %in% common_degs) %>%
  rename(FDR_1 = FDR, logFC_1 = logFC)

common2 <- degs_nawm_wm2 %>% filter(gene %in% common_degs) %>%
  rename(FDR_2 = FDR, logFC_2 = logFC)

# Merge both sets of results by gene
common_degs_data <- inner_join(common1, common2, by = "gene") %>%
  mutate(Category = "Common DEGs")


##############

# Combine all categories
all_degs <- bind_rows(common_degs_data, unique_nawm_wm1_data, unique_nawm_wm2_data)

# Save to CSV file
write.csv(all_degs, "DEGs_with_FDR_FC_Categories_CORRECTED.csv", row.names = FALSE)



