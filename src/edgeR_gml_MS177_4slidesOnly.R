## mamba activate R4.4 (R4.3)

## WM_line_Lesion_1 : Iron Positive
## WM_line_Lesion_2 : Iron Negative

#f <- "RawCounts_MS177-6-5Only_ROIbyGene_AnnotAdded.txt"
f <- "RawCounts_MS177-6-5_4slides_ROIbyGene_AnnotAdded.txt"
#f <- "RawCounts_MS177-6-5_4slides_ROIbyGene_AnnotAdded.txt"
tt=read.table(f,header=T,skip=2); colnames(tt)=gsub(" ","_",colnames(tt)); row.names(tt)=tt$gene;
lines <- readLines(f, n = 3)
g=strsplit(gsub(" ","_",lines[1]),"\t")[[1]]; g=g[2:length(g)]
pa=strsplit(gsub(" ","_",lines[2]),"\t")[[1]]; pa=pa[2:length(pa)]
gene <- tt$gene;
x=aggregate(. ~ gene, data = tt, FUN = mean); row.names(x)=x$gene; x=x[,-1];

library(edgeR)
library(limma)
g=as.factor(g);
y = DGEList(counts=x,group=g)
keep = filterByExpr(y, group=g)
y = normLibSizes(y)
# Load necessary libraries
library(edgeR)
library(dplyr)  # or you can use library(magrittr)

# Now the rest of your code should work as expected with %>%
g=as.factor(g);
# Create DGEList object
y <- DGEList(counts = as.matrix(x), group = g)

# Filter lowly expressed genes
keep <- filterByExpr(y, group = g)
y <- y[keep, , keep.lib.sizes=FALSE]

# Normalize the library sizes
y <- calcNormFactors(y)

# Design matrix based on group
design <- model.matrix(~0 + g)

y <- estimateDisp(y, design)

# Now fit the model
fit <- glmQLFit(y, design)

# Function to perform contrast and extract significant genes
get_significant_genes <- function(contrast_name) {
  contrast <- makeContrasts(contrasts = contrast_name, levels = design)
  lrt <- glmLRT(fit, contrast = contrast)
  topTags(lrt, adjust.method = "fdr", n = Inf)$table %>%
    as.data.frame() %>%
    #filter(FDR < 0.05 & abs(logFC) >= log(1.5))
    filter(FDR < 0.05)
    #filter( abs(logFC) >= log(0))
}

# Define all pairwise comparisons
comparisons <- list(
  "gWM_line_Lesion_1 - gWM_line_Lesion_2",
  "gNAWM - gWM_line_Lesion_1",
  "gNAWM - gWM_line_Lesion_2"
#  "gWML_Lesion_1 - gWML_Lesion_2",
#  "gNAWM - gWML_Lesion_1",
#  "gNAWM - gWML_Lesion_2",
#  "gWM_line_Lesion_1 - gWML_Lesion_1",
#  "gWM_line_Lesion_1 - gWML_Lesion_2",
#  "gWM_line_Lesion_2 - gWML_Lesion_1",
#  "gWM_line_Lesion_2 - gWML_Lesion_2"
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
 
 
# Output combined results
write.table(as.matrix(combined_results), file = "NAWMnLine_pairwise_comparisons_results_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
# Optional: Generate correlation heatmap of group averages
# Compute group means
group_means <- t(sapply(levels(g), function(group) {
  group_samples <- which(g == group)
  rowMeans(cpm(y)[, group_samples])
}))
group_means <-t(group_means)
colnames(group_means) <- levels(g)
# Compute correlations
group_correlations <- cor(group_means)
# Plot heatmap
library(pheatmap)
pheatmap(group_correlations,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "complete",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Correlation Heatmap of Groups")
# Load necessary library
color_palette <- colorRampPalette(c("white", "yellow", "red"))(100)
pheatmap(group_correlations,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "complete",
         color = color_palette,  # Custom palette from white to red
         na_col = "grey",  # Set the color for NA values to grey
         breaks = seq(0, 1, length.out = 101),  # Scale from 0 to 1
         main = "Correlation Heatmap of Groups")
 
 
 ########### PCA Plot ###########
 
 # Perform PCA on normalized expression values (CPM for each sample)
y <- calcNormFactors(y)  # Ensure normalization factors are applied
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
  scale_color_brewer(palette = "Set1")

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
  scale_color_brewer(palette = "Set1")



#### modify colors #####
# Load necessary library
library(ggplot2)

# Define color mapping
#group_colors <- c("NAWM" = "green", "WM_line_Lesion_1" = "blue", "WM_line_Lesion_2" = "red")
group_colors <- c("NAWM" = "#1B9E77", "WM_line_Lesion_1" = "#377EB8", "WM_line_Lesion_2" = "#E41A1C")

# Plot PCA using ggplot2 with manually assigned colors
ggplot(pca_data_selected, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(size = 4) +
#  stat_ellipse(type = "t", linetype = "dashed", size = 1) +  # Confidence ellipses
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
  scale_color_brewer(palette = "Set1")



#### modify colors #####
# Load necessary library
library(ggplot2)

pca_data_selected <- data.frame(
  PC1 = pca_result_selected$x[, 1],
  PC2 = pca_result_selected$x[, 2],
  Sample = colnames(cpm_data_selected),
  Group = g[selected_samples]  # Corresponding group labels
)

# Define color mapping
#group_colors <- c("NAWM" = "green", "WM_line_Lesion_1" = "blue", "WM_line_Lesion_2" = "red")
group_colors <- c("NAWM" = "#1B9E77", "WML_Lesion_1" = "#377EB8", "WML_Lesion_2" = "#E41A1C")

# Plot PCA using ggplot2 with manually assigned colors
ggplot(pca_data_selected, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(size = 4) +
#  stat_ellipse(type = "t", linetype = "dashed", size = 1) +  # Confidence ellipses
#  geom_text(vjust = 1.5, size = 3) +  # Optionally label points with sample names
  labs(title = "PCA of Selected Groups (NAWM, WML_Lesion_1, WML_Lesion_2)", x = "PC1", y = "PC2") +
  theme_minimal() +
  scale_color_manual(values = group_colors)




########################################################################
#####       Split Heatmap With Grouped Columns

########################################################################
library(pheatmap)
library(dplyr)

y <- calcNormFactors(y)  # Ensure normalization factors are applied
cpm_data <- cpm(y, log = TRUE)  # Use log-transformed CPM values for PCA





# Extract DEGs and add gene names as a column
degs_nawm_wm1 <- get_significant_genes("gNAWM - gWM_line_Lesion_1") %>%
  mutate(gene = rownames(.))

degs_nawm_wm2 <- get_significant_genes("gNAWM - gWM_line_Lesion_2") %>%
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
#expr_subset <- x[, selected_samples[ordered_samples]]
expr_subset <- cpm_data[, selected_samples[ordered_samples]]

# Extract expression for DEGs and order rows by category
ordered_genes <- c(common_degs, unique_nawm_wm1, unique_nawm_wm2)
heatmap_data <- expr_subset[ordered_genes, ]

# Scale expression values row-wise
heatmap_data_scaled <- t(scale(t(heatmap_data)))

# Create row annotation for gene category
gene_category <- c(rep("Common DEGs", length(common_degs)),
                   rep("Unique NAWM vs WM_line_Lesion_1", length(unique_nawm_wm1)),
                   rep("Unique NAWM vs WM_line_Lesion_2", length(unique_nawm_wm2)))
gene_annotation <- data.frame(Category = gene_category)
rownames(gene_annotation) <- ordered_genes

# Create column annotation for sample groups
column_annotation <- data.frame(Group = g[selected_samples][ordered_samples])
rownames(column_annotation) <- colnames(heatmap_data_scaled)

# Generate heatmap with both row and column splits
pheatmap( #t(scale(t(log(heatmap_data_scaled + 10)))),
		  heatmap_data_scaled,
         cluster_rows = FALSE,  # No gene clustering
         cluster_cols = FALSE,  # No sample clustering
         annotation_row = gene_annotation, 
         annotation_col = column_annotation,  
         gaps_row = c(length(common_degs), length(common_degs) + length(unique_nawm_wm1)),  # Row splits
         gaps_col = c(sum(g[selected_samples][ordered_samples] == "NAWM"), 
                      sum(g[selected_samples][ordered_samples] %in% c("NAWM", "WM_line_Lesion_1"))),  # Column splits
         scale = "row",
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100)) 
 
 
####################################################
############
############            Clustering within Gene Categories (common DEGs, mutual exclusive DEGs)
###########
####################################################

library(pheatmap)
library(dplyr)

y <- calcNormFactors(y)  # Ensure normalization factors are applied
cpm_data <- cpm(y, log = TRUE)  # Use log-transformed CPM values for PCA

# Extract DEGs and add gene names as a column
degs_nawm_wm1 <- get_significant_genes("gNAWM - gWM_line_Lesion_1") %>%
  mutate(gene = rownames(.))

degs_nawm_wm2 <- get_significant_genes("gNAWM - gWM_line_Lesion_2") %>%
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
         scale = "row",
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100))

################ End of Heatmap of Categori clsutering ####################

# Extract FDR and logFC for each gene category
degs_nawm_wm1 <- get_significant_genes("gNAWM - gWM_line_Lesion_1") %>%
  mutate(gene = rownames(.), Category = "Unique NAWM vs WM_line_Lesion_1")

degs_nawm_wm2 <- get_significant_genes("gNAWM - gWM_line_Lesion_2") %>%
  mutate(gene = rownames(.), Category = "Unique NAWM vs WM_line_Lesion_2")

common_degs_data <- degs_nawm_wm1 %>%
  filter(gene %in% common_degs) %>%
  mutate(Category = "Common DEGs")

unique_nawm_wm1_data <- degs_nawm_wm1 %>%
  filter(gene %in% unique_nawm_wm1)

unique_nawm_wm2_data <- degs_nawm_wm2 %>%
  filter(gene %in% unique_nawm_wm2)

# Combine all categories
all_degs <- bind_rows(common_degs_data, unique_nawm_wm1_data, unique_nawm_wm2_data)

# Save to CSV file
write.csv(all_degs, "DEGs_with_FDR_FC_Categories.csv", row.names = FALSE)

