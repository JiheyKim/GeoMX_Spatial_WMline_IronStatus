###### mamba activate R4.4 #######

library(readr)
library(dplyr)
library(edgeR)

# Step 1: Read your expression matrix (assuming tab-separated values)
expr <- read.delim("5Only_4slidesOnly_ROIbyGene_CORRECTED_AnnotAdded.txt", header = FALSE, stringsAsFactors = FALSE)


# Step 2: Assign the first three rows as metadata
region <- as.character(expr[1, -1])
slide <- as.character(expr[2, -1])
sampleID <- as.character(expr[3, -1])
gene_id <- expr[3:nrow(expr), 1]
data <- expr[3:nrow(expr), -1]


# Clean up expression data
colnames(data) <- paste0(region, "-", sampleID, "_", 1:ncol(data))  # Make unique sample names
rownames(data) <- gene_id

# Convert expression to numeric
data <- as.data.frame(sapply(data, as.numeric))
rownames(data) <- gene_id

# Step 3: Load your gene list
# From IPA 

gene_list <- readLines("IL-12_DEGs.txt")
#gene_list <- readLines("CellularResponseToHeatStress_DEGs.txt")
#gene_list <- readLines("ComplementCascade_DEGs.txt")
#gene_list <- readLines("S100Family_DEGs.txt")
#gene_list <- readLines("Neuroinflammation_DEGs.txt")
#gene_list <- readLines("MacrophageAlternative_DEGs.txt")

# Step 4: Filter genes
filtered_data <- data[rownames(data) %in% gene_list, ]

# Step 5: Reorder/select columns by region
region_map <- data.frame(
  sample = colnames(data),
  region = sapply(colnames(data), function(x) strsplit(x, "-")[[1]][1])
)

keep_regions <- c("NAWM", "WM_line_Lesion_1", "WM_line_Lesion_2")
selected_samples <- region_map$sample[region_map$region %in% keep_regions]

# Final table
final_matrix <- filtered_data[, selected_samples]

# Step 6: Create a DGEList object for edgeR
# We assume that each sample in selected_samples corresponds to a unique column in the expression data
# You may need to define your group information if it's not included (in this case, it's not in the code)

# Define a dummy group factor for the samples (you can modify this to match your actual groupings)
group <- factor(rep("group1", ncol(final_matrix)))  # Adjust as per your experimental design

# Create DGEList object
y <- DGEList(counts = as.matrix(final_matrix), group = group)

# Step 7: Normalize library sizes using calcNormFactors
y <- calcNormFactors(y)

# After normalization, you can access the normalized counts using y$counts
normalized_counts <- y$counts

# Optional: Log-transform the normalized counts (if desired)
log_normalized_data <- log1p(normalized_counts)  # log(x + 1) transformation

# Optional: Save the normalized output
write.csv(log_normalized_data, "IL12_normalized_expression.csv")
#write.csv(log_normalized_data, "CellularResponseToHeatStress_normalized_expression.csv")
#write.csv(log_normalized_data, "ComplementCascade_normalized_expression.csv")
#write.csv(log_normalized_data, "S100Family_normalized_expression.csv")
#write.csv(log_normalized_data, "Neuroinflammation_normalized_expression.csv")
#write.csv(log_normalized_data, "MacrophageAlternative_normalized_expression.csv")



############### Heatmap ################


library(pheatmap)

# Load the filtered expression matrix
expr <- read.csv("IL12_normalized_expression.csv", row.names = 1)
#expr <- read.csv("CellularResponseToHeatStress_normalized_expression.csv", row.names = 1)
#expr <- read.csv("ComplementCascade_normalized_expression.csv", row.names = 1)
#expr <- read.csv("S100Family_normalized_expression.csv", row.names = 1)
#expr <- read.csv("Neuroinflammation_normalized_expression.csv", row.names = 1)
#expr <- read.csv("MacrophageAlternative_normalized_expression.csv", row.names = 1)

# Extract region names from column names (everything before the first dot)
region_for_expr <- sub("\\..*", "", colnames(expr))

# Create annotation dataframe
annotation_col <- data.frame(Region = region_for_expr)
rownames(annotation_col) <- colnames(expr)

# Define desired order of regions
region_order <- c("NAWM", "WM_line_Lesion_1", "WM_line_Lesion_2")

# Reorder columns by region
ordered_cols <- rownames(annotation_col)[order(factor(annotation_col$Region, levels = region_order))]
expr_ordered <- expr[, ordered_cols]
annotation_ordered <- annotation_col[ordered_cols, , drop = FALSE]

# Optional: Z-score scaling
scaled_expr <- t(scale(t(expr_ordered)))

# Plot heatmap
pheatmap(
  scaled_expr,
  cluster_rows = TRUE,
  cluster_cols = FALSE,  # disable column clustering to preserve order
  annotation_col = annotation_ordered,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  fontsize_row = 10,
  fontsize_col = 10,
  main = "Heatmap of IL-12 Signaling"
#  main = "Heatmap of Cellular Response To Heat Stress"
#  main = "Heatmap of Complement Cascade"
#  main = "Heatmap of S100 Family Signaling"
#  main = "Heatmap of Neuroinflammation Signaling"
#  main = "Heatmap of Macrophage Alternative Signaling"
)
