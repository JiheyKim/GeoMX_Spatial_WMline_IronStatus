library(edgeR)
library(dplyr)
library(ComplexHeatmap)

# Prepare the DGEList and filter
y = DGEList(counts=x, group=g)
keep = filterByExpr(y, group=g)
y = normLibSizes(y)
m = cpm(y, normalized.lib.sizes=TRUE)

# Contrast for gWM_line_Lesion_1 vs. gWM_line_Lesion_2
contrast <- makeContrasts(gWM_line_Lesion_1 - gNAWM, levels = design)
lrt <- glmLRT(fit, contrast = contrast)


# Extract all genes (no filtering on FDR or logFC)
all_genes <- topTags(lrt, adjust.method = "fdr", n = Inf)$table %>%
  as.data.frame()

# Save all genes
write.table(all_genes, file="WMline_IronPos_vs_NAWM_AllGenes.txt", quote=F, sep="\t")



# Extract significant genes
significant_genes <- topTags(lrt, adjust.method = "fdr", n = Inf)$table %>%
  as.data.frame() %>%
#  filter(FDR < 0.05 & abs(logFC) >= log(2))
  filter(FDR < 0.05 )

# Get the gene names of significantly differentially expressed genes
significant_gene_names <- row.names(significant_genes)

# Extract the expression data for the significant genes
significant_expr_data <- m[significant_gene_names, ]

# Subset the data to include only samples from the two groups of interest
group_samples <- which(g %in% c("WM_line_Lesion_1", "NAWM"))
expr_data_subset <- significant_expr_data[, group_samples]

# Ensure that column names of expr_data_subset match the sample names from y
sample_names <- rownames(y$samples)  # Ensure this matches your sample naming structure
colnames(expr_data_subset) <- sample_names[group_samples]

# Check for NA, NaN, or Inf values and handle them
expr_data_subset[is.na(expr_data_subset)] <- 0
expr_data_subset[is.nan(expr_data_subset)] <- 0
expr_data_subset[is.infinite(expr_data_subset)] <- 0

# Generate the heatmap
output_file <- "WMline_IronPos_vs_NAWM_FDR005_wGeneName.jpg"

# Create a JPG file to save the heatmap
jpeg(output_file, width = 800, height = 600)  # Adjust width and height as needed

# Generate the heatmap with split by groups
Heatmap(t(scale(t(log(expr_data_subset + 10)))), 
        column_split = g[group_samples],  # Split by group
        show_row_names = TRUE, 
        show_column_names = TRUE,
        row_names_gp = gpar(fontsize = 5))  # Adjust fontsize value as needed

# Turn off the device
dev.off()
