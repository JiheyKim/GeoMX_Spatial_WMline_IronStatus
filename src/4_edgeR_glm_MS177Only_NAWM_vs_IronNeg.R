
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

# Contrast for gNAWM vs. gWM_line_Lesion_2
contrast <- makeContrasts(WM_line_Lesion_2 - NAWM, levels = design)
qlf <- glmQLFTest(fit, contrast = contrast)

# Extract significant genes
significant_genes <- topTags(qlf, adjust.method = "fdr", n = Inf)$table %>%
  as.data.frame() %>%
#  filter(FDR < 0.05 & abs(logFC) >= log(2))
  filter(FDR < 0.05 )

# Get the gene names of significantly differentially expressed genes
significant_gene_names <- row.names(significant_genes)

# Extract the expression data for the significant genes
significant_expr_data <- m[significant_gene_names, ]

# Subset the data to include only samples from the two groups of interest
group_samples <- which(g %in% c("WM_line_Lesion_2", "NAWM"))

# Ensure that group_samples are ordered so WM_line_Lesion_2 comes first
group_samples <- group_samples[order(g[group_samples], decreasing = FALSE)]

# Subset the expression data based on the reordered group samples
expr_data_subset <- significant_expr_data[, group_samples]

# Ensure that column names of expr_data_subset match the sample names from y
sample_names <- rownames(y$samples)  # Ensure this matches your sample naming structure
colnames(expr_data_subset) <- sample_names[group_samples]

# Check for NA, NaN, or Inf values and handle them
expr_data_subset[is.na(expr_data_subset)] <- 0
expr_data_subset[is.nan(expr_data_subset)] <- 0
expr_data_subset[is.infinite(expr_data_subset)] <- 0

# Generate the heatmap
output_file <- "WMline_NAWM_vs_WMline_IronNeg_FDR005_wGeneName_CORRECTED.jpg"

# Create a JPG file to save the heatmap
jpeg(output_file, width = 800, height = 600)  # Adjust width and height as needed

# Generate the heatmap with split by groups
Heatmap(t(scale(t(log(expr_data_subset + 1))))[, ncol(expr_data_subset):1], 
        column_split = g[group_samples],  # Split by group
        show_row_names = TRUE, 
        show_column_names = TRUE,
        row_names_gp = gpar(fontsize = 5))  # Adjust fontsize value as needed

# Turn off the device
dev.off()

