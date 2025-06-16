# Load libraries
library(dplyr)
library(readr)

# === Load Ferroptosis gene lists ===
ferr_act <- readLines("FerrAct_List.txt")
ferr_inhb <- readLines("FerrInhb_List.txt")
ferr_deg <- readLines("Sig_Ferr_GeneNames_NAWM_WMlines.txt")
ferr_tony <- readLines("TonyFerrList.txt")
ferr_gpt <- readLines("FerrGPT_List.txt")
# === Load expression matrix ===
# Skip the first 3 metadata rows, and read gene counts
expr_raw <- read.delim("../RawCounts_MS177-6-5Only_4slidesOnly_ROIbyGene_CORRECTED_AnnotAdded_4heatmap.txt", skip = 3, check.names = FALSE)
rownames(expr_raw) <- expr_raw$gene
expr_raw <- expr_raw[, -1]  # Remove gene name column

# === Extract rows corresponding to ferroptosis genes ===
ferr_genes <- unique(c(ferr_act, ferr_inhb))
expr_ferr <- expr_raw[rownames(expr_raw) %in% ferr_genes, ]

# === (Optional) Split into activators and inhibitors ===
expr_ferr_act <- expr_raw[rownames(expr_raw) %in% ferr_act, ]
expr_ferr_inhb <- expr_raw[rownames(expr_raw) %in% ferr_inhb, ]
expr_ferr_deg <- expr_raw[rownames(expr_raw) %in% ferr_deg, ]
expr_ferr_tony <- expr_raw[rownames(expr_raw) %in% ferr_tony, ]
expr_ferr_gpt <- expr_raw[rownames(expr_raw) %in% ferr_gpt, ]

# === Save outputs ===
write.table(expr_ferr, "Ferroptosis_All.txt", sep = "\t", quote = FALSE)
write.table(expr_ferr_act, "Ferroptosis_Activators.txt", sep = "\t", quote = FALSE)
write.table(expr_ferr_inhb, "Ferroptosis_Inhibitors.txt", sep = "\t", quote = FALSE)
write.table(expr_ferr_deg, "Ferroptosis_DEGs.txt", sep = "\t", quote = FALSE)
write.table(expr_ferr_tony, "Ferroptosis_Tony.txt", sep = "\t", quote = FALSE)
write.table(expr_ferr_gpt, "Ferroptosis_GPT.txt", sep = "\t", quote = FALSE)


########## Heatmap 

# Load packages
library(pheatmap)
library(RColorBrewer)
library(readr)

# === Load data ===
#ferr_counts <- read.delim("Ferroptosis_All.txt", row.names = 1, check.names = FALSE)
#ferr_counts <- read.delim("Ferroptosis_Activators.txt", row.names = 1, check.names = FALSE)
#ferr_counts <- read.delim("Ferroptosis_Inhibitors.txt", row.names = 1, check.names = FALSE)
ferr_counts <- read.delim("Ferroptosis_DEGs.txt", row.names = 1, check.names = FALSE)
ferr_counts <- read.delim("Ferroptosis_Tony.txt", row.names = 1, check.names = FALSE)
ferr_counts <- read.delim("Ferroptosis_GPT.txt", row.names = 1, check.names = FALSE)

# === Log2 transform with pseudo-count ===
log_expr <- log2(ferr_counts + 1)

# === Optional: row scaling for heatmap (z-score per gene) ===
log_expr_scaled <- t(scale(t(log_expr)))  # rows = genes

# === Extract Region information from metadata ===
meta_lines <- readLines("../RawCounts_MS177-6-5Only_4slidesOnly_ROIbyGene_CORRECTED_AnnotAdded_4heatmap.txt", n = 3)
region_labels <- strsplit(meta_lines[1], "\t")[[1]][-1]  # remove "Region" label

# Check alignment
stopifnot(length(region_labels) == ncol(expr_raw))  # Ensure the region matches number of samples

# Create annotation data frame
sample_anno <- data.frame(Region = region_labels)
rownames(sample_anno) <- colnames(expr_raw)  # match by column names

# Define colors (adjust as needed)
anno_colors <- list(
  Region = c(
    NAWM = "#A6CEE3",
    WM_line_Lesion = "#1F78B4",
    WML_Lesion = "#B2DF8A"
  )
)


# Draw heatmap
pheatmap(log_expr_scaled,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         annotation_col = sample_anno["Region"],
         annotation_colors = anno_colors["Region"],
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         main = "Ferroptosis Gene Expression (log2 normalized)",
         show_rownames = TRUE,
         show_colnames = FALSE)



########### ADD SubRegiosns ########

meta_lines <- readLines("../RawCounts_MS177-6-5Only_4slidesOnly_ROIbyGene_CORRECTED_AnnotAdded_4heatmap.txt", n = 3)

region_labels <- strsplit(meta_lines[1], "\t")[[1]][-1]      # 1st line: Region
subregion_labels <- strsplit(meta_lines[2], "\t")[[1]][-1]   # 2nd line: Subregion (edit as needed)

sample_anno <- data.frame(
  Subregion = subregion_labels,
  Region = region_labels
)
rownames(sample_anno) <- colnames(expr_raw)

anno_colors <- list(
  Subregion = c(
    NAWM = "#FDBF6F",
    WM_line_Lesion_IronPos = "#FF7F00",
    WM_line_Lesion_IronNeg = "#CAB2D6",
    WML_Lesion_IronPos = "#6A3D9A",
    WML_Lesion_IronNeg = "#B15928"
  ),
    Region = c(
    NAWM = "#A6CEE3",
    WM_line_Lesion = "#1F78B4",
    WML_Lesion = "#B2DF8A"
  )
)

pheatmap(log_expr_scaled,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         annotation_col = sample_anno,
         annotation_colors = anno_colors,
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         main = "Ferroptosis Gene Expression (log2 normalized)",
         show_rownames = TRUE,
         show_colnames = FALSE)



########## Only NAWM and WM lines ###########

# Subset to desired subregions
keep_subregions <- c("NAWM", "WM_line_Lesion_IronPos", "WM_line_Lesion_IronNeg")
subset_samples <- rownames(sample_anno)[sample_anno$Subregion %in% keep_subregions]

# Subset expression matrix and annotation
log_expr_subset <- log_expr_scaled[, subset_samples]
sample_anno_subset <- sample_anno[subset_samples, , drop = FALSE]

# Subset annotation colors
anno_colors_subset <- list(
  Region = anno_colors$Region[unique(sample_anno_subset$Region)],
  Subregion = anno_colors$Subregion[unique(sample_anno_subset$Subregion)]
)

pheatmap(log_expr_subset,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         annotation_col = sample_anno_subset,
         annotation_colors = anno_colors_subset,
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         main = "Ferroptosis Gene Expression (Selected Subregions)",
         show_rownames = TRUE,
         show_colnames = FALSE)


########## Only  WM lines ###########



# Subset to desired subregions
keep_subregions2 <- c("WM_line_Lesion_IronPos", "WM_line_Lesion_IronNeg")
subset_samples <- rownames(sample_anno)[sample_anno$Subregion %in% keep_subregions2]

# Subset expression matrix and annotation
log_expr_subset <- log_expr_scaled[, subset_samples]
sample_anno_subset <- sample_anno[subset_samples, , drop = FALSE]

# Subset annotation colors
anno_colors_subset <- list(
  Region = anno_colors$Region[unique(sample_anno_subset$Region)],
  Subregion = anno_colors$Subregion[unique(sample_anno_subset$Subregion)]
)

pheatmap(log_expr_subset,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         annotation_col = sample_anno_subset,
         annotation_colors = anno_colors_subset,
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         main = "Ferroptosis Gene Expression (Selected Subregions)",
         show_rownames = TRUE,
         show_colnames = FALSE)





