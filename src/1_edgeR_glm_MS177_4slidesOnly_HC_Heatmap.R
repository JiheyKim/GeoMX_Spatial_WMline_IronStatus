library(edgeR)
library(limma)
library(magrittr)

# Read the data
#f <- "RawCounts_LineSamplesOnly_RegionROI2byGene_PatientAdded.txt"
#f <- "RawCounts_MS177-6-5_4slides_ROIbyGene_AnnotAdded_sorted_SubRegionAdded.txt"
f <- "RawCounts_MS177-6-5_4slides_ROIbyGene_AnnotAdded_sorted_SubRegionAdded_4Heatmap.txt"
tt = read.table(f, header=T, skip=3)
colnames(tt) = gsub(" ", "_", colnames(tt))
row.names(tt) = tt$gene
lines <- readLines(f, n = 3)
g = strsplit(gsub(" ", "_", lines[1]), "\t")[[1]]; g = g[2:length(g)]
subg = strsplit(gsub(" ", "_", lines[2]), "\t")[[1]]; subg = subg[2:length(subg)]
pa = strsplit(gsub(" ", "_", lines[3]), "\t")[[1]]; pa = pa[2:length(pa)]
gene <- tt$gene
x = aggregate(. ~ gene, data = tt, FUN = mean)
row.names(x) = x$gene; x = x[, -1]

# Create DGEList object
g = as.factor(g)
y <- DGEList(counts = as.matrix(x), group = g)

# Filter lowly expressed genes
keep <- filterByExpr(y, group = g)
y <- y[keep, , keep.lib.sizes=FALSE]

# Normalize the library sizes
y <- calcNormFactors(y)

# Estimate dispersion
y <- estimateDisp(y, design = model.matrix(~0 + g))

# Design matrix based on group
design <- model.matrix(~0 + g)

# Fit the model using quasi-likelihood F-test
fit <- glmQLFit(y, design)

# Check the fit
head(fit$coefficients)

library(pheatmap)
library(ComplexHeatmap)

# Prepare normalized expression matrix
y = DGEList(counts = x, group = g)
keep = filterByExpr(y, group = g)
y = y[keep, , keep.lib.sizes=FALSE]
y = calcNormFactors(y)
m = cpm(y, normalized.lib.sizes = TRUE)

# Create annotation data frame
annotation_col <- data.frame(
  Slide = pa,
  SubRegion = subg,
  Region = g
)
row.names(annotation_col) <- colnames(m)

# Generate the heatmap
#### Bright Colors 
pheatmap(
  t(scale(t(log(m + 10)))),
  scale = "row", 
  show_rownames = FALSE,
  show_colnames = FALSE,
  border_color = NA,
  column_split = g,
  clustering_method = "average",
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  color = colorRampPalette(c("blue", "white", "red"))(100),  # Bright colors
  annotation_col = annotation_col,
  use_raster = FALSE
)

