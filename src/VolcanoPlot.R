# Install the ggrepel package if you don't have it already
# install.packages("ggrepel")

library(ggplot2)
library(ggrepel)  # Load the ggrepel package

################################################################
####
####           Iron Pos vs Iron Neg       
####
################################################################


# Read the data
rd <- read.table("WML_IronPos_vs_IronNeg_4Volcano.txt", sep='\t', header=TRUE)
#rd <- read.table("WMline_IronPos_vs_IronNeg_4Volcano.txt", sep='\t', header=TRUE)


de <- rd
de$diffexpressed <- "NO change"
de$diffexpressed[de$logFC > 0 & de$FDR < 0.05] <- "Iron+"
de$diffexpressed[de$logFC < -0 & de$FDR < 0.05] <- "Iron-"

# Define custom colors
mycolors <- c("Iron+" = "blue", "Iron-" = "red", "NO change" = "black")

# Create the base plot
p <- ggplot(data=de, aes(x=logFC, y=-log10(FDR), col=diffexpressed)) + 
  geom_point() + 
  theme_minimal() + 
  scale_colour_manual(values = mycolors) +
  scale_x_reverse() +  # <<== Mirror image here
  xlab("logFC (mirrored: Iron+ on right, Iron- on left)")

# Highlight top genes
top_genes <- de[order(de$FDR),][1:20,]

# Add gene labels
p3 <- p + 
  geom_text_repel(data=top_genes, aes(x=logFC, y=-log10(FDR), label=gene), 
                  size=3, color="black", box.padding=0.35, point.padding=0.5, 
                  segment.color="black")

# Display the mirrored plot
p3


############### top 10 genes each side #########

# Select the top 10 most significant genes for Iron+ and Iron-
top_iron_pos <- de[de$diffexpressed == "Iron+", ][order(de$FDR), ][1:10, ]
top_iron_neg <- de[de$diffexpressed == "Iron-", ][order(de$FDR), ][1:10, ]

# Combine the selected genes
top_genes <- rbind(top_iron_pos, top_iron_neg)

# Add gene names to the plot using geom_text_repel to avoid overlaps
p3 <- p + 
  geom_text_repel(data=top_genes, aes(x=logFC, y=-log10(FDR), label=gene), 
                  size=3, color="black", box.padding=0.35, point.padding=0.5, 
                  segment.color="black")

# Display the plot
p3



################################################################
####
####           Iron Pos vs NAWM      
####
################################################################

library(ggplot2)
library(ggrepel)  # Load the ggrepel package

# Read the data

rd <- read.table("WMline_IronPos_vs_NAWM_4Volcano.txt", sep='\t', header=TRUE)
#rd <- read.table("WML_IronPos_vs_NAWM_4Volcano.txt", sep='\t', header=TRUE)

de <- rd
de$diffexpressed <- "NO change"
de$diffexpressed[de$logFC > 0 & de$FDR < 0.05] <- "NAWM"
de$diffexpressed[de$logFC < -0 & de$FDR < 0.05] <- "Iron+"

# Define custom colors
mycolors <- c("NAWM" = "blue", "Iron+" = "red", "NO change" = "black")

# Create the base plot
p <- ggplot(data=de, aes(x=logFC, y=-log10(FDR), col=diffexpressed)) + 
  geom_point() + 
  theme_minimal() + 
  scale_colour_manual(values = mycolors) +
  scale_x_reverse() +  # <<== Mirror image here
  xlab("logFC (mirrored: Iron+ on right, NAWM on left)")

# Highlight top genes
top_genes <- de[order(de$FDR),][1:20,]

# Add gene labels
p3 <- p + 
  geom_text_repel(data=top_genes, aes(x=logFC, y=-log10(FDR), label=gene), 
                  size=3, color="black", box.padding=0.35, point.padding=0.5, 
                  segment.color="black")

# Display the mirrored plot
p3


############### top 10 genes each side #########

# Select the top 10 most significant genes for Iron+ and Iron-
top_iron_pos <- de[de$diffexpressed == "Iron+", ][order(de$FDR), ][1:10, ]
top_nawm <- de[de$diffexpressed == "NAWM", ][order(de$FDR), ][1:10, ]

# Combine the selected genes
top_genes <- rbind(top_iron_pos, top_nawm)

# Add gene names to the plot using geom_text_repel to avoid overlaps
p3 <- p + 
  geom_text_repel(data=top_genes, aes(x=logFC, y=-log10(FDR), label=gene), 
                  size=3, color="black", box.padding=0.35, point.padding=0.5, 
                  segment.color="black")

# Display the plot
p3


################################################################
####
####           Iron Neg vs NAWM      
####
################################################################

library(ggplot2)
library(ggrepel)  # Load the ggrepel package

# Read the data

rd <- read.table("WMline_IronNeg_vs_NAWM_4Volcano.txt", sep='\t', header=TRUE)
#rd <- read.table("WML_IronNeg_vs_NAWM_4Volcano.txt", sep='\t', header=TRUE)

de <- rd
de$diffexpressed <- "NO change"
de$diffexpressed[de$logFC > 0 & de$FDR < 0.05] <- "NAWM"
de$diffexpressed[de$logFC < -0 & de$FDR < 0.05] <- "Iron-"

# Define custom colors
mycolors <- c("NAWM" = "blue", "Iron-" = "red", "NO change" = "black")

# Create the base plot
p <- ggplot(data=de, aes(x=logFC, y=-log10(FDR), col=diffexpressed)) + 
  geom_point() + 
  theme_minimal() + 
  scale_colour_manual(values = mycolors) +
  scale_x_reverse() +  # <<== Mirror image here
  xlab("logFC (mirrored: Iron- on right, NAWM on left)")

# Highlight top genes
top_genes <- de[order(de$FDR),][1:20,]

# Add gene labels
p3 <- p + 
  geom_text_repel(data=top_genes, aes(x=logFC, y=-log10(FDR), label=gene), 
                  size=3, color="black", box.padding=0.35, point.padding=0.5, 
                  segment.color="black")

# Display the mirrored plot
p3


############### top 10 genes each side #########

# Select the top 10 most significant genes for Iron+ and Iron-
top_iron_neg <- de[de$diffexpressed == "Iron-", ][order(de$FDR), ][1:10, ]
top_nawm <- de[de$diffexpressed == "NAWM", ][order(de$FDR), ][1:10, ]

# Combine the selected genes
top_genes <- rbind(top_iron_neg, top_nawm)

# Add gene names to the plot using geom_text_repel to avoid overlaps
p3 <- p + 
  geom_text_repel(data=top_genes, aes(x=logFC, y=-log10(FDR), label=gene), 
                  size=3, color="black", box.padding=0.35, point.padding=0.5, 
                  segment.color="black")

# Display the plot
p3

