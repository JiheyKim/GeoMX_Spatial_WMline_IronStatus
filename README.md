# Spatial Transcriptomics of Chronic-Active Lesions in Multiple Sclerosis

This project investigates the molecular mechanisms underlying chronic active white matter lesions in multiple sclerosis (MS), with a particular focus on the role of iron-positive and iron-negative microglia/macrophages. MS is a chronic demyelinating disease of the central nervous system, characterized by a transition from a relapsing-remitting phase to a progressive stage—where existing therapies are often ineffective.

## Overview

Using NanoString GeoMx Digital Spatial Profiling (DSP), we performed spatial transcriptomic analysis on postmortem MS brain tissue. Transcript profiles were obtained from:

- Iron-positive lesion rims  
- Iron-negative lesion rims  
- Lesion cores  
- Surrounding normal-appearing white matter (NAWM)

Our analysis uncovered region- and iron-status-specific transcriptional signatures linked to:

- Immune activation
- IL-12 signaling pathways
- Oxidative and cellular stress responses

This repository provides the data processing scripts, analysis code, and documentation to ensure reproducibility and transparency. The ultimate goal is to better understand the localized inflammatory environments driving lesion expansion and disease progression in MS.

---

## Repository Contents

- `data/` – Raw and processed expression count tables
- `scripts/` – All R scripts for preprocessing, normalization, QC, clustering, and visualization
- `figures/` – Output plots including UMAP, PCA, heatmaps, and hierarchical clustering
- `results/` – Differential expression results and summary tables
- `README.txt` – Project overview and setup instructions (this file)

---

## Requirements

### R Version

- **Raw data reading & QC**: R 4.4
- **Analysis & modeling**: R 4.3

### Required R Packages

#### For R 4.4 (QC and preprocessing)

```r
NanoStringNCTools
GeomxTools
GeoMxWorkflows
knitr
dplyr
ggforce
ggplot2
scales
ggfortify
reshape2
cowplot
pheatmap

#### For R 4.3 (statistical analysis and modeling)

```r
edgeR
limma
dplyr
magrittr
pheatmap
ComplexHeatmap
ggplot2
tibble
randomForest
