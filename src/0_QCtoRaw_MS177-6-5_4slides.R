library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)


datadir <- file.path("../AllCombined_WM_GM_out")
DCCFiles <- dir(file.path(datadir), pattern = ".dcc$",
                full.names = TRUE, recursive = TRUE)
PKCFiles <- dir(file.path(datadir, "pkcs"), pattern = ".pkc$",
                                full.names = TRUE, recursive = TRUE)


PKCFiles <- dir(file.path(datadir), pattern = ".pkc$",
                                full.names = TRUE, recursive = TRUE)

SampleAnnotationFile <-
#    dir(file.path(datadir), pattern = "^Annotation_TEST1.xlsx$", full.names = TRUE, recursive = TRUE)
#    dir(file.path(datadir), pattern = "^Annotation_All_Combined_CII_v1_NAremoved_TC.xlsx$", full.names = TRUE, recursive = TRUE)
#    dir(file.path(datadir), pattern = "^Annotation_All_Combined_CII_v1_NAremoved_TC_MonoWaveRemoved.xlsx$", full.names = TRUE, recursive = TRUE)
#    dir(file.path(datadir), pattern = "^Annotation_MS177-6-5Only_5slides.xlsx$", full.names = TRUE, recursive = TRUE)
    dir(file.path(datadir), pattern = "^Annotation_MS177-6-5Only_4slides.xlsx$", full.names = TRUE, recursive = TRUE)
#   dir(file.path(datadir), pattern = "^Annotation_TEST1.xlsx$", full.names = TRUE, recursive = TRUE)

demoData <-
    readNanoStringGeoMxSet(dccFiles = DCCFiles,
                           pkcFiles = PKCFiles,
                           phenoDataFile = SampleAnnotationFile,
                           phenoDataSheet = "Template",
                           phenoDataDccColName = "Sample_ID",
                           protocolDataColNames = c("aoi", "roi"),
                           experimentDataColNames = c("panel"))


library(knitr)
pkcs <- annotation(demoData)
modules <- gsub(".pkc", "", pkcs)
kable(data.frame(PKCs = pkcs, modules = modules))

library(dplyr)
library(ggforce)

# select the annotations we want to show, use `` to surround column names with
# spaces or special symbols
count_mat <- count(pData(demoData), `slide name`, class, segment, Region, SubRegion, Group1, Group2, Group3, Group4, Patient)

test_gr <- gather_set_data(count_mat, 1:4)



test_gr$x <- factor(test_gr$x,
                    levels = c("class", "slide name", "segment", "Region", "SubRegion", "Group1", "Group2", "Group3", "Group4", "Patient"))


ggplot(test_gr, aes(x, id = id, split = y, value = n)) +
    geom_parallel_sets(aes(fill = Oligodendrocytes), alpha = 0.5, axis.width = 0.1) +
    geom_parallel_sets_axes(axis.width = 0.2) +
    geom_parallel_sets_labels(color = "white", size = 5) +
    theme_classic(base_size = 17) + 
    theme(legend.position = "bottom",
          axis.ticks.y = element_blank(),
          axis.line = element_blank(),
          axis.text.y = element_blank()) +
    scale_y_continuous(expand = expansion(0)) + 
    scale_x_discrete(expand = expansion(0)) +
    labs(x = "", y = "") +
    annotate(geom = "segment", x = 4.25, xend = 4.25,
             y = 20, yend = 120, lwd = 2) +
    annotate(geom = "text", x = 4.19, y = 70, angle = 90, size = 5,
             hjust = 0.5, label = "100 segments")
             
             
             
# Shift counts to one             
demoData <- shiftCountsOne(demoData, useDALogic = TRUE)




################################################################

########             Segment QC                     ###########

################################################################




# Default QC cutoffs are commented in () adjacent to the respective parameters
# study-specific values were selected after visualizing the QC results in more
# detail below
QC_params <-
    list(minSegmentReads = 1000, # Minimum number of reads (1000)
         percentTrimmed = 80,    # Minimum % of reads trimmed (80%)
         percentStitched = 80,   # Minimum % of reads stitched (80%)
         percentAligned = 75,    # Minimum % of reads aligned (80%)
         percentSaturation = 30, # Minimum sequencing saturation (50%)
         minNegativeCount = 1,   # Minimum negative control counts (10)
         maxNTCCount = 9000,     # Maximum counts observed in NTC well (1000)
         minNuclei = 20,         # Minimum # of nuclei estimated (100)
         minArea = 1000)         # Minimum segment area (5000)
demoData <-
    setSegmentQCFlags(demoData, 
                      qcCutoffs = QC_params)        

# Collate QC Results
QCResults <- protocolData(demoData)[["QCFlags"]]
flag_columns <- colnames(QCResults)
QC_Summary <- data.frame(Pass = colSums(!QCResults[, flag_columns]),
                         Warning = colSums(QCResults[, flag_columns]))
QCResults$QCStatus <- apply(QCResults, 1L, function(x) {
    ifelse(sum(x) == 0L, "PASS", "WARNING")
})
QC_Summary["TOTAL FLAGS", ] <-
    c(sum(QCResults[, "QCStatus"] == "PASS"),
      sum(QCResults[, "QCStatus"] == "WARNING"))


library(ggplot2)

col_by <- "segment"

# Graphical summaries of QC statistics plot function
QC_histogram <- function(assay_data = NULL,
                         annotation = NULL,
                         fill_by = NULL,
                         thr = NULL,
                         scale_trans = NULL) {
    plt <- ggplot(assay_data,
                  aes_string(x = paste0("unlist(`", annotation, "`)"),
                             fill = fill_by)) +
        geom_histogram(bins = 50) +
        geom_vline(xintercept = thr, lty = "dashed", color = "black") +
        theme_bw() + guides(fill = "none") +
        facet_wrap(as.formula(paste("~", fill_by)), nrow = 4) +
        labs(x = annotation, y = "Segments, #", title = annotation)
    if(!is.null(scale_trans)) {
        plt <- plt +
            scale_x_continuous(trans = scale_trans)
    }
    plt
}

QC_histogram(sData(demoData), "Trimmed (%)", col_by, 80)

QC_histogram(sData(demoData), "Stitched (%)", col_by, 80)

QC_histogram(sData(demoData), "Aligned (%)", col_by, 75)

QC_histogram(sData(demoData), "Saturated (%)", col_by, 50) +
    labs(title = "Sequencing Saturation (%)",
         x = "Sequencing Saturation (%)")

#QC_histogram(sData(demoData), "area", col_by, 1000, scale_trans = "log10")



# calculate the negative geometric means for each module
negativeGeoMeans <- 
    esBy(negativeControlSubset(demoData), 
         GROUP = "Module", 
         FUN = function(x) { 
             assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs") 
         }) 
protocolData(demoData)[["NegGeoMean"]] <- negativeGeoMeans

# explicitly copy the Negative geoMeans from sData to pData
negCols <- paste0("NegGeoMean_", modules)
pData(demoData)[, negCols] <- sData(demoData)[["NegGeoMean"]]
for(ann in negCols) {
    plt <- QC_histogram(pData(demoData), ann, col_by, 2, scale_trans = "log10")
   # print(plt)
}


### Finally we plot all of the QC Summary information in a table.

kable(QC_Summary, caption = "QC Summary Table for each Segment")

# Table: QC Summary Table for each Segment
# 
# |              | Pass| Warning|
# |:-------------|----:|-------:|
# |LowReads      |  680|       1|
# |LowTrimmed    |  681|       0|
# |LowStitched   |  675|       6|
# |LowAligned    |  675|       6|
# |LowSaturation |  679|       2|
# |LowNegatives  |  681|       0|
# |TOTAL FLAGS   |  672|       9|



####### CHECK WARNING ROIs by JIHYE #########

demoData@protocolData@data[demoData@protocolData@data$QCFlags$LowReads == "TRUE",]
#DSP-1001660019352-D-E02.dcc  # Type 1_WM line 
demoData@protocolData@data[demoData@protocolData@data$QCFlags$LowStitched == "TRUE",]
#DSP-1001660018524-F-D08.dcc # NAWM
#DSP-1001660018524-F-E05.dcc  # GML
#DSP-1001660018524-F-E07.dcc  # GML
#DSP-1012340086401-C-A01.dcc   # NAWM
#DSP-1012340093321-B-A11.dcc  # GM line_type 4 
#DSP-1012340093321-B-D07.dcc  # GM line_type 3
demoData@protocolData@data[demoData@protocolData@data$QCFlags$LowAligned == "TRUE",]
#DSP-1001660018524-F-D08.dcc # NAWM
#DSP-1001660018524-F-E07.dcc  # GML
#DSP-1001660019352-D-E11.dcc  # Type 1_WM line 
#DSP-1012340086401-C-A01.dcc # NAWM
#DSP-1012340093321-B-A11.dcc  # GM line_type 4 
#DSP-1012340093321-B-D07.dcc  # GM line_type 3
demoData@protocolData@data[demoData@protocolData@data$QCFlags$LowSaturation == "TRUE",]
#DSP-1001660019352-D-E02.dcc  # Type 1_WM line 
#DSP-1012340086401-C-A08.dcc # WM line



demoData <- demoData[, QCResults$QCStatus == "PASS"]

# Subsetting our dataset has removed samples which did not pass QC
dim(demoData)



################################################################

########              Probe QC                        ###########

################################################################



# Generally keep the qcCutoffs parameters unchanged. Set removeLocalOutliers to 
# FALSE if you do not want to remove local outliers
demoData <- setBioProbeQCFlags(demoData, 
                               qcCutoffs = list(minProbeRatio = 0.1,
                                                percentFailGrubbs = 20), 
                               removeLocalOutliers = TRUE)

ProbeQCResults <- fData(demoData)[["QCFlags"]]

# Define QC table for Probe QC
qc_df <- data.frame(Passed = sum(rowSums(ProbeQCResults[, -1]) == 0),
                    Global = sum(ProbeQCResults$GlobalGrubbsOutlier),
                    Local = sum(rowSums(ProbeQCResults[, -2:-1]) > 0
                                & !ProbeQCResults$GlobalGrubbsOutlier))


#Subset object to exclude all that did not pass Ratio & Global testing
ProbeQCPassed <- 
    subset(demoData, 
           fData(demoData)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
               fData(demoData)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)
dim(ProbeQCPassed)
demoData <- ProbeQCPassed 



################################################################

########        Create Gene-level Count Data         ###########

################################################################



# Check how many unique targets the object has
length(unique(featureData(demoData)[["TargetName"]]))

target_demoData <- aggregateCounts(demoData)
dim(target_demoData)


exprs(target_demoData)[1:5, 1:2]

################################################################

########        Limit of Quantification               ###########

################################################################


# Define LOQ SD threshold and minimum value
cutoff <- 2
minLOQ <- 2

# Calculate LOQ per module tested
LOQ <- data.frame(row.names = colnames(target_demoData))
for(module in modules) {
    vars <- paste0(c("NegGeoMean_", "NegGeoSD_"),
                   module)
    if(all(vars[1:2] %in% colnames(pData(target_demoData)))) {
        LOQ[, module] <-
            pmax(minLOQ,
                 pData(target_demoData)[, vars[1]] * 
                     pData(target_demoData)[, vars[2]] ^ cutoff)
    }
}
pData(target_demoData)$LOQ <- LOQ


################################################################

########               Filtering                      ###########

################################################################

LOQ_Mat <- c()
for(module in modules) {
    ind <- fData(target_demoData)$Module == module
    Mat_i <- t(esApply(target_demoData[ind, ], MARGIN = 1,
                       FUN = function(x) {
                           x > LOQ[, module]
                       }))
    LOQ_Mat <- rbind(LOQ_Mat, Mat_i)
}
# ensure ordering since this is stored outside of the geomxSet
LOQ_Mat <- LOQ_Mat[fData(target_demoData)$TargetName, ]



################################################################

########            Segment Gene Detection            ###########

################################################################

# Save detection rate information to pheno data
pData(target_demoData)$GenesDetected <- 
    colSums(LOQ_Mat, na.rm = TRUE)
pData(target_demoData)$GeneDetectionRate <-
    pData(target_demoData)$GenesDetected / nrow(target_demoData)

# Determine detection thresholds: 1%, 5%, 10%, 15%, >15%
pData(target_demoData)$DetectionThreshold <- 
    cut(pData(target_demoData)$GeneDetectionRate,
        breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
        labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%"))



# stacked bar plot of different cut points (1%, 5%, 10%, 15%)
ggplot(pData(target_demoData),
       aes(x = DetectionThreshold)) +
    geom_bar(aes(fill = Region)) +
    #geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
    geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +
    theme_bw() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(x = "Gene Detection Rate",
         y = "Segment, #",
         fill = "Region")
         
# cut percent genes detected at 1, 5, 10, 15
kable(table(pData(target_demoData)$DetectionThreshold,
            pData(target_demoData)$class))       

            
##########################################################################              
##############   
##############    !!!   Number of ROIs that we will use  !!!!    ######## 
##############       
##########################################################################

## Default              
#target_demoData <-
#    target_demoData[, pData(target_demoData)$GeneDetectionRate >= .1]

#dim(target_demoData[, pData(target_demoData)$GeneDetectionRate >= .1])


# We will use every genes
#dim(target_demoData[, pData(target_demoData)$GeneDetectionRate >= .01])

#target_demoData <-
#    target_demoData[, pData(target_demoData)$GeneDetectionRate >= .01]
#dim(target_demoData)

##
target_demoData <-
    target_demoData[, pData(target_demoData)$GeneDetectionRate >= 0.0]
dim(target_demoData)
         
         
################################################################

########            Gene Detection Rate               ###########

################################################################
         
         
         
library(scales) # for percent

# Calculate detection rate:
LOQ_Mat <- LOQ_Mat[, colnames(target_demoData)]
fData(target_demoData)$DetectedSegments <- rowSums(LOQ_Mat, na.rm = TRUE)
fData(target_demoData)$DetectionRate <-
    fData(target_demoData)$DetectedSegments / nrow(pData(target_demoData))

# Gene of interest detection table
goi <- c("PDCD1", "CD274", "IFNG", "CD8A", "CD68", "EPCAM",
         "KRT18", "NPHS1", "NPHS2", "CALB1", "CLDN8")
goi_df <- data.frame(
    Gene = goi,
    Number = fData(target_demoData)[goi, "DetectedSegments"],
    DetectionRate = percent(fData(target_demoData)[goi, "DetectionRate"]))         
         
         
         
################################################################

########            Gene filtering                    ###########

################################################################
         
   # Plot detection rate:
plot_detect <- data.frame(Freq = c(1, 5, 10, 20, 30, 50))
plot_detect$Number <-
    unlist(lapply(c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5),
                  function(x) {sum(fData(target_demoData)$DetectionRate >= x)}))
plot_detect$Rate <- plot_detect$Number / nrow(fData(target_demoData))
rownames(plot_detect) <- plot_detect$Freq

ggplot(plot_detect, aes(x = as.factor(Freq), y = Rate, fill = Rate)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = formatC(Number, format = "d", big.mark = ",")),
              vjust = 1.6, color = "black", size = 4) +
    scale_fill_gradient2(low = "orange2", mid = "lightblue",
                         high = "dodgerblue3", midpoint = 0.65,
                         limits = c(0,1),
                         labels = scales::percent) +
    theme_bw() +
    scale_y_continuous(labels = scales::percent, limits = c(0,1),
                       expand = expansion(mult = c(0, 0))) +
    labs(x = "% of Segments",
         y = "Genes Detected, % of Panel > LOQ")      
         

negativeProbefData <- subset(fData(target_demoData), CodeClass == "Negative")
neg_probes <- unique(negativeProbefData$TargetName)

##########################################################################              
##############   
##############    !!!   Number of Genes that we will use  !!!!    ######## 
##############       
##########################################################################

### Default

#target_demoData <- 
#    target_demoData[fData(target_demoData)$DetectionRate >= 0.1 |
#                        fData(target_demoData)$TargetName %in% neg_probes, ]
#dim(target_demoData)



dim(target_demoData[fData(target_demoData)$DetectionRate >= 0.01 |
                        fData(target_demoData)$TargetName %in% neg_probes, ])
                        
                        
### We will use > 0.01 more genes
target_demoData <- 
    target_demoData[fData(target_demoData)$DetectionRate >= 0.01 |
                        fData(target_demoData)$TargetName %in% neg_probes, ]
dim(target_demoData)

write.table(exprs(target_demoData), file="RawCounts_MS177-6-5Only_4slidesOnly_ROIbyGene_CORRECTED.txt", quote=F, sep="\t")


