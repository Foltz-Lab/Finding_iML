#########################################################################################
####  Here, we generate Annotated Volcano plots of differentially accessible regions in the ML Subsets ATAC-Seq data. 
####  The DARs were called using DiffBind, and peaks were annotated using ChipSeeker. 
#########################################################################################

library(tidyverse)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(EnhancedVolcano)

# read in the csv files
print("Loading in the data...")
contrast7 <- read.csv("/storage1/fs1/jennifer.a.foltz/Active/Matt/JT048_JT051/MLsubsets_DiffBind/EnhancedVolcano/AnnotatedVolcano/20230829_contrast7_peakAnno_vol.csv")
contrast8 <- read.csv("/storage1/fs1/jennifer.a.foltz/Active/Matt/JT048_JT051/MLsubsets_DiffBind/EnhancedVolcano/AnnotatedVolcano/20230829_contrast8_peakAnno_vol.csv")
contrast9 <- read.csv("/storage1/fs1/jennifer.a.foltz/Active/Matt/JT048_JT051/MLsubsets_DiffBind/EnhancedVolcano/AnnotatedVolcano/20230829_contrast9_peakAnno_vol.csv")
contrast10 <- read.csv("/storage1/fs1/jennifer.a.foltz/Active/Matt/JT048_JT051/MLsubsets_DiffBind/EnhancedVolcano/AnnotatedVolcano/20230829_contrast10_peakAnno_vol.csv")
contrast11 <- read.csv("/storage1/fs1/jennifer.a.foltz/Active/Matt/JT048_JT051/MLsubsets_DiffBind/EnhancedVolcano/AnnotatedVolcano/20230829_contrast11_peakAnno_vol.csv")
contrast12 <- read.csv("/storage1/fs1/jennifer.a.foltz/Active/Matt/JT048_JT051/MLsubsets_DiffBind/EnhancedVolcano/AnnotatedVolcano/20230829_contrast12_peakAnno_vol.csv")

# ─── Helper: build a label vector for all DARs ───────────────────────────────
# Plots ALL peaks, but only labels the most-significant promoter peak per gene
# (and only for genes listed in select_genes).
make_promoter_labels <- function(dat, select_genes) {
  labels <- rep("", nrow(dat))                           # default: no label
  
  # Subset to promoter peaks belonging to genes of interest
  is_promoter <- grepl("Promoter", dat$annotation)
  is_selected <- !is.na(dat$SYMBOL) & dat$SYMBOL %in% select_genes
  promoter_sel <- dat[is_promoter & is_selected, ]
  
  # Keep only the single most-significant (lowest FDR) peak per gene
  promoter_sel <- promoter_sel[order(promoter_sel$FDR), ]
  promoter_sel <- promoter_sel[!duplicated(promoter_sel$SYMBOL), ]
  
  # Stamp those row-names with the gene symbol in the full label vector
  labels[rownames(dat) %in% rownames(promoter_sel)] <- promoter_sel$SYMBOL
  return(labels)
}

##Generate volcano plots##

keycontrast7 <- ifelse(contrast7$Fold > 0.5, 'cyan4', ifelse(contrast7$Fold < -0.5, '#008000', 'grey45'))
keycontrast7[is.na(keycontrast7)] <- 'grey45'
names(keycontrast7)[keycontrast7 == 'cyan4'] <- 'Up in eML'
names(keycontrast7)[keycontrast7 == 'grey45'] <- 'Not Significant'
names(keycontrast7)[keycontrast7 == '#008000'] <- 'Up in effcNK'

# Make a labels vector so each gene is only labeled once
best_peak_indices <- sapply(split(1:nrow(contrast7), contrast7$SYMBOL), function(idx) {
  idx[which.min(contrast7$FDR[idx])]
})


plot_labels <- as.character(contrast7$SYMBOL)

plot_labels[!(1:length(plot_labels) %in% best_peak_indices)] <- ""

# Generate Volcano Plot

EV_7 <- EnhancedVolcano(toptable = contrast7, lab = plot_labels, pointSize = 12, title = "eML vs. effcNK", titleLabSize = 100, axisLabSize = 100, xlim = c(-1.7, 1.7), ylim = c(0, 17), x = 'Fold', y = 'FDR', pCutoff = 0.05, FCcutoff = 0.5, cutoffLineType = 'dashed', cutoffLineWidth = 2.0, widthConnectors = 1.0, drawConnectors = TRUE, labSize = 30, labFace = "italic", boxedLabels = TRUE, colCustom = keycontrast7, legendLabSize = 50, legendIconSize = 25, colConnectors = 'black', colAlpha = 0.4, vlineWidth = 3, selectLab = c("KLF2", "LTB", "TBX21", "CD160", "TGFBR3L", "SPON2", "S1PR5", "CCL4L2", "FOXP1", "IFRD1", "ZFR", "SEMA4D", "IL12RB1", "CISH", "NR4A1", "ICAM2", "CRTAM", "TGIF1", "TGFBR2", "SEMA6", "IL12RB2", "CD96", "CD44", "KLRC1", "TGFBR3", "TIGIT", "IL18R1", "TCF7", "IL7R", "ENTPD1"), ylab = bquote(~-Log[10]~adjusted~italic(FDR)))

tiff("AnnotatedVolcano_contrast7.tiff", width = 36, height = 26, res = 100, units = "in")
print(EV_7)

dev.off()
print("Completed Contrast 7 Volcano Plot.")

select8 <- c("RAPGEF1", "KLF2", "LTB", "TBX21", "CD160", "TGFBR3L", "SPON2", "S1PR5", "CCL4L2", 
            "FOXP1", "IFRD1", "SEMA4D", "IL12RB1", "CISH", "NR4A1", "ICAM2", "CRTAM", "TGIF1", "TGFBR2", "SEMA6", 
            "IL12RB2", "CD96", "CD44", "KLRC1", "TGFBR3", "TIGIT", "IL18R1", "TCF7", "IL7R", "ENTPD1")

#keycontrast8 <- ifelse(contrast8$Fold > 0.5, 'cyan4', ifelse(contrast8$Fold < -0.5, '#551A8B', 'grey45'))
keycontrast8 <- case_when(
  contrast8$FDR < 0.05 & contrast8$Fold > 0.5  ~ 'cyan4',
  contrast8$FDR < 0.05 & contrast8$Fold < -0.5 ~ '#551A8B',
  TRUE                                         ~ 'grey45'
)
keycontrast8[is.na(keycontrast8)] <- 'grey45'
names(keycontrast8)[keycontrast8 == 'cyan4'] <- 'Up in eML'
names(keycontrast8)[keycontrast8 == 'grey45'] <- 'Not Significant'
names(keycontrast8)[keycontrast8 == '#551A8B'] <- 'Up in cNK'

# Generate Volcano Plot

EV_8 <- EnhancedVolcano(toptable = contrast8, lab = make_promoter_labels(contrast8, select8), pointSize = 12, title = "eML vs cNK", titleLabSize = 100, axisLabSize = 100, 
xlim = c(-3.0, 3.0), ylim = c(0, 18), x = 'Fold', y = 'FDR', pCutoff = 0.05, FCcutoff = 0.5, cutoffLineType = 'dashed', cutoffLineWidth = 2.0, 
widthConnectors = 1.0, drawConnectors = TRUE, labSize = 30, labFace = "italic", boxedLabels = TRUE, colCustom = keycontrast8, legendLabSize = 50, 
legendIconSize = 25, colConnectors = 'black', colAlpha = 0.4, vlineWidth = 3, maxoverlapsConnectors = Inf, ylab = bquote(~-Log[10]~adjusted~italic(FDR)))

tiff("AnnotatedVolcano_contrast8v2.tiff", width = 40, height = 30, res = 100, units = "in") 
print(EV_8)

dev.off()
print("Completed Contrast 8 Volcano Plot.")

# invert to flip iML and eML on plot
contrast9$Fold_inverted <- -contrast9$Fold

select9 <- c("IL7R", "SEMA6A", "ZFR", "FOXP1", "IFRD1", "SEMA4D")

#keycontrast9 <- ifelse(contrast9$Fold_inverted > 0.5, '#008000', ifelse(contrast9$Fold_inverted < -0.5, 'cyan4', 'grey45'))
keycontrast9 <- case_when(
  contrast9$FDR < 0.05 & contrast9$Fold_inverted > 0.5  ~ '#008000',
  contrast9$FDR < 0.05 & contrast9$Fold_inverted < -0.5 ~ 'cyan4',
  TRUE                                                  ~ 'grey45'
)
keycontrast9[is.na(keycontrast9)] <- 'grey45'
names(keycontrast9)[keycontrast9 == '#008000'] <- 'Up in iML'
names(keycontrast9)[keycontrast9 == 'grey45'] <- 'Not Significant'
names(keycontrast9)[keycontrast9 == 'cyan4'] <- 'Up in eML'

# Generate Volcano Plot

EV_9 <- EnhancedVolcano(toptable = contrast9, lab = make_promoter_labels(contrast9, select9), pointSize = 12, title = "iML vs eML", titleLabSize = 100, axisLabSize = 75, 
xlim = c(-3.0, 3.0), ylim = c(0, 18), x = 'Fold_inverted', y = 'FDR', pCutoff = 0.05, FCcutoff = 0.5, cutoffLineType = 'dashed', cutoffLineWidth = 2.0, 
widthConnectors = 1.0, drawConnectors = TRUE, boxedLabels = TRUE, labSize = 20, labFace = "italic", colCustom = keycontrast9, legendLabSize = 50, 
legendIconSize = 25, colConnectors = 'black', colAlpha = 0.4, vlineWidth = 3, maxoverlapsConnectors = Inf, ylab = bquote(~-Log[10]~adjusted~italic(FDR)))

tiff("AnnotatedVolcano_contrast9v3.tiff", width = 40, height = 30, res = 100, units = "in")
print(EV_9)

dev.off()
print("Completed Contrast 9 Volcano Plot.")

keycontrast10 <- ifelse(contrast10$Fold > 0.5, 'cyan4', ifelse(contrast10$Fold < -0.5, '#008000', 'grey45'))
keycontrast10[is.na(keycontrast10)] <- 'grey45'
names(keycontrast10)[keycontrast10 == 'cyan4'] <- 'Up in effcNK'
names(keycontrast10)[keycontrast10 == 'grey45'] <- 'Not Significant'
names(keycontrast10)[keycontrast10 == '#008000'] <- 'Up in cNK'

# Make a labels vector so each gene is only labeled once
best_peak_indices <- sapply(split(1:nrow(contrast10), contrast10$SYMBOL), function(idx) {
  idx[which.min(contrast10$FDR[idx])]
})

plot_labels <- as.character(contrast10$SYMBOL)

plot_labels[!(1:length(plot_labels) %in% best_peak_indices)] <- ""

# Generate Volcano Plot

EV_10 <- EnhancedVolcano(toptable = contrast10, lab = plot_labels, pointSize = 12, title = "effcNK vs cNK", titleLabSize = 100, axisLabSize = 100, xlim = c(-1.7, 1.7), ylim = c(0, 17), x = 'Fold', y = 'FDR', pCutoff = 0.05, FCcutoff = 0.5, cutoffLineType = 'dashed', cutoffLineWidth = 2.0, widthConnectors = 1.0, drawConnectors = TRUE, labSize = 20, labFace = "italic", colCustom = keycontrast10, legendLabSize = 50, legendIconSize = 25, colConnectors = 'black', colAlpha = 0.4, vlineWidth = 3, boxedLabels = TRUE, selectLab = c("KLF2", "LTB", "TBX21", "CD160", "TGFBR3L", "SPON2", "S1PR5", "CCL4L2", "FOXP1", "IFRD1", "ZFR", "SEMA4D", "IL12RB1", "CISH", "NR4A1", "ICAM2", "CRTAM", "TGIF1", "TGFBR2", "SEMA6", "IL12RB2", "CD96", "CD44", "KLRC1", "TGFBR3", "TIGIT", "IL18R1", "TCF7", "IL7R", "ENTPD1"), ylab = bquote(~-Log[10]~adjusted~italic(FDR)))

tiff("AnnotatedVolcano_contrast10.tiff", width = 36, height = 26, res = 100, units = "in")
print(EV_10)

dev.off()
print("Completed Contrast 10 Volcano Plot.")

keycontrast11 <- ifelse(contrast11$Fold > 0.5, 'cyan4', ifelse(contrast11$Fold < -0.5, '#008000', 'grey45'))
keycontrast11[is.na(keycontrast11)] <- 'grey45'
names(keycontrast11)[keycontrast11 == 'cyan4'] <- 'Up in effcNK'
names(keycontrast11)[keycontrast11 == 'grey45'] <- 'Not Significant'
names(keycontrast11)[keycontrast11 == '#008000'] <- 'Up in iML'

# Make a labels vector so each gene is only labeled once
best_peak_indices <- sapply(split(1:nrow(contrast11), contrast11$SYMBOL), function(idx) {
  idx[which.min(contrast11$FDR[idx])]
})

plot_labels <- as.character(contrast11$SYMBOL)

plot_labels[!(1:length(plot_labels) %in% best_peak_indices)] <- ""

# Generate Volcano Plot

EV_11 <- EnhancedVolcano(toptable = contrast11, lab = plot_labels, pointSize = 12, title = "effcNK vs iML", titleLabSize = 100, axisLabSize = 100, xlim = c(-1.7, 1.7), ylim = c(0, 17), x = 'Fold', y = 'FDR', pCutoff = 0.05, FCcutoff = 0.5, cutoffLineType = 'dashed', cutoffLineWidth = 2.0, widthConnectors = 1.0, drawConnectors = TRUE, labSize = 20, labFace = "italic", colCustom = keycontrast11, legendLabSize = 50, legendIconSize = 25, colConnectors = 'black', colAlpha = 0.4, vlineWidth = 3, boxedLabels = TRUE, selectLab = c("KLF2", "LTB", "TBX21", "CD160", "TGFBR3L", "SPON2", "S1PR5", "CCL4L2", "FOXP1", "IFRD1", "ZFR", "SEMA4D", "IL12RB1", "CISH", "NR4A1", "ICAM2", "CRTAM", "TGIF1", "TGFBR2", "SEMA6", "IL12RB2", "CD96", "CD44", "KLRC1", "TGFBR3", "TIGIT", "IL18R1", "TCF7", "IL7R", "ENTPD1"), ylab = bquote(~-Log[10]~adjusted~italic(FDR)))

tiff("AnnotatedVolcano_contrast11.tiff", width = 36, height = 26, res = 100, units = "in")
print(EV_11)

dev.off()
print("Completed Contrast 11 Volcano Plot.")

select12 <- c("KLF2", "CRTAM", "KLRC1", "S1PR5", "CCL4L2", "SPON2", "LTB",
              "TBX21", "ENTPD1", "IL7R", "ICAM2", "TCF7", "IL18R1",
              "IL12RB2", "TIGIT", "CD96", "CD44", "TGIF1", "TGFBR2", "TGFBR3",
              "CD160", "TGFBR3L", "CISH", "IL12RB1", "NR4A1", "RAPGEF1")

#keycontrast12 <- ifelse(contrast12$Fold > 0.5, '#008000', ifelse(contrast12$Fold < -0.5, '#551A8B', 'grey45'))
keycontrast12 <- case_when(
  contrast12$FDR < 0.05 & contrast12$Fold > 0.5  ~ '#008000',
  contrast12$FDR < 0.05 & contrast12$Fold < -0.5 ~ '#551A8B',
  TRUE                                           ~ 'grey45'
)
keycontrast12[is.na(keycontrast12)] <- 'grey45'
names(keycontrast12)[keycontrast12 == '#008000'] <- 'Up in iML'
names(keycontrast12)[keycontrast12 == 'grey45'] <- 'Not Significant'
names(keycontrast12)[keycontrast12 == '#551A8B'] <- 'Up in cNK'


# Generate Volcano Plot

EV_12 <- EnhancedVolcano(toptable = contrast12, lab = make_promoter_labels(contrast12, select12), pointSize = 12, title = " iML vs cNK", titleLabSize = 100, axisLabSize = 100, 
xlim = c(-3.0, 3.0), ylim = c(0, 18), x = 'Fold', y = 'FDR', pCutoff = 0.05, FCcutoff = 0.5, cutoffLineType = 'dashed', cutoffLineWidth = 2.0, 
widthConnectors = 1.0, drawConnectors = TRUE, labSize = 30, labFace = "italic", colCustom = keycontrast12, legendLabSize = 50, legendIconSize = 25, 
colConnectors = 'black', colAlpha = 0.4, vlineWidth = 3, boxedLabels = TRUE, maxoverlapsConnectors = Inf, ylab = bquote(~-Log[10]~adjusted~italic(FDR)))

tiff("AnnotatedVolcano_contrast12v2.tiff", width = 60, height = 45, res = 100, units = "in") 
print(EV_12)

dev.off()
print("Completed Contrast 12 Volcano Plot.")
