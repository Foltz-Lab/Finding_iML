# Propeller tool- to test statistically significant differences 
# between iML and cNK cell types across difference source types
# normal vs tumor, PBMC vs other sources

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("speckle")
# Load necessary library

library(speckle)
library(limma)
library(dplyr)

###################################################################################################################################
# Load the data


pan_cancer <- read.csv("/Volumes/jennifer.a.foltz/Active/Veda/pan_cancer/adata_combined_malm_ding_Tcell_filtered_fan_included/output/Malm_Ding_combined_selected_sources_Marsilea_ciml_included.csv")
unique(pan_cancer$source)
pan_cancer <- pan_cancer %>%
  filter(classification != "unclassified") %>%  # Remove 'unclassified'
  mutate(classification = case_when(
    classification %in% c("CD56bright", "CD56dim") ~ "cNK",
    classification %in% c("ML1", "ML2", "ML_transition") ~ "iML",
    TRUE ~ classification
  ))

# Plot cell type proportions
plotCellTypeProps(clusters = pan_cancer$classification, sample = pan_cancer$patient)

####################################################################################################################################################################

##### Lung #####

# Subset for lung
lung <- subset(pan_cancer, source %in% c("lung_normal", "lung_tumor"))
unique(lung$source)

# Perform propeller testing
propeller_results_lung <- propeller(
  clusters = lung$classification,
  sample = lung$patient,
  group = lung$source,transform="asin"
)

# Print the results
print(propeller_results_lung)

# Plot cell type proportions
plotCellTypeProps(clusters = lung$classification, sample = lung$patient)

#####Breast####

breast <- subset(pan_cancer, source %in% c("breast_normal", "breast_tumor"))
unique(breast$source)
table(breast$classification)

# Perform propeller testing
propeller_results_breast <- propeller(
  clusters = breast$classification,
  sample = breast$patient,
  group = breast$source, transform="asin"
)


# Print the results
print(propeller_results_breast)

# Plot cell type proportions
plotCellTypeProps(clusters = breast$classification, sample = breast$patient)

####Pancreas####

# Subset for pancreas
pancreas <- subset(pan_cancer, source %in% c("pancreas_normal", "pancreas_tumor"))
unique(pancreas$source)

# Perform propeller testing
propeller_results_pancreas <- propeller(
  clusters = pancreas$classification,
  sample = pancreas$patient,
  group = pancreas$source, transform="asin"
)

# Print the results
print(propeller_results_pancreas)

# Plot cell type proportions
plotCellTypeProps(clusters = pancreas$classification, sample = pancreas$patient)

####Prostate####

# Subset for prostate
prostate <- subset(pan_cancer, source %in% c("prostate_normal", "prostate_tumor"))
unique(prostate$source)

# Perform propeller testing
propeller_results_prostate <- propeller(
  clusters = prostate$classification,
  sample = prostate$patient,
  group = prostate$source, transform="asin"
)

# Print the results
print(propeller_results_prostate)

# Plot cell type proportions
plotCellTypeProps(clusters = prostate$classification, sample = prostate$patient)

#### Skin #####

# Subset for skin
skin <- subset(pan_cancer, source %in% c("skin_normal", "melanoma"))
unique(skin$source)

# Perform propeller testing
propeller_results_skin <- propeller(
  clusters = skin$classification,
  sample = skin$patient,
  group = skin$source, transform="asin"
)

# Print the results
print(propeller_results_skin)

# Plot cell type proportions
plotCellTypeProps(clusters = skin$classification, sample = skin$patient)

unique(pan_cancer$source)

####################################################################################################################################################################

#PLOT PROPORTION OF NK SUBSETS ACROSS SAMPLES
# Extract the p-values from the data frame
p_values <- propeller_results_lung$P.Value

cell_types <- rownames(propeller_results_lung)

bar_colors <- rainbow(length(cell_types))  

barplot(-log10(p_values), 
        names.arg = cell_types, 
        las = 1,                # las=1 for horizontal axis labels
        col = bar_colors,       # Use the color vector for different bar colors
        main = "P-values for NK Cells subsets across samples- lung source type", 
        xlab = "Cell Types", 
        ylab = "-log10(P-value)", 
        cex.names = 0.6
)

####################################################################################################################################################################

# ASSOCIATION of NK cell subsets across nk_type

# Define significance thresholds
p_threshold <- 0.05
fdr_threshold <- 0.05

# Check significance based on p-value and FDR
propeller_results_lung$p_significance <- ifelse(propeller_results_lung$P.Value < p_threshold, "Significant", "Not Significant")
propeller_results_lung$fdr_significance <- ifelse(propeller_results_lung$FDR < fdr_threshold, "Significant", "Not Significant")

# Display results
propeller_results_lung

####################################################################################################################################################################

####PBMC and lung_normal#####

PBMC_lung_normal <- subset(pan_cancer, source %in% c("PBMC", "lung_normal"))
unique(PBMC_lung_normal$source)

# Perform propeller testing
propeller_results_PBMC_lung_normal <- propeller(
  clusters = PBMC_lung_normal$classification,
  sample = PBMC_lung_normal$patient,
  group = PBMC_lung_normal$source, transform="asin"
)

# Print the results
print(propeller_results_PBMC_lung_normal)


#####PBMC and lung tumor########
PBMC_lung_tumor<- subset(pan_cancer, source %in% c("PBMC", "lung_tumor"))
unique(PBMC_lung_tumor$source)

# Perform propeller testing
propeller_results_PBMC_lung_tumor <- propeller(
  clusters = PBMC_lung_tumor$classification,
  sample = PBMC_lung_tumor$patient,
  group = PBMC_lung_tumor$source, transform="asin"
)

# Print the results
print(propeller_results_PBMC_lung_tumor)

#####PBMC and breast tumor########
PBMC_breast_tumor <- subset(pan_cancer, source %in% c("PBMC", "breast_tumor"))
unique(PBMC_breast_tumor$source)

# Perform propeller testing
propeller_results_PBMC_breast_tumor <- propeller(
  clusters = PBMC_breast_tumor$classification,
  sample = PBMC_breast_tumor$patient,
  group = PBMC_breast_tumor$source, transform="asin"
)

# Print the results
print(propeller_results_PBMC_breast_tumor)



#####PBMC and breast_normal########
PBMC_breast_normal <- subset(pan_cancer, source %in% c("PBMC", "breast_normal"))
unique(PBMC_breast_normal$source)

# Perform propeller testing
propeller_results_PBMC_breast_normal <- propeller(
  clusters = PBMC_breast_normal$classification,
  sample = PBMC_breast_normal$patient,
  group = PBMC_breast_normal$source, transform="asin"
)

# Print the results
print(propeller_results_PBMC_breast_normal)


#####PBMC and GBM########

PBMC_GBM <- subset(pan_cancer, source %in% c("PBMC", "glioblastoma"))
unique(PBMC_GBM$source)

# Perform propeller testing
propeller_results_PBMC_GBM <- propeller(
  clusters = PBMC_GBM$classification,
  sample = PBMC_GBM$patient,
  group = PBMC_GBM$source, transform="asin"
)

# Print the results
print(propeller_results_PBMC_GBM)

#####PBMC and Sarcoma########
PBMC_S <- subset(pan_cancer, source %in% c("PBMC", "sarcoma"))
unique(PBMC_S$source)

# Perform propeller testing
propeller_results_PBMC_S <- propeller(
  clusters = PBMC_S$classification,
  sample = PBMC_S$patient,
  group = PBMC_S$source, transform="asin"
)

# Print the results
print(propeller_results_PBMC_S)

#####PBMC and cervical tumor########
PBMC_cervical <- subset(pan_cancer, source %in% c("PBMC", "cervical_tumor"))
unique(PBMC_cervical$source)

# Perform propeller testing
propeller_results_PBMC_cervical <- propeller(
  clusters = PBMC_cervical$classification,
  sample = PBMC_cervical$patient,
  group = PBMC_cervical$source, transform="asin"
)

# Print the results
print(propeller_results_PBMC_cervical)

#####PBMC and HNSCC ########
PBMC_HNSCC <- subset(pan_cancer, source %in% c("PBMC", "HNSCC"))
unique(PBMC_HNSCC$source)

# Perform propeller testing
propeller_results_PBMC_HNSCC <- propeller(
  clusters = PBMC_HNSCC$classification,
  sample = PBMC_HNSCC$patient,
  group = PBMC_HNSCC$source, transform="asin"
)

# Print the results
print(propeller_results_PBMC_HNSCC)

#####PBMC and multiple_myeloma ########
PBMC_multiple_myeloma <- subset(pan_cancer, source %in% c("PBMC", "multiple_myeloma"))
unique(PBMC_multiple_myeloma$source)

# Perform propeller testing
propeller_results_PBMC_multiple_myeloma <- propeller(
  clusters = PBMC_multiple_myeloma$classification,
  sample = PBMC_multiple_myeloma$patient,
  group = PBMC_multiple_myeloma$source, transform="asin"
)

# Print the results
print(propeller_results_PBMC_multiple_myeloma)

#####PBMC and renal_tumor ########
PBMC_renal_tumor <- subset(pan_cancer, source %in% c("PBMC", "renal_tumor"))
unique(PBMC_renal_tumor$source)

# Perform propeller testing
propeller_results_PBMC_renal_tumor <- propeller(
  clusters = PBMC_renal_tumor$classification,
  sample = PBMC_renal_tumor$patient,
  group = PBMC_renal_tumor$source, transform="asin"
)

# Print the results
print(propeller_results_PBMC_renal_tumor)

#####PBMC and colorectal_tumor ########
PBMC_colorectal_tumor <- subset(pan_cancer, source %in% c("PBMC", "colorectal_tumor"))
unique(PBMC_colorectal_tumor$source)

# Perform propeller testing
propeller_results_PBMC_colorectal_tumor <- propeller(
  clusters = PBMC_colorectal_tumor$classification,
  sample = PBMC_colorectal_tumor$patient,
  group = PBMC_colorectal_tumor$source, transform="asin"
)

# Print the results
print(propeller_results_PBMC_colorectal_tumor)

