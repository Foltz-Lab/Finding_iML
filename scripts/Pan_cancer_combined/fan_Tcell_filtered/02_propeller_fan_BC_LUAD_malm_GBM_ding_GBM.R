# Propeller tool- to test statistically significant differences 
# between iML and cNK cell types across difference source types
# GBM vs BrM_BC vs BrM_LC(LUAD subsetted)

library(speckle)
library(limma)
library(dplyr)

csv<-read.csv("/Volumes/jennifer.a.foltz/Active/Veda/pan_cancer/adata_combined_malm_ding_Tcell_filtered_fan_included/output/fan_BC_LUAD_malm_GBM_ding_GBM.csv")
csv %>%
  group_by(source) %>%
  summarise(unique_patients = n_distinct(patient))

##### GBM vs BrM_BC #####

# Subset for lung
GBM_BC <- subset(csv, source %in% c("glioblastoma", "BC"))
unique(GBM_BC$source)

# Recode classifications
GBM_BC <- GBM_BC %>%
  filter(classification != "unclassified") %>%  # Remove 'unclassified'
  mutate(classification = case_when(
    classification %in% c("CD56bright", "CD56dim") ~ "cNK",
    classification %in% c("ML1", "ML2", "ML_transition") ~ "eML",
    TRUE ~ classification
  ))

# Perform propeller testing
propeller_results_GBM_BC<- propeller(
  clusters = GBM_BC$classification,
  sample = GBM_BC$patient,
  group = GBM_BC$source,transform="asin"
)

# Print the results
print(propeller_results_GBM_BC)
unique(GBM_BC$source)
# Plot cell type proportions
plotCellTypeProps(clusters = GBM_BC$classification, sample = GBM_BC$patient)

##### GBM vs BrM_LC #####

# Subset for lung
GBM_LC <- subset(csv, source %in% c("glioblastoma", "LC"))
unique(GBM_LC$source)

# Recode classifications
GBM_LC <- GBM_LC %>%
  filter(classification != "unclassified") %>%  # Remove 'unclassified'
  mutate(classification = case_when(
    classification %in% c("CD56bright", "CD56dim") ~ "cNK",
    classification %in% c("ML1", "ML2", "ML_transition") ~ "eML",
    TRUE ~ classification))

# Perform propeller testing
propeller_results_GBM_LC <- propeller(
  clusters = GBM_LC$classification,
  sample = GBM_LC$patient,
  group = GBM_LC$source,
  transform = "asin")

# Print the results
print(propeller_results_GBM_LC)

# Plot cell type proportions
plotCellTypeProps(clusters = GBM_LC$classification, sample = GBM_LC$patient)


##### BrM_BC vs BrM_LC #####

# Subset for BC and LC
BC_LC <- subset(csv, source %in% c("BC", "LC"))
unique(BC_LC$source)

# Recode classifications
BC_LC <- BC_LC %>%
  filter(classification != "unclassified") %>%  # Remove 'unclassified'
  mutate(classification = case_when(
    classification %in% c("CD56bright", "CD56dim") ~ "cNK",
    classification %in% c("ML1", "ML2", "ML_transition") ~ "eML",
    TRUE ~ classification
  ))

# Perform propeller testing
propeller_results_BC_LC <- propeller(
  clusters = BC_LC$classification,
  sample = BC_LC$patient,
  group = BC_LC$source,
  transform = "asin"
)

# Print the results
print(propeller_results_BC_LC)

# Plot cell type proportions
plotCellTypeProps(clusters = BC_LC$classification, sample = BC_LC$patient)

