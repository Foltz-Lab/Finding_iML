# Propeller tool- to test statistically significant differences 
# between iML and cNK cell types across difference source types
# malm primary LUAD vs fan BrM_LUAD

library(speckle)
library(limma)
library(dplyr)

csv2<-read.csv("/Volumes/jennifer.a.foltz/Active/Veda/pan_cancer/adata_combined_malm_ding_Tcell_filtered_fan_included/output/fan_BrM_LUAD_vs_malm_prim_LUAD_matched_pts.csv")

##### #lung_tumor vs BrM_NSCLC-LUAD #####

# Subset for lung_tumor and NSCLC-LUAD
lung_NSCLC <- subset(csv2, source %in% c("lung_tumor", "NSCLC-LUAD"))
table(lung_NSCLC$source)

# Recode classifications
lung_NSCLC <- lung_NSCLC %>%
  filter(classification != "unclassified") %>%  # Remove 'unclassified'
  mutate(classification = case_when(
    classification %in% c("CD56bright", "CD56dim") ~ "cNK",
    classification %in% c("ML1", "ML2", "ML_transition") ~ "eML",
    TRUE ~ classification
  ))

# Perform propeller testing
propeller_results_lung_NSCLC <- propeller(
  clusters = lung_NSCLC$classification,
  sample = lung_NSCLC$patient,
  group = lung_NSCLC$source,
  transform = "asin"
)

# Print the results
print(propeller_results_lung_NSCLC)

# Plot cell type proportions
plotCellTypeProps(clusters = lung_NSCLC$classification, sample = lung_NSCLC$patient)

