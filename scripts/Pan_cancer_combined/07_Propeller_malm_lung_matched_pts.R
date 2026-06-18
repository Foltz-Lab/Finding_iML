# Propeller tool- to test statistically significant differences - iML vs cNK
# matched tumor and adjacent lung normal
# https://rdrr.io/github/Oshlack/speckle/man/propeller.ttest.html

library(speckle)
library(limma)
library(dplyr)
library(ggplot2)

pan_cancer <- read.csv("/Volumes/jennifer.a.foltz/Active/Veda/pan_cancer/profiling/updated_batch_pan_cancer_v4_dataset_v2_PtID_v2_extra_datasets_removal/pan_cancer_eMLclassified_adata_batch_dataset_v2_after_extra_datasets_removed_classification_allGenes.csv")

##### Lung matched patients #####
lung <- subset(pan_cancer, source %in% c("lung_normal", "lung_tumor"))
multi_source_patients <- lung %>%
  distinct(patient_ID_v2, source) %>%             # Remove duplicate rows
  group_by(patient_ID_v2) %>%
  summarise(n_sources = n_distinct(source)) %>%   # Count unique sources per patient
  filter(n_sources > 1) %>%
  pull(patient_ID_v2)                             # Get patient IDs
# filtering patients with multiple sources and having tumor adjacent lung_normal
lung_f <- lung %>%
  filter((patient_ID_v2 %in% multi_source_patients))

table(lung_f["source"])
lung_f %>%
  group_by(source) %>%
  summarise(unique_patients = n_distinct(patient_ID_v2))

# Subset for lung
unique(lung_f$source)

# Recode classifications
lung_f <- lung_f %>%
  filter(classification != "unclassified") %>%  # Remove 'unclassified'
  mutate(classification = case_when(
    classification %in% c("CD56bright", "CD56dim") ~ "cNK",
    classification %in% c("ML1", "ML2", "ML_transition") ~ "eML",
    TRUE ~ classification
  ))

table(lung_f["source"])

# Perform propeller testing- resulted in warning
propeller_results_lung_f <- propeller(
  clusters = lung_f$classification,
  sample = lung_f$patient_ID_v2,
  group = lung_f$source, transform = "asin"
)

#### paired samples ####
#May be it is due to paired samples from same patient, need to try by paired samples way from propeller

lung_f$sample_id <- paste(lung_f$patient_ID_v2, lung_f$source, sep = "_")
metadata <- lung_f %>%
  distinct(sample_id, patient_ID_v2, source) %>%
  rename(
    sample = sample_id,
    pair = patient_ID_v2,
    group = source
  )
metadata$group <- factor(metadata$group, levels = c("lung_normal", "lung_tumor"))
metadata$pair <- factor(metadata$pair)
design <- model.matrix(~ 0 + group + pair, data = metadata)
colnames(design) <- make.names(colnames(design))
mycontr <- makeContrasts(grouplung_tumor - grouplung_normal, levels = design)

#props <- getTransformedProps(clusters = lung_f$classification, sample = lung_f$sample_id)
props <- getTransformedProps(
  clusters = lung_f$classification,
  sample = lung_f$sample_id, transform = "asin"
)
props <- lapply(props, function(x) x[, metadata$sample])
propeller_paired_ttest <- propeller.ttest(
  prop.list = props,
  design = design,
  contrasts = mycontr,
  robust = TRUE,
  trend = FALSE,
  sort = TRUE)

# Define significance thresholds
p_threshold <- 0.05
fdr_threshold <- 0.05

# Check significance based on p-value and FDR
propeller_paired_ttest$p_significance <- ifelse(propeller_paired_ttest$P.Value < p_threshold, "Significant", "Not Significant")
propeller_paired_ttest$fdr_significance <- ifelse(propeller_paired_ttest$FDR < fdr_threshold, "Significant", "Not Significant")
# Display results
propeller_paired_ttest

barplot(props$Proportions, ,legend=FALSE, 
        ylab="Proportions", xlab="samples",main = "Transformed proportions for NK Cells subsets across samples", xaxt="n") 


