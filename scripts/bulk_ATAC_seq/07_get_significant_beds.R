#########################################################################################
####  In this script, we filter our differential peaks bed files on an |FDR| > 0.5. These files only include 
####  peaks with a p-value < 0.05, so we do not filter on p-value here. The output of this script will be used
####  with deepTools to generate heatmaps and profiles of significant DARs.
#########################################################################################

library(tidyverse)

# 1. Read the data
eML_contrast8 <- read_csv("20230829_contrast8_eML_SIG.csv")

eML_sig8 <- eML_contrast8[eML_contrast8$FDR < 0.05,]

# 2. Clean and select proper BED columns
diffpeaksbed_eMLsig8 <- eML_sig8 %>%
  filter(!is.na(seqnames), !is.na(start), !is.na(end)) %>%
  select(seqnames, start, end, strand)

# 3. Write out for deepTools
write.table(diffpeaksbed_eMLsig8, 
            file="eML_contrast8_diffpeaks_significant.bed", 
            sep="\t", 
            quote=FALSE, 
            row.names=FALSE,
            col.names=FALSE)



# 1. Read the data
eML_contrast9 <- read_csv("20230829_contrast9_eML_SIG.csv")

eML_sig9 <- eML_contrast9[eML_contrast9$FDR < 0.05,]

# 2. Clean and select proper BED columns
diffpeaksbed_eMLsig9 <- eML_sig9 %>%
  filter(!is.na(seqnames), !is.na(start), !is.na(end)) %>%
  select(seqnames, start, end, strand)

# 3. Write out for deepTools
write.table(diffpeaksbed_eMLsig9,
            file="eML_contrast9_diffpeaks_significant.bed",
            sep="\t",
            quote=FALSE,
            row.names=FALSE,
            col.names=FALSE)


# 1. Read the data
iML_contrast9 <- read_csv("20230829_contrast9_iML_SIG.csv")

iML_sig9 <- iML_contrast9[iML_contrast9$FDR < 0.05,]

# 2. Clean and select proper BED columns
diffpeaksbed_iMLsig9 <- iML_sig9 %>%
  filter(!is.na(seqnames), !is.na(start), !is.na(end)) %>%
  select(seqnames, start, end, strand)

# 3. Write out for deepTools
write.table(diffpeaksbed_iMLsig9,
            file="iML_contrast9_diffpeaks_significant.bed",
            sep="\t",
            quote=FALSE,
            row.names=FALSE,
            col.names=FALSE)

# 1. Read the data
iML_contrast12 <- read_csv("20230829_contrast12_iML_SIG.csv")

iML_sig12 <- iML_contrast12[iML_contrast12$FDR < 0.05,]

# 2. Clean and select proper BED columns
diffpeaksbed_iMLsig12 <- iML_sig12 %>%
  filter(!is.na(seqnames), !is.na(start), !is.na(end)) %>%
  select(seqnames, start, end, strand)

# 3. Write out for deepTools
write.table(diffpeaksbed_iMLsig12,
            file="iML_contrast12_diffpeaks_significant.bed",
            sep="\t",
            quote=FALSE,
            row.names=FALSE,
            col.names=FALSE)

# 1. Read the data
cNK_contrast8 <- read_csv("20230829_contrast8_cNK_SIG.csv")

cNK_sig8 <- cNK_contrast8[cNK_contrast8$FDR < 0.05,]

# 2. Clean and select proper BED columns
diffpeaksbed_cNKsig8 <- cNK_sig8 %>%
  filter(!is.na(seqnames), !is.na(start), !is.na(end)) %>%
  select(seqnames, start, end, strand)

# 3. Write out for deepTools
write.table(diffpeaksbed_cNKsig8,
            file="cNK_contrast8_diffpeaks_significant.bed",
            sep="\t",
            quote=FALSE,
            row.names=FALSE,
            col.names=FALSE)


# 1. Read the data
cNK_contrast12 <- read_csv("20230829_contrast12_cNK_SIG.csv")

cNK_sig12 <- cNK_contrast12[cNK_contrast12$FDR < 0.05,]

# 2. Clean and select proper BED columns
diffpeaksbed_cNKsig12 <- cNK_sig12 %>%
  filter(!is.na(seqnames), !is.na(start), !is.na(end)) %>%
  select(seqnames, start, end, strand)

# 3. Write out for deepTools
write.table(diffpeaksbed_cNKsig12,
            file="cNK_contrast12_diffpeaks_significant.bed",
            sep="\t",
            quote=FALSE,
            row.names=FALSE,
            col.names=FALSE)
