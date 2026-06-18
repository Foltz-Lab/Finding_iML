# Added updated metadata to Ding et al., pan-cancer data
# Extracted RNA counts, metadata and Umap coordinates from seurat object 

library(Seurat)

meta <- read.table("/storage1/fs1/jennifer.a.foltz/Active/Veda/pan_cancer/Ding_lab/RNA/Old_data/PanImmune_integrated_allRNA_by_chemistry_v8.0_df_T_cells_NK_no_doublets_metadata_v8.7_03262024_from_sublineage_objects_CLEAN.tsv",
                   sep="\t", header=TRUE, row.names = 1)
object <- readRDS("/storage1/fs1/jennifer.a.foltz/Active/Veda/pan_cancer/Ding_lab/RNA/Old_data/PanImmune_integrated_allRNA_by_cancer_chemistry_v8.0_df_T_cells_NK_no_doublets_protein_coding_mito_removed.rds")
# Ensure row names of meta match the column names of the Seurat object
if (!all(rownames(meta) == colnames(object))) {
  stop("Row names of meta do not match column names of the Seurat object. Fix this before proceeding.")
}
setdiff(rownames(meta), colnames(object))  # Returns items in meta but not in Seurat object
setdiff(colnames(object), rownames(meta))

head(meta)
#rownames(meta) <- meta[, 1]
head(rownames(meta))

# Completely replace the metadata with the meta data frame
object@meta.data <- meta

# Verify the replacement
head(object@meta.data)

library(Seurat)
saveRDS(object, "PanImmune_allRNA_cell_type_v8_7_updated_meta_CLEAN.rds")
object<-readRDS("PanImmune_allRNA_cell_type_v8_7_updated_meta_CLEAN.rds")
str(object)
head(object@meta.data)

# Desired cell types
cell_types <- c("NK CD56dim", "trNK CD56bright ZNF683", "NK CD56bright KLRC1", "trNK CD56bright KLRF1")
library(ggplot2)
# Subset Seurat object
object_subset <- subset(object, subset = cell_type_v8.7_rna %in% cell_types)

# Check the subset
print(object_subset)
head(object_subset@meta.data)

unique(object_subset$cell_type_v8.7_rna)
saveRDS(object_subset, file='PanImmune_allRNA_cell_type_v8_7_updated_meta_CLEAN_NKcells.rds')

object<-readRDS('PanImmune_allRNA_cell_type_v8_7_updated_meta_CLEAN_NKcells.rds')
head(object@meta.data)

umap_coords <- as.data.frame(Embeddings(object, reduction = "umap"))
head(umap_coords)

rna_counts <- as.matrix(GetAssayData(object, assay = "RNA", slot = "counts"))
rna_counts[1:5, 1:5] # Display a small portion of the matrix

metadata <- object@meta.data
head(metadata)

# Extract RNA counts
rna_counts <- GetAssayData(object, assay = "RNA", slot = "counts")

# Check if all counts are zero
all_zero <- all(rna_counts == 0)
print(all_zero)  # TRUE if all counts are zero, FALSE otherwise

# Extract metadata
metadata <- object@meta.data

# Save metadata to CSV
write.csv(metadata, file = "meta_CLEAN_nkcells_RNA.csv", row.names = TRUE)

# Extract RNA counts
rna_counts <- as.matrix(GetAssayData(object, assay = "RNA", slot = "counts"))

# Save RNA counts to CSV
write.csv(as.data.frame(rna_counts), file = "rna_counts_CLEAN_nkcells_RNA.csv", row.names = TRUE)

# Extract UMAP coordinates
umap_coords <- as.data.frame(Embeddings(object, reduction = "umap"))

# Save UMAP coordinates to CSV
write.csv(umap_coords, file = "umap_coordinates_CLEAN_nkcells_RNA.csv", row.names = TRUE)
