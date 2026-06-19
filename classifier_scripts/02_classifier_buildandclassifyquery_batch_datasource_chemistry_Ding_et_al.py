#Implementation of Machine learning classifier on Ding et al., data with Data.source_Chemistry as batch 
#Docker used - evelyns2000/foltz_tools:classifier2024 
#classifier script - https://www.science.org/doi/abs/10.1126/sciimmunol.adk4893

# input your variables below:
patient = "nk_ding"
protein = False
batch = 'Data.source_Chemistry' # be sure to specify batch for appropriate read in of vae

# the ref variables are the same for all queries
ref_model = "totalvi_vae_reference_model_withclassifiers" 
ref_adata = "fig4_reference_adata_train_posttotalvi.h5ad"

# import libraries required
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

import scanpy as sc
import os
import anndata
import pandas as pd
from matplotlib import rcParams
import numpy as np
from scvi.model import TOTALVI

sc.settings.verbosity = 3
sc.logging.print_header()
sc.settings.set_figure_params(dpi=200, frameon=False, figsize=(9,3), facecolor = 'white', color_map = 'magma', format = 'png')

# load RNA counts, metadata, ADT and umap files 
RNApath = f"rna_counts_CLEAN_nkcells_RNA.csv"
metapath =  f"meta_CLEAN_nkcells_RNA_DataSource_Chemistry.csv"
ADTpath = f"ADTcounts.csv"
umappath = "/storage1/fs1/jennifer.a.foltz/Active/Veda/pan_cancer/Ding_lab/RNA/updated_batch_datasource_chemistry/umap_coordinates_CLEAN_nkcells_RNA.csv"

# where adata is query, and ref is reFerence anndata file
adata = sc.read_csv(RNApath)
adata = adata.transpose()

# inspect the barcodes
print('RNA Indexes Before')
print(adata.obs.index)

adata.obs.index = adata.obs.index.str.replace('\.', '-', regex = True)
print('RNA Indexes After')
print(adata.obs.index)
print()

# read and inspect metadata
meta = pd.read_csv(metapath, index_col = 0)
adata.obs = meta
print('Meta Indexes')
print(meta.index)
print()

ref = sc.read_h5ad(ref_adata)
if protein == True:
# Creating the reference adata object from csv's exported from seurat, adding protein counts data ==========
    protein_adata = sc.read_csv(ADTpath)
    protein_adata = protein_adata.transpose()

# inspect:
    print('Protein Indexes Before')
    print(protein_adata.obs.index)
# if needed, uncomment next two lines
    protein_adata.obs.index = protein_adata.obs.index.str.replace('X', '') # note default value of regex changing in other versions, current default is True
    protein_adata.obs.index = protein_adata.obs.index.str.replace('\.', '-', regex = True)
    print('Protein Indexes After')
    print(protein_adata.obs.index)
    print()

    protein_adata.obs = meta
    adata.obsm["protein_expression"] = protein_adata.to_df()

# now confirm equality in 2 ways!
# first is sum where True is = 1 so if all equal should be the total length of the two indexes being compared
    print('Confirm the indexes are equal')
    print(sum(protein_adata.obs.index == adata.obs.index))
    print(sum(meta.index == protein_adata.obs.index))
# and option 2 through pandas:
    print(protein_adata.obs.index.equals(adata.obs.index))
    print(protein_adata.obs.index.equals(meta.index))
    adata.obsm["protein_expression"].columns = adata.obsm["protein_expression"].columns.str.replace('-TotalSeqC', 'ADT', regex = True) # make names for proteins the same
    proteins_to_check = ['IgG2aADT', 'IgG2bADT', 'IgG1ADT', 'PD-1ADT', 'CD8ADT', 'KIR2DL1-S1-S3-S5ADT', 'KIR2DL2-3ADT',
       'KIR3DL1ADT', 'KIR2DL5ADT'] # remove isotype controls, proteins only present in some batches, and KIR proteins

# see what proteins are present before removal
    print(adata.obsm["protein_expression"].columns)
    print(len(adata.obsm["protein_expression"].columns))

    adata.obsm["protein_expression"] = adata.obsm["protein_expression"][adata.obsm["protein_expression"].columns.difference(proteins_to_check)]

# verify proteins were removed
    print(adata.obsm["protein_expression"].columns)
    print(len(adata.obsm["protein_expression"].columns))
elif protein == False:
    # if no protein data:
# put matrix of zeros for protein expression (considered missing)
    pro_exp = ref.obsm["protein_expression"]
    data = np.zeros((adata.n_obs, pro_exp.shape[1]))
    adata.obsm["protein_expression"] = pd.DataFrame(
        columns=pro_exp.columns, index=adata.obs_names, data=data
        )
print(adata.obsm["protein_expression"]) # visualize query
print(ref.obsm["protein_expression"]) # visualize ref
print(sum(adata.obsm["protein_expression"].columns == ref.obsm["protein_expression"].columns)) # confirm protein features are the same
# below for all patient samples regardless of whether they have ADT data
print(sum(meta.index == adata.obs.index))
print(adata.obs.index.equals(meta.index))
adata.obs['batch'] = adata.obs[batch]

adata.raw = adata

# copy raw counts to counts layer
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

adata.raw = adata # move normalized data to adata raw

adata.obs["celltype.l2"] = "Unknown"

# reorganize query proteins, missing proteins become all 0
for p in ref.obsm["protein_expression"].columns:
    if p not in adata.obsm["protein_expression"].columns:
        adata.obsm["protein_expression"][p] = 0.0

# ensure columns are in same order
adata.obsm["protein_expression"] = adata.obsm["protein_expression"].loc[
    :, ref.obsm["protein_expression"].columns]

ref.obs["dataset_name"] = "Reference"
adata.obs["dataset_name"] = "Query"

#adata.write(f'{patient}_prepped.h5ad')

# begin classification
vae = TOTALVI.load(ref_model, ref)

TOTALVI.prepare_query_anndata(adata, reference_model = vae)
vae_q = TOTALVI.load_query_data(adata, vae)

vae_q.train(
     max_epochs=150,
     plan_kwargs=dict(weight_decay=0.0, scale_adversarial_loss=0.0)
) 

adata.obsm["X_totalvi_scarches"] = vae_q.get_latent_representation(adata)

predictionsbbc = vae_q.latent_space_classifer_bbc_.predict(adata.obsm["X_totalvi_scarches"])
#optional BRF classifier inclusion
predictionsbrf = vae_q.latent_space_classifer_brf_.predict(adata.obsm["X_totalvi_scarches"])

probsbbc = vae_q.latent_space_classifer_bbc_.predict_proba(adata.obsm["X_totalvi_scarches"])
probsbrf = vae_q.latent_space_classifer_brf_.predict_proba(adata.obsm["X_totalvi_scarches"])

# check indexes 
df_bbc = pd.DataFrame(probsbbc, columns = vae_q.latent_space_classifer_bbc_.classes_, index = adata.obs_names) # add names to probabilities array and convert to df
df_brf = pd.DataFrame(probsbrf, columns = vae_q.latent_space_classifer_brf_.classes_, index = adata.obs_names) # add names to probabilities array and convert to df

df_bbc.to_csv(f'{patient}_probabilitiesBBCoutput.csv') # may need to to do this after pd.dataframe
df_brf.to_csv(f'{patient}_probabilitiesBRFoutput.csv')

dfi_bbc = df_bbc.loc[adata.obs_names] # IMPORTANT organize df so that it is in the same order as the adata object! since the barcodes are not checked internally
dfi_brf = df_brf.loc[adata.obs_names] # IMPORTANT organize df so that it is in the same order as the adata object! since the barcodes are not checked internally

# now add the metadata
adata.obs["CD56brightBBCprob"] = dfi_bbc["CD56bright"]
adata.obs["CD56dimBBCprob"] = dfi_bbc["CD56dim"]
adata.obs["eML1BBCprob"] = dfi_bbc["ML1"]
adata.obs["eML2BBCprob"] = dfi_bbc["ML2"]

adata.obs["CD56brightBRFprob"] = dfi_brf["CD56bright"]
adata.obs["CD56dimBRFprob"] = dfi_brf["CD56dim"]
adata.obs["eML1BRFprob"] = dfi_brf["ML1"]
adata.obs["eML2BRFprob"] = dfi_brf["ML2"]

adata.obs["predictionsBBC"] = predictionsbbc
adata.obs["predictionsBRF"] = predictionsbrf

print(pd.crosstab(adata.obs['predictionsBBC'],adata.obs['Sample_RNA'], normalize = 'columns'))
print(pd.crosstab(adata.obs['predictionsBRF'],adata.obs['Sample_RNA'], normalize = 'columns'))

predsbytime = pd.crosstab(adata.obs['predictionsBBC'],adata.obs['Sample_RNA'], normalize = 'columns') # now stash the data for later in an easy to use format
predsbytime.to_csv(f'{patient}_cumulative_probabilitiesBBCbytime.csv') 


umap = pd.read_csv(umappath, index_col = 0)
adata = adata[umap.index] 
umap = umap.to_numpy() # convert to array
adata.obsm["X_umap"] = umap

sc.pl.umap(adata, color = ['CD56brightBBCprob', 'CD56dimBBCprob', 'eML1BBCprob', 'eML2BBCprob'], cmap = 'viridis', save = f'_{patient}_probsofBBCclassifier.png')
sc.pl.umap(adata, color = ['CD56brightBRFprob', 'CD56dimBRFprob', 'eML1BRFprob', 'eML2BRFprob'], cmap = 'viridis', save = f'_{patient}_probsofBRFclassifier.png') 
sc.pl.umap(adata, color = ['predictionsBBC', 'predictionsBRF'], save = f'_{patient}_predictionsofclassifier.png') 

# Check data types in adata.obs
print(adata.obs.dtypes)

# Identify columns with non-string objects
for col in adata.obs.columns:
    if not all(isinstance(x, str) or pd.isnull(x) for x in adata.obs[col]):
        print(f"Column '{col}' has non-string values.")

# Repeat for other identified columns
for col in adata.obs.columns:
    if not all(isinstance(x, str) or pd.isnull(x) for x in adata.obs[col]):
        adata.obs[col] = adata.obs[col].astype(str)

# Check if all columns in obs are valid
for col in adata.obs.columns:
    if not all(isinstance(x, str) or pd.isnull(x) for x in adata.obs[col]):
        raise ValueError(f"Column '{col}' still contains invalid values!")

adata.write_h5ad(f'{patient}_eMLclassified_adata_batch_DataSource_Chemistry.h5ad')
vae_q.save(f'{patient}withtotalvi_vae_model_withclassifiers_batch_DataSource_Chemistry', overwrite=True)
