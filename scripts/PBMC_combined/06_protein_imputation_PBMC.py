#Protein imputation - PBMC subsetted malmberg data with batch as dataset_v2 unique variables
#Docker used - evelyns2000/foltz_tools:classifier2024 

#!/usr/bin/env python
# coding: utf-8
#https://docs.scvi-tools.org/en/1.0.0/tutorials/notebooks/scarches_scvi_tools.html#impute-protein-data-for-query-and-visualize

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

import scanpy as sc
import pandas as pd
from scvi.model import TOTALVI
import psutil
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages  # To save plots into a single PDF file

# Function to check system resource usage
def check_resource_usage():
    print(f"Memory usage: {psutil.virtual_memory().percent}%")
    print(f"CPU usage: {psutil.cpu_percent(interval=1)}%")

# Step 1: Load the original AnnData object
adata = sc.read_h5ad('/storage1/fs1/jennifer.a.foltz/Active/Veda/pan_cancer/profiling/updated_batch_pan_cancer_v4_dataset_v2_PtID_v2_extra_datasets_removal/pan_cancer_eMLclassified_adata_batch_dataset_v2_after_extra_datasets_removed_classification.h5ad')


# Step 2: Filter to PBMC only
tissue_data = adata[adata.obs['source'] == 'PBMC'].copy()

if tissue_data.n_obs == 0:
    raise ValueError("No PBMC samples found in the dataset.")

# Step 3: Get unique values from dataset_v2 for transform_batch
transform_batch_values = tissue_data.obs['dataset_v2'].unique().tolist()
print(f"Transform batch values from dataset_v2: {transform_batch_values}")

# Step 4: Load the TOTALVI model
try:
    vae_q = TOTALVI.load("/storage1/fs1/jennifer.a.foltz/Active/Veda/pan_cancer/profiling/updated_batch_pan_cancer_v4_dataset_v2_PtID_v2_extra_datasets_removal/classifier_output/pan_cancerwithtotalvi_vae_model_withclassifiers_batch_dataset_v2_after_extra_datasets_removed", tissue_data)
except Exception as e:
    raise RuntimeError(f"Error loading model for PBMC: {e}")

check_resource_usage()

# Step 5: Impute protein expression
_, imputed_proteins = vae_q.get_normalized_expression(
    tissue_data,
    n_samples=25,
    return_mean=True,
    transform_batch=transform_batch_values,
)

print(f"Imputed proteins for PBMC: {imputed_proteins.shape[1]}")

# Step 6: Concatenate imputed proteins with observations
tissue_data.obs = pd.concat([tissue_data.obs, imputed_proteins], axis=1)

tissue_data.obs['classification'] = tissue_data.obs['classification'].replace({
    'ML1': 'iML',
    'ML2': 'iML',
    'ML_transition': 'iML'
})
tissue_data = tissue_data[tissue_data.obs["classification"] != "unclassified"].copy()
tissue_data.obs['classification'] = tissue_data.obs['classification'].astype('category')

# Step 7: Define custom color palette
custom_colors = [ '#1f77b4' ,'green',"#CD96CD"]

# Get unique classifications and map colors to them
classifications = tissue_data.obs['classification'].unique().tolist()
color_palette = {cls: custom_colors[i % len(custom_colors)] for i, cls in enumerate(classifications)}

print(f"Color mapping: {color_palette}")

# Step 8: Generate and save all plots to PDF
pdf_filename = '/storage1/fs1/jennifer.a.foltz/Active/Veda/pan_cancer/adata_combined_malm_ding_Tcell_filtered_fan_included/PBMC_combined_Malm_ding_Tcell_filtered/output/PBMC_protein_expression_plots.pdf'
print("Generating violin plots...")
with PdfPages(pdf_filename) as pdf:
    for protein in imputed_proteins.columns:
        fig, ax = plt.subplots(figsize=(10, 6))
        sns.violinplot(
            x='classification', y=protein,
            data=tissue_data.obs,
            ax=ax,
            inner='box',
            scale='width',
            cut=0,
            palette=color_palette
        )
        ax.set_title(f'PBMC: Violin plot of {protein} expression')
        ax.set_ylabel(f'{protein} expression')
        ax.set_xlabel('Classification')
        plt.xticks(ha='right')
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

print(f'All violin plots for PBMC saved in {pdf_filename}')

check_resource_usage()

# Step 10: Save specific proteins as individual PNG images
proteins_of_interest = ['CD117ADT', 'CD57ADT', 'NKG2AADT']
image_output_dir = '/storage1/fs1/jennifer.a.foltz/Active/Veda/pan_cancer/adata_combined_malm_ding_Tcell_filtered_fan_included/PBMC_combined_Malm_ding_Tcell_filtered/output/'

print("Saving selected protein plots as PNG images...")
for protein in proteins_of_interest:
    if protein not in imputed_proteins.columns:
        print(f"Warning: {protein} not found in imputed proteins, skipping.")
        continue

    fig, ax = plt.subplots(figsize=(8, 6))
    sns.violinplot(
        x='classification', y=protein,
        data=tissue_data.obs,
        ax=ax,
        inner='box',
        scale='width',
        cut=0,
        palette=color_palette,
    )
    ax.set_title(f'PBMC: Violin plot of {protein} expression' , fontsize= 18)
    ax.set_ylabel(f'{protein} expression', fontsize= 18)
    ax.set_xlabel('Classification', fontsize= 18)
    plt.xticks(fontsize= 18)
    plt.tight_layout()
    fig.savefig(f'{image_output_dir}{protein}_violin_plot.png', dpi=600, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved {protein}_violin_plot.png')

print("Done saving individual protein images.")