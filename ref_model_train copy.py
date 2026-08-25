import scanpy as sc
import numpy as np
import pandas as pd
import cell2location
import torch
torch.set_num_threads(8) #pytorch is for deep learning, cell2location calls it during model training. default pytorch may use 1 cpu core so we tell it to use 8

#tutorial: https://cell2location.readthedocs.io/en/latest/notebooks/cell2location_tutorial.html

print("Step1", flush=True)
#1 Read data
adata_ref = sc.read('/expanse/lustre/projects/csd670/akao1/Osteosarcoma/osteosarcoma_reference.h5ad')

print("Step2", flush=True)
#2 Filter data
from cell2location.utils.filtering import filter_genes
selected = filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
adata_ref = adata_ref[:, selected].copy()
#cell_count_cutoff = min number of cells expressing that gene
#cell_percentage = min percent of cells expressing that gene
#nonz_mean_cutoff = min mean expression of gene
#this step uses defaults from the tutorial

print("Step3", flush=True)
#3 prepare anndata for the regression model ("registers data, tells cell2location where to find data")
cell2location.models.RegressionModel.setup_anndata(
                        adata=adata_ref,
                        # different batch/patient, so gene expression difference treated as batch effect
                        batch_key='sample',
                        # cell type, covariate/colname used for constructing signatures
                        labels_key='cell2loc_labels')

print("Step4", flush=True)
#4 create regression model
from cell2location.models import RegressionModel
mod = RegressionModel(adata_ref)

print("Step5", flush=True)
#5 train model
mod.train(max_epochs=250)
#saving the model data at designated directory
mod.save("/expanse/lustre/projects/csd670/akao1/Osteosarcoma/ref_model", overwrite=True)
print("model saved in directory", flush=True)

print("Step6", flush=True)
#6 exported estimated cell abundance (gene expression levels for each cell) as posterior distribution from bayesian model
#adds on model data to adata_ref object
adata_ref = mod.export_posterior( #batch size means step run in increments to decrease crashing
    adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500}
)

print("Step7", flush=True)
#7 extracting genes x cell type dataframe for downstream use
# export estimated expression in each cluster
inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
            for i in adata_ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_ref.uns['mod']['factor_names'] #list of cell names as colnames

inf_aver.to_csv("/expanse/lustre/projects/csd670/akao1/Osteosarcoma/osteosarcoma_inf_aver.csv")
print("downstream-use dataframe saved in directory", flush=True)
print("script finished", flush=True)

