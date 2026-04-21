import os
import scprep
import phate
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import scprep
from PIL import Image
from anndata import AnnData
from pyslingshot import Slingshot
import anndata as ad

%load_ext autoreload
%autoreload 2
%matplotlib inline

projDir = '/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/FL/'

path_to_RNA_counts = projDir + 'FL_Cartilage_RNA_counts.csv'
path_to_metadata = projDir + 'FL_Cartilage_metadata.csv'
path_to_harmony = projDir + 'FL_Cartilage_harmony.csv'

counts = pd.read_csv(path_to_RNA_counts, index_col=0)
counts = scprep.filter.filter_rare_genes(counts, min_cells=10)
counts = scprep.normalize.library_size_normalize(counts)
counts = scprep.transform.sqrt(counts)

metadata =  pd.read_csv(path_to_metadata, index_col = 0)

harmony = pd.read_csv(path_to_harmony, index_col = 0)

phate_operator = phate.PHATE(n_jobs=-1, knn = 3, n_pca=50, decay = 50, t = 50, gamma = 0)
Y_phate = phate_operator.fit_transform(harmony)

scprep.plot.scatter2d(Y_phate, figsize=(12,8), c = metadata['cellcluster'],
                      ticks=False, label_prefix="PHATE", fontsize=10)
scprep.plot.scatter2d(Y_phate, figsize=(12,8), c = metadata['stage'],
                      ticks=False, label_prefix="PHATE", fontsize=10)
scprep.plot.scatter2d(Y_phate, figsize=(12,8), c = metadata['genotype'],
                      ticks=False, label_prefix="PHATE", fontsize=10)

adata = AnnData(X=counts, obs=metadata)
adata.obsm["X_phate"] = Y_phate

sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.pl.highly_variable_genes(adata)
sc.tl.pca(adata)
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)
sc.pp.neighbors(adata, use_rep='X_phate')
sc.tl.leiden(adata, resolution=1)
sc.pl.embedding(adata, basis='phate', color='leiden')
adata.write("/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/FL/FL_Cartilage.h5ad")

projDir = '/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/HL/'

path_to_RNA_counts = projDir + 'HL_Cartilage_RNA_counts.csv'
path_to_metadata = projDir + 'HL_Cartilage_metadata.csv'
path_to_harmony = projDir + 'HL_Cartilage_harmony.csv'

counts = pd.read_csv(path_to_RNA_counts, index_col=0)
counts = scprep.filter.filter_rare_genes(counts, min_cells=10)
counts = scprep.normalize.library_size_normalize(counts)
counts = scprep.transform.sqrt(counts)

metadata =  pd.read_csv(path_to_metadata, index_col = 0)

harmony = pd.read_csv(path_to_harmony, index_col = 0)

phate_operator = phate.PHATE(n_jobs=-1, knn = 3, n_pca=50, decay = 50, t = 50, gamma = 0)
Y_phate = phate_operator.fit_transform(harmony)

scprep.plot.scatter2d(Y_phate, figsize=(12,8), c = metadata['cellcluster'],
                      ticks=False, label_prefix="PHATE", fontsize=10)
scprep.plot.scatter2d(Y_phate, figsize=(12,8), c = metadata['stage'],
                      ticks=False, label_prefix="PHATE", fontsize=10)
scprep.plot.scatter2d(Y_phate, figsize=(12,8), c = metadata['genotype'],
                      ticks=False, label_prefix="PHATE", fontsize=10)

adata = AnnData(X=counts, obs=metadata)
adata.obsm["X_phate"] = Y_phate

sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.pl.highly_variable_genes(adata)
sc.tl.pca(adata)
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)
sc.pp.neighbors(adata, use_rep='X_phate')
sc.tl.leiden(adata, resolution=1)
sc.pl.embedding(adata, basis='phate', color='leiden')

adata.write("/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/HL/HL_Cartilage.h5ad")

projDir = '/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/PA/'

path_to_RNA_counts = projDir + 'PA_Osteoblasts_RNA_counts.csv'
path_to_metadata = projDir + 'PA_Osteoblasts_metadata.csv'
path_to_harmony = projDir + 'PA_Osteoblasts_harmony.csv'

counts = pd.read_csv(path_to_RNA_counts, index_col=0)
counts = scprep.filter.filter_rare_genes(counts, min_cells=10)
counts = scprep.normalize.library_size_normalize(counts)
counts = scprep.transform.sqrt(counts)

metadata =  pd.read_csv(path_to_metadata, index_col = 0)

harmony = pd.read_csv(path_to_harmony, index_col = 0)

phate_operator = phate.PHATE(n_jobs=-1, knn = 3, n_pca=50, decay = 50, t = 50, gamma = 0)
Y_phate = phate_operator.fit_transform(harmony)

scprep.plot.scatter2d(Y_phate, figsize=(12,8), c = metadata['cellcluster'],
                      ticks=False, label_prefix="PHATE", fontsize=10)
scprep.plot.scatter2d(Y_phate, figsize=(12,8), c = metadata['stage'],
                      ticks=False, label_prefix="PHATE", fontsize=10)
scprep.plot.scatter2d(Y_phate, figsize=(12,8), c = metadata['genotype'],
                      ticks=False, label_prefix="PHATE", fontsize=10)
scprep.plot.scatter2d(Y_phate, figsize=(12,8), c = metadata['tissue'],
                      ticks=False, label_prefix="PHATE", fontsize=10)

adata = AnnData(X=counts, obs=metadata)
adata.obsm["X_phate"] = Y_phate

sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.pl.highly_variable_genes(adata)
sc.tl.pca(adata)
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)
sc.pp.neighbors(adata, use_rep='X_phate')
sc.tl.leiden(adata, resolution=1)
sc.pl.embedding(adata, basis='phate', color='leiden')

adata.write("/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/PA/PA_Osteoblasts.h5ad")

