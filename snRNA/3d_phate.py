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
import scvelo as scv

%load_ext autoreload
%autoreload 2
%matplotlib inline

FL_harmony = pd.read_csv("/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/Phate/FL/FL_harmony.csv", index_col = 0)
HL_harmony = pd.read_csv("/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/Phate/HL/HL_harmony.csv", index_col = 0)
PA_harmony = pd.read_csv("/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/Phate/PA/PA_harmony.csv", index_col = 0)

FL_metadata =  pd.read_csv("/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/Phate/FL/FL_metadate.csv",, index_col = 0)
HL_metadata =  pd.read_csv("/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/Phate/HL/HL_metadate.csv",, index_col = 0)
PA_metadata =  pd.read_csv("/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/Phate/PA/PA_metadate.csv",, index_col = 0)

phate_operator = phate.PHATE(n_jobs=-1, knn = 3, n_pca=100, decay = 15, t = 35, gamma = 0)
phate_operator.set_params(n_components=3)
Y_phate = phate_operator.fit_transform(FL_harmony)

scprep.plot.scatter3d(Y_phate, figsize=(12,8), c = FL_metadata['cellcluster'],
                      ticks=False, label_prefix="PHATE", fontsize=10)
scprep.plot.scatter3d(Y_phate, figsize=(12,8), c = FL_metadata['genotype'],
                      ticks=False, label_prefix="PHATE", fontsize=10)
scprep.plot.scatter3d(Y_phate, figsize=(12,8), c = FL_metadata['stage'],
                      ticks=False, label_prefix="PHATE", fontsize=10)

phate_operator = phate.PHATE(n_jobs=-1, knn = 3, n_pca=100, decay = 15, t = 35, gamma = 0)
phate_operator.set_params(n_components=3)
Y_phate = phate_operator.fit_transform(HL_harmony)

scprep.plot.scatter3d(Y_phate, figsize=(12,8), c = HL_metadata['cellcluster'],
                      ticks=False, label_prefix="PHATE", fontsize=10)
scprep.plot.scatter3d(Y_phate, figsize=(12,8), c = HL_metadata['genotype'],
                      ticks=False, label_prefix="PHATE", fontsize=10)
scprep.plot.scatter3d(Y_phate, figsize=(12,8), c = HL_metadata['stage'],
                      ticks=False, label_prefix="PHATE", fontsize=10)

phate_operator = phate.PHATE(n_jobs=-1, knn = 3, n_pca=100, decay = 15, t = 35, gamma = 0)
phate_operator.set_params(n_components=3)
Y_phate = phate_operator.fit_transform(PA_harmony)

scprep.plot.scatter3d(Y_phate, figsize=(12,8), c = PA_metadata['cellcluster'],
                      ticks=False, label_prefix="PHATE", fontsize=10)
scprep.plot.scatter3d(Y_phate, figsize=(12,8), c = PA_metadata['genotype'],
                      ticks=False, label_prefix="PHATE", fontsize=10)
scprep.plot.scatter3d(Y_phate, figsize=(12,8), c = pA_metadata['stage'],
                      ticks=False, label_prefix="PHATE", fontsize=10)
scprep.plot.scatter3d(Y_phate, figsize=(12,8), c = pA_metadata['tissue'],
                      ticks=False, label_prefix="PHATE", fontsize=10)
