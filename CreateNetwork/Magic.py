import magic
import scprep
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

%matplotlib inline

data = scprep.io.load_csv('/home/yj345/palmer_scratch/snRNAseq/AllCells/counts/all_RNA_counts.csv')

magic_op = magic.MAGIC()
data_magic = magic_op.fit_transform(data, genes='all_genes')

data_magic.to_csv("/home/yj345/palmer_scratch/snRNAseq/AllCells/counts/all_RNA_counts_magic.csv")
