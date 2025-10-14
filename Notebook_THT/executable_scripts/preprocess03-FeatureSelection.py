#!/usr/bin/env python
# coding: utf-8

import sys
import os
import subprocess
import tempfile
import json

import matplotlib.pyplot as plt
import numpy as np
np.float_ = np.float64
import scanpy as sc
import seaborn as sns

sc.settings.verbosity = 0
sc.settings.set_figure_params(
    dpi=80,
    facecolor="white",
    frameon=False,
)

# ## 挑选高变基因

input_fpath = "RData/merge_normalization.h5ad"
adata = sc.read(input_fpath)

hvg_n=int(sys.argv[1])
sc.pp.highly_variable_genes(adata, layer="log1p_norm",n_top_genes=hvg_n,batch_key="sample")

output_fpath = "RData/hvgn"+sys.argv[1]+"_feature_selection.h5ad"
adata.write(output_fpath)

print("脚本执行完成")
