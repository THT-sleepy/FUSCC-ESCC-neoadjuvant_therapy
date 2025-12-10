#!/usr/bin/env python
# coding: utf-8

# ## 修改路径为conda路径

# In[1]:

import sys
import os
import harmonypy as hm


import scanpy as sc

sc.settings.verbosity = 0
sc.settings.set_figure_params(
    dpi=80,
    facecolor="white",
    frameon=False
)
sc.settings.figdir = "plots"
sc.settings.autoshow = False

# In[9]:

input_fpath="RData/hvgn"+sys.argv[1]+"_feature_selection.h5ad"
adata = sc.read(
    filename=input_fpath
)


adata.X = adata.layers["log1p_norm"]


# setting highly variable as highly deviant to use scanpy 'use_highly_variable' argument in sc.pp.pca
#adata.var["highly_variable"] = adata.var["highly_deviant"]
sc.pp.pca(adata, svd_solver="arpack", use_highly_variable=True)

#pcaplt_path = "hvg"+ sys.argv[1] + "pre04_pca.pdf"
#sc.pl.pca_scatter(adata, color="total_counts",save=pcaplt_path)

#sc.tl.tsne(adata, use_rep="X_pca")

#tsneplt_path = "hvg"+ sys.argv[1] +"pre04_tsne.pdf"
#sc.pl.tsne(adata, color="total_counts",save=tsneplt_path)

sc.pp.neighbors(adata) #默认PC=50
sc.tl.umap(adata)


# In[16]:

umapplt_path = "hvg"+ sys.argv[1] +"umap.pdf"
sc.pl.umap(adata, color="total_counts",save=umapplt_path)

umapplt_path = "hvg"+ sys.argv[1] +"umap_annoQC.pdf"
sc.pl.umap(
    adata,
    color=["n_genes_by_counts", "pct_counts_mt", "scDblFinder_score", "scDblFinder_class"],
    save=umapplt_path
)
def assign_group(sample_name):
    s = str(sample_name).lower()
    if 'pbmc' in s:
        return 'pbmc_op'
    elif 'pre' in s:
        return 'pre'
    else:
        return 'op'

adata.obs['sample_group'] = adata.obs['sample'].apply(assign_group)

umapplt_path = "hvg" +sys.argv[1] + "umap_annoBatch.pdf"
sc.pl.umap(
    adata,
    color=["sample_group"],
    save=umapplt_path
)
output_fpath="RData/hvg"+ sys.argv[1] +"dimensionality_reduction.h5ad"
adata.write(output_fpath)

print("Over!")
