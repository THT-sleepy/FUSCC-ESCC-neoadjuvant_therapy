#!/usr/bin/env python
# coding: utf-8

# ## 修改路径为conda路径

# In[1]:

import sys
import os



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

input_fpath="RData/"+"mad"+sys.argv[1]+"_"+"hvgn"+sys.argv[2]+"_feature_selection.h5ad"
adata = sc.read(
    filename=input_fpath
)


# In[10]:


adata.X = adata.layers["log1p_norm"]


# In[11]:


# setting highly variable as highly deviant to use scanpy 'use_highly_variable' argument in sc.pp.pca
adata.var["highly_variable"] = adata.var["highly_deviant"]
sc.pp.pca(adata, svd_solver="arpack", use_highly_variable=True)


# In[12]:

pcaplt_path = "mad"+sys.argv[1] + "hvg"+ sys.argv[2] + "pre04_pca.pdf"
sc.pl.pca_scatter(adata, color="total_counts",save=pcaplt_path)


# In[13]:


sc.tl.tsne(adata, use_rep="X_pca")


# In[14]:

tsneplt_path = "mad"+sys.argv[1] + "hvg"+ sys.argv[2] +"pre04_tsne.pdf"
sc.pl.tsne(adata, color="total_counts",save=tsneplt_path)


# In[15]:


sc.pp.neighbors(adata)
sc.tl.umap(adata)


# In[16]:

umapplt_path1 = "mad"+sys.argv[1] + "hvg"+ sys.argv[2] +"pre04_umap.pdf"
sc.pl.umap(adata, color="total_counts",save=umapplt_path1)


# In[17]:

umapplt_path2 = "mad"+sys.argv[1] + "hvg"+ sys.argv[2] +"pre04_umap_annoQC.pdf"
sc.pl.umap(
    adata,
    color=["total_counts", "pct_counts_mt", "scDblFinder_score", "scDblFinder_class"],
    save=umapplt_path2
)
##基因数目分布直方图
sc.pp.calculate_qc_metrics(adata, inplace=True)
import matplotlib.pyplot as plt
n_genes_per_cell = adata.obs['n_genes_by_counts']
plt.figure(figsize=(10, 6))
plt.hist(n_genes_per_cell, bins=50, color='skyblue', edgecolor='black')
plt.title('Distribution of Number of Genes per Cell', fontsize=16)
plt.xlabel('Number of Genes', fontsize=12)
plt.ylabel('Number of Cells', fontsize=12)
plt.grid(axis='y', linestyle='--', alpha=0.7)
barplot_fpath = "plots/" + "mad" + sys.argv[1] + "n_genes_distribution.pdf"
plt.savefig(barplot_fpath, dpi=300)
plt.close()


# In[18]:

output_fpath="RData/"+ "mad"+sys.argv[1]+"_"+"hvg"+ sys.argv[2] +"_dimensionality_reduction.h5ad"
adata.write(output_fpath)

