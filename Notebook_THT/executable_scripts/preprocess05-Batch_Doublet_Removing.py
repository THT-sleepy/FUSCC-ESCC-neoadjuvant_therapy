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

input_fpath="RData/hvg"+ sys.argv[1] +"dimensionality_reduction.h5ad"
adata = sc.read(
    filename=input_fpath
) #X是log1p_norm

#去掉双细胞
adata = adata[adata.obs.scDblFinder_class == 'singlet'].copy()

#重新normalize和feature selection
adata.X = adata.layers["soupX_counts"]

##normalize
scales_counts = sc.pp.normalize_total(adata, target_sum=None, inplace=False)
# log1p transform
adata.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)

##feature selection
hvg_n=int(sys.argv[1])
sc.pp.highly_variable_genes(adata, layer="log1p_norm",n_top_genes=hvg_n,batch_key="sample")

adata.X = adata.layers["log1p_norm"]

##pca
sc.pp.pca(adata, svd_solver="arpack", use_highly_variable=True)
plt_path="hvg"+sys.argv[1]+".png"
sc.pl.pca_variance_ratio(adata, log=True,save=plt_path)

#min dist=0.05
adata1 = adata.copy()
sc.external.pp.bbknn(adata1, batch_key="sample",n_pcs=30)  # running bbknn 1.6.0
sc.tl.umap(adata1,min_dist=0.05)
umapplt_path = "hvg" +sys.argv[1] + "PC30"+"_bbknn_batch_corrected" + "_doubletremoved_umapmd0.05.pdf"
sc.pl.umap(
    adata1,
    color=["sample_group","n_genes_by_counts", "pct_counts_mt","scDblFinder_class"],
    save=umapplt_path
)
output_fpath="RData/hvg"+ sys.argv[1] +"_PC30"+"_batchdoubletremoved_umapmd0.05.h5ad"
adata1.write(output_fpath)

#min dist=0.1
adata2 = adata.copy()
sc.external.pp.bbknn(adata2, batch_key="sample",n_pcs=30)  # running bbknn 1.6.0
sc.tl.umap(adata2,min_dist=0.1)
umapplt_path = "hvg" +sys.argv[1] + "PC30"+"_bbknn_batch_corrected" + "_doubletremoved_umapmd0.1.pdf"
sc.pl.umap(
    adata2,
    color=["sample_group","n_genes_by_counts", "pct_counts_mt","scDblFinder_class"],
    save=umapplt_path
)
output_fpath="RData/hvg"+ sys.argv[1] +"_PC30"+"_batchdoubletremoved_umapmd0.1.h5ad"
adata2.write(output_fpath)

#min dist=0.2
adata3 = adata.copy()
sc.external.pp.bbknn(adata3, batch_key="sample",n_pcs=30)  # running bbknn 1.6.0
sc.tl.umap(adata3,min_dist=0.2)
umapplt_path = "hvg" +sys.argv[1] + "PC30"+"_bbknn_batch_corrected" + "_doubletremoved_umapmd0.2.pdf"
sc.pl.umap(
    adata3,
    color=["sample_group","n_genes_by_counts", "pct_counts_mt","scDblFinder_class"],
    save=umapplt_path
)
output_fpath="RData/hvg"+ sys.argv[1] +"_PC30"+"_batchdoubletremoved_umapmd0.2.h5ad"
adata3.write(output_fpath)

#min dist=0.4
adata4 = adata.copy()
sc.external.pp.bbknn(adata4, batch_key="sample",n_pcs=30)  # running bbknn 1.6.0
sc.tl.umap(adata4,min_dist=0.4)
umapplt_path = "hvg" +sys.argv[1] + "PC30"+"_bbknn_batch_corrected" + "_doubletremoved_umapmd0.4.pdf"
sc.pl.umap(
    adata4,
    color=["sample_group","n_genes_by_counts", "pct_counts_mt","scDblFinder_class"],
    save=umapplt_path
)
output_fpath="RData/hvg"+ sys.argv[1] +"_PC30"+"_batchdoubletremoved_umapmd0.4.h5ad"
adata4.write(output_fpath)



print("Over!")
