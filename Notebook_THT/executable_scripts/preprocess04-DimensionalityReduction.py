#!/usr/bin/env python
# coding: utf-8

# ## 修改路径为conda路径

# In[1]:


import os


# In[2]:


print(os.environ['PATH'])


# In[3]:


conda_path = "/data1/liyi/zhaoyue/miniconda3/envs/sc_preprocess/bin:/data1/liyi/zhaoyue/miniconda3/condabin:/data1/liyi/softwares/STAR-2.7.11b/source:/data1/liyi/softwares/RSEM-master:/home/liyi301/sratoolkit.3.2.1-ubuntu64/bin:/home/liyi301/.local/bin:/home/liyi301/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/usr/libexec/gcc/x86_64-redhat-linux/8:/home/liyi301/.local/bin:/home/liyi301/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin"


# In[4]:


os.environ['PATH'] = conda_path + os.environ['PATH']


# In[5]:


print(os.environ['PATH'])


# ## 导入包及设置参数

# In[6]:


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


adata = sc.read(
    filename="RData/feature_selection.h5ad"
)


# In[10]:


adata.X = adata.layers["log1p_norm"]


# In[11]:


# setting highly variable as highly deviant to use scanpy 'use_highly_variable' argument in sc.pp.pca
adata.var["highly_variable"] = adata.var["highly_deviant"]
sc.pp.pca(adata, svd_solver="arpack", use_highly_variable=True)


# In[12]:


sc.pl.pca_scatter(adata, color="total_counts",save="pre04_pca.pdf")


# In[13]:


sc.tl.tsne(adata, use_rep="X_pca") ##默认是50个PC


# In[14]:


sc.pl.tsne(adata, color="total_counts",save="pre04_tsne.pdf")


# In[15]:


sc.pp.neighbors(adata) ##默认是50个PC,15个neighbors
sc.tl.umap(adata)


# In[16]:


sc.pl.umap(adata, color="total_counts",save="pre04_umap.pdf")


# In[17]:


sc.pl.umap(
    adata,
    color=["total_counts", "pct_counts_mt", "scDblFinder_score", "scDblFinder_class"],
    save="pre04_umap_annoQC.pdf"
)


# In[18]:


adata.write("RData/dimensionality_reduction.h5ad")
