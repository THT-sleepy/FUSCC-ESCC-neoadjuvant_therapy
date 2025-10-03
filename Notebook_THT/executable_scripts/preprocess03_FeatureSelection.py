#!/usr/bin/env python
# coding: utf-8

# ## 修改conda路径

# In[1]:


import os
print(os.environ['PATH'])


# In[2]:


conda_path = "/data1/liyi/zhaoyue/miniconda3/envs/sc_preprocess/bin:/data1/liyi/zhaoyue/miniconda3/condabin:/data1/liyi/softwares/STAR-2.7.11b/source:/data1/liyi/softwares/RSEM-master:/home/liyi301/sratoolkit.3.2.1-ubuntu64/bin:/home/liyi301/.local/bin:/home/liyi301/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/usr/libexec/gcc/x86_64-redhat-linux/8:/home/liyi301/.local/bin:/home/liyi301/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin"


# In[3]:


os.environ['PATH'] = conda_path + os.environ['PATH']


# In[4]:


print(os.environ['PATH'])


# ## 导入包

# In[6]:


import logging

import matplotlib.pyplot as plt
import numpy as np
np.float_ = np.float64
import anndata2ri
import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro
import scanpy as sc
import seaborn as sns

sc.settings.verbosity = 0
sc.settings.set_figure_params(
    dpi=80,
    facecolor="white",
    frameon=False,
)

rcb.logger.setLevel(logging.ERROR)
ro.pandas2ri.activate()
anndata2ri.activate()

get_ipython().run_line_magic('load_ext', 'rpy2.ipython')


# ## 挑选高变基因

# In[7]:


get_ipython().run_cell_magic('R', '', 'library(scry)\nlibrary(zellkonverter)\n')


# In[11]:


adata = sc.read("RData/normalization.h5ad")


# In[12]:


get_ipython().run_cell_magic('R', '', 'adata <- readH5AD("RData/normalization.h5ad",use_hdf5=TRUE)\nsce = devianceFeatureSelection(adata, assay="X")\n')


# In[13]:


binomial_deviance = ro.r("rowData(sce)$binomial_deviance").T


# In[14]:


idx = binomial_deviance.argsort()[-4000:]
mask = np.zeros(adata.var_names.shape, dtype=bool)
mask[idx] = True

adata.var["highly_deviant"] = mask
adata.var["binomial_deviance"] = binomial_deviance


# In[15]:


sc.pp.highly_variable_genes(adata, layer="log1p_norm")


# In[16]:


adata.write("RData/feature_selection.h5ad")

