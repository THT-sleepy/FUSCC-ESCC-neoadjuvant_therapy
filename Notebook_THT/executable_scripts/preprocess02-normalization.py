#!/usr/bin/env python
# coding: utf-8

# ## 修改路径为conda路径

# In[1]:


import os
print(os.environ['PATH'])


# In[2]:


conda_path = "/data1/liyi/zhaoyue/miniconda3/envs/sc_preprocess/bin:/data1/liyi/zhaoyue/miniconda3/condabin:/data1/liyi/softwares/STAR-2.7.11b/source:/data1/liyi/softwares/RSEM-master:/home/liyi301/sratoolkit.3.2.1-ubuntu64/bin:/home/liyi301/.local/bin:/home/liyi301/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/usr/libexec/gcc/x86_64-redhat-linux/8:/home/liyi301/.local/bin:/home/liyi301/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin"


# In[3]:


os.environ['PATH'] = conda_path + os.environ['PATH']


# In[4]:


print(os.environ['PATH'])


# ## 导入包,设置scanpy参数

# In[5]:


import pandas as pd
import ast
import logging
import numpy as np
np.float_ = np.float64
import anndata2ri
import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro
import scanpy as sc
import seaborn as sns
from matplotlib import pyplot as plt
from scipy.sparse import issparse

sc.settings.verbosity = 0
sc.settings.set_figure_params(
    dpi=80,
    facecolor="white",
    # color_map="YlGnBu",
    frameon=False,
)

rcb.logger.setLevel(logging.ERROR)
ro.pandas2ri.activate()
anndata2ri.activate()

get_ipython().run_line_magic('load_ext', 'rpy2.ipython')


# ## 先把前面的质控信息整理一下输出成excel

# In[6]:


data_root = "input/final_cellranger_data"


# In[7]:


samples = []


# In[8]:


for sample_id in os.listdir(data_root):
    # 样本文件夹的完整路径
    sample_dir = os.path.join(data_root, sample_id)
    
    # 指定质控后h5ad文件输出路径
    quanlity_control_h5_path = os.path.join(sample_dir,"quanlity_control.h5ad")
    # 指定质控信息输出路径
    quanlity_control_info_path = os.path.join(sample_dir,"quanlity_control_info.txt")
    # 构建样本信息字典
    sample_info = {
                "sample_id": sample_id,
                "quanlity_control_h5_path" : quanlity_control_h5_path,
                "quanlity_control_info_path": quanlity_control_info_path,
            }
    # 添加到样本列表
    samples.append(sample_info)



# In[10]:


qc_info = []
for sample_dict in samples:
    try:
        sample_id = sample_dict['sample_id']
        file_path = sample_dict['quanlity_control_info_path']
        # 读取文件内容
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read().strip()
        # 转换为元组
        data_tuple = ast.literal_eval(content)
        #添加到qc_info列表中
        qc_info.append(data_tuple)
        print(f"已读取样本: {sample_id}")
    except Exception as e:
        print(f"读取样本 {sample_id} 失败: {str(e)}")


# In[11]:


qc_info


# ## 把质控好的矩阵合并到一起

# In[12]:


adatas = {}
for sample_dict in samples:
    try:
        sample_id = sample_dict['sample_id']
        file_path = sample_dict['quanlity_control_h5_path']
        # 读取h5ad文件
        sample_adata = sc.read_h5ad(file_path)
        adatas[sample_id] = sample_adata
        print(f"已读取样本: {sample_id}")
    except Exception as e:
        print(f"读取样本 {sample_id} 失败: {str(e)}")


# In[13]:


adata = sc.concat(adatas,label="sample")


# In[14]:


adata


# In[15]:


## 确保细胞名唯一
adata.obs_names_make_unique()


# In[16]:


##进行合并后var属性会丢失我们需要自己添加回来
all_var = [x.var for x in adatas.values()]
all_var = pd.concat(all_var,join="outer")
all_var = all_var[~all_var.index.duplicated()]
adata.var = all_var.loc[adata.var_names]


# In[17]:


adata


# ## 标化

# ###  1 Shifted logarithm

# In[18]:


scales_counts = sc.pp.normalize_total(adata, target_sum=None, inplace=False)
# log1p transform
adata.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)

adata.write("RData/normalization.h5ad")

