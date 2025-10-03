#!/usr/bin/env python
# coding: utf-8

# ## 修改路径为conda路径

# In[1]:


import os


# In[2]:


print(os.environ['PATH'])


# In[3]:


conda_path ="/data1/liyi/zhaoyue/miniconda3/envs/sc_preprocess/bin:/data1/liyi/zhaoyue/miniconda3/condabin:/data1/liyi/softwares/STAR-2.7.11b/source:/data1/liyi/softwares/RSEM-master:/home/liyi301/sratoolkit.3.2.1-ubuntu64/bin:/home/liyi301/.local/bin:/home/liyi301/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/usr/libexec/gcc/x86_64-redhat-linux/8:/home/liyi301/.local/bin:/home/liyi301/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin"


# In[4]:


os.environ['PATH'] = conda_path + os.environ['PATH']


# In[5]:


print(os.environ['PATH'])


# ## 导入库

# In[6]:


import numpy as np
np.float_ = np.float64
import scanpy as sc
import seaborn as sns
import pandas as pd
from scipy.stats import median_abs_deviation

sc.settings.verbosity = 0
sc.settings.set_figure_params(
    dpi=80,
    facecolor="white",
    frameon=False,
)


# In[7]:


data_root = "input/final_cellranger_data"


# In[57]:


samples = []


# In[58]:


for sample_id in os.listdir(data_root):
    # 样本文件夹的完整路径
    sample_dir = os.path.join(data_root, sample_id)
    
    # 构建filtered特征矩阵文件夹路径
    filtered_matrix_dir = os.path.join(sample_dir, "sample_filtered_feature_bc_matrix")
    # 构建raw 特征矩阵文件路径
    raw_matrix_path = os.path.join(sample_dir, "raw_feature_bc_matrix.h5")
    # 构建cluster信息文件路径
    cluster_path = os.path.join(sample_dir,"clusters.csv")
    # 指定质控后h5ad文件输出路径
    quanlity_control_h5_path = os.path.join(sample_dir,"quanlity_control.h5ad")
    # 指定质控信息输出路径
    quanlity_control_info_path = os.path.join(sample_dir,"quanlity_control_info.txt")
    # 构建样本信息字典
    sample_info = {
                "sample_id": sample_id,
                "filtered_matrix_dir_path": filtered_matrix_dir,
                "raw_matrix_h5_path": raw_matrix_path,
                "cluster_path": cluster_path,
                "quanlity_control_h5_path" : quanlity_control_h5_path,
                "quanlity_control_info_path": quanlity_control_info_path,
            }
    # 添加到样本列表
    samples.append(sample_info)


# In[59]:


sample = samples[n]


# In[60]:


sample


# In[10]:


samples


# In[11]:


print(f"找到 {len(samples)} 个有效样本")


# ### 统计下每个样本的细胞数,污染比例,soupx过滤后细胞数,去除低质量细胞后细胞数，双细胞数目过滤双细胞后剩余细胞数

# In[14]:


sample_qc_info = {
                "sample_id" : sample["sample_id"],
                "cellranger_filtered_cell_count": "",
                "contamination_fraction_bysoupx": "",
                "cell_count_after_soupx": "",
                "cell_count_after_lowquanlity_filter": "",
                "double_cell_numbers":"",
            }


# ### 读取filtered count matrix

# In[15]:


adata = sc.read_10x_mtx(
            sample["filtered_matrix_dir_path"],
            var_names='gene_symbols'
        )
# 确保基因名唯一
adata.var_names_make_unique()
sample_qc_info["cellranger_filtered_cell_count"] = adata.n_obs


# ## 质控步骤1 去掉cell-free RNA-soupX

# In[16]:


import logging


# In[17]:


import anndata2ri


# In[18]:


import rpy2.rinterface_lib.callbacks as rcb


# In[19]:


import rpy2.robjects as ro


# In[20]:


rcb.logger.setLevel(logging.ERROR)


# In[21]:


ro.pandas2ri.activate()


# In[22]:


anndata2ri.activate()


# In[23]:


get_ipython().run_line_magic('load_ext', 'rpy2.ipython')


# In[24]:


get_ipython().run_cell_magic('R', '', 'library(SoupX)\n')


# In[25]:


### SoupX需要使用的信息包括简单的cluster信息,cellranger的raw count以及filtered count
# 来自cellranger 的couster信息
soupx_groups = pd.read_csv(sample["cluster_path"], encoding='utf-8')
soupx_groups = soupx_groups.set_index("Barcode")["Cluster"]


#来自cellranger 的filtered count matrix
cells = adata.obs_names
genes = adata.var_names
data = adata.X.T ##因为SoupX要的矩阵是features x barcodes所以要转置一下

#来自cellranger 的raw count matrix
adata_raw = sc.read_10x_h5(
    sample["raw_matrix_h5_path"]
)
adata_raw.var_names_make_unique()
data_tod = adata_raw.X.T ###tod指的是table of droplets也就是cellranger的原始矩阵


# In[26]:


get_ipython().run_cell_magic('R', '-i data -i data_tod -i genes -i cells -i soupx_groups -o out -o contamination_fraction', '\n# specify row and column names of data\nrownames(data) = genes\ncolnames(data) = cells\n# ensure correct sparse format for table of counts and table of droplets\ndata <- as(data, "sparseMatrix")\ndata_tod <- as(data_tod, "sparseMatrix")\n\n# Generate SoupChannel Object for SoupX \nsc = SoupChannel(data_tod, data, calcSoupProfile = FALSE)\n\n# Add extra meta data to the SoupChannel object\nsoupProf = data.frame(row.names = rownames(data), est = rowSums(data)/sum(data), counts = rowSums(data))\nsc = setSoupProfile(sc, soupProf)\n\n# Set cluster information in SoupChannel\nsc = setClusters(sc, soupx_groups)\n\n# Estimate contamination fraction\nsc  = autoEstCont(sc, doPlot=FALSE)\ncontamination_fraction = sc$fit$rhoEst\n# Infer corrected table of counts and rount to integer\nout = adjustCounts(sc, roundToInt = TRUE)\n')


# In[27]:


adata.layers["counts"] = adata.X
adata.layers["soupX_counts"] = out.T
adata.X = adata.layers["soupX_counts"]


# In[28]:


sample_qc_info["contamination_fraction_bysoupx"] = contamination_fraction.item()
sample_qc_info["cell_count_after_soupx"] = adata.n_obs


# In[29]:


sample_qc_info


# ## 质控步骤2 去掉低质量细胞

# In[30]:


# mitochondrial genes
adata.var["mt"] = adata.var_names.str.startswith("MT-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes.
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")


# In[31]:


sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True
)
adata


# In[32]:


def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier


# In[33]:


adata.obs["outlier"] = (
    is_outlier(adata, "log1p_total_counts", 5)
    | is_outlier(adata, "log1p_n_genes_by_counts", 5)
    | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
)
adata.obs.outlier.value_counts()


# In[34]:


adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (
    adata.obs["pct_counts_mt"] > 8
)
adata.obs.mt_outlier.value_counts()


# In[35]:


adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()


# In[36]:


sample_qc_info["cell_count_after_lowquanlity_filter"] = adata.n_obs


# In[37]:


print(f"Total number of genes: {adata.n_vars}")

# Min 20 cells - filters out 0 count genes
sc.pp.filter_genes(adata, min_cells=20)
print(f"Number of genes after cell filter: {adata.n_vars}")


# ## 质控步骤3 标记双细胞

# In[38]:


get_ipython().run_cell_magic('R', '', 'library(Seurat)\nlibrary(scater)\nlibrary(scDblFinder)\nlibrary(BiocParallel)\n')


# In[39]:


data_mat = adata.X.T


# In[40]:


get_ipython().run_cell_magic('R', '-i data_mat -o doublet_score -o doublet_class', '\nset.seed(123)\nsce = scDblFinder(\nSingleCellExperiment(\n    list(counts=data_mat),\n) \n)\ndoublet_score = sce$scDblFinder.score\ndoublet_class = sce$scDblFinder.class\n')


# In[41]:


adata.obs["scDblFinder_score"] = doublet_score
adata.obs["scDblFinder_class"] = doublet_class
adata.obs.scDblFinder_class.value_counts()


# In[42]:


sample_qc_info["double_cell_numbers"] = int((adata.obs['scDblFinder_class'] == 'doublet').sum())
sample_qc_info


# In[61]:


adata.write(sample["quanlity_control_h5_path"])
# 保存到文本文件
with open(sample["quanlity_control_info_path"], "w") as f:
    f.write(str(sample_qc_info))  # 将元组转为字符串写入

