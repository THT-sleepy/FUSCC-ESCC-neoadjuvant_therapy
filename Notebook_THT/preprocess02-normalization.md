## 修改路径为conda路径


```python
import os
print(os.environ['PATH'])
```

    /usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/snap/bin



```python
conda_path = "/home/data/t190513/miniconda3/envs/sc_preprocess/bin:/home/data/t190513/miniconda3/condabin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/snap/bin:/home/data/t190513/miniconda3/bin"
```


```python
os.environ['PATH'] = conda_path + os.environ['PATH']
```


```python
print(os.environ['PATH'])
```

    /home/data/t190513/miniconda3/envs/sc_preprocess/bin:/home/data/t190513/miniconda3/condabin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/snap/bin:/home/data/t190513/miniconda3/bin/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/snap/bin


## 导入包,设置scanpy参数


```python
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

%load_ext rpy2.ipython
```

## 先把前面的质控信息整理一下输出成excel


```python
data_root = "input/final_cellranger_data"
```


```python
samples = []
```


```python
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
```


```python
samples = samples[0:2]
```


```python
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
```

    已读取样本: XZM_preC1_2285332
    已读取样本: YXY_preC2_2282374



```python
qc_info
```




    [{'sample_id': 'XZM_preC1_2285332',
      'cellranger_filtered_cell_count': 11162,
      'contamination_fraction_bysoupx': 0.01,
      'cell_count_after_soupx': 11162,
      'cell_count_after_lowquanlity_filter': 5182,
      'double_cell_numbers': 389},
     {'sample_id': 'YXY_preC2_2282374',
      'cellranger_filtered_cell_count': 12216,
      'contamination_fraction_bysoupx': 0.019,
      'cell_count_after_soupx': 12216,
      'cell_count_after_lowquanlity_filter': 8061,
      'double_cell_numbers': 742}]



## 把质控好的矩阵合并到一起


```python
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
```

    已读取样本: XZM_preC1_2285332
    已读取样本: YXY_preC2_2282374



```python
adata = sc.concat(adatas,label="sample")
```

    /home/data/t190513/miniconda3/envs/sc_preprocess/lib/python3.12/site-packages/anndata/_core/anndata.py:1791: UserWarning: Observation names are not unique. To make them unique, call `.obs_names_make_unique`.
      utils.warn_names_duplicates("obs")



```python
adata
```




    AnnData object with n_obs × n_vars = 13243 × 18289
        obs: 'n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'pct_counts_in_top_20_genes', 'total_counts_mt', 'log1p_total_counts_mt', 'pct_counts_mt', 'total_counts_ribo', 'log1p_total_counts_ribo', 'pct_counts_ribo', 'total_counts_hb', 'log1p_total_counts_hb', 'pct_counts_hb', 'outlier', 'mt_outlier', 'scDblFinder_score', 'scDblFinder_class', 'sample'
        layers: 'counts', 'soupX_counts'




```python
## 确保细胞名唯一
adata.obs_names_make_unique()
```


```python
##进行合并后var属性会丢失我们需要自己添加回来
all_var = [x.var for x in adatas.values()]
all_var = pd.concat(all_var,join="outer")
all_var = all_var[~all_var.index.duplicated()]
adata.var = all_var.loc[adata.var_names]
```


```python
adata
```




    AnnData object with n_obs × n_vars = 13243 × 18289
        obs: 'n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'pct_counts_in_top_20_genes', 'total_counts_mt', 'log1p_total_counts_mt', 'pct_counts_mt', 'total_counts_ribo', 'log1p_total_counts_ribo', 'pct_counts_ribo', 'total_counts_hb', 'log1p_total_counts_hb', 'pct_counts_hb', 'outlier', 'mt_outlier', 'scDblFinder_score', 'scDblFinder_class', 'sample'
        var: 'gene_ids', 'feature_types', 'mt', 'ribo', 'hb', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts', 'n_cells'
        layers: 'counts', 'soupX_counts'



## 标化

### 1 Shifted logarithm


```python
scales_counts = sc.pp.normalize_total(adata, target_sum=None, inplace=False)
# log1p transform
adata.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)
```

### 2 scran


```python
from scipy.sparse import csr_matrix
```


```r
%%R
library(scran)
library(BiocParallel)
```


```python
#这一步差不多要1h,可以搞点其它事情
# Preliminary clustering for differentiated normalisation
adata_pp = adata.copy()
sc.pp.normalize_total(adata_pp)
sc.pp.log1p(adata_pp)
sc.pp.pca(adata_pp, n_comps=15)
sc.pp.neighbors(adata_pp)
sc.tl.leiden(adata_pp, key_added="groups")
```

    /tmp/ipykernel_1914047/2223024553.py:8: FutureWarning: In the future, the default backend for leiden will be igraph instead of leidenalg.
    
     To achieve the future defaults please pass: flavor="igraph" and n_iterations=2.  directed must also be False to work with igraph's implementation.
      sc.tl.leiden(adata_pp, key_added="groups")



```python
data_mat = adata_pp.X.T
# convert to CSC if possible. See https://github.com/MarioniLab/scran/issues/70
if issparse(data_mat):
    if data_mat.nnz > 2**31 - 1:
        data_mat = data_mat.tocoo()
    else:
        data_mat = data_mat.tocsc()
ro.globalenv["data_mat"] = data_mat
ro.globalenv["input_groups"] = adata_pp.obs["groups"]
```


```python
del adata_pp
```


```r
%%R -o size_factors

size_factors = sizeFactors(
    computeSumFactors(
        SingleCellExperiment(
            list(counts=data_mat)), 
            clusters = input_groups,
            min.mean = 0.1,
            BPPARAM = SerialParam()
    )
)
```


```python
adata.obs["size_factors"] = size_factors
scran = adata.X / adata.obs["size_factors"].values[:, None]
adata.layers["scran_normalization"] = csr_matrix(sc.pp.log1p(scran.tocsr()))
```

### Analytic Pearson residuals


```python
analytic_pearson = sc.experimental.pp.normalize_pearson_residuals(adata, inplace=False)
adata.layers["analytic_pearson_residuals"] = csr_matrix(analytic_pearson["X"])
```


```python
adata.write("RData/merge_normalization.h5ad")
```
