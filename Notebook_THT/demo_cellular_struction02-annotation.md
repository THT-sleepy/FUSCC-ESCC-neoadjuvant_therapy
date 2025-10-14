## 修改路径为conda路径


```python
import os
```


```python
print(os.environ['PATH'])
```

    /usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/snap/bin



```python
conda_path = "/home/data/t190513/miniconda3/envs/sc_annotation/bin:/home/data/t190513/miniconda3/condabin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/snap/bin:/home/data/t190513/miniconda3/bin"
```


```python
os.environ['PATH'] = conda_path + os.environ['PATH']
```


```python
print(os.environ['PATH'])
```

    /home/data/t190513/miniconda3/envs/sc_annotation/bin:/home/data/t190513/miniconda3/condabin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/snap/bin:/home/data/t190513/miniconda3/bin/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/snap/bin



```python
import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning)

from numba.core.errors import NumbaDeprecationWarning

warnings.simplefilter("ignore", category=NumbaDeprecationWarning)
```

## 加载需要的模块


```python
import urllib.request
from pathlib import Path

import celltypist
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scarches as sca
import seaborn as sns
from celltypist import models
from scipy.sparse import csr_matrix
```

    In order to use the mouse gastrulation seqFISH datsets, please install squidpy (see https://github.com/scverse/squidpy).
    In order to use sagenet models, please install pytorch geometric (see https://pytorch-geometric.readthedocs.io) and 
     captum (see https://github.com/pytorch/captum).
    mvTCR is not installed. To use mvTCR models, please install it first using "pip install mvtcr"
    multigrate is not installed. To use multigrate models, please install it first using "pip install multigrate".



```python
warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)
```

### 设置绘图参数


```python
sc.set_figure_params(figsize=(5, 5))
```


```python
adata = sc.read(
    filename="RData/mergedemo_clustered.h5ad"
)
```

## 第一轮大类注释


```python
## 主要是分出上皮，免疫，基质(纤维和内皮)三大类细胞
## 上皮的标志是EPCAM+(EPCAM)
## 免疫细胞的标志是CD45+(PTPRC)
## 纤维细胞的标志是CD10+(MME)
## 内皮细胞的标志是CD31+(PECAM1)
```

### 设置LEVEL1 marker gene


```python
marker_genes = {
    "Epithelial Cell": ["EPCAM", "KRT5", "KRT13", "KRT14", "KRT18", "SFN"],
    "Immune Cell": ["PTPRC"],
    "Fibro Cell" : ["MME", "COL1A2", "DCN", "FBLN1",
                    "COL1A1","COL3A1","FN1","COL6A1"],
    "Endo Cell" : ["PECAM1"]
   , 
}
```

### 仅保留我们的数据中表达的基因，这样绘图的时候不会有错


```python
marker_genes_in_data = {}
for ct, markers in marker_genes.items():
    markers_found = []
    for marker in markers:
        if marker in adata.var.index:
            markers_found.append(marker)
    marker_genes_in_data[ct] = markers_found
```


```python
marker_genes_in_data
```




    {'Epithelial Cell': ['KRT18'],
     'Immune Cell': ['PTPRC'],
     'Fibro Cell': ['FN1', 'COL6A1'],
     'Endo Cell': ['PECAM1']}



## 手动注释


```python
level1_cts = [
    "Epithelial Cell",
    "Immune Cell",
    "Fibro Cell",
    "Endo Cell",
] 
```


```python
marker_genes_in_data
```




    {'Epithelial Cell': ['KRT18'],
     'Immune Cell': ['PTPRC'],
     'Fibro Cell': ['FN1', 'COL6A1'],
     'Endo Cell': ['PECAM1']}




```python
for ct in level1_cts:
    print(f"{ct.upper()}:")  # print cell subtype name
    sc.pl.umap(
        adata,
        color=marker_genes_in_data[ct],
        vmin=0,
        vmax="p99",  # set vmax to the 99th percentile of the gene count instead of the maximum, to prevent outliers from making expression in other cells invisible. Note that this can cause problems for extremely lowly expressed genes.
        sort_order=False,  # do not plot highest expression on top, to not get a biased view of the mean expression among cells
        frameon=False,
        cmap="Reds",  # or choose another color map e.g. from here: https://matplotlib.org/stable/tutorials/colors/colormaps.html
    )
    print("\n\n\n")  # print white space for legibility
```

    EPITHELIAL CELL:



    
![png](output_23_1.png)
    


    
    
    
    
    IMMUNE CELL:



    
![png](output_23_3.png)
    


    
    
    
    
    FIBRO CELL:



    
![png](output_23_5.png)
    


    
    
    
    
    ENDO CELL:



    
![png](output_23_7.png)
    


    
    
    
    



```python
adata
```




    AnnData object with n_obs × n_vars = 28911 × 15674
        obs: 'n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'pct_counts_in_top_20_genes', 'total_counts_mt', 'log1p_total_counts_mt', 'pct_counts_mt', 'total_counts_ribo', 'log1p_total_counts_ribo', 'pct_counts_ribo', 'total_counts_hb', 'log1p_total_counts_hb', 'pct_counts_hb', 'outlier', 'mt_outlier', 'scDblFinder_score', 'scDblFinder_class', 'sample', 'leiden_res0_1', 'leiden_res0_25', 'leiden_res0_5', 'leiden_res1'
        var: 'gene_ids', 'feature_types', 'mt', 'ribo', 'hb', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts', 'n_cells', 'highly_variable', 'means', 'dispersions', 'dispersions_norm', 'highly_variable_nbatches', 'highly_variable_intersection'
        uns: 'hvg', 'leiden_res0_1', 'leiden_res0_1_colors', 'leiden_res0_25', 'leiden_res0_25_colors', 'leiden_res0_5', 'leiden_res0_5_colors', 'leiden_res1', 'leiden_res1_colors', 'neighbors', 'pca', 'sample_colors', 'scDblFinder_class_colors', 'tsne', 'umap'
        obsm: 'X_pca', 'X_pca_harmony', 'X_tsne', 'X_umap'
        varm: 'PCs'
        layers: 'counts', 'log1p_norm', 'soupX_counts'
        obsp: 'connectivities', 'distances'




```python
?sc.tl.umap
```


    [0;31mSignature:[0m
    [0msc[0m[0;34m.[0m[0mtl[0m[0;34m.[0m[0mumap[0m[0;34m([0m[0;34m[0m
    [0;34m[0m    [0madata[0m[0;34m:[0m [0;34m'AnnData'[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0;34m*[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0mmin_dist[0m[0;34m:[0m [0;34m'float'[0m [0;34m=[0m [0;36m0.5[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0mspread[0m[0;34m:[0m [0;34m'float'[0m [0;34m=[0m [0;36m1.0[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0mn_components[0m[0;34m:[0m [0;34m'int'[0m [0;34m=[0m [0;36m2[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0mmaxiter[0m[0;34m:[0m [0;34m'int | None'[0m [0;34m=[0m [0;32mNone[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0malpha[0m[0;34m:[0m [0;34m'float'[0m [0;34m=[0m [0;36m1.0[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0mgamma[0m[0;34m:[0m [0;34m'float'[0m [0;34m=[0m [0;36m1.0[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0mnegative_sample_rate[0m[0;34m:[0m [0;34m'int'[0m [0;34m=[0m [0;36m5[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0minit_pos[0m[0;34m:[0m [0;34m'_InitPos | np.ndarray | None'[0m [0;34m=[0m [0;34m'spectral'[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0mrandom_state[0m[0;34m:[0m [0;34m'_LegacyRandom'[0m [0;34m=[0m [0;36m0[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0ma[0m[0;34m:[0m [0;34m'float | None'[0m [0;34m=[0m [0;32mNone[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0mb[0m[0;34m:[0m [0;34m'float | None'[0m [0;34m=[0m [0;32mNone[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0mmethod[0m[0;34m:[0m [0;34m"Literal['umap', 'rapids']"[0m [0;34m=[0m [0;34m'umap'[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0mkey_added[0m[0;34m:[0m [0;34m'str | None'[0m [0;34m=[0m [0;32mNone[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0mneighbors_key[0m[0;34m:[0m [0;34m'str'[0m [0;34m=[0m [0;34m'neighbors'[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0mcopy[0m[0;34m:[0m [0;34m'bool'[0m [0;34m=[0m [0;32mFalse[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m[0;34m)[0m [0;34m->[0m [0;34m'AnnData | None'[0m[0;34m[0m[0;34m[0m[0m
    [0;31mDocstring:[0m
    Embed the neighborhood graph using UMAP :cite:p:`McInnes2018`.
    
    UMAP (Uniform Manifold Approximation and Projection) is a manifold learning
    technique suitable for visualizing high-dimensional data. Besides tending to
    be faster than tSNE, it optimizes the embedding such that it best reflects
    the topology of the data, which we represent throughout Scanpy using a
    neighborhood graph. tSNE, by contrast, optimizes the distribution of
    nearest-neighbor distances in the embedding such that these best match the
    distribution of distances in the high-dimensional space.
    We use the implementation of umap-learn_ :cite:p:`McInnes2018`.
    For a few comparisons of UMAP with tSNE, see :cite:t:`Becht2018`.
    
    .. _umap-learn: https://github.com/lmcinnes/umap
    
    Parameters
    ----------
    adata : 'AnnData'
        Annotated data matrix.
    min_dist : 'float', optional (default: 0.5)
        The effective minimum distance between embedded points. Smaller values
        will result in a more clustered/clumped embedding where nearby points on
        the manifold are drawn closer together, while larger values will result
        on a more even dispersal of points. The value should be set relative to
        the ``spread`` value, which determines the scale at which embedded
        points will be spread out. The default of in the `umap-learn` package is
        0.1.
    spread : 'float', optional (default: 1.0)
        The effective scale of embedded points. In combination with `min_dist`
        this determines how clustered/clumped the embedded points are.
    n_components : 'int', optional (default: 2)
        The number of dimensions of the embedding.
    maxiter : 'int | None', optional (default: None)
        The number of iterations (epochs) of the optimization. Called `n_epochs`
        in the original UMAP.
    alpha : 'float', optional (default: 1.0)
        The initial learning rate for the embedding optimization.
    gamma : 'float', optional (default: 1.0)
        Weighting applied to negative samples in low dimensional embedding
        optimization. Values higher than one will result in greater weight
        being given to negative samples.
    negative_sample_rate : 'int', optional (default: 5)
        The number of negative edge/1-simplex samples to use per positive
        edge/1-simplex sample in optimizing the low dimensional embedding.
    init_pos : '_InitPos | np.ndarray | None', optional (default: 'spectral')
        How to initialize the low dimensional embedding. Called `init` in the
        original UMAP. Options are:
    
        * Any key for `adata.obsm`.
        * 'paga': positions from :func:`~scanpy.pl.paga`.
        * 'spectral': use a spectral embedding of the graph.
        * 'random': assign initial embedding positions at random.
        * A numpy array of initial embedding positions.
    random_state : '_LegacyRandom', optional (default: 0)
        If `int`, `random_state` is the seed used by the random number generator;
        If `RandomState` or `Generator`, `random_state` is the random number generator;
        If `None`, the random number generator is the `RandomState` instance used
        by `np.random`.
    a : 'float | None', optional (default: None)
        More specific parameters controlling the embedding. If `None` these
        values are set automatically as determined by `min_dist` and
        `spread`.
    b : 'float | None', optional (default: None)
        More specific parameters controlling the embedding. If `None` these
        values are set automatically as determined by `min_dist` and
        `spread`.
    method : "Literal['umap', 'rapids']", optional (default: 'umap')
        Chosen implementation.
    
        ``'umap'``
            Umap’s simplical set embedding.
        ``'rapids'``
            GPU accelerated implementation.
    
            .. deprecated:: 1.10.0
                Use :func:`rapids_singlecell.tl.umap` instead.
    key_added : 'str | None', optional (default: None)
        If not specified, the embedding is stored as
        :attr:`~anndata.AnnData.obsm`\ `['X_umap']` and the the parameters in
        :attr:`~anndata.AnnData.uns`\ `['umap']`.
        If specified, the embedding is stored as
        :attr:`~anndata.AnnData.obsm`\ ``[key_added]`` and the the parameters in
        :attr:`~anndata.AnnData.uns`\ ``[key_added]``.
    neighbors_key : 'str', optional (default: 'neighbors')
        Umap looks in
        :attr:`~anndata.AnnData.uns`\ ``[neighbors_key]`` for neighbors settings and
        :attr:`~anndata.AnnData.obsp`\ ``[.uns[neighbors_key]['connectivities_key']]`` for connectivities.
    copy : 'bool', optional (default: False)
        Return a copy instead of writing to adata.
    
    Returns
    -------
    Returns `None` if `copy=False`, else returns an `AnnData` object. Sets the following fields:
    
    `adata.obsm['X_umap' | key_added]` : :class:`numpy.ndarray` (dtype `float`)
        UMAP coordinates of data.
    `adata.uns['umap' | key_added]` : :class:`dict`
        UMAP parameters.
    [0;31mFile:[0m      ~/miniconda3/envs/sc_annotation/lib/python3.10/site-packages/scanpy/tools/_umap.py
    [0;31mType:[0m      function



```python
sc.pl.umap(adata, color=["leiden_res0_1","scDblFinder_class"],legend_loc="on data")
```


    
![png](output_26_0.png)
    


### 用点图来看一下


```python
markers = {
    ct: [m for m in ct_markers if m in adata.var.index]
    for ct, ct_markers in marker_genes.items()
    if ct in level1_cts
}
```


```python
sc.pl.dotplot(
    adata,
    groupby="leiden_res0_1",
    var_names=markers,
    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
)
```


    
![png](output_29_0.png)
    



```python
cl_annotation = {
    "0": "Epithelial Cell",
    "1": "Immune Cell",
    "2": "doublet (PTPRC+PECAM1)",
    "3": "Immune Cell",
    "4": "Fibro Cell",
    "5": "Immune Cell",
    "6": "Endo Cell",
    "7": "doublet (PTPRC+PECAM1)",
}
```


```python
adata.obs["manual_celltype_level1"] = adata.obs.leiden_res0_1.map(cl_annotation)
```


```python
sc.pl.umap(adata, color=["manual_celltype_level1"])
```


    
![png](output_32_0.png)
    



```python
sc.tl.rank_genes_groups(
    adata, groupby="leiden_res0_1", method="wilcoxon", key_added="dea_leiden_res0_1"
)
```

### 查看每个cluster特异性表达的top5基因


```python
sc.pl.rank_genes_groups_dotplot(
    adata, groupby="leiden_res0_1", standard_scale="var", n_genes=5, key="dea_leiden_res0_1"
)
```


    
![png](output_35_0.png)
    


### 有些基因在多个cluster高表达，过滤一下这些基因


```python
sc.tl.filter_rank_genes_groups(
    adata,
    min_in_group_fraction=0.2,
    max_out_group_fraction=0.2,
    key="dea_leiden_res0_1",
    key_added="dea_leiden_res0_1_filtered",
)
```


```python
sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby="leiden_res0_1",
    standard_scale="var",
    n_genes=5,
    key="dea_leiden_res0_1_filtered",
)
```


    
![png](output_38_0.png)
    


## 可以再详细看每个cluster


```python
adata.layers
```




    Layers with keys: counts, log1p_norm, soupX_counts




```python
adata
```




    AnnData object with n_obs × n_vars = 28911 × 15674
        obs: 'n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'pct_counts_in_top_20_genes', 'total_counts_mt', 'log1p_total_counts_mt', 'pct_counts_mt', 'total_counts_ribo', 'log1p_total_counts_ribo', 'pct_counts_ribo', 'total_counts_hb', 'log1p_total_counts_hb', 'pct_counts_hb', 'outlier', 'mt_outlier', 'scDblFinder_score', 'scDblFinder_class', 'sample', 'leiden_res0_1', 'leiden_res0_25', 'leiden_res0_5', 'leiden_res1', 'manual_celltype_level1', 'celltypist_cell_label_coarse', 'celltypist_conf_score_coarse', 'celltypist_cell_label_fine', 'celltypist_conf_score_fine'
        var: 'gene_ids', 'feature_types', 'mt', 'ribo', 'hb', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts', 'n_cells', 'highly_variable', 'means', 'dispersions', 'dispersions_norm', 'highly_variable_nbatches', 'highly_variable_intersection'
        uns: 'hvg', 'leiden_res0_1', 'leiden_res0_1_colors', 'leiden_res0_25', 'leiden_res0_25_colors', 'leiden_res0_5', 'leiden_res0_5_colors', 'leiden_res1', 'leiden_res1_colors', 'neighbors', 'pca', 'sample_colors', 'scDblFinder_class_colors', 'tsne', 'umap', 'manual_celltype_level1_colors', 'dea_leiden_res0_1', 'dendrogram_leiden_res0_1', 'leiden_res0_1_filtered', 'dea_leiden_res0_1_filtered', 'celltypist_cell_label_coarse_colors', 'celltypist_cell_label_fine_colors'
        obsm: 'X_pca', 'X_pca_harmony', 'X_tsne', 'X_umap'
        varm: 'PCs'
        layers: 'counts', 'log1p_norm', 'soupX_counts'
        obsp: 'connectivities', 'distances'




```python
if 'leiden_res0_1' in adata.uns:
    names_data = adata.uns['leiden_res0_1'].get('names', None)
    print("当前 names 字段的数据类型:", type(names_data))
    print("当前 names 字段的前 10 个值（如果是数组/列表）:", names_data[:10] if hasattr(names_data, '__getitem__') else names_data)
    print("是否包含非字符串值:", any(not isinstance(x, str) for x in names_data) if isinstance(names_data, (list, np.ndarray)) else "非列表/数组类型")
else:
    print("uns 中没有 leiden_res0_1 字段，请检查聚类结果的存储路径")
```

    当前 names 字段的数据类型: <class 'numpy.recarray'>
    当前 names 字段的前 10 个值（如果是数组/列表）: [(nan, 'NKG7', 'AIF1', 'IL7R', 'CALD1', nan, nan, 'LST1')
     (nan, nan, 'SPI1', 'LTB', 'SPARC', nan, 'PECAM1', 'AIF1')
     ('PERP', nan, 'FCER1G', 'CD3E', 'COL6A2', 'HLA-DRA', 'SPARC', nan)
     ('S100A2', nan, 'TYROBP', nan, nan, 'CD79A', 'ENG', 'LILRB2')
     (nan, 'CST7', nan, nan, 'FN1', nan, nan, 'SERPINA1')
     (nan, 'CCL5', nan, 'TCF7', 'IGFBP4', nan, 'A2M', nan)
     (nan, nan, 'SERPINA1', nan, nan, nan, 'CRIP2', 'FCER1G')
     ('JUP', 'GZMA', nan, nan, nan, 'HLA-DQB1', nan, 'PELATON')
     ('CD9', 'HCST', 'CD14', 'SELL', 'SERPING1', 'HLA-DQA1', 'RNASE1', 'MS4A7')
     ('CTTN', 'PRF1', 'CTSS', 'CD69', 'MYL9', 'LTB', nan, 'LRRC25')]
    是否包含非字符串值: True



```python
import numpy as np

# 1. 提取 recarray 并查看其结构化 dtype（确认字段名，通常是 '0','1'...'7'）
leiden_group = adata.uns['leiden_res0_1']
names_recarray = leiden_group['names']
print("recarray 原始 dtype:", names_recarray.dtype)  # 会显示类似 [('0','O'), ('1','O'), ..., ('7','O')]

# 2. 提取所有字段名（如 '0','1','2'...'7'）
field_names = names_recarray.dtype.names  # 自动获取所有字段名，避免手动写死
print("recarray 字段名:", field_names)  # 输出: ('0', '1', '2', '3', '4', '5', '6', '7')

# 3. 逐个字段处理：替换 nan 为字符串，再组合成普通二维数组
fixed_columns = []
for field in field_names:
    # 提取当前字段的列数据（每个字段对应一列）
    col_data = names_recarray[field]  # 每列是一个 numpy 数组（元素可能是 object 类型，含 nan）
    
    # 处理该列的 nan：将 nan 转为空字符串 ""，其他元素转为字符串
    # 注意：object 类型数组中的 nan 需用 np.isnan 先判断（但需避免非数值类型报错）
    def safe_replace_nan(x):
        if isinstance(x, float) and np.isnan(x):
            return ""  #  nan 替换为空字符串
        else:
            return str(x)  # 其他元素转为字符串
    
    # 对列中每个元素应用处理函数
    fixed_col = np.array([safe_replace_nan(x) for x in col_data], dtype=str)
    fixed_columns.append(fixed_col)

# 4. 将处理后的列组合成普通二维字符串数组（形状：[样本数, 8]）
# 转置是因为 fixed_columns 是按列存储，需转为行×列的二维数组
names_fixed = np.column_stack(fixed_columns).astype(str)
print("修复后数组 dtype:", names_fixed.dtype)  # 应显示 <Uxx（字符串 dtype）
print("修复后数组形状:", names_fixed.shape)    # 应显示 (n_samples, 8)

# 5. 覆盖原有的 recarray 字段，用普通字符串数组替代
leiden_group['names'] = names_fixed

# 6. 验证修复结果（可选）
print("\n修复后前 3 行数据:")
print(names_fixed[:3])
print("是否存在非字符串值:", any(not isinstance(x, str) for row in names_fixed for x in row))  # 应输出 False
```

    recarray 原始 dtype: (numpy.record, [('0', 'O'), ('1', 'O'), ('2', 'O'), ('3', 'O'), ('4', 'O'), ('5', 'O'), ('6', 'O'), ('7', 'O')])
    recarray 字段名: ('0', '1', '2', '3', '4', '5', '6', '7')
    修复后数组 dtype: <U15
    修复后数组形状: (15674, 8)
    
    修复后前 3 行数据:
    [['' 'NKG7' 'AIF1' 'IL7R' 'CALD1' '' '' 'LST1']
     ['' '' 'SPI1' 'LTB' 'SPARC' '' 'PECAM1' 'AIF1']
     ['PERP' '' 'FCER1G' 'CD3E' 'COL6A2' 'HLA-DRA' 'SPARC' '']]
    是否存在非字符串值: False



```python
import numpy as np

# 1. 提取 recarray 并查看其结构化 dtype（确认字段名，通常是 '0','1'...'7'）
leiden_group = adata.uns['leiden_res0_1_filtered']
names_recarray = leiden_group['names']
print("recarray 原始 dtype:", names_recarray.dtype)  # 会显示类似 [('0','O'), ('1','O'), ..., ('7','O')]

# 2. 提取所有字段名（如 '0','1','2'...'7'）
field_names = names_recarray.dtype.names  # 自动获取所有字段名，避免手动写死
print("recarray 字段名:", field_names)  # 输出: ('0', '1', '2', '3', '4', '5', '6', '7')

# 3. 逐个字段处理：替换 nan 为字符串，再组合成普通二维数组
fixed_columns = []
for field in field_names:
    # 提取当前字段的列数据（每个字段对应一列）
    col_data = names_recarray[field]  # 每列是一个 numpy 数组（元素可能是 object 类型，含 nan）
    
    # 处理该列的 nan：将 nan 转为空字符串 ""，其他元素转为字符串
    # 注意：object 类型数组中的 nan 需用 np.isnan 先判断（但需避免非数值类型报错）
    def safe_replace_nan(x):
        if isinstance(x, float) and np.isnan(x):
            return ""  #  nan 替换为空字符串
        else:
            return str(x)  # 其他元素转为字符串
    
    # 对列中每个元素应用处理函数
    fixed_col = np.array([safe_replace_nan(x) for x in col_data], dtype=str)
    fixed_columns.append(fixed_col)

# 4. 将处理后的列组合成普通二维字符串数组（形状：[样本数, 8]）
# 转置是因为 fixed_columns 是按列存储，需转为行×列的二维数组
names_fixed = np.column_stack(fixed_columns).astype(str)
print("修复后数组 dtype:", names_fixed.dtype)  # 应显示 <Uxx（字符串 dtype）
print("修复后数组形状:", names_fixed.shape)    # 应显示 (n_samples, 8)

# 5. 覆盖原有的 recarray 字段，用普通字符串数组替代
leiden_group['names'] = names_fixed

# 6. 验证修复结果（可选）
print("\n修复后前 3 行数据:")
print(names_fixed[:3])
print("是否存在非字符串值:", any(not isinstance(x, str) for row in names_fixed for x in row))  # 应输出 False
```

    recarray 原始 dtype: (numpy.record, [('0', 'O'), ('1', 'O'), ('2', 'O'), ('3', 'O'), ('4', 'O'), ('5', 'O'), ('6', 'O'), ('7', 'O')])
    recarray 字段名: ('0', '1', '2', '3', '4', '5', '6', '7')
    修复后数组 dtype: <U15
    修复后数组形状: (15674, 8)
    
    修复后前 3 行数据:
    [['' 'NKG7' 'AIF1' 'IL7R' 'CALD1' '' '' 'LST1']
     ['' '' 'SPI1' 'LTB' 'SPARC' '' 'PECAM1' 'AIF1']
     ['PERP' '' 'FCER1G' 'CD3E' 'COL6A2' 'HLA-DRA' 'SPARC' '']]
    是否存在非字符串值: False



```python
import numpy as np

# 1. 提取 recarray 并查看其结构化 dtype（确认字段名，通常是 '0','1'...'7'）
leiden_group = adata.uns['dea_leiden_res0_1_filtered']
names_recarray = leiden_group['names']
print("recarray 原始 dtype:", names_recarray.dtype)  # 会显示类似 [('0','O'), ('1','O'), ..., ('7','O')]

# 2. 提取所有字段名（如 '0','1','2'...'7'）
field_names = names_recarray.dtype.names  # 自动获取所有字段名，避免手动写死
print("recarray 字段名:", field_names)  # 输出: ('0', '1', '2', '3', '4', '5', '6', '7')

# 3. 逐个字段处理：替换 nan 为字符串，再组合成普通二维数组
fixed_columns = []
for field in field_names:
    # 提取当前字段的列数据（每个字段对应一列）
    col_data = names_recarray[field]  # 每列是一个 numpy 数组（元素可能是 object 类型，含 nan）
    
    # 处理该列的 nan：将 nan 转为空字符串 ""，其他元素转为字符串
    # 注意：object 类型数组中的 nan 需用 np.isnan 先判断（但需避免非数值类型报错）
    def safe_replace_nan(x):
        if isinstance(x, float) and np.isnan(x):
            return ""  #  nan 替换为空字符串
        else:
            return str(x)  # 其他元素转为字符串
    
    # 对列中每个元素应用处理函数
    fixed_col = np.array([safe_replace_nan(x) for x in col_data], dtype=str)
    fixed_columns.append(fixed_col)

# 4. 将处理后的列组合成普通二维字符串数组（形状：[样本数, 8]）
# 转置是因为 fixed_columns 是按列存储，需转为行×列的二维数组
names_fixed = np.column_stack(fixed_columns).astype(str)
print("修复后数组 dtype:", names_fixed.dtype)  # 应显示 <Uxx（字符串 dtype）
print("修复后数组形状:", names_fixed.shape)    # 应显示 (n_samples, 8)

# 5. 覆盖原有的 recarray 字段，用普通字符串数组替代
leiden_group['names'] = names_fixed

# 6. 验证修复结果（可选）
print("\n修复后前 3 行数据:")
print(names_fixed[:3])
print("是否存在非字符串值:", any(not isinstance(x, str) for row in names_fixed for x in row))  # 应输出 False
```

    recarray 原始 dtype: (numpy.record, [('0', 'O'), ('1', 'O'), ('2', 'O'), ('3', 'O'), ('4', 'O'), ('5', 'O'), ('6', 'O'), ('7', 'O')])
    recarray 字段名: ('0', '1', '2', '3', '4', '5', '6', '7')
    修复后数组 dtype: <U15
    修复后数组形状: (15674, 8)
    
    修复后前 3 行数据:
    [['' 'NKG7' 'AIF1' 'IL7R' 'CALD1' '' '' 'LST1']
     ['' '' 'SPI1' 'LTB' 'SPARC' '' 'PECAM1' 'AIF1']
     ['PERP' '' 'FCER1G' 'CD3E' 'COL6A2' 'HLA-DRA' 'SPARC' '']]
    是否存在非字符串值: False



```python
adata.write("RData/mergedemo_annotation.h5ad")
```

## 自动注释-celltypist

### celltypist需要标化每个细胞到10000counts


```python
adata_celltypist = adata.copy()  # make a copy of our adata
adata_celltypist.X = adata.layers["soupX_counts"]  # set adata.X to raw counts
sc.pp.normalize_total(
    adata_celltypist, target_sum=10**4
)  # normalize to 10,000 counts per cell
sc.pp.log1p(adata_celltypist)  # log-transform
# make .X dense instead of sparse, for compatibility with celltypist:
adata_celltypist.X = adata_celltypist.X.toarray()
```


```python
models.download_models(
    force_update=True, model=["Immune_All_Low.pkl", "Immune_All_High.pkl"]
)
```

    📜 Retrieving model list from server https://celltypist.cog.sanger.ac.uk/models/models.json
    📚 Total models in list: 56
    📂 Storing models in /home/data/t190513/.celltypist/data/models
    💾 Total models to download: 2
    💾 Downloading model [1/2]: Immune_All_Low.pkl
    💾 Downloading model [2/2]: Immune_All_High.pkl



```python
model_low = models.Model.load(model="Immune_All_Low.pkl")
model_high = models.Model.load(model="Immune_All_High.pkl")
```


```python
model_high.cell_types
```




    array(['B cells', 'B-cell lineage', 'Cycling cells', 'DC', 'DC precursor',
           'Double-negative thymocytes', 'Double-positive thymocytes', 'ETP',
           'Early MK', 'Endothelial cells', 'Epithelial cells',
           'Erythrocytes', 'Erythroid', 'Fibroblasts', 'Granulocytes',
           'HSC/MPP', 'ILC', 'ILC precursor', 'MNP', 'Macrophages',
           'Mast cells', 'Megakaryocyte precursor',
           'Megakaryocytes/platelets', 'Mono-mac', 'Monocyte precursor',
           'Monocytes', 'Myelocytes', 'Plasma cells', 'Promyelocytes',
           'T cells', 'pDC', 'pDC precursor'], dtype=object)




```python
model_low.cell_types
```




    array(['Age-associated B cells', 'Alveolar macrophages', 'B cells',
           'CD16+ NK cells', 'CD16- NK cells', 'CD8a/a', 'CD8a/b(entry)',
           'CMP', 'CRTAM+ gamma-delta T cells', 'Classical monocytes',
           'Cycling B cells', 'Cycling DCs', 'Cycling NK cells',
           'Cycling T cells', 'Cycling gamma-delta T cells',
           'Cycling monocytes', 'DC', 'DC precursor', 'DC1', 'DC2', 'DC3',
           'Double-negative thymocytes', 'Double-positive thymocytes', 'ELP',
           'ETP', 'Early MK', 'Early erythroid', 'Early lymphoid/T lymphoid',
           'Endothelial cells', 'Epithelial cells', 'Erythrocytes',
           'Erythrophagocytic macrophages', 'Fibroblasts',
           'Follicular B cells', 'Follicular helper T cells', 'GMP',
           'Germinal center B cells', 'Granulocytes', 'HSC/MPP',
           'Hofbauer cells', 'ILC', 'ILC precursor', 'ILC1', 'ILC2', 'ILC3',
           'Intermediate macrophages', 'Intestinal macrophages',
           'Kidney-resident macrophages', 'Kupffer cells',
           'Large pre-B cells', 'Late erythroid', 'MAIT cells', 'MEMP', 'MNP',
           'Macrophages', 'Mast cells', 'Megakaryocyte precursor',
           'Megakaryocyte-erythroid-mast cell progenitor',
           'Megakaryocytes/platelets', 'Memory B cells',
           'Memory CD4+ cytotoxic T cells', 'Mid erythroid', 'Migratory DCs',
           'Mono-mac', 'Monocyte precursor', 'Monocytes', 'Myelocytes',
           'NK cells', 'NKT cells', 'Naive B cells',
           'Neutrophil-myeloid progenitor', 'Neutrophils',
           'Non-classical monocytes', 'Plasma cells', 'Plasmablasts',
           'Pre-pro-B cells', 'Pro-B cells',
           'Proliferative germinal center B cells', 'Promyelocytes',
           'Regulatory T cells', 'Small pre-B cells', 'T(agonist)',
           'Tcm/Naive cytotoxic T cells', 'Tcm/Naive helper T cells',
           'Tem/Effector helper T cells', 'Tem/Effector helper T cells PD1+',
           'Tem/Temra cytotoxic T cells', 'Tem/Trm cytotoxic T cells',
           'Transitional B cells', 'Transitional DC', 'Transitional NK',
           'Treg(diff)', 'Trm cytotoxic T cells', 'Type 1 helper T cells',
           'Type 17 helper T cells', 'gamma-delta T cells', 'pDC',
           'pDC precursor'], dtype=object)




```python
predictions_high = celltypist.annotate(
    adata_celltypist, model=model_high, majority_voting=True
)
```

    🔬 Input data has 28911 cells and 15674 genes
    🔗 Matching reference genes in the model
    🧬 4450 features used for prediction
    ⚖️ Scaling input data
    🖋️ Predicting labels
    ✅ Prediction done!
    👀 Detected a neighborhood graph in the input object, will run over-clustering on the basis of it
    ⛓️ Over-clustering input data with resolution set to 15
    🗳️ Majority voting the predictions
    ✅ Majority voting done!



```python
predictions_high_adata = predictions_high.to_adata()
```


```python
adata.obs["celltypist_cell_label_coarse"] = predictions_high_adata.obs.loc[
    adata.obs.index, "majority_voting"
]
adata.obs["celltypist_conf_score_coarse"] = predictions_high_adata.obs.loc[
    adata.obs.index, "conf_score"
]
```


```python
predictions_low = celltypist.annotate(
    adata_celltypist, model=model_low, majority_voting=True
)
```

    🔬 Input data has 28911 cells and 15674 genes
    🔗 Matching reference genes in the model
    🧬 4450 features used for prediction
    ⚖️ Scaling input data
    🖋️ Predicting labels
    ✅ Prediction done!
    👀 Detected a neighborhood graph in the input object, will run over-clustering on the basis of it
    ⛓️ Over-clustering input data with resolution set to 15
    🗳️ Majority voting the predictions
    ✅ Majority voting done!



```python
predictions_low_adata = predictions_low.to_adata()
```


```python
adata.obs["celltypist_cell_label_fine"] = predictions_low_adata.obs.loc[
    adata.obs.index, "majority_voting"
]
adata.obs["celltypist_conf_score_fine"] = predictions_low_adata.obs.loc[
    adata.obs.index, "conf_score"
]
```


```python
sc.pl.umap(
    adata,
    color=["celltypist_cell_label_coarse", "celltypist_conf_score_coarse"],
    frameon=False,
    sort_order=False,
    wspace=1,
)
```


    
![png](output_60_0.png)
    



```python
sc.pl.umap(
    adata,
    color=["celltypist_cell_label_fine", "celltypist_conf_score_fine","manual_celltype_level1"],
    frameon=False,
    sort_order=False,
    wspace=1,
)
```


    
![png](output_61_0.png)
    



```python

```


```python
sc.pl.dendrogram(adata, groupby="celltypist_cell_label_fine")
```

    WARNING: dendrogram data not found (using key=dendrogram_celltypist_cell_label_fine). Running `sc.tl.dendrogram` with default parameters. For fine tuning it is recommended to run `sc.tl.dendrogram` independently.



    
![png](output_63_1.png)
    





    <Axes: >




```python
pd.crosstab(adata.obs.leiden_2, adata.obs.celltypist_cell_label_fine).loc[
    "12", :
].sort_values(ascending=False)
```




    celltypist_cell_label_fine
    Tcm/Naive helper T cells       313
    Tem/Trm cytotoxic T cells       15
    Tem/Temra cytotoxic T cells      7
    MAIT cells                       1
    Alveolar macrophages             0
    CD16+ NK cells                   0
    Classical monocytes              0
    DC2                              0
    HSC/MPP                          0
    Memory B cells                   0
    Mid erythroid                    0
    NK cells                         0
    Naive B cells                    0
    Non-classical monocytes          0
    Pro-B cells                      0
    Small pre-B cells                0
    Name: 12, dtype: int64



## 迁移注释-scArches

### scArches需要没有标化过的counts


```python
adata_to_map = adata.copy()
for layer in list(adata_to_map.layers.keys()):
    if layer != "counts":
        del adata_to_map.layers[layer]
adata_to_map.X = adata_to_map.layers["counts"]
```


```python
reference_model_features = pd.read_csv(
    "https://figshare.com/ndownloader/files/41436645", index_col=0
)
```

### 需要把我们的基因顺序弄成和reference数据集一致


```python
adata_to_map.var["gene_names"] = adata_to_map.var.index
adata_to_map.var.set_index("gene_ids", inplace=True)
```


```python
reference_model_features["gene_names"] = reference_model_features.index
reference_model_features.set_index("gene_ids", inplace=True)
```

### 看一下是否所有需要的基因我们都有了


```python
print("Total number of genes needed for mapping:", reference_model_features.shape[0])
```

    Total number of genes needed for mapping: 4000



```python
print(
    "Number of genes found in query dataset:",
    adata_to_map.var.index.isin(reference_model_features.index).sum(),
)
```

    Number of genes found in query dataset: 3998


### 把缺失的基因的count设为0添加到Anndata对象里


```python
missing_genes = [
    gene_id
    for gene_id in reference_model_features.index
    if gene_id not in adata_to_map.var.index
]
```


```python
missing_gene_adata = sc.AnnData(
    X=csr_matrix(np.zeros(shape=(adata.n_obs, len(missing_genes))), dtype="float32"),
    obs=adata.obs.iloc[:, :1],
    var=reference_model_features.loc[missing_genes, :],
)
missing_gene_adata.layers["counts"] = missing_gene_adata.X
```


```python
if "PCs" in adata_to_map.varm.keys():
    del adata_to_map.varm["PCs"]
```


```python
adata_to_map_augmented = sc.concat(
    [adata_to_map, missing_gene_adata],
    axis=1,
    join="outer",
    index_unique=None,
    merge="unique",
)
```


```python
adata_to_map_augmented = adata_to_map_augmented[
    :, reference_model_features.index
].copy()
```


```python
(adata_to_map_augmented.var.index == reference_model_features.index).all()
```




    True



### 现在可以把index改回gene_name了


```python
adata_to_map_augmented.var["gene_ids"] = adata_to_map_augmented.var.index
adata_to_map_augmented.var.set_index("gene_names", inplace=True)
```


```python
adata_to_map_augmented.obs.batch.unique()
```




    ['12']
    Categories (1, object): ['12']




```python
# loading model.pt from figshare
if not Path("./reference_model").exists():
    Path("./reference_model").mkdir()
elif not Path("./reference_model/model.pt").exists():
    urllib.request.urlretrieve(
        "https://figshare.com/ndownloader/files/41436648",
        filename="reference_model/model.pt",
    )
```


```python
scarches_model = sca.models.SCVI.load_query_data(
    adata=adata_to_map_augmented,
    reference_model="./reference_model",
    freeze_dropout=True,
)
```

    [34mINFO    [0m File .[35m/reference_model/[0m[95mmodel.pt[0m already downloaded                                                        


    INFO:2025-10-02 16:23:51,985:jax._src.xla_bridge:925: Unable to initialize backend 'rocm': module 'jaxlib.xla_extension' has no attribute 'GpuAllocatorConfig'
    Unable to initialize backend 'rocm': module 'jaxlib.xla_extension' has no attribute 'GpuAllocatorConfig'
    INFO:2025-10-02 16:23:52,041:jax._src.xla_bridge:925: Unable to initialize backend 'tpu': INTERNAL: Failed to open libtpu.so: libtpu.so: cannot open shared object file: No such file or directory
    Unable to initialize backend 'tpu': INTERNAL: Failed to open libtpu.so: libtpu.so: cannot open shared object file: No such file or directory


### train query data


```python
scarches_model.train(max_epochs=500, plan_kwargs={"weight_decay": 0.0})
```

    INFO: GPU available: False, used: False
    GPU available: False, used: False
    INFO: TPU available: False, using: 0 TPU cores
    TPU available: False, using: 0 TPU cores
    INFO: IPU available: False, using: 0 IPUs
    IPU available: False, using: 0 IPUs
    INFO: HPU available: False, using: 0 HPUs
    HPU available: False, using: 0 HPUs


    Epoch 500/500: 100%|██████████| 500/500 [14:45<00:00,  1.80s/it, v_num=1, train_loss_step=997, train_loss_epoch=1.04e+3]    

    INFO: `Trainer.fit` stopped: `max_epochs=500` reached.
    `Trainer.fit` stopped: `max_epochs=500` reached.


    Epoch 500/500: 100%|██████████| 500/500 [14:45<00:00,  1.77s/it, v_num=1, train_loss_step=997, train_loss_epoch=1.04e+3]



```python
adata.obsm["X_scVI"] = scarches_model.get_latent_representation()
```


```python
sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.umap(adata)
```


```python
sc.pl.umap(
    adata,
    color=["IGHD", "IGHM", "PRDM1"],
    vmin=0,
    vmax="p99",  # set vmax to the 99th percentile of the gene count instead of the maximum, to prevent outliers from making expression in other cells invisible. Note that this can cause problems for extremely lowly expressed genes.
    sort_order=False,  # do not plot highest expression on top, to not get a biased view of the mean expression among cells
    frameon=False,
    cmap="Reds",  # or choose another color map e.g. from here: https://matplotlib.org/stable/tutorials/colors/colormaps.html
)
```


    
![png](output_91_0.png)
    



```python
ref_emb = sc.read(
    filename="RData/reference_embedding.h5ad"
)
```


```python
ref_emb
```




    AnnData object with n_obs × n_vars = 86332 × 10
        obs: 'donor', 'batch', 'site', 'cell_type'
        uns: 'neighbors', 'umap'
        obsm: 'X_umap'
        obsp: 'connectivities', 'distances'




```python
ref_emb.obs["reference_or_query"] = "reference"
```


```python
adata_emb = sc.AnnData(X=adata.obsm["X_scVI"], obs=adata.obs)
```


```python
adata_emb.obs["reference_or_query"] = "query"
```


```python
emb_ref_query = sc.concat(
    [ref_emb, adata_emb],
    axis=0,
    join="outer",
    index_unique=None,
    merge="unique",
)
```


```python
sc.pp.neighbors(emb_ref_query)
sc.tl.umap(emb_ref_query)
```


```python
sc.pl.umap(
    emb_ref_query,
    color=["reference_or_query"],
    sort_order=False,
    frameon=False,
)
```

    ... storing 'reference_or_query' as categorical



    
![png](output_99_1.png)
    



```python
sc.set_figure_params(figsize=(8, 8))
```


```python
sc.pl.umap(
    emb_ref_query,
    color=["cell_type"],
    sort_order=False,
    frameon=False,
    legend_loc="on data",
    legend_fontsize=10,
    na_color="black",
)
```


    
![png](output_101_0.png)
    



```python
knn_transformer = sca.utils.knn.weighted_knn_trainer(
    train_adata=ref_emb,
    train_adata_emb="X",  # location of our joint embedding
    n_neighbors=15,
)
```

    Weighted KNN with n_neighbors = 15 ... 


```python
labels, uncert = sca.utils.knn.weighted_knn_transfer(
    query_adata=adata_emb,
    query_adata_emb="X",  # location of our embedding, query_adata.X in this case
    label_keys="cell_type",  # (start of) obs column name(s) for which to transfer labels
    knn_model=knn_transformer,
    ref_adata_obs=ref_emb.obs,
)
```


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    Cell In[1], line 1
    ----> 1 labels, uncert = sca.utils.knn.weighted_knn_transfer(
          2     query_adata=adata_emb,
          3     query_adata_emb="X",  # location of our embedding, query_adata.X in this case
          4     label_keys="cell_type",  # (start of) obs column name(s) for which to transfer labels
          5     knn_model=knn_transformer,
          6     ref_adata_obs=ref_emb.obs,
          7 )


    NameError: name 'sca' is not defined



```python
adata_emb.obs["transf_cell_type"] = labels.loc[adata_emb.obs.index, "cell_type"]
adata_emb.obs["transf_cell_type_unc"] = uncert.loc[adata_emb.obs.index, "cell_type"]
```


```python
adata.obs.loc[adata_emb.obs.index, "transf_cell_type"] = adata_emb.obs[
    "transf_cell_type"
]
adata.obs.loc[adata_emb.obs.index, "transf_cell_type_unc"] = adata_emb.obs[
    "transf_cell_type_unc"
]
```


```python
sc.set_figure_params(figsize=(5, 5))
```


```python
sc.pl.umap(adata, color="transf_cell_type", frameon=False)
```


```python
sc.pl.umap(adata, color="transf_cell_type_unc", frameon=False)
```


```python
fig, ax = plt.subplots(figsize=(8, 3))
ct_order = (
    adata.obs.groupby("transf_cell_type")
    .agg({"transf_cell_type_unc": "median"})
    .sort_values(by="transf_cell_type_unc", ascending=False)
)
sns.boxplot(
    adata.obs,
    x="transf_cell_type",
    y="transf_cell_type_unc",
    color="grey",
    ax=ax,
    order=ct_order.index,
)
ax.tick_params(rotation=90, axis="x")
```


```python
adata.obs["transf_cell_type_certain"] = adata.obs.transf_cell_type.tolist()
adata.obs.loc[adata.obs.transf_cell_type_unc > 0.2, "transf_cell_type_certain"] = (
    "Unknown"
)
```


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    Cell In[3], line 1
    ----> 1 adata.obs["transf_cell_type_certain"] = adata.obs.transf_cell_type.tolist()
          2 adata.obs.loc[adata.obs.transf_cell_type_unc > 0.2, "transf_cell_type_certain"] = (
          3     "Unknown"
          4 )


    NameError: name 'adata' is not defined



```python
sc.pl.umap(adata, color="transf_cell_type_certain", frameon=False)
```


```python
sc.pl.umap(adata, color="transf_cell_type_certain", groups="Unknown")
```


```python
cell_types_to_check = [
    "CD14+ Mono",
    "cDC2",
    "NK",
    "B1 B",
    "CD4+ T activated",
    "T naive",
    "MK/E prog",
]
```


```python
sc.pl.dotplot(
    adata,
    var_names={
        ct: marker_genes_in_data[ct] for ct in cell_types_to_check
    },  # gene names grouped by cell type in a dictionary
    groupby="transf_cell_type_certain",
    standard_scale="var",  # normalize gene scores from 0 to 1
)
```


```python
sc.pl.umap(
    adata, color=["transf_cell_type_unc", "transf_cell_type_certain"], frameon=False
)
```


```python
import scanpy as sc
```


```python
?sc.concat
```


    [0;31mSignature:[0m
    [0msc[0m[0;34m.[0m[0mconcat[0m[0;34m([0m[0;34m[0m
    [0;34m[0m    [0madatas[0m[0;34m:[0m [0;34m'Collection[AnnData] | typing.Mapping[str, AnnData]'[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0;34m*[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0maxis[0m[0;34m:[0m [0;34m'Literal[0, 1]'[0m [0;34m=[0m [0;36m0[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0mjoin[0m[0;34m:[0m [0;34m"Literal['inner', 'outer']"[0m [0;34m=[0m [0;34m'inner'[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0mmerge[0m[0;34m:[0m [0;34m'StrategiesLiteral | Callable | None'[0m [0;34m=[0m [0;32mNone[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0muns_merge[0m[0;34m:[0m [0;34m'StrategiesLiteral | Callable | None'[0m [0;34m=[0m [0;32mNone[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0mlabel[0m[0;34m:[0m [0;34m'str | None'[0m [0;34m=[0m [0;32mNone[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0mkeys[0m[0;34m:[0m [0;34m'Collection | None'[0m [0;34m=[0m [0;32mNone[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0mindex_unique[0m[0;34m:[0m [0;34m'str | None'[0m [0;34m=[0m [0;32mNone[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0mfill_value[0m[0;34m:[0m [0;34m'Any | None'[0m [0;34m=[0m [0;32mNone[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0mpairwise[0m[0;34m:[0m [0;34m'bool'[0m [0;34m=[0m [0;32mFalse[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m[0;34m)[0m [0;34m->[0m [0;34m'AnnData'[0m[0;34m[0m[0;34m[0m[0m
    [0;31mDocstring:[0m
    Concatenates AnnData objects along an axis.
    
    See the :doc:`concatenation <../concatenation>` section in the docs for a more in-depth description.
    
    Params
    ------
    adatas
        The objects to be concatenated. If a Mapping is passed, keys are used for the `keys`
        argument and values are concatenated.
    axis
        Which axis to concatenate along.
    join
        How to align values when concatenating. If "outer", the union of the other axis
        is taken. If "inner", the intersection. See :doc:`concatenation <../concatenation>`
        for more.
    merge
        How elements not aligned to the axis being concatenated along are selected.
        Currently implemented strategies include:
    
        * `None`: No elements are kept.
        * `"same"`: Elements that are the same in each of the objects.
        * `"unique"`: Elements for which there is only one possible value.
        * `"first"`: The first element seen at each from each position.
        * `"only"`: Elements that show up in only one of the objects.
    uns_merge
        How the elements of `.uns` are selected. Uses the same set of strategies as
        the `merge` argument, except applied recursively.
    label
        Column in axis annotation (i.e. `.obs` or `.var`) to place batch information in.
        If it's None, no column is added.
    keys
        Names for each object being added. These values are used for column values for
        `label` or appended to the index if `index_unique` is not `None`. Defaults to
        incrementing integer labels.
    index_unique
        Whether to make the index unique by using the keys. If provided, this
        is the delimiter between "{orig_idx}{index_unique}{key}". When `None`,
        the original indices are kept.
    fill_value
        When `join="outer"`, this is the value that will be used to fill the introduced
        indices. By default, sparse arrays are padded with zeros, while dense arrays and
        DataFrames are padded with missing values.
    pairwise
        Whether pairwise elements along the concatenated dimension should be included.
        This is False by default, since the resulting arrays are often not meaningful.
    
    Notes
    -----
    
    .. warning::
    
        If you use `join='outer'` this fills 0s for sparse data when
        variables are absent in a batch. Use this with care. Dense data is
        filled with `NaN`.
    
    Examples
    --------
    
    Preparing example objects
    
    >>> import anndata as ad, pandas as pd, numpy as np
    >>> from scipy import sparse
    >>> a = ad.AnnData(
    ...     X=sparse.csr_matrix(np.array([[0, 1], [2, 3]])),
    ...     obs=pd.DataFrame({"group": ["a", "b"]}, index=["s1", "s2"]),
    ...     var=pd.DataFrame(index=["var1", "var2"]),
    ...     varm={
    ...         "ones": np.ones((2, 5)),
    ...         "rand": np.random.randn(2, 3),
    ...         "zeros": np.zeros((2, 5)),
    ...     },
    ...     uns={"a": 1, "b": 2, "c": {"c.a": 3, "c.b": 4}},
    ... )
    >>> b = ad.AnnData(
    ...     X=sparse.csr_matrix(np.array([[4, 5, 6], [7, 8, 9]])),
    ...     obs=pd.DataFrame(
    ...         {"group": ["b", "c"], "measure": [1.2, 4.3]}, index=["s3", "s4"]
    ...     ),
    ...     var=pd.DataFrame(index=["var1", "var2", "var3"]),
    ...     varm={"ones": np.ones((3, 5)), "rand": np.random.randn(3, 5)},
    ...     uns={"a": 1, "b": 3, "c": {"c.b": 4}},
    ... )
    >>> c = ad.AnnData(
    ...     X=sparse.csr_matrix(np.array([[10, 11], [12, 13]])),
    ...     obs=pd.DataFrame({"group": ["a", "b"]}, index=["s1", "s2"]),
    ...     var=pd.DataFrame(index=["var3", "var4"]),
    ...     uns={"a": 1, "b": 4, "c": {"c.a": 3, "c.b": 4, "c.c": 5}},
    ... )
    
    Concatenating along different axes
    
    >>> ad.concat([a, b]).to_df()
        var1  var2
    s1     0     1
    s2     2     3
    s3     4     5
    s4     7     8
    >>> ad.concat([a, c], axis=1).to_df()
        var1  var2  var3  var4
    s1     0     1    10    11
    s2     2     3    12    13
    
    Inner and outer joins
    
    >>> inner = ad.concat([a, b])  # Joining on intersection of variables
    >>> inner
    AnnData object with n_obs × n_vars = 4 × 2
        obs: 'group'
    >>> (inner.obs_names, inner.var_names)  # doctest: +NORMALIZE_WHITESPACE
    (Index(['s1', 's2', 's3', 's4'], dtype='object'),
    Index(['var1', 'var2'], dtype='object'))
    >>> outer = ad.concat([a, b], join="outer")  # Joining on union of variables
    >>> outer
    AnnData object with n_obs × n_vars = 4 × 3
        obs: 'group', 'measure'
    >>> outer.var_names
    Index(['var1', 'var2', 'var3'], dtype='object')
    >>> outer.to_df()  # Sparse arrays are padded with zeroes by default
        var1  var2  var3
    s1     0     1     0
    s2     2     3     0
    s3     4     5     6
    s4     7     8     9
    
    Keeping track of source objects
    
    >>> ad.concat({"a": a, "b": b}, label="batch").obs
       group batch
    s1     a     a
    s2     b     a
    s3     b     b
    s4     c     b
    >>> ad.concat([a, b], label="batch", keys=["a", "b"]).obs  # Equivalent to previous
       group batch
    s1     a     a
    s2     b     a
    s3     b     b
    s4     c     b
    >>> ad.concat({"a": a, "b": b}, index_unique="-").obs
         group
    s1-a     a
    s2-a     b
    s3-b     b
    s4-b     c
    
    Combining values not aligned to axis of concatenation
    
    >>> ad.concat([a, b], merge="same")
    AnnData object with n_obs × n_vars = 4 × 2
        obs: 'group'
        varm: 'ones'
    >>> ad.concat([a, b], merge="unique")
    AnnData object with n_obs × n_vars = 4 × 2
        obs: 'group'
        varm: 'ones', 'zeros'
    >>> ad.concat([a, b], merge="first")
    AnnData object with n_obs × n_vars = 4 × 2
        obs: 'group'
        varm: 'ones', 'rand', 'zeros'
    >>> ad.concat([a, b], merge="only")
    AnnData object with n_obs × n_vars = 4 × 2
        obs: 'group'
        varm: 'zeros'
    
    The same merge strategies can be used for elements in `.uns`
    
    >>> dict(ad.concat([a, b, c], uns_merge="same").uns)
    {'a': 1, 'c': {'c.b': 4}}
    >>> dict(ad.concat([a, b, c], uns_merge="unique").uns)
    {'a': 1, 'c': {'c.a': 3, 'c.b': 4, 'c.c': 5}}
    >>> dict(ad.concat([a, b, c], uns_merge="only").uns)
    {'c': {'c.c': 5}}
    >>> dict(ad.concat([a, b, c], uns_merge="first").uns)
    {'a': 1, 'b': 2, 'c': {'c.a': 3, 'c.b': 4, 'c.c': 5}}
    [0;31mFile:[0m      ~/miniconda3/envs/sc_annotation/lib/python3.10/site-packages/anndata/_core/merge.py
    [0;31mType:[0m      function

