# Ttr-like细胞亚群的相关分析       (What)

* Dec 18, 2025                                 (When)
* liyi /data1/liyi/zhaoyue/FUSCC-ESCC-neoadjuvant_thrapy           (Where)
* 看一下Ttr-like细胞具体包含哪些亚群 (Why)


### 导入模块
```python
import scanpy as sc
import pandas as pd
import numpy as np

sc.settings.set_figure_params(
    dpi=300,
    facecolor="white",
    frameon=False
)
sc.settings.figdir = "plots"

adata = sc.read("RData/1128_final_escc120.h5ad")
```


#### 将clone_id以及是否是Ttr细胞的信息整理进来
```
### 整合clone.id
df_clone_id = pd.read_csv("input/clone_id_1216.csv", sep=",")
clone_map = dict(zip(df_clone_id["cell_id"], df_clone_id["clone.id"]))
adata.obs["clone_id"] = adata.obs.index.map(clone_map)

### 整合Ttr
df_Ttr = pd.read_csv("input/Ttr_cellids.csv", sep=",")
ttr_map = dict(zip(df_Ttr["cell_id"], df_Ttr["Ttr"]))
adata.obs["Ttr"] = adata.obs.index.map(ttr_map)
```


























#### 对Ttr细胞再做聚类和注释
```
adata_sub = adata[adata.obs["Ttr"]==True].copy()

### 1.去掉细胞数小于10个的样本
sample_counts = adata_sub.obs['sample'].value_counts()
keep_samples = sample_counts[sample_counts >= 10].index
adata_sub = adata_sub[adata_sub.obs['sample'].isin(keep_samples)].copy()

# 2. 重新normalize和feature selection
adata_sub.X = adata_sub.layers["soupX_counts"]  # 用soupX校正后的counts
scales_counts = sc.pp.normalize_total(adata_sub, target_sum=None, inplace=False)
adata_sub.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)  # log1p转换
sc.pp.highly_variable_genes(  # 筛选高变基因
       adata_sub, layer="log1p_norm",batch_key="sample",n_top_genes=2000
)
adata_sub.layers["regress"] = adata_sub.layers["log1p_norm"]
adata_sub.X = adata_sub.layers["regress"]  # 将log1p后的数据设为默认X

#3 regressout mt ribo
sc.pp.regress_out(adata_sub, ['pct_counts_mt','pct_counts_ribo'])

#4 PCA降维
sc.pp.pca(adata_sub, svd_solver="arpack", use_highly_variable=True)

#5. 循环运行不同min_dist的UMAP+Leiden聚类（0.1/0.2/0.4）
    for min_dist in [0.2,0.4]:
        adata_md = adata_sub.copy()

        # BBKNN批次校正 + UMAP降维
        sc.external.pp.bbknn(adata_md, batch_key="sample", n_pcs=30)
        sc.tl.umap(adata_md, min_dist=min_dist,random_state=4)

        # 多分辨率Leiden聚类
        sc.tl.leiden(adata_md, key_added="leiden_res0_6", resolution=0.6)
        sc.tl.leiden(adata_md, key_added="leiden_res0_8", resolution=0.8)
        sc.tl.leiden(adata_md, key_added="leiden_res1", resolution=1.0)

        # 保存UMAP聚类图
        plot_path = f"Ttr_umapmd{min_dist}_cluster.png"
        sc.pl.umap(
            adata_md,
            color=["leiden_res0_6", "leiden_res0_8", "leiden_res1"],
            legend_loc="on data",
            save=plot_path,
            show=False  # 避免弹出图片窗口，仅保存
        )

        # 保存AnnData对象
        output_fpath = f"RData/Ttr_umapmd{min_dist}.h5ad"
        adata_md.X = adata_md.layers["log1p_norm"]
        adata_md.obs["Ttr"] = adata_md.obs["Ttr"].astype(str)
        adata_md.write(output_fpath)
```
### 没啥区别,可以就用md0.4+res0.6
<img src="..\figures\umapTtr_umapmd0.2_cluster.png">
<img src="..\figures\umapTtr_umapmd0.4_cluster.png">
