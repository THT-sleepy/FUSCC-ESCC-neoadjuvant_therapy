# 把注释好的亚群整合好后做一下处理         (What)

* Nov 26, 2025                                 (When)
* liyi /data1/liyi/zhaoyue/FUSCC-ESCC-neoadjuvant_thrapy           (Where)
*  no why(Why)

### 有一些细胞会乱跑需要去掉,另外还要去掉血液来源的非血液细胞
```
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import pandas as pd

sc.settings.set_figure_params(
    dpi=300,
    facecolor="white",
    frameon=True
)
sc.settings.figdir = "plots"

adata = sc.read("RData/merge2ndanno.h5ad")

fig, ax = plt.subplots(1, 1)
sc.pl.umap(adata, color='major_celltype', title="1st round anno ct umap", ax=ax, show=False)
ax.set_xscale('linear')
ax.set_yscale('linear')
ax.xaxis.set_major_locator(MaxNLocator(nbins=8))  # x轴最多显示8个刻度（可按需调整）
ax.yaxis.set_major_locator(MaxNLocator(nbins=8))  # y轴最多显示8个刻度
plt.savefig("plots/umapmergedannotated_withaxis.png", dpi=300, bbox_inches="tight")
plt.show()
```
<img src="..\figures\umapmergedannotated_withaxis.png">

```
umap1 = adata.obsm['X_umap'][:, 0]  # 提取UMAP1坐标
umap2 = adata.obsm['X_umap'][:, 1]  # 提取UMAP2坐标

##B
mask = (umap1 >= 0) & (umap1 <= 20) & (umap2 >= 0) & (umap2 <= 6) & (adata.obs['major_celltype'] == "B cell")
adata_B = adata[mask, :].copy()

##Endo
mask = (umap1 >= 15) & (umap1 <= 25) & (umap2 >= 4) & (umap2 <= 8) & (adata.obs['major_celltype'] == "Endothelial cell")
adata_Endo = adata[mask, :].copy()

##上皮
mask = (umap1 >= 3) & (umap1 <= 15) & (umap2 >= 6) & (umap2 <= 20) & (adata.obs['major_celltype'] == "Epithelial cell")
adata_Epi = adata[mask, :].copy()

##纤维
mask = (umap1 >= 15) & (umap1 <= 25) & (umap2 >= 8) & (umap2 <= 16) & (adata.obs['major_celltype'] == "Fibrocyte")
adata_Fib = adata[mask, :].copy()

##髓系
mask = (umap1 >= -10) & (umap1 <= 5) & (umap2 >= -5) & (umap2 <= 12) & (adata.obs['major_celltype'] == "Myeloid cell")
adata_Myeloid = adata[mask, :].copy()

##T
mask = (umap1 >= 0) & (umap1 <= 15) & (umap2 >= -12) & (umap2 <= 0) & (adata.obs['major_celltype'] == "T&ILC cell")
adata_T = adata[mask, :].copy()

adatas ={
"B":adata_B,
"Endo":adata_Endo,
"Epi":adata_Epi,
"Fib":adata_Fib,
"Myeloid":adata_Myeloid,
"T":adata_T,
}
adata_filtered = sc.concat(adatas,join="outer")
all_var = [x.var for x in adatas.values()]
all_var = pd.concat(all_var,join="outer")
all_var = all_var[~all_var.index.duplicated()]
adata_filtered.var = all_var.loc[adata_filtered.var_names]
sc.pl.umap(adata_filtered,color="major_celltype",save="merged_filtered_1126.png")
adata_filtered.write("RData/2nd_round_anno_umapfiltered.h5ad")
```
<img src="..\figures\umapmerged_filtered_1126.png">



# 2开始分析
<img src="..\figures\cnv_umap1126_cnv_umap.png">
<img src="..\figures\heatmap1126_clustered_cnv_heatmap.png">
<img src="..\figures\heatmap1126_cnv_heatmap.png">
<img src="..\figures\umap1126_cnv_scanpyumap.png">

### 从cnv均值出发,看下每个cnv cluser的cnv score 均值
```
adata = sc.read("RData/infercnv.h5ad")
adata.obs.groupby("cnv_leiden")["cnv_score"].mean()
```
<img src="..\figures\1128cnv_cluster_cnvscore均值.png"><br>
0,2,3,4,6,7,8,11,13,19都是非上皮,其cnv均值最大为第3群0.003856
其它都是上皮,cnv均值最小的是第1群0.002344
所以用cnv均值来区分良恶性感觉不太靠谱，感觉是通过聚类找不到正常上皮的次选

通过聚类可以知道，9和30两群上皮和非上皮的cnv是最接近的，可以把这两群认为是正常上皮，
不过第9群其实仔细看也像肿瘤上皮，看起来是没有完美的办法，不过这种方法是infercnv推荐的，
那就这么弄吧

```
adata.obs["cnv_status"] = "tumor"
adata.obs.loc[adata.obs["cnv_leiden"].isin(["0", "2", "3", "4", "6", "7", "8", "11", "13","19","9","30"]), "cnv_status"] = (
    "normal"
)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={"wspace": 0.5})
cnv.pl.umap(adata, color="cnv_status", ax=ax1, show=False)
sc.pl.umap(adata, color="cnv_status", ax=ax2,save="1128umaptumor_nomal.png")
```
<img src="..\figures\umap1128umaptumor_nomal.png">

最后添加上上皮的注释,另外还要去掉血液来源的非血液细胞

```
adata.obs["minor_celltype"] = adata.obs["minor_celltype"].cat.add_categories(["c51_Epi_Normal", "c52_Epi_Tumor"])

cond_tumor = (adata.obs["major_celltype"] == "Epithelial cell") & (adata.obs["cnv_status"] == "tumor")
adata.obs.loc[cond_tumor, "minor_celltype"] = "c52_Epi_Tumor"

cond_normal = (adata.obs["major_celltype"] == "Epithelial cell") & (adata.obs["cnv_status"] == "normal")
adata.obs.loc[cond_normal, "minor_celltype"] = "c51_Epi_Normal"

#删除血液来源的非血液细胞
cond_pbmc = adata.obs["sample"].str.contains("pbmc", case=False, na=False)
cond_celltype = adata.obs["major_celltype"].isin(["Epithelial cell", "Endothelial cell", "Fibrocyte"])
cond_drop = cond_pbmc & cond_celltype

# 步骤2：保留“非需删除”的细胞（~表示取反）
adata_filtered = adata[~cond_drop, :].copy()

#删除SZR_op_2279307
adata_filtered = adata_filtered[~adata_filtered.obs['sample']=="SZR_op_2279307", :].copy()

#修改c10的注释(从Teff修改为Tem/Teff)
adata_filtered.obs['minor_celltype'] = adata_filtered.obs['minor_celltype'].replace(
    'c10_CD8_Teff_GZMK',
    'c10_CD8_Tem/Teff_GZMK'
)
adata_filtered.write("RData/1128_final_escc120.h5ad")
```

### 后面还把有一个基因数特别少的样本给去掉了 2279307_post_tumor
