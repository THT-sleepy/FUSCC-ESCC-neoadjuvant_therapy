# inferCNV过程及结果记录         (What)

* Oct 26, 2025                                 (When)
* liyi /data1/liyi/zhaoyue/FUSCC-ESCC-neoadjuvant_thrapy           (Where)
* 通过inferCNV区分正常和恶性上皮 (Why)

## 1 在进行分析前,需要对亚群注释后整合的adata去一下乱跑的细胞
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
sc.pl.umap(adata, color='first_round_celltype', title="1st round anno ct umap", ax=ax, show=False)
ax.set_xscale('linear')
ax.set_yscale('linear')
ax.xaxis.set_major_locator(MaxNLocator(nbins=8))  # x轴最多显示8个刻度（可按需调整）
ax.yaxis.set_major_locator(MaxNLocator(nbins=8))  # y轴最多显示8个刻度
plt.savefig("plots/umap_2ndanno_withaxis.png", dpi=300, bbox_inches="tight")
plt.show()
```
<img src="..\figures\umap_2ndanno_withaxis.png">



```
umap1 = adata.obsm['X_umap'][:, 0]  # 提取UMAP1坐标
umap2 = adata.obsm['X_umap'][:, 1]  # 提取UMAP2坐标

##B
mask = (umap1 >= -4) & (umap1 <= 12) & (umap2 >= -8) & (umap2 <= 4) & (adata.obs['first_round_celltype'] == "B")
adata_B = adata[mask, :].copy()

##Endo
mask = (umap1 >= 0) & (umap1 <= 8) & (umap2 >= 8) & (umap2 <= 12) & (adata.obs['first_round_celltype'] == "Endo")
adata_Endo = adata[mask, :].copy()

##上皮
mask = (umap1 >= 8) & (umap1 <= 20) & (umap2 >= 4) & (umap2 <= 20) & (adata.obs['first_round_celltype'] == "Epi")
adata_Epi = adata[mask, :].copy()

##纤维
mask = (umap1 >= -4) & (umap1 <= 8) & (umap2 >= 12) & (umap2 <= 20) & (adata.obs['first_round_celltype'] == "Fib")
adata_Fib = adata[mask, :].copy()

##髓系
mask = (umap1 >= -8) & (umap1 <= 8) & (umap2 >= -4) & (umap2 <= 12) & (adata.obs['first_round_celltype'] == "Myeloid")
adata_Myeloid = adata[mask, :].copy()

##T
mask = (umap1 >= 4) & (umap1 <= 20) & (umap2 >= -12) & (umap2 <= 0) & (adata.obs['first_round_celltype'] == "T")
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
sc.pl.umap(adata_filtered,color="first_round_celltype")
adata_filtered.write("RData/2nd_round_anno_umapfiltered.h5ad")
```

# 2开始分析
<img src="..\figures\cnv_umapcnv_umap.png">
<img src="..\figures\heatmapclustered_cnv_heatmap.png">
<img src="..\figures\heatmapcnv_heatmap.png">
<img src="..\figures\umapcnv_scanpyumap.png">
<img src="..\figures\umaptumor_nomal.png">

```
最后添加上注释吧epi全是肿瘤细胞(顺便把前面写漏了的补上来)
cluster_mapping = {
    # B细胞相关
    'Bn_FCER2': 'c16_Bn_FCER2',
    'Bgc_HMGB2': 'c17_Bgc_HMGB2',
    'Bm_TNFRSF13B': 'c18_Bm_TNFRSF13B',
    'Bp_MZB1': 'c19_Bp_MZB1',

    # 成纤维细胞相关
    'Fib_CFD': 'c28_Fib_CFD',
    'Fib_MMP1': 'c29_Fib_MMP1',
    'Fib_CST1': 'c30_Fib_CST1',
    'Fib_PTGDS': 'c31_Fib_PTGDS',
    'apFib_CD74': 'c32_apFib_CD74',

    # 周细胞/平滑肌
    'Pericyte_RGS5': 'c33_Pericyte_RGS5',
    'Pericyte_MYH11': 'c34_Pericyte_MYH11',
    'SMC_TAGLN': 'c35_SMC_TAGLN',
    'Epi': 'c41_Epi_Tumor'
}
df['final_anno_celltype'] = df['final_anno_celltype'].replace(cluster_mapping)





```
