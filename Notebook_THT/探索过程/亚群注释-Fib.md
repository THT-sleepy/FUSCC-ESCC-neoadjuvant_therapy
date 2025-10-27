# 纤维细胞注释过程及结果记录         (What)

* Oct 23, 2025                                 (When)
* liyi /data1/liyi/zhaoyue/FUSCC-ESCC-neoadjuvant_thrapy           (Where)
* 单细胞需要注释后才能进行后续分析 (Why)

## 1 再次检查有无混杂的其它细胞系如上皮的marker表达
*聚类情况
<img src="..\figures\umapFib_umapmd0.4_cluster.png">
```
import scanpy as sc

sc.settings.set_figure_params(
    dpi=300,
    facecolor="white",
    frameon=False
)
sc.settings.figdir = "plots"

adata = sc.read("RData/Fib__umapmd0.4.h5ad")

marker_genes = {
"Immune":["PTPRC"],
"Fib": ["DCN","COL1A1","COL1A2","FN1","COL3A1","COL6A1"],
"Endo" : ["VWF","PECAM1","ENG","CDH5"],
"Epi" : ["EPCAM","SFN","KRT5","KRT14"],
"T":["CD3D","CD2","CD3E","CD3G"],
"B":["CD79A","JCHAIN","CD19","MZB1"],
"Myeloid":["LYZ","CD14","CD68",
           "HPGDS","TPSAB1",
          "GCA"],
}
marker_genes_in_data = {}
for ct, markers in marker_genes.items():
    markers_found = []
    for marker in markers:
        if marker in adata.var.index:
            markers_found.append(marker)
    marker_genes_in_data[ct] = markers_found
sc.pl.umap(adata,color=marker_genes_in_data["T"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Fibsubset_T.png")
sc.pl.umap(adata,color=marker_genes_in_data["Endo"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Fibsubset_Endo.png")
sc.pl.umap(adata,color=marker_genes_in_data["Epi"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Fibsubset_Epi.png")
sc.pl.umap(adata,color=marker_genes_in_data["B"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Fibsubset_B.png")
sc.pl.umap(adata,color=marker_genes_in_data["Myeloid"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Fibsubset_Myeloid.png")

sc.pl.dotplot(
    adata,
    groupby="leiden_res0_6",
    var_names=marker_genes_in_data,
    standard_scale="var",
    save="Fibsubset_大类marker.png"
)
```
<img src="..\figures\umapFibsubset_T.png">
看起来没有混的T细胞<br>

<img src="..\figures\umapFibsubset_Endo.png">
好像也没有混的内皮细胞,ENG其实算不上特别好的marker

<img src="..\figures\umapFibsubset_Epi.png">
有点上皮但没有成群

<img src="..\figures\umapFibsubset_B.png">
B细胞也没有混

<img src="..\figures\umapFibsubset_Myeloid.png">
也没有混髓系的

* 点图
<img src="..\figures\dotplot_Fibsubset_大类marker.png">
除开ENG其它都没啥

## 1.5 先看看celltypist的结果
<img src="..\figures\umapcelltypist__Fib.png">
## 2 用umap优化marker，确定好分辨率
```
marker_genes = {
    "Pericyte": ["RGS5", "MCAM", "ACTA2", "MYH11", "PDGFRB", "GJA4", "CD36", "NOTCH3"],
    "VSMC": ["TAGLN", "ACTA2", "ACTG2", "MYH11", "CNN1", "HHIP"],
    "NMF": ["SLPI", "PI16", "CLU"],
    "NAF": ["CLU", "IGF1", "C7", "APOD"],
    "CAF": ["APOD", "CXCL1", "CXCL6", "CCL2", "CCL11", "CXCL14", "CXCL5", "CXCL8", "CSF3", "MMP1", "MMP3", "MMP11"],
    "myofibroblast": ["ACAT2", "COL1A2", "PDGFRB", "HOPX", "IGFBP5", "TIMP1", "MMP11", "COL10A1", "POSTN", "LRRC15", "SFRP4", "SFRP2", "COMP", "MMP1", "COL7A1", "WNT5A", "ISG15", "IL7R"],
    "mCAF": ["POSTN", "FN1", "LUM", "DCN", "VCAN", "COL5A1", "COL5A2", "COL6A3"],
    "Inflammatory fibroblasts": ["FBLN1", "IGF1", "CXCL1", "IGFBP6", "SLPI", "SAA1", "C3", "C7", "CEBPD", "CLU", "CTGF", "HGF", "HSPA6", "DNAJB1", "MYC", "AFT4"],
    "ap-fibroblasts": ["CD74", "HLA-DRA", "HLA-DRB1", "CXCL12"],
    "Lip-fibroblasts": ["APOPA2", "FABP1", "FABP4", "FRZB", "CFD", "APOD"],
    "proliferative fibroblasts": ["MKI67", "TOP2A", "RBP1", "STAR", "STMN1"],
    "CTNNB1+ fibroblasts": ["WSB1", "DDX17", "CTNNB1"],
    "progenitor-like fibroblasts": ["MFAP5", "PI16", "CD34", "APOD", "FGF7", "COL15A1"],
}
marker_genes_in_data = {}
for ct, markers in marker_genes.items():
    markers_found = []
    for marker in markers:
        if marker in adata.var.index:
            markers_found.append(marker)
    marker_genes_in_data[ct] = markers_found
sc.pl.umap(adata,color=marker_genes_in_data["Pericyte"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Fibresubset_pericyte.png")
sc.pl.umap(adata,color=marker_genes_in_data["VSMC"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Fibresubset_vsmc.png")
sc.pl.umap(adata,color=marker_genes_in_data["NMF"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Fibresubset_NMF.png")
sc.pl.umap(adata,color=marker_genes_in_data["NAF"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Fibresubset_NAF.png")
sc.pl.umap(adata,color=marker_genes_in_data["CAF"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Fibresubset_CAF.png")
sc.pl.umap(adata,color=marker_genes_in_data["myofibroblast"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Fibresubset_myofibroblast.png")
sc.pl.umap(adata,color=marker_genes_in_data["mCAF"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Fibresubset_mCAF.png")
sc.pl.umap(adata,color=marker_genes_in_data["Inflammatory fibroblasts"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Fibresubset_iF.png")
sc.pl.umap(adata,color=marker_genes_in_data["ap-fibroblasts"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Fibresubset_ap-fibroblasts.png")
sc.pl.umap(adata,color=marker_genes_in_data["Lip-fibroblasts"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Fibresubset_Lip-fibroblasts.png")
sc.pl.umap(adata,color=marker_genes_in_data["proliferative fibroblasts"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="proliferative fibroblasts.png")
sc.pl.umap(adata,color=marker_genes_in_data["CTNNB1+ fibroblasts"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Fibresubset_ctnnb1fibroblasts.png")
sc.pl.umap(adata,color=marker_genes_in_data["progenitor-like fibroblasts"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Fibresubset_profibroblasts.png")
```

* pericyte
<img src="..\figures\umapFib_umapmd0.4_cluster.png">
<img src="..\figures\umapFibresubset_pericyte.png">
就RGS5+GJA4就挺好的,对应右边这一大团

* vsmc
<img src="..\figures\umapFib_umapmd0.4_cluster.png">
<img src="..\figures\umapFibresubset_vsmc.png">
不要HHIP和MYH11,看起来像是res0_8/1的5群

* NMF
<img src="..\figures\umapFib_umapmd0.4_cluster.png">
<img src="..\figures\umapFibresubset_nmf.png">
就要这三个marker吧，看起来像是res0_8的第0群

* NAF
<img src="..\figures\umapFib_umapmd0.4_cluster.png">
<img src="..\figures\umapFibresubset_naf.png">
感觉和NMF是同一块区域,后面再看看吧

* CAF
<img src="..\figures\umapFib_umapmd0.4_cluster.png">
<img src="..\figures\umapFibresubset_caf.png">
CXCL5,CXCL8,,CSF3,MMP1,MMP3有点像res1的第9群

* myofibroblast
<img src="..\figures\umapFib_umapmd0.4_cluster.png">
<img src="..\figures\umapFibresubset_myofibroblast.png">
感觉这些marker都不太特异，先不管这一群

* mCAF
<img src="..\figures\umapFib_umapmd0.4_cluster.png">
<img src="..\figures\umapFibresubset_mCAF.png">
这些marker也是不太特异

* Inflammatory fibroblasts
<img src="..\figures\umapFib_umapmd0.4_cluster.png">
<img src="..\figures\umapFibresubset_iF.png">
这些marker也不是很特异

* ap-fibroblasts
<img src="..\figures\umapFib_umapmd0.4_cluster.png">
<img src="..\figures\umapFibresubset_ap-fibroblasts.png">
应该是下面这个小细胞团,保留前3个marker

* Lip-fibroblasts
<img src="..\figures\umapFib_umapmd0.4_cluster.png">
<img src="..\figures\umapFibresubset_Lip-fibroblasts.png">
也是并不特异

* proliferative fibroblasts
<img src="..\figures\umapFib_umapmd0.4_cluster.png">
<img src="..\figures\umapFibresubset_profibroblasts.png">
不特异

* CTNNB1+ fibroblasts
<img src="..\figures\umapFib_umapmd0.4_cluster.png">
<img src="..\figures\umapFibresubset_ctnnb1fibroblasts.png">
不特异

* progenitor-like fibroblasts
<img src="..\figures\umapFib_umapmd0.4_cluster.png">
<img src="..\figures\umapFibresubset_profibroblasts.png">
不特异

总结一下就是右边是周细胞,下面的小细胞团是ap-fibroyte,左边的左下角是VSMC,其它是分不清功能的Fibrocyte
,左边的可以就用特征基因区分吧,res选0.8

## 3 DEG找marker

```
sc.tl.rank_genes_groups(
 adata, groupby="leiden_res0_8", method="wilcoxon", key_added="dea_leiden_res0_8"
)
sc.tl.filter_rank_genes_groups(
  adata,
  min_in_group_fraction=0.2,
  max_out_group_fraction=0.2,
  key="dea_leiden_res0_8",
  key_added="dea_leiden_res0_8_filtered",
 )
sc.pl.rank_genes_groups_dotplot(
  adata,
  groupby="leiden_res0_8",
  standard_scale="var",
  n_genes=10,
  key="dea_leiden_res0_8_filtered",
  save="Myeloid_top10gene.png"
)
```
<img src="..\figures\umapFib_umapmd0.4_cluster.png">
<img src="..\figures\dotplot_Myeloid_top10gene.png">


```
adata = sc.read("RData/Fib__umapmd0.4.h5ad")
cl_annotation = {
    "0": "Fib_CFD",
    "1": "Fib_MMP1",
    "2": "Pericyte_RGS5",
    "3": "Fib_CST1",
    "4": "Fib_PTGDS",
    "5":"SMC_TAGLN",
    "6": "Pericyte_MYH11",
    "7": "apFib_CD74",
}
adata.obs["final_anno_celltype"]=adata.obs.leiden_res0_8.map(cl_annotation)
```

## final umap和点图
```
sc.pl.umap(adata,color="final_anno_celltype",save="Fibsubset_anno.png")


marker_genes = {
"Fib_CFD":["CFD","ABI3BP","C7"],
"Fib_MMP1":["MMP1","CHI3L1","CCN4",],
"Pericyte_RGS5":["RGS5","GJA4","NDUFA4L2"],
"Fib_CST1":["CST1","SFRP4"],
"Fib_PTGDS":["PTGDS","CXCL14","COL6A5"],
"SMC":["TAGLN", "ACTA2", "ACTG2","CNN1","ALDH1B1"],
"Pericyte_MYH11":["RGS5","GJA4","MYH11"],
"apFib_CD74":["CD74","HLA-DRA", "HLA-DRB1"]
}
sc.pl.dotplot(
    adata,
    groupby="final_anno_celltype",
    var_names=marker_genes,
    standard_scale="var",
    save="Fib_finalanno.png"
)

```
* final umap
<img src="..\figures\umapFibsubset_anno.png">
* dotplot
<img src="..\figures\dotplot_Fib_finalanno.png">
