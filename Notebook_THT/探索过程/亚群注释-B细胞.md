# B细胞注释过程及结果记录         (What)

* Oct 20, 2025                                 (When)
* liyi /data1/liyi/zhaoyue/FUSCC-ESCC-neoadjuvant_thrapy           (Where)
* 单细胞需要注释后才能进行后续分析 (Why)

## 1 再次检查有无混杂的其它细胞系如上皮的marker表达
*聚类情况
<img src="..\figures\umapB_umapmd0.4_cluster.png">
```
import scanpy as sc

sc.settings.set_figure_params(
    dpi=300,
    facecolor="white",
    frameon=False
)
sc.settings.figdir = "plots"

adata = sc.read("RData/B__umapmd0.4.h5ad")

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
"NK":["GZMA","NCAM1","NKG7","KLRC1","IFNG","GZMB","KLRD1"],
}
marker_genes_in_data = {}
for ct, markers in marker_genes.items():
    markers_found = []
    for marker in markers:
        if marker in adata.var.index:
            markers_found.append(marker)
    marker_genes_in_data[ct] = markers_found
sc.pl.umap(adata,color=marker_genes_in_data["Fib"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Bsubset_Fib.png")
sc.pl.umap(adata,color=marker_genes_in_data["Endo"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Bsubset_Endo.png")
sc.pl.umap(adata,color=marker_genes_in_data["Epi"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Bsubset_Epi.png")
sc.pl.umap(adata,color=marker_genes_in_data["T"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Bsubset_T.png")
sc.pl.umap(adata,color=marker_genes_in_data["Myeloid"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Bsubset_Myeloid.png")
sc.pl.umap(adata,color=marker_genes_in_data["NK"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Bsubset_NK.png")
sc.pl.dotplot(
    adata,
    groupby="leiden_res0_6",
    var_names=marker_genes_in_data,
    standard_scale="var",
    save="Bsubset_大类marker.png"
)
```
<img src="..\figures\umapBsubset_Fib.png">
看起来没有混的纤维细胞<br>

<img src="..\figures\umapBsubset_Endo.png">
不知道是不是混了,PECAM1不太特异

<img src="..\figures\umapBsubset_Epi.png">
res0_8的c4是B+Epi要去掉

<img src="..\figures\umapBsubset_T.png">
T细胞没有混

<img src="..\figures\umapBsubset_Myeloid.png">
res0_8的cl13是B+中性粒要去掉
```


* 点图
<img src="..\figures\dotplot_Bsubset_大类marker.png">

### 去掉双细胞群
```
adata = adata[adata.obs.leiden_res0_8 != '4' ].copy()
adata = adata[adata.obs.leiden_res0_8 != '13' ].copy()
```

### 2 开始注释
```
marker_genes = {
"Bn":["TCL1A",
      "YBX3",
      "FCER2",
      "IL4R",
      "MS4A1",
      "FCRL1",
      "CD72",
      "BACH2",
      "IGHD",
      "IGHM"] ,
"Bm":["CD27",
      "IGHG1",
      "AIM2",
      "TNFRSF13B",
      "CRIP2",
      "ITGB1"],
"GCB":["BCL6",
       "RGS13",
       "AICDA",
       "IL21R",
       "MKI67",
       "STMN1",
       "HMGB2",
       "TOP2A",
       "TCL1A",
       "SUGCT",
       "NEIL1",
       "NEF2B",
       "LMO2",
       "GMDS",
       "MARCKSL1"],
"activated B":["CD69","CD83"],
"Transitional B": ["MME", "CD38", "CD24", "ACSM3", "MSI2"],
"Plasma cells": ["MZB1",
                "HSP90B1",
                "FNDC3B",
                "PRDM1",
                "IGKC",
                "JCHAIN",
                "CD38",
                "IRF4",
                "XBP1"],
"Plasmablast": ["XBP1",
                "RF4",
                "PRDM1",
                "PAX5",
                "MZB1",
                "JCHAIN",
                "STMN1",
                "HMGB2",
                "MKI67",
                "SDC1"],
}
marker_genes_in_data = {}
for ct, markers in marker_genes.items():
    markers_found = []
    for marker in markers:
        if marker in adata.var.index:
            markers_found.append(marker)
    marker_genes_in_data[ct] = markers_found
sc.pl.umap(adata,color=marker_genes_in_data["Bn"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Bsubset_BnMarkers.png")
sc.pl.umap(adata,color=marker_genes_in_data["Bm"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Bsubset_BmMarkers.png")
sc.pl.umap(adata,color=marker_genes_in_data["GCB"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Bsubset_GCBMarkers.png")
sc.pl.umap(adata,color=marker_genes_in_data["activated B"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Bsubset_ACBMarkers.png")
sc.pl.umap(adata,color=marker_genes_in_data["Transitional B"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Bsubset_Transitional_BMarkers.png")
sc.pl.umap(adata,color=marker_genes_in_data["Plasma cells"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Bsubset_PlasmMarkers.png")
sc.pl.umap(adata,color=marker_genes_in_data["Plasmablast"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Bsubset_Plasmblast.png")
```
* Bn markers
<img src="..\figures\umapBsubset_BnMarkers.png">
都不赖
* Bm markers
<img src="..\figures\umapBsubset_BmMarkers.png">
IGHG1不对(生信加油站真坑啊),其它也都还行
* GCB markers
<img src="..\figures\umapBsubset_GCBMarkers.png">
RGS13,IL21R,MKI67,STMN1,HMGB2,TOP2A,TCL1A取这几个
* ACB markers
<img src="..\figures\umapBsubset_ACBMarkers.png">
感觉不需要这两个marker
* Transitional B markers
<img src="..\figures\umapBsubset_Transitional_BMarkers.png">
这些marker也不需要,里面应该是没有Transitional B
* Plasma cells markers
<img src="..\figures\umapBsubset_PlasmMarkers.png">
IRF4不要,其它都还行
* Plasmablast markers
<img src="..\figures\umapBsubset_Plasmblast.png">
和Plasm重了的可以不要

### 优化marker后再绘制点图
```
marker_genes_in_data = {
"Bn":["TCL1A",
      "YBX3",
      "FCER2",
      "IL4R",
      "MS4A1",
      "FCRL1",
      "CD72",
      "BACH2",
      "IGHD",
      "IGHM"] ,
"Bm":["CD27",
      "AIM2",
      "TNFRSF13B",
      "CRIP2",
      "ITGB1"],
"GCB":["RGS13",
       "IL21R",
       "MKI67",
       "STMN1",
       "HMGB2",
       "TOP2A",
       "TCL1A",],
"Plasma cells": ["MZB1",
                "HSP90B1",
                "FNDC3B",
                "PRDM1",
                "IGKC",
                "JCHAIN",
                "CD38",
                "XBP1"],
"Plasmablast": ["XBP1",
                "PRDM1",
                "PAX5",
                "MZB1",
                "JCHAIN",
                "STMN1",
                "HMGB2",
                "MKI67",
                "SDC1"],
}
sc.pl.dotplot(
    adata,
    groupby="leiden_res0_8",
    var_names=marker_genes_in_data,
    standard_scale="var",
    save="Bresubset.png"
)
```
<img src="..\figures\umapB_umapmd0.4_cluster.png">
<img src="..\figures\dotplot_Bresubset.png">

结合umap和点图进行注释
```
cl_annotation = {
    "0": "Bm",
    "1": "Bn",
    "2": "B Plasma",
    "3": "B Plasma",
    "5": "B Plasma",
    "6": "B Plasma",
    "7": "B Plasma",
    "8": "B Plasma",
    "9": "B Plasma",
    "10": "B Plasma",
    "12": "B Plasma",
    "14":"Bn",
    "15": "B Plasma",
    "16": "Bgc"
}
```
再添上特征基因
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
  save="B_top10gene.png"
)
adata = sc.read("RData/B__umapmd0.4.h5ad")
adata = adata[adata.obs.leiden_res0_8 != '4' ].copy()
adata = adata[adata.obs.leiden_res0_8 != '13' ].copy()
cl_annotation = {
    "0": "Bm_TNFRSF13B",
    "1": "Bn_FCER2",
    "2": "Bp_MZB1",
    "3": "Bp_MZB1",
    "5": "Bp_MZB1",
    "6": "Bp_MZB1",
    "7": "Bp_MZB1",
    "8": "Bp_MZB1",
    "9": "Bp_MZB1",
    "10": "Bp_MZB1",
    "11":"Bp_MZB1",
    "12": "Bp_MZB1",
    "14":"Bn_FCER2",
    "15": "Bp_MZB1",
    "16": "Bgc_HMGB2"
}
adata.obs["final_anno_celltype"]=adata.obs.leiden_res0_8.map(cl_annotation)
###umap
sc.pl.umap(adata,color="final_anno_celltype",save="Bsubset_anno.png")
###dotplot
marker_genes_in_data = {
"Bn_FCER2":[
      "TCL1A",
      "YBX3",
      "FCER2",
      "IL4R",
      "MS4A1",
      "FCRL1",
      "CD72",
      "BACH2",
      "IGHD",
      "IGHM"] ,
"Bm_TNFRSF13B":["CD27",
      "AIM2",
      "TNFRSF13B",
      "CRIP2",
      "ITGB1"],
"Bgc_HMGB2":["RGS13",
       "IL21R",
       "MKI67",
       "STMN1",
       "HMGB2",
       "TOP2A",
       "TCL1A",],
"Bp_MZB1": ["MZB1",
                "HSP90B1",
                "FNDC3B",
                "PRDM1",
                "IGKC",
                "JCHAIN",
                "CD38",
                "XBP1"],
}
sc.pl.dotplot(
    adata,
    groupby="final_anno_celltype",
    var_names=marker_genes_in_data,
    standard_scale="var",
    save="Bsubset.png"
)
adata.write("RData/B_final_anno.h5ad")
```
<img src="..\figures\dotplot_B_top10gene.png">

* 最终umap
<img src="..\figures\umapBsubset_anno.png">
* 最终dotplot
<img src="..\figures\dotplot_Bsubset.png">

## 最后用celltypist验证一下结果
<img src="..\figures\umapcelltypist_coarse_B.png">
<img src="..\figures\umapcelltypist_fine_B.png">
可以的
