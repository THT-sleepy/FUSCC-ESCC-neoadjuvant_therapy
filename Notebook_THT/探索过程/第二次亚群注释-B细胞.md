# B细胞注释过程及结果记录         (What)

* Nov 17, 2025                                 (When)
* liyi /data1/liyi/zhaoyue/FUSCC-ESCC-neoadjuvant_thrapy           (Where)
* 单细胞需要注释后才能进行后续分析 (Why)

### 先看看聚类和celltypist注释情况
<img src="..\figures\umapB cell_umapmd0.4_regressmtribo_cluster.png">
<img src="..\figures\umap1117_celltypist_coarse_B.png">
<img src="..\figures\umap1117_celltypist_fine_B.png">

B细胞除了注释外，有一个问题是为什么浆细胞聚类出来这么奇怪？

## 1 把min_dist调小一点试一下？


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
