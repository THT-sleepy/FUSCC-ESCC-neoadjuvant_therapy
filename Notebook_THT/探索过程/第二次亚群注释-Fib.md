# 纤维细胞注释过程及结果记录         (What)

* Nov 22, 2025                                 (When)
* liyi /data1/liyi/zhaoyue/FUSCC-ESCC-neoadjuvant_thrapy           (Where)
* 单细胞需要注释后才能进行后续分析 (Why)

#### 1 先看看聚类和celltypist注释情况
<img src="..\figures\umapFibrocyte_umapmd0.4_regressmtribo_cluster.png">
<img src="..\figures\umap1122_celltypist__Fib.png">
居然有一群Glia细胞

#### 2 大类marker排除双细胞
```
import scanpy as sc
sc.settings.set_figure_params(
    dpi=300,
    facecolor="white",
    frameon=False
)
sc.settings.figdir = "plots"

adata = sc.read("RData/Fibrocyte__umapmd0.4.h5ad")

marker_genes = ["PTPRC","CD3E","CD79A","COL1A1",
"VWF","PECAM1","EPCAM","SFN","KRT5","LYZ","GCA","TPSAB1"]
sc.pl.umap(adata,color=marker_genes,vmin=0,vmax="p99",sort_order=False,cmap="Blues",save="1122_121_Fib_majormarker.png")
```
<img src="..\figures\umap1122_121_Fib_majormarker.png">
还是比较干净的

#### 3 筛marker
```
marker_genes = {
    "Pericyte": ["RGS5", "MCAM", "ACTA2", "MYH11", "PDGFRB", "CD36", "NOTCH3"],
    "VSMC": ["TAGLN", "ACTA2", "ACTG2", "MYH11", "CNN1", "HHIP"],
    "iFib": ["CXCL1", "CXCL2", "IL6", "CEBPD", "CLU", "CTGF", "HGF", "HSPA6", "DNAJB1", "MYC", "AFT4"],
    "apFib": ["CXCL12", "CD74", "HLA-DRB1", "HLA-DRA"],
    "myoFib": ["HOPX", "IGFBP5", "TIMP1", "MMP11", "COL10A1", "POSTN", "LRRC15", "SFRP4", "SFRP2", "COMP", "MMP1", "COL7A1", "WNT5A", "ISG15", "IL7R"],
    "prfFib": ["RBP1", "STAR", "STMN1"],
    "proFib": ["MFAP5", "PI16", "CD34", "APOD", "FGF7", "COL15A1"],
    "Glia": ["SOX10", "S100B", "PLP", "GFAP"]
}
marker_genes_in_data = {}
for ct, markers in marker_genes.items():
    markers_found = []
    for marker in markers:
        if marker in adata.var.index:
            markers_found.append(marker)
    marker_genes_in_data[ct] = markers_found

sc.pl.umap(adata,color=marker_genes_in_data["Pericyte"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="1122_121Fib_Pericytemarkers.png")
sc.pl.umap(adata,color=marker_genes_in_data["VSMC"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="1122_121Fib_VSMCmarkers.png")
sc.pl.umap(adata,color=marker_genes_in_data["iFib"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="1122_121Fib_iFibmarkers.png")
sc.pl.umap(adata,color=marker_genes_in_data["apFib"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="1122_121Fib_apFibmarkers.png")
sc.pl.umap(adata,color=marker_genes_in_data["myoFib"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="1122_121Fib_myoFibmarkers.png")
sc.pl.umap(adata,color=marker_genes_in_data["prfFib"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="1122_121Fib_prfFibmarkers.png")
sc.pl.umap(adata,color=marker_genes_in_data["proFib"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="1122_121Fib_proFibmarkers.png")
sc.pl.umap(adata,color=marker_genes_in_data["Glia"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="1122_121Fib_Gliamarkers.png")
```
<img src="..\figures\umapFibrocyte_umapmd0.4_regressmtribo_cluster.png">
右边那一团应该是周细胞
<img src="..\figures\umap1122_121Fib_Pericytemarkers.png">
左上是VSMC？
<img src="..\figures\umap1122_121Fib_VSMCmarkers.png">
不太特异？
<img src="..\figures\umap1122_121Fib_iFibmarkers.png">
右下角这一小搓确实应该是apFib
<img src="..\figures\umap1122_121Fib_apFibmarkers.png">
不太能看出来是不是有这一群
<img src="..\figures\umap1122_121Fib_myoFibmarkers.png">
没有prfFib
<img src="..\figures\umap1122_121Fib_prfFibmarkers.png">
？左下这一群好像有点特异性，第0群是proFib
<img src="..\figures\umap1122_121Fib_proFibmarkers.png">
只有一个S100B不能说明问题呀
<img src="..\figures\umap1122_121Fib_Gliamarkers.png">

res0.6,0是proFib,4和5是周细胞，6是apFib,其它不知道

####4 看下DEG
```
sc.tl.rank_genes_groups(
 adata, groupby="leiden_res0_6", method="wilcoxon", key_added="dea_leiden_res0_6"
)
sc.tl.filter_rank_genes_groups(
  adata,
  min_in_group_fraction=0.2,
  max_out_group_fraction=0.2,
  key="dea_leiden_res0_6",
  key_added="dea_leiden_res0_6_filtered",
 )
sc.pl.rank_genes_groups_dotplot(
  adata,
  groupby="leiden_res0_6",
  standard_scale="var",
  n_genes=10,
  key="dea_leiden_res0_6_filtered",
  save="regressmtribo_121_res0_6_Fib_top10gene.png"
)
```
<img src="..\figures\dotplot_regressmtribo_121_res0_6_Fib_top10gene.png">

####5 添加注释
0:proFib,"MFAP5", "PI16", "CD34";CD34<br>;
1:Fib, "COL11A1","COL10A1","SDC1","PLPP4";COL11A1<br>
2:Fib,"MMP1","CHI3L1","CXCL6","IL24";MMP1<br>
3:Fib,"PTGDS","CXCL14","COL6A5","ITGA8";COL6A5<br>
4:Pericyte:"RGS5", "MCAM", "ACTA2","CCDC102B","TRPC6","ARHGDIB";RGS5<br>
5:Pericyte:"MYH11","LMOD1","BCAM";MYH11<>
6:apFib,"CD74", "HLA-DRB1", "HLA-DRA";HLA-DRB1<br>

```
adata = sc.read("RData/Fibrocyte__umapmd0.4.h5ad")
cl_annotation = {
        "0": "c47_proFib_CD34",
        "1": "c44_Fib_COL11A1",
        "2": "c45_Fib_MMP1",
        "3": "c46_Fib_COL6A5",
        "4": "c49_PC_RGS5",
        "5": "c50_PC_MYH11",
        "6": "c48_apFib_HLA-DRB1",
    }
adata.obs["minor_celltype"]=adata.obs.leiden_res0_6.map(cl_annotation)
adata.write("RData/Fib_annotated.h5ad")
```

#### 6 检查每个cluster患者数
```
adata.obs['Patient_ID'] = adata.obs['sample'].str.split('_').str[-1]
adata.obs.groupby('leiden_res0_6')['Patient_ID'].nunique()
```
<img src="..\figures\Fib_患者亚群数.png">

#### 7 umap sample_site
```
sc.settings.set_figure_params(
    dpi=300,
    facecolor="white",
    frameon=True
)
def determine_sample_site(group_name):
    # 检查 group_name 是否是字符串并且包含 'pbmc'
    if isinstance(group_name, str) and 'pbmc' in group_name.lower():
        return 'blood'
    else:
        return 'tumor'

# 使用 .apply() 将函数应用到 'sample_group' 列的每一个元素
adata.obs['sample_site'] = adata.obs['sample_group'].apply(determine_sample_site)
palette = {
     "blood":"#3A6FA1",
     "tumor":"#BA3931",
}
sc.pl.umap(adata,color="sample_site",palette=palette,save="Fib_samplesite.png")
```
<img src="..\figures\umapFib_samplesite.png">

#### dotplot
```
marker_genes = [
"COL11A1","COL10A1","SDC1","PLPP4",
"MMP1","CHI3L1","CXCL6","IL24",
"PTGDS","CXCL14","COL6A5","ITGA8",
"MFAP5","PI16","CD34",
"CD74","HLA-DRB1","HLA-DRA",
"RGS5", "MCAM", "ACTA2","CCDC102B","TRPC6","ARHGDIB",
"MYH11","LMOD1","BCAM",
]
from matplotlib.colors import LinearSegmentedColormap
custom_cmap = LinearSegmentedColormap.from_list("custom", ["#0A3D7C", "#F6F6F6", "#810426"])
dp = sc.pl.dotplot(
        adata,
        groupby="minor_celltype",
        categories_order = [
             "c44_Fib_COL11A1",
             "c45_Fib_MMP1",
             "c46_Fib_COL6A5",
             "c47_proFib_CD34",
             "c48_apFib_HLA-DRB1",
             "c49_PC_RGS5",
             "c50_PC_MYH11",],
        var_names=marker_genes,
        standard_scale="var",
        return_fig=True,
        var_group_labels="       ",
        var_group_positions =[(0,3),(4,7),(8,11),
        (12,14),(15,17),(18,23),(24,26)],
        cmap = custom_cmap
    )
dp.style(dot_edge_color='black', dot_edge_lw=1).savefig("plots/1125_dotplot_Fib.png")
```
<img src="..\figures\1125_dotplot_Fib.png">

后面让傅师兄再弄10个颜色吧，豆包弄的颜色不太对
### umap clt
```
sc.settings.set_figure_params(
    dpi=300,
    facecolor="white",
    frameon=True
)
palette = {
  "c43": ,
  "c44": ,
  "c45": ,
  "c46": ,
  "c47": ,
  "c48": ,
  "c49": ,
}
cl_annotation = {
  "0": "c46",
  "1": "c43",
  "2": "c44",
  "3": "c45",
  "4": "c48",
  "5": "c49",
  "6": "c47",
    }
adata.obs["minor_clt"]=adata.obs.leiden_res1.map(cl_annotation)
sc.pl.umap(adata,color="minor_clt",legend_loc="on data",title="DC cells\n(n=7798)",palette=palette,save="1122_DC.png")
```
<img src="..\figures\umap1122_DC.png">
