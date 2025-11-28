# B细胞注释过程及结果记录         (What)

* Nov 17, 2025                                 (When)
* liyi /data1/liyi/zhaoyue/FUSCC-ESCC-neoadjuvant_thrapy           (Where)
* 单细胞需要注释后才能进行后续分析 (Why)

### 先看看聚类和celltypist注释情况
<img src="..\figures\umapB cell_umapmd1.2_regressmtribo_cluster.png">
<img src="..\figures\umap1117_celltypist_coarse_B.png">
<img src="..\figures\umap1117_celltypist_fine_B.png">

B细胞除了注释外，有一个问题是为什么浆细胞聚类出来这么奇怪？

### 先看看大类marker
```
import scanpy as sc

sc.settings.set_figure_params(
    dpi=300,
    facecolor="white",
    frameon=False
)
sc.settings.figdir = "plots"

adata = sc.read("RData/B cell__umapmd1.2.h5ad")
marker_genes = ["PTPRC","CD3E","CD79A","COL1A1",
"VWF","PECAM1","EPCAM","SFN","KRT5","LYZ","GCA","TPSAB1"]
sc.pl.umap(adata,color=marker_genes,vmin=0,vmax="p99",sort_order=False,cmap="Blues",save="121_B_majormarker.png")

```
<img src="..\figures\umap121_B_majormarker.png">
顶上的像是B混了上皮，下面的res0_8的14群有可能是混了中性粒细胞

### 确认下是不是混了中性粒细胞
```
marker_genes = ["CSF3R","FCGR3B","CXCR2"]
sc.pl.umap(adata,color=marker_genes,vmin=0,vmax="p99",sort_order=False,cmap="Blues",save="121_B_neutrolphilmarker.png")
```
确实混了

### 去掉混的细胞吧
```
adata = adata[adata.obs.leiden_res0_8 != '8'].copy()
adata = adata[adata.obs.leiden_res0_8 != '14'].copy()
```

### 筛一下经典marker
```
marker_genes = {
  "Bn": ["MS4A1", "IGHD", "FCER2", "TCL1A", "IL4R", "CD72", "BACH2", "IGHM", "YBX3"],
  "Bm": ["MS4A1", "CD27", "IGHM", "IGHG", "IGHA", "TNFRSF13B", "CRIP2", "ITGB1"],
  "Bgc": ["MS4A1", "BCL6", "TCL1A", "IL4R", "BACH2", "RGS13", "AICDA", "IL21R", "MKI67", "STMN1", "HMGB2", "TOP2A", "SUGCT", "NEIL1", "NEF2B", "LMO2", "GMDS", "MARCKSL1", "ACTG1", "ACTB"],
  "Plasma cell": ["CD38", "SDC1", "TNFRSF17", "MUM1", "MZB1", "XBP1", "PRDM1", "IRF4", "JCHAIN", "SSR4"],
  "Plasmablast": ["MS4A1", "CD27", "CD38", "MKI67", "MUM1", "MZB1", "XBP1", "JCHAIN", "STMN1", "HMGB2"]}
marker_genes_in_data = {}
for ct, markers in marker_genes.items():
      markers_found = []
      for marker in markers:
          if marker in adata.var.index:
              markers_found.append(marker)
      marker_genes_in_data[ct] = markers_found
sc.pl.umap(adata,color=marker_genes_in_data["Bn"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="121B_Bnmarkers.png")
sc.pl.umap(adata,color=marker_genes_in_data["Bm"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="121B_Bmmarkers.png")
sc.pl.umap(adata,color=marker_genes_in_data["Bgc"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="121B_Bgcmarkers.png")
sc.pl.umap(adata,color=marker_genes_in_data["Plasma cell"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="121B_plasmamarkers.png")
sc.pl.umap(adata,color=marker_genes_in_data["Plasmablast"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="121B_plasmablastmarkers.png")
```
"MS4A1","TCL1A","FCER2","IGHD","IGHM","BACH2",Bn对应res0.8的第1群
<img src="..\figures\umap121B_Bnmarkers.png">
TNFRSF13B,Bm对应0和2群
<img src="..\figures\umap121B_Bmmarkers.png">
之前能分出Bgc，现在反而不行了
<img src="..\figures\umap121B_Bgcmarkers.png">
显然右边的都是浆细胞
<img src="..\figures\umap121B_plasmamarkers.png">
分不出浆母细胞
<img src="..\figures\umap121B_plasmablastmarkers.png">

### 再结合DEG挑一下marker
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
  save="1118_B_top10gene.png"
)

```
<img src="..\figures\dotplot_1118_B_top10gene.png">

### 添加最终注释
0:Bm markergene是"TUBA1A",
1:Bn markergene是"MS4A1","TCL1A","FCER2","IGHD","IGHM","BACH2"
2:Bm markergenes是"CHORDC1","PPP1R16B","CIITA","TNFRSF13B,"CXCR5",
其它:Plasma cell markergenes是"MZB1","XBP1"


```
adata = sc.read("RData/B cell__umapmd1.2.h5ad")
adata = adata[adata.obs.leiden_res0_8 != '8'].copy()
adata = adata[adata.obs.leiden_res0_8 != '14'].copy()
cl_annotation = {
        "0": "c21_Bm_TUBA1A",
        "1": "c19_Bn_TCL1A",
        "2": "c20_Bm_TNFRSF13B",
        "3": "c22_Bp_MZB1",
        "4": "c22_Bp_MZB1",
        "5": "c22_Bp_MZB1",
        "6": "c22_Bp_MZB1",
        "7": "c22_Bp_MZB1",
        "9": "c22_Bp_MZB1",
        "10": "c22_Bp_MZB1",
        "11": "c22_Bp_MZB1",
        "12": "c22_Bp_MZB1",
        "13": "c22_Bp_MZB1",
        "15": "c22_Bp_MZB1"
    }
adata.obs["minor_celltype"]=adata.obs.leiden_res0_8.map(cl_annotation)
adata.write("RData/B_annotated.h5ad")
```
### umap clt

```
sc.settings.set_figure_params(
    dpi=300,
    facecolor="white",
    frameon=True
)
palette = {
     "c19":"#5BA0D4",
     "c20":"#5AA452",
     "c21":"#F29696",
     "c22":"#D8432E"
}
cl_annotation = {
        "0": "c21",
        "1": "c19",
        "2": "c20",
        "3": "c22",
        "4": "c22",
        "5": "c22",
        "6": "c22",
        "7": "c22",
        "9": "c22",
        "10": "c22",
        "11": "c22",
        "12": "c22",
        "13": "c22",
        "15": "c22"
    }
adata.obs["minor_clt"]=adata.obs.leiden_res0_8.map(cl_annotation)
sc.pl.umap(adata,color="minor_clt",legend_loc="on data",title="B cells\n(n=43183)",palette=palette,save="1125_B.png")
```
<img src="..\figures\umapB.png">
<img src="..\figures\umap1125_B.png">

### umap sample_site
```
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
sc.pl.umap(adata,color="sample_site",palette=palette,save="B_samplesite.png")
```
<img src="..\figures\umapB_samplesite.png">

### 检查每个cluster患者数
```
adata.obs['Patient_ID'] = adata.obs['sample'].str.split('_').str[-1]
patient_count_by_celltype = adata.obs.groupby('minor_celltype')['Patient_ID'].nunique()
```
<img src="..\figures\B_亚群患者数.png">


### dotplot
```
marker_genes = ["MS4A1","TCL1A","FCER2","IGHD","IGHM","BACH2",
"CHORDC1","PPP1R16B","CIITA","TNFRSF13B","CXCR5",
"TUBA1A",
"MZB1","XBP1"
]

from matplotlib.colors import LinearSegmentedColormap
custom_cmap = LinearSegmentedColormap.from_list("custom", ["#0A3D7C", "#F6F6F6", "#810426"])
dp = sc.pl.dotplot(
        adata,
        groupby="minor_celltype",
        categories_order = [
          "c19_Bn_TCL1A",
          "c20_Bm_TNFRSF13B",
          "c21_Bm_TUBA1A",
          "c22_Bp_MZB1",],
        var_names=marker_genes,
        standard_scale="var",
        return_fig=True,
        var_group_labels="    ",
        var_group_positions =[(0,5),(6,10),(11,11),(12,13)],
        cmap = custom_cmap
    )
dp.style(dot_edge_color='black', dot_edge_lw=1).savefig("plots/1125_dotplot_B.png")

```
<img src="..\figures\1125_dotplot_B.png">
