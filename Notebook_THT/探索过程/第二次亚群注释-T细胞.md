# T细胞注释过程及结果记录         (What)

* Nov 11, 2025                                 (When)
* liyi /data1/liyi/zhaoyue/FUSCC-ESCC-neoadjuvant_thrapy           (Where)
* 单细胞需要注释后才能进行后续分析 (Why)

## 1 再次检查有无混杂的其它细胞系如上皮的marker表达
* 聚类情况
<img src="..\figures\umapT&ILC cell_umapmd0.4_cluster.png">

```
import scanpy as sc

sc.settings.set_figure_params(
    dpi=300,
    facecolor="white",
    frameon=False
)
sc.settings.figdir = "plots"

adata = sc.read("RData/T&ILC cell__umapmd0.4.h5ad")

marker_genes = ["PTPRC","CD3E","CD79A","COL1A1",
"VWF","PECAM1","EPCAM","SFN","KRT5","LYZ","GCA","TPSAB1"]
sc.pl.umap(adata,color=marker_genes,vmin=0,vmax="p99",sort_order=False,cmap="Blues",save="121_T_majormarker.png")
```
<img src="..\figures\umap121_T_majormarker.png">
res0_8的第12群显然是T混了上皮,可以直接去掉。第10群有可能是混了中性粒，需要再看看中性粒的表达情况
```
marker_genes = ["CSF3R","FCGR3B","CXCR2"]
sc.pl.umap(adata,color=marker_genes,vmin=0,vmax="p99",sort_order=False,cmap="Blues",save="121_T_neutrolphilmarker.png")
```
<img src="..\figures\umap121_T_neutrolphilmarker.png">
确实是混了中性粒细胞<br>


去掉这两群
```
adata = adata[adata.obs.leiden_res0_8 != '10'].copy()
adata = adata[adata.obs.leiden_res0_8 != '12'].copy()
```

## 2 把ILC,CD4T,CD8T分出来再做亚群聚类
```
marker_genes = {
"T" :["CD4","CD8A","CD8B","CD3D","CD3E","CD3G"],
"NK":["FCGR3A","KLRC1","NCAM1"],
"Tprf":["MKI67","TOP2A","STMN1"],
}
marker_genes_in_data = {}
for ct, markers in marker_genes.items():
    markers_found = []
    for marker in markers:
        if marker in adata.var.index:
            markers_found.append(marker)
    marker_genes_in_data[ct] = markers_found
sc.pl.umap(adata,color=marker_genes_in_data["T"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="121Tsubset_T.png")
sc.pl.umap(adata,color=marker_genes_in_data["NK"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="121Tsubset_NK.png")
sc.pl.umap(adata,color=marker_genes_in_data["Tprf"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="121Tsubset_prf.png")
```
<img src="..\figures\umap121Tsubset_T.png">
<img src="..\figures\umap121Tsubset_NK.png">
<img src="..\figures\umap121Tsubset_prf.png">
<img src="..\figures\umapT&ILC cell_umapmd0.4_cluster.png">

res0_8是合适的
ct0,3,5,7是CD4+
ct11是Tprf
ct2,4,6,8是CD8+
ct1,9是ILC/NKT<br>
用点图再看看
```
sc.pl.dotplot(
    adata,
    groupby="leiden_res0_8",
    var_names=marker_genes_in_data,
    standard_scale="var",
    save="Tsubset_4_8_NKmarkers.png"
)
```
<img src="..\figures\dotplot_Tsubset_4_8_NKmarkers.png">

```
cl_annotation = {
    "0": "CD4+T",
    "1": "NKT_ILC",
    "2": "CD8+T",
    "3": "CD4+T",
    "4": "CD8+T",
    "5": "CD4+T",
    "6": "CD8+T",
    "7": "CD4+T",
    "8": "CD8+T",
    "9": "NKT_ILC",
    "11":"Tprf"
}
adata.obs["TILC_celltype"]=adata.obs.leiden_res0_8.map(cl_annotation)
adata.write("RData/T_1st_round_anno.h5ad")
###Tprf可以直接注释保存
adata_prf = adata[adata.obs.TILC_celltype == 'Tprf'].copy()
adata_prf.obs["final_anno_celltype"] = "Tprf"
adata_prf.write("RData/Tprf_finalanno.h5ad")
```

### 3 再跑亚群聚类
```
import scanpy as sc
import pandas as pd

sc.settings.verbosity = 0
sc.settings.set_figure_params(
    dpi=300,
    facecolor="white",
    frameon=False
)
sc.settings.figdir = "plots"
sc.settings.autoshow = False

adata = sc.read("RData/T_1st_round_anno.h5ad") #X是log1p_norm

def process_celltype(adata, celltype_name, save_dir="RData/", plot_prefix=""):
    """
    批量处理特定细胞类型的单细胞数据（归一化→PCA→BBKNN→UMAP→Leiden聚类）

    参数：
    adata: 原始AnnData对象（需包含first_round_celltype列和soupX_counts层）
    celltype_name: 目标细胞类型（如"T"、"B"、"Epi"）
    save_dir: h5ad文件保存目录（默认RData/）
    plot_prefix: 图片文件名前缀（默认空，可加项目标识）
    """
    # 1. 提取目标细胞类型子集
    adata_sub = adata[adata.obs.TILC_celltype == celltype_name].copy()
    print(f"开始处理 {celltype_name} 细胞，共 {adata_sub.n_obs} 个细胞")

    # 1.5 去掉细胞数小于10个的样本
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
    adata_sub.X = adata_sub.layers["log1p_norm"]  # 将log1p后的数据设为默认X

    # 3. PCA降维
    sc.pp.pca(adata_sub, svd_solver="arpack", use_highly_variable=True)

    # 4. 循环运行不同min_dist的UMAP+Leiden聚类（0.4）
    for min_dist in [0.4]:
        # 复制数据避免相互干扰
        adata_md = adata_sub.copy()

        # BBKNN批次校正 + UMAP降维
        sc.external.pp.bbknn(adata_md, batch_key="sample", n_pcs=30)
        sc.tl.umap(adata_md, min_dist=min_dist)

        # 多分辨率Leiden聚类
        sc.tl.leiden(adata_md, key_added="leiden_res0_6", resolution=0.6)
        sc.tl.leiden(adata_md, key_added="leiden_res0_8", resolution=0.8)
        sc.tl.leiden(adata_md, key_added="leiden_res1", resolution=1.0)

        # 保存UMAP聚类图
        plot_path = f"{plot_prefix}{celltype_name}_umapmd{min_dist}_cluster.png"
        sc.pl.umap(
            adata_md,
            color=["leiden_res0_6", "leiden_res0_8", "leiden_res1"],
            legend_loc="on data",
            save=plot_path,
            show=False  # 避免弹出图片窗口，仅保存
        )

        # 保存AnnData对象
        output_fpath = f"{save_dir}{celltype_name}__umapmd{min_dist}.h5ad"
        adata_md.write(output_fpath)
        print(f"  {celltype_name} (min_dist={min_dist}) 已保存：{output_fpath}")

    print(f"{celltype_name} 细胞处理完成！\n")
# 假设你的原始AnnData对象已加载为 'adata'
# 1. 定义要处理的细胞类型列表（从T开始，按需求调整顺序）
celltypes_to_process = ["CD4+T", "CD8+T", "NKT_ILC"]

# 2. 批量运行（可自定义保存目录和图片前缀）
for celltype in celltypes_to_process:
    process_celltype(
        adata=adata,
        celltype_name=celltype,
        save_dir="RData/",  # h5ad文件保存到RData目录
        plot_prefix=""  # 图片名前缀（如scRNA_T_umapmd0.1_cluster.png）
    )
```

## final anno CD4T
T细胞的注释比较特殊,不同细胞群会有相同的marker表达，比如Tcm是在Tn基础上表达一定的杀伤marker


### 先看看celltypist的结果
<img src="..\figures\umap121celltypist_coarse_CD4T.png">
<img src="..\figures\umap121celltypist_fine_CD4T.png">

### 再看看聚类情况
<img src="..\figures\umapCD4+T_umapmd0.4_cluster.png">

### 筛marker
```
import scanpy as sc
sc.settings.set_figure_params(
    dpi=300,
    facecolor="white",
    frameon=False
)
sc.settings.figdir = "plots"
adata = sc.read("RData/CD4+T__umapmd0.4.h5ad")
marker_genes = {
    "Tn": ["CCR7", "LEF1", "SELL", "TCF7", "MAL", "TXK", "CD27", "CD28", "S1PR1", "IL7R", "CD45RA", "CD38", "KLF2", "CHMP1B", "CYCS", "CD55", "AREG"],
    "Temra": ["CD45RA", "GZMH", "GZMB", "TBX21", "ASCL2", "CX3CR1", "KLRG1", "FCGR3A", "FGFBP2", "PRF1", "EOMES", "S1PR1", "S1PR5", "NKG7", "GNLY", "CTSW"],
    "Tcm": ["CCR7", "CCR4", "CXCR3", "CXCR4", "ANXA1", "GPR183", "LMNA", "CDKN1A", "GPR183", "TNFSF9", "GZMK", "SELL", "IL7R", "CD27", "CD28", "PRF1", "GZMA", "CCL5", "S1PR1", "PTGER2", "ICAM2", "ANXA2", "RGS1", "TCF7", "CD69"],
    "Tem": ["GNLY", "FOS", "JUN", "IL7R", "CXCR3", "CXCR4", "CXCR5", "GZMK", "GZMH", "GZMA", "CD74", "CD28", "CCL4", "CCL4L2", "CD44", "NR4A1", "TNFSF9"],
    "Trm": ["ZNF683", "CD52", "HOPX", "ID2", "CXCR6", "XCL1", "XCL2", "KLRC1", "GZMH", "GZMA", "CCR6", "HSPA1A", "HSPA1B", "CAPG", "CD6", "MYADM", "RORA", "NR4A1", "NR4A2", "NR4A3", "CD69", "ITGAE", "KLRB1", "PTGER4", "IL7R"],
    "Tex": ["PDCD1", "HAVCR2", "LAG3", "LAYN", "TIGIT", "ENTPD1", "CXCL13", "HLA-DR", "CTLA4", "TNFRSF18", "TOX", "IFNG", "MIR155HG", "TNFRSF9", "ITGAE"],
    "Treg": ["FOXP3"],
    "Th17": ["RORA", "RORC", "CCR6", "IL23R", "IL22", "IL17A", "IL17F", "IL26", "FURIN", "CTSH", "KLRB1", "CAPG", "ITGAE"],
    "Th1-like_Tfh": ["TOX","TOX2","IL21","CXCL13","GNG4","CD200","ZBED2","CCL3","CCL4","IFNG","GZMB","LAG3","HAVCR2","CTLA4","CXCR6","CXCR3","BHLHE40","GZMB","PDCD1","ICOS","IGFLR1","ITGAE","BCL6","ICA1","IL6ST","MAGEH1","BTLA","CD200"],
}
marker_genes_in_data = {}
for ct, markers in marker_genes.items():
    markers_found = []
    for marker in markers:
        if marker in adata.var.index:
            markers_found.append(marker)
    marker_genes_in_data[ct] = markers_found
sc.pl.umap(adata,color=marker_genes_in_data["Tn"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Tresubset_Tnmarkers.png")
sc.pl.umap(adata,color=marker_genes_in_data["Temra"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Tresubset_Temramarkers.png")
sc.pl.umap(adata,color=marker_genes_in_data["Tcm"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Tsubset_Tcmmarkers.png")
sc.pl.umap(adata,color=marker_genes_in_data["Tem"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Tsubset_Temmarkers.png")
sc.pl.umap(adata,color=marker_genes_in_data["Trm"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Tsubset_Trmmarkers.png")
sc.pl.umap(adata,color=marker_genes_in_data["Tex"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Tsubset_Texmarkers.png")
sc.pl.umap(adata,color=marker_genes_in_data["Treg"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Tsubset_Tregmarkers.png")
sc.pl.umap(adata,color=marker_genes_in_data["Th17"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Tsubset_Th17markers.png")
sc.pl.umap(adata,color=marker_genes_in_data["Th1-like_Tfh"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Tsubset_Th1-like_Tfhmarkers.png")
```
<img src="..\figures\umapTresubset_Tnmarkers.png">
CD38,CHMP1B,CYCS没有特异性
<img src="..\figures\umapTresubset_Temramarkers.png">
可能没有这一群
<img src="..\figures\umapTsubset_Tcmmarkers.png">
ANXA1,IL7R,SELL,GZMK这几个看起来比较对
<img src="..\figures\umapTsubset_Temmarkers.png">
FOS,JUN
<img src="..\figures\umapTsubset_Trmmarkers.png">
感觉是没有这一群
<img src="..\figures\umapTsubset_Texmarkers.png">
PDCD1,HAVCR2,LAG3,LAYN,TIGIT,ENTPD1,CXCL13,CTLA4,TNFRSF18,TNFRSF9,TOX
<img src="..\figures\umapTsubset_Tregmarkers.png">
<img src="..\figures\umapTsubset_Th17markers.png">
感觉没有这一群
<img src="..\figures\umapTsubset_Th1-like_Tfhmarkers.png">

### 单纯用umap很难分得清各个cluster到底是什么细胞,用点图尝试下
```
marker_genes = {
    "Tn": ["CCR7", "LEF1", "SELL", "TCF7", "MAL",],
    "Tcm": ["ANXA1","GZMK","IL7R",],
    "Tem": ["FOS", "JUN",],
    "Trm": ["ZNF683"],
    "Tex": ["PDCD1", "HAVCR2", "LAG3", "LAYN", "TIGIT", "ENTPD1", "CXCL13", "CTLA4", "TNFRSF18", "TOX",],
    "Treg": ["FOXP3"],
}
sc.pl.dotplot(
        adata,
        groupby="leiden_res1",
        var_names=marker_genes,
        standard_scale="var",
        save="CD4亚群.png"
    )
```
如果是选res=1的话，明确的cluster是0,3是Tem,1,6,8是Treg,2,7是Tcm,4,5是Tn,9可能是Tem
<img src="..\figures\dotplot_CD4亚群.png">



### 利用差异基因找marker
```
sc.tl.rank_genes_groups(
 adata, groupby="leiden_res1", method="wilcoxon", key_added="dea_leiden_res1"
)
sc.tl.filter_rank_genes_groups(
  adata,
  min_in_group_fraction=0.2,
  max_out_group_fraction=0.2,
  key="dea_leiden_res1",
  key_added="dea_leiden_res1_filtered",
 )
sc.pl.rank_genes_groups_dotplot(
  adata,
  groupby="leiden_res1",
  standard_scale="var",
  n_genes=10,
  key="dea_leiden_res1_filtered",
  save="CD4T_top10gene.png"
)

```
画完这个图发现第9群其实是T混了上皮
0群:NR4A1,FOSB,UBE2S
1,6群事实上分不开:TNFRSF18,FOXP3,TIGIT
2,7群事实上也分不开
3群好像也有问题

把res降到0.8重新弄下
```
marker_genes = {
    "Tn": ["CCR7", "LEF1", "SELL", "TCF7", "MAL",],
    "Tcm": ["ANXA1","GZMK","IL7R",],
    "Tem": ["FOS", "JUN",],
    "Trm": ["ZNF683"],
    "Tex": ["PDCD1", "HAVCR2", "LAG3", "LAYN", "TIGIT", "ENTPD1", "CXCL13", "CTLA4", "TNFRSF18", "TOX",],
    "Treg": ["FOXP3"],
}
sc.pl.dotplot(
        adata,
        groupby="leiden_res0_8",
        var_names=marker_genes,
        standard_scale="var",
        save="res0_8CD4亚群.png"
    )
```
<img src="..\figures\dotplot_res0_8CD4亚群.png">

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
  save="CD4T_top10gene.png"
)
```
<img src="..\figures\dotplot_CD4T_top10gene.png">

0群:Tcm,特征基因为ANXA1,TIMP1,PTGER2
1群:Tn,特征基因为TCF7,SELL,CCR7,TXK
2群：双细胞要去掉
3群:Treg,特征基因为FOXP3,TNFRSF18,TNFRSF4
4群:Tem,特征基因为NR4A1,FOS,JUN,
5群:Treg,特征基因为SHMT2,HLA-DPB1,HLA-DRB1
6群：双细胞要去掉

绘制放到S1最后的图
```
adata = sc.read("RData/CD4+T__umapmd0.4.h5ad")
adata = adata[adata.obs.leiden_res0_8 != '2'].copy()
adata = adata[adata.obs.leiden_res0_8 != '6'].copy()
cl_annotation = {
    "0": "c02_CD4_Tcm_ANXA1",
    "1": "c01_CD4_Tn_TCF7",
    "3": "c04_CD4_Treg_FOXP3",
    "4": "c03_CD4_Tem_NR4A1",
    "5": "c05_CD4_Treg_SHMT2",
}
adata.obs["minor_celltype"]=adata.obs.leiden_res0_8.map(cl_annotation)
adata.write("RData/CD4T_annotated.h5ad")
```

#### umap
```
sc.settings.set_figure_params(
    dpi=300,
    facecolor="white",
    frameon=True
)
palette = {
     "c01":"#4279AD",
     "c02":"#AFCCD3",
     "c03":"#A25936",
     "c04":"#6171A8",
     "c05":"#BD5BA5",
}
cl_annotation = {
    "0": "c02",
    "1": "c01",
    "3": "c04",
    "4": "c03",
    "5": "c05",
}
adata.obs["minor_clt"]=adata.obs.leiden_res0_8.map(cl_annotation)
sc.pl.umap(adata,color="minor_clt",legend_loc="on data",title="CD4+ T cells",palette=palette,save="CD4T.png")

```
<img src="..\figures\umapCD4T.png">

#### dotplot
```
marker_genes = ["TCF7","SELL","CCR7","TXK",
                "ANXA1","TIMP1","PTGER2",
                "NR4A1","FOS","JUN",
                "FOXP3","TNFRSF18","TNFRSF4",
                "SHMT2","HLA-DPB1","HLA-DRB1"]
from matplotlib.colors import LinearSegmentedColormap
custom_cmap = LinearSegmentedColormap.from_list("custom", ["#0A3D7C", "#F6F6F6", "#810426"])
dp = sc.pl.dotplot(
        adata,
        groupby="minor_celltype",
        categories_order = ["c01_CD4_Tn_TCF7",
                            "c02_CD4_Tcm_ANXA1",
                            "c03_CD4_Tem_NR4A1",
                            "c04_CD4_Treg_FOXP3",
                            "c05_CD4_Treg_SHMT2",],
        var_names=marker_genes,
        standard_scale="var",
        return_fig=True,
        var_group_labels="     ",
        var_group_positions = [(0,3),(4,6),(7,9),
        (10,12),(13,15)],
        cmap = custom_cmap
    )
dp.style(dot_edge_color='black', dot_edge_lw=1).savefig("plots/dotplot_CD4T.png")
dp.add_totals().style(dot_edge_color='black', dot_edge_lw=1).savefig("plots/dotplot_CD4T_withtotal.png")

```
<img src="..\figures\dotplot_CD4T.png">
<img src="..\figures\dotplot_CD4T_withtotal.png">
* 每个亚群的患者数目
<img src="..\figures\CD4_亚群患者数.png">

### 再尝试下res0.1
前面发现如果是选res=1的话，明确的cluster是0是Tem,1,6,8是Treg,2,7是Tcm,4,5是Tn,3和9是T混了上皮


```
sc.tl.rank_genes_groups(
 adata, groupby="leiden_res1", method="wilcoxon", key_added="dea_leiden_res1"
)
sc.tl.filter_rank_genes_groups(
  adata,
  min_in_group_fraction=0.2,
  max_out_group_fraction=0.2,
  key="dea_leiden_res1",
  key_added="dea_leiden_res1_filtered",
 )
sc.pl.rank_genes_groups_dotplot(
  adata,
  groupby="leiden_res1",
  standard_scale="var",
  n_genes=10,
  key="dea_leiden_res1_filtered",
  save="res1_CD4T_top10gene.png"
)
```
<img src="..\figures\dotplot_CD4亚群.png">
<img src="..\figures\umapCD4+T_umapmd0.4_cluster.png">
<img src="..\figures\dotplot_res1_CD4T_top10gene.png">
0群：Tem，特征基因是NR4A1,FOS,JUN
1群：Treg，特征基因是FOXP3,TNFRSF18
6群：Treg，特征基因是FOXP3,TNFRSF18,FOSB,NR4A3
8群：Treg，特征基因是FOXP3,SHMT2,HLA-DPB1,HLA-DRB1
2群:Tcm,特征基因是ANXA1,TIMP1
7群:Tcm，特征基因是ANXA1,TIMP1,SESN3
4群:Tn,特征基因是CCR7,SELL,TCF,LEF1,FOXJ3
5群:Tn,特征基因是CCR7,SELL,TCF,LEF1


第1群和第6群经典marker和top10基因都分不开，单独做一个DEG看下
```
mask = adata.obs['leiden_res1'].isin(['1', '6'])
adata_filtered = adata[mask, :]
sc.tl.rank_genes_groups(
 adata_filtered, groupby="leiden_res1", method="wilcoxon", key_added="dea_leiden_res1"
)
sc.tl.filter_rank_genes_groups(
  adata_filtered,
  min_in_group_fraction=0.2,
  max_out_group_fraction=0.2,
  key="dea_leiden_res1",
  key_added="dea_leiden_res1_filtered",
 )
sc.pl.rank_genes_groups_dotplot(
  adata_filtered,
  groupby="leiden_res1",
  standard_scale="var",
  n_genes=10,
  key="dea_leiden_res1_filtered",
  save="Treg_CD4T_top10gene.png"
)
```
<img src="..\figures\dotplot_Treg_CD4T_top10gene.png">

第2群和第7群也是的
```
mask = adata.obs['leiden_res1'].isin(['2','7'])
adata_filtered = adata[mask, :]
sc.tl.rank_genes_groups(
 adata_filtered, groupby="leiden_res1", method="wilcoxon", key_added="dea_leiden_res1"
)
sc.tl.filter_rank_genes_groups(
  adata_filtered,
  min_in_group_fraction=0.2,
  max_out_group_fraction=0.2,
  key="dea_leiden_res1",
  key_added="dea_leiden_res1_filtered",
 )
sc.pl.rank_genes_groups_dotplot(
  adata_filtered,
  groupby="leiden_res1",
  standard_scale="var",
  n_genes=10,
  key="dea_leiden_res1_filtered",
  save="Tcm_CD4T_top10gene.png"
)
```
第4群和第5群也是的
```
mask = adata.obs['leiden_res1'].isin(['4','5'])
adata_filtered = adata[mask, :]
sc.tl.rank_genes_groups(
 adata_filtered, groupby="leiden_res1", method="wilcoxon", key_added="dea_leiden_res1"
)
sc.tl.filter_rank_genes_groups(
  adata_filtered,
  min_in_group_fraction=0.2,
  max_out_group_fraction=0.2,
  key="dea_leiden_res1",
  key_added="dea_leiden_res1_filtered",
 )
sc.pl.rank_genes_groups_dotplot(
  adata_filtered,
  groupby="leiden_res1",
  standard_scale="var",
  n_genes=10,
  key="dea_leiden_res1_filtered",
  save="Tn_CD4T_top10gene.png"
)
```
<img src="..\figures\dotplot_Tcm_CD4T_top10gene.png">
<img src="..\figures\dotplot_Tn_CD4T_top10gene.png">


```
adata = sc.read("RData/CD4+T__umapmd0.4.h5ad")
adata = adata[adata.obs.leiden_res1 != '3'].copy()
adata = adata[adata.obs.leiden_res1 != '9'].copy()
cl_annotation = {
    "0": "c02_CD4_Tcm_ANXA1",
    "1": "c01_CD4_Tn_TCF7",
    "3": "c04_CD4_Treg_FOXP3",
    "4": "c03_CD4_Tem_NR4A1",
    "5": "c05_CD4_Treg_SHMT2",
}
adata.obs["minor_celltype"]=adata.obs.leiden_res0_8.map(cl_annotation)
```
<img src="..\figures\CD4_res1亚群患者数.png">

### dotplot
```
marker_genes = ["TCF7","SELL","CCR7","TXK",
                "ANXA1","TIMP1","PTGER2",
                "NR4A1","FOS","JUN",
                "FOXP3","TNFRSF18","TNFRSF4",
                "SHMT2","HLA-DPB1","HLA-DRB1"]
from matplotlib.colors import LinearSegmentedColormap
custom_cmap = LinearSegmentedColormap.from_list("custom", ["#0A3D7C", "#F6F6F6", "#810426"])
dp = sc.pl.dotplot(
        adata,
        groupby="minor_celltype",
        categories_order = ["c01_CD4_Tn_TCF7",
                            "c02_CD4_Tcm_ANXA1",
                            "c03_CD4_Tem_NR4A1",
                            "c04_CD4_Treg_FOXP3",
                            "c05_CD4_Treg_SHMT2",],
        var_names=marker_genes,
        standard_scale="var",
        return_fig=True,
        var_group_labels="     ",
        var_group_positions = [(0,3),(4,6),(7,9),
        (10,12),(13,15)],
        cmap = custom_cmap
    )
dp.style(dot_edge_color='black', dot_edge_lw=1).savefig("plots/dotplot_CD4T.png")
dp.add_totals().style(dot_edge_color='black', dot_edge_lw=1).savefig("plots/dotplot_CD4T_withtotal.png")
```

事实上有些群就是分不开，没必要再细分了，就res0.8是对的

## final anno CD8T
### 先看看celltypist的结果
<img src="..\figures\umapcelltypist_coarse_CD8T.png">
<img src="..\figures\umapcelltypist_fine_CD8T.png">

### 再看看聚类情况
<img src="..\figures\umapCD8+T_umapmd0.4_cluster.png">

### 筛marker
```
import scanpy as sc
sc.settings.set_figure_params(
    dpi=300,
    facecolor="white",
    frameon=False
)
sc.settings.figdir = "plots"
adata = sc.read("RData/CD8+T__umapmd0.4.h5ad")
marker_genes = {
    "Tn": ["CCR7", "LEF1", "SELL", "TCF7", "MAL", "TXK", "CD27", "CD28", "S1PR1", "IL7R", "CD45RA", "CD38", "KLF2", "CHMP1B", "CYCS", "CD55", "AREG"],
    "Temra": ["CD45RA", "GZMH", "GZMB", "TBX21", "ASCL2", "CX3CR1", "KLRG1", "FCGR3A", "FGFBP2", "PRF1", "EOMES", "S1PR1", "S1PR5", "NKG7", "GNLY", "CTSW"],
    "Tcm": ["CCR7", "CCR4", "CXCR3", "CXCR4", "ANXA1", "GPR183", "LMNA", "CDKN1A", "GPR183", "TNFSF9", "GZMK", "SELL", "IL7R", "CD27", "CD28", "PRF1", "GZMA", "CCL5", "S1PR1", "PTGER2", "ICAM2", "ANXA2", "RGS1", "TCF7", "CD69"],
    "Tem": ["GNLY", "FOS", "JUN", "IL7R", "CXCR3", "CXCR4", "CXCR5", "GZMK", "GZMH", "GZMA", "CD74", "CD28", "CCL4", "CCL4L2", "CD44", "NR4A1", "TNFSF9"],
    "Trm": ["ZNF683", "CD52", "HOPX", "ID2", "CXCR6", "XCL1", "XCL2", "KLRC1", "GZMH", "GZMA", "CCR6", "HSPA1A", "HSPA1B", "CAPG", "CD6", "MYADM", "RORA", "NR4A1", "NR4A2", "NR4A3", "CD69", "ITGAE", "KLRB1", "PTGER4", "IL7R"],
    "Tex": ["PDCD1", "HAVCR2", "LAG3", "LAYN", "TIGIT", "ENTPD1", "CXCL13", "HLA-DR", "CTLA4", "TNFRSF18", "TOX", "IFNG", "MIR155HG", "TNFRSF9", "ITGAE"],
    "Tc17":["SLC4A10","KLRB1","TMIGD2","RORA","RORC","ZBTB16","IL26","IL17A","IL23R"],
    "MAIT":["SLC4A10","KLRB1","TMIGD2","IL7R","NCR3","LTB","CCL20"],
    "T_NK-like":["KLRG1","KLRD1","TYROBP",
                "KIR2DL1","KIR2DL3","KIR3DL1",
                "KIR3DL2","CD160","EOMES","TXK",
                "KLRC1","KIR2DL4"]
}
marker_genes_in_data = {}
for ct, markers in marker_genes.items():
    markers_found = []
    for marker in markers:
        if marker in adata.var.index:
            markers_found.append(marker)
    marker_genes_in_data[ct] = markers_found
sc.pl.umap(adata,color=marker_genes_in_data["Tn"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="CD8Tresubset_Tnmarkers.png")
sc.pl.umap(adata,color=marker_genes_in_data["Temra"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="CD8Tresubset_Temramarkers.png")
sc.pl.umap(adata,color=marker_genes_in_data["Tcm"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="CD8Tresubset_Tcmmarkers.png")
sc.pl.umap(adata,color=marker_genes_in_data["Tem"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="CD8Tresubset_Temmarkers.png")
sc.pl.umap(adata,color=marker_genes_in_data["Trm"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="CD8Tresubset_Trmmarkers.png")
sc.pl.umap(adata,color=marker_genes_in_data["Tex"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="CD8Tresubset_Texmarkers.png")
sc.pl.umap(adata,color=marker_genes_in_data["Tc17"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="CD8Tresubset_Tc17markers.png")
sc.pl.umap(adata,color=marker_genes_in_data["MAIT"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="CD8Tresubset_MAITmarkers.png")
sc.pl.umap(adata,color=marker_genes_in_data["T_NK-like"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="CD8Tresubset_T_NKliemarkers.png")
```
"CCR7", "LEF1", "SELL", "TCF7"
<img src="..\figures\umapCD8Tresubset_Tnmarkers.png">
"CX3CR1","KLRG1","TBX21","FGFBP2","GNLY","GZMH","GZMB","ASCL2","S1PR1","S1PR5"
<img src="..\figures\umapCD8Tresubset_Temramarkers.png">
感觉是没有这一群
<img src="..\figures\umapCD8Tresubset_Tcmmarkers.png">
"FOS","JUN","CXCR4"
<img src="..\figures\umapCD8Tresubset_Temmarkers.png">
"CD69",可以注意下HSPA1A和HSPA1B
<img src="..\figures\umapCD8Tresubset_Trmmarkers.png">
"PDCD1","HAVCR2","LAG3","LAYN","TIGIT","ENTPD1","CXCL13","CTLA4","TNFRSF18"
<img src="..\figures\umapCD8Tresubset_Texmarkers.png">
"SLC4A10","KLRB1","RORC","ZBTB16"
<img src="..\figures\umapCD8Tresubset_Tc17markers.png">
"SLC4A10","KLRB1","NCR3"
<img src="..\figures\umapCD8Tresubset_MAITmarkers.png">
应该是没有这一群
<img src="..\figures\umapCD8Tresubset_T_NKliemarkers.png">



这样看来分辨率1都小了一点，试试1.2，1.4,1.6,1.8
<img src="..\figures\umapCD8+T_umapmd0.4_cluster.png">
<img src="..\figures\umapCD8_cluster_res1_2to1_8.png">
我靠升高居然也分不出MAIT这一群吗？
先用res1分析看看吧

### 只用umap很难分得出各个亚群，结合点图再看看
```
marker_genes = {
    "Tn": ["CCR7", "LEF1", "SELL", "TCF7", "MAL","IL7R"],
    "Temra": ["CX3CR1","KLRG1","TBX21","FGFBP2","GNLY","GZMH","GZMB","ASCL2","S1PR1","S1PR5"],
    "Tm" :["CXCR4","ZFP36L2"],
    "Tem": ["FOS", "JUN"],
    "Trm": ["ZNF683","CXCR6","CD69"],
    "Tex": ["PDCD1", "HAVCR2", "LAG3", "LAYN", "TIGIT", "ENTPD1", "CXCL13", "CTLA4", "TNFRSF18", "TOX",],
    "MAIT": ["SLC4A10","KLRB1","RORC","ZBTB16","NCR3"],
}
sc.pl.dotplot(
        adata,
        groupby="leiden_res1",
        var_names=marker_genes,
        standard_scale="var",
        save="121_CD8亚群.png"
    )
```
<img src="..\figures\dotplot_121_CD8亚群.png">
第0群应该是Trm/Tem
第1群应该是Temra
第2群应该是Trm/Tem
第3群像是Tex
第4群像是Tex
第5群像是Tm,但没法归类到Tcm/Tem/Trm
第6群感觉是Tn混了MAIT

看看DEG的结果
```
sc.tl.rank_genes_groups(
 adata, groupby="leiden_res1", method="wilcoxon", key_added="dea_leiden_res1"
)
sc.tl.filter_rank_genes_groups(
  adata,
  min_in_group_fraction=0.2,
  max_out_group_fraction=0.2,
  key="dea_leiden_res1",
  key_added="dea_leiden_res1_filtered",
 )
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
sc.pl.rank_genes_groups_dotplot(
  adata,
  groupby="leiden_res1",
  standard_scale="var",
  n_genes=10,
  key="dea_leiden_res1_filtered",
  save="121_res1_CD8T_top10gene.png"
)
```
<img src="..\figures\dotplot_121_res1_CD8T_top10gene.png">
第3群是混了上皮细胞要去掉
其它的还是分的开的，试下把res升高到1.4
```
sc.pl.dotplot(
        adata,
        groupby="leiden_res1_4",
        var_names=marker_genes,
        standard_scale="var",
        save="res1_4_121_CD8亚群.png"
    )
```
<img src="..\figures\dotplot_res1_4_121_CD8亚群.png">

```
sc.tl.rank_genes_groups(
 adata, groupby="leiden_res1_4", method="wilcoxon", key_added="dea_leiden_res1_4"
)
sc.tl.filter_rank_genes_groups(
  adata,
  min_in_group_fraction=0.2,
  max_out_group_fraction=0.2,
  key="dea_leiden_res1_4",
  key_added="dea_leiden_res1_4_filtered",
 )
sc.pl.rank_genes_groups_dotplot(
  adata,
  groupby="leiden_res1_4",
  standard_scale="var",
  n_genes=10,
  key="dea_leiden_res1_4_filtered",
  save="121_res1_4_CD8T_top10gene.png"
)
```
<img src="..\figures\dotplot_121_res1_4_CD8T_top10gene.png">
1.4在去双细胞方面要优于1，把第10群去了应该就行
0,4,8群是分不开的;3，6，9群是分不开的，其它勉强是分得开的
看下res1_4的患者数是不是可以接受
可以接受<br>
<img src="..\figures\CD8_res1_4亚群患者数.png"><br>
分辨率再往上，0群分开两群，没有必要

0,4,8群：Temra,特征基因是CX3CR1,KLRG1?
1群: Tem/Trm,归为Tem，理由是Trm应该耗竭marker也要阳性，特征基因是FOS,JUN,FOSB,ERG1
2群: 和第一群类似，Tem，特征基因是FOS,JUN,ITM2C,TNFRSF9
3,6,9群：Tex，特征基因是CXCL13,PDCD1,LAG3
5群: Tem,特征基因是FOS,JUN,NR4A3,CRTAM
7群:Tn,特征基因是TCF7,LEF1,IL7R
11群Tem，特征基因是FOS,JUN,MT1X,MT1E,MT1F

### 把CD8T中可能混的上皮细胞利用打分给去掉了,再来重新看下
现在的聚类情况
<img src="..\figures\umapCD8+T_umapmd0.4filterepi_cluster.png">
之前的聚类情况
<img src="..\figures\umapCD8+T_umapmd0.4_cluster.png">
<img src="..\figures\umapCD8_cluster_res1_2to1_8.png">

感觉结果是要好一些，感觉之前的聚类有一些小cluster就是受到上皮基因的影响

```
import scanpy as sc
sc.settings.set_figure_params(
    dpi=300,
    facecolor="white",
    frameon=False
)
sc.settings.figdir = "plots"
adata = sc.read("RData/CD8+T__umapmd0.4_filterepi.h5ad")
marker_genes = {
    "Tn": ["CCR7", "LEF1", "SELL", "TCF7", "MAL","IL7R"],
    "Temra": ["CX3CR1","KLRG1","TBX21","FGFBP2","GNLY","GZMH","GZMB","ASCL2","S1PR1","S1PR5"],
    "Tm" :["CXCR4","ZFP36L2"],
    "Tem": ["FOS", "JUN"],
    "Trm": ["ZNF683","CXCR6","CD69"],
    "Tex": ["PDCD1", "HAVCR2", "LAG3", "LAYN", "TIGIT", "ENTPD1", "CXCL13", "CTLA4", "TNFRSF18", "TOX",],
    "MAIT": ["SLC4A10","KLRB1","RORC","ZBTB16","NCR3"],
}
sc.pl.dotplot(
        adata,
        groupby="leiden_res1_4",
        var_names=marker_genes,
        standard_scale="var",
        save="epifilter_121_CD8亚群.png"
    )

```
<img src="..\figures\dotplot_epifilter_121_CD8亚群.png">
0,3,7,8:Tem
1,2:Temra
4,9,10:Tex
5:Tm?
6:Tn-like

### 看下DEG
```
sc.tl.rank_genes_groups(
 adata, groupby="leiden_res1_4", method="wilcoxon", key_added="dea_leiden_res1_4"
)
sc.tl.filter_rank_genes_groups(
  adata,
  min_in_group_fraction=0.2,
  max_out_group_fraction=0.2,
  key="dea_leiden_res1_4",
  key_added="dea_leiden_res1_4_filtered",
 )
sc.pl.rank_genes_groups_dotplot(
  adata,
  groupby="leiden_res1_4",
  standard_scale="var",
  n_genes=10,
  key="dea_leiden_res1_4_filtered",
  save="epifilter121_res1_4_CD8T_top10gene.png"
)
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False,
save="121_res1_4_CD8T_top25gene.png",key="dea_leiden_res1_4",)
```
<img src="..\figures\dotplot_epifilter121_res1_4_CD8T_top10gene.png">

把相同功能亚群做DEG看下是否有必要细分
#### Temra
```
mask = adata.obs['leiden_res1_4'].isin(['1','2'])
adata_filtered = adata[mask, :]
sc.tl.rank_genes_groups(
 adata_filtered, groupby="leiden_res1_4", method="wilcoxon", key_added="dea_leiden_res1_4"
)
sc.tl.filter_rank_genes_groups(
  adata_filtered,
  min_in_group_fraction=0.2,
  max_out_group_fraction=0.2,
  key="dea_leiden_res1_4",
  key_added="dea_leiden_res1_4_filtered",
 )
sc.pl.rank_genes_groups_dotplot(
  adata_filtered,
  groupby="leiden_res1_4",
  standard_scale="var",
  n_genes=10,
  key="dea_leiden_res1_4_filtered",
  save="Temra_epifilter121_res1_4_CD8T_top10gene.png"
)
```

<img src="..\figures\dotplot_Temra_epifilter121_res1_4_CD8T_top10gene.png">
虽然只有3个差异基因,但这三个基因都还挺重要，分或不分感觉都可以

#### Tem
```
mask = adata.obs['leiden_res1_4'].isin(['0','3','7','8'])
adata_filtered = adata[mask, :]
sc.tl.rank_genes_groups(
 adata_filtered, groupby="leiden_res1_4", method="wilcoxon", key_added="dea_leiden_res1_4"
)
sc.tl.filter_rank_genes_groups(
  adata_filtered,
  min_in_group_fraction=0.2,
  max_out_group_fraction=0.2,
  key="dea_leiden_res1_4",
  key_added="dea_leiden_res1_4_filtered",
 )
sc.pl.rank_genes_groups_dotplot(
  adata_filtered,
  groupby="leiden_res1_4",
  standard_scale="var",
  n_genes=10,
  key="dea_leiden_res1_4_filtered",
  save="Tem_epifilter121_res1_4_CD8T_top10gene.png"
)
```
<img src="..\figures\dotplot_Tem_epifilter121_res1_4_CD8T_top10gene.png">
3和8其实分不太开,0和7和另外的是分得开的

#### Tex
```
mask = adata.obs['leiden_res1_4'].isin(['4','9','10'])
adata_filtered = adata[mask, :]
sc.tl.rank_genes_groups(
 adata_filtered, groupby="leiden_res1_4", method="wilcoxon", key_added="dea_leiden_res1_4"
)
sc.tl.filter_rank_genes_groups(
  adata_filtered,
  min_in_group_fraction=0.2,
  max_out_group_fraction=0.2,
  key="dea_leiden_res1_4",
  key_added="dea_leiden_res1_4_filtered",
 )
sc.pl.rank_genes_groups_dotplot(
  adata_filtered,
  groupby="leiden_res1_4",
  standard_scale="var",
  n_genes=10,
  key="dea_leiden_res1_4_filtered",
  save="Tex_epifilter121_res1_4_CD8T_top10gene.png"
)
```
<img src="..\figures\dotplot_Tex_epifilter121_res1_4_CD8T_top10gene.png">
第10群能分开

看下来是res1/1.2要合适一点
我勒个豆这CD8T的注释有点曲折
```
marker_genes = {
    "Tn": ["CCR7", "LEF1", "SELL", "TCF7", "MAL","IL7R"],
    "Temra": ["CX3CR1","KLRG1","TBX21","FGFBP2","GNLY","GZMH","GZMB","ASCL2","S1PR1","S1PR5"],
    "Tm" :["CXCR4","ZFP36L2"],
    "Tem": ["FOS", "JUN"],
    "Trm": ["ZNF683","CXCR6","CD69"],
    "Tex": ["PDCD1", "HAVCR2", "LAG3", "LAYN", "TIGIT", "ENTPD1", "CXCL13", "CTLA4", "TNFRSF18", "TOX",],
    "MAIT": ["SLC4A10","KLRB1","RORC","ZBTB16","NCR3"],
}
sc.pl.dotplot(
        adata,
        groupby="leiden_res1_2",
        var_names=marker_genes,
        standard_scale="var",
        save="res1_2_epifilter_121_CD8亚群.png"
    )
```
<img src="..\figures\umapCD8+T_umapmd0.4filterepi_cluster.png">
<img src="..\figures\dotplot_res1_2_epifilter_121_CD8亚群.png">
0:Temra
1,2,7:Tem/Trm(如果是Tumor来源，大概率是Trm)
4:T_KLRG1?
3,5:Tex
6:Tn-like

### DEG
```
sc.tl.rank_genes_groups(
 adata, groupby="leiden_res1_2", method="wilcoxon", key_added="dea_leiden_res1_2"
)
sc.tl.filter_rank_genes_groups(
  adata,
  min_in_group_fraction=0.2,
  max_out_group_fraction=0.2,
  key="dea_leiden_res1_2",
  key_added="dea_leiden_res1_2_filtered",
 )
sc.pl.rank_genes_groups_dotplot(
  adata,
  groupby="leiden_res1_2",
  standard_scale="var",
  n_genes=10,
  key="dea_leiden_res1_2_filtered",
  save="epifilter121_res1_2_CD8T_top10gene.png"
)
```
<img src="..\figures\dotplot_epifilter121_res1_2_CD8T_top10gene.png">
#### Tem/Trm
```
mask = adata.obs['leiden_res1_2'].isin(['1','2','7'])
adata_filtered = adata[mask, :]
sc.tl.rank_genes_groups(
 adata_filtered, groupby="leiden_res1_2", method="wilcoxon", key_added="dea_leiden_res1_2"
)
sc.tl.filter_rank_genes_groups(
  adata_filtered,
  min_in_group_fraction=0.2,
  max_out_group_fraction=0.2,
  key="dea_leiden_res1_2",
  key_added="dea_leiden_res1_2_filtered",
 )
sc.pl.rank_genes_groups_dotplot(
  adata_filtered,
  groupby="leiden_res1_2",
  standard_scale="var",
  n_genes=10,
  key="dea_leiden_res1_2_filtered",
  save="Tem_epifilter121_res1_2_CD8T_top10gene.png"
)
```
<img src="..\figures\dotplot_Tem_epifilter121_res1_2_CD8T_top10gene.png">

#### Tex
```
mask = adata.obs['leiden_res1_2'].isin(['3','5'])
adata_filtered = adata[mask, :]
sc.tl.rank_genes_groups(
 adata_filtered, groupby="leiden_res1_2", method="wilcoxon", key_added="dea_leiden_res1_2"
)
sc.tl.filter_rank_genes_groups(
  adata_filtered,
  min_in_group_fraction=0.2,
  max_out_group_fraction=0.2,
  key="dea_leiden_res1_2",
  key_added="dea_leiden_res1_2_filtered",
 )
sc.pl.rank_genes_groups_dotplot(
  adata_filtered,
  groupby="leiden_res1_2",
  standard_scale="var",
  n_genes=10,
  key="dea_leiden_res1_2",
  save="Tex_epifilter121_res1_2_CD8T_top10gene.png"
)
```
<img src="..\figures\dotplot_Tex_epifilter121_res1_2_CD8T_top10gene.png">
第5群和第3群比，有更高的GNLY,CCL4,NKG7基因表达，杀伤分子表达要更多一些,cytotoxity打分应该会更高

### 看下第6群里面的MAIT和Tn到底可不可以分开
```
marker_genes = ["CCR7", "LEF1", "SELL", "TCF7","IL7R","SLC4A10","KLRB1","NCR3"]
sc.pl.umap(adata_filtered,color=marker_genes,vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="CD8T_clt6_TnMAITmarkers.png")
```
看起来是分得开的呀，为啥leiden就是分不开呢
<img src="..\figures\umapCD8T_clt6_TnMAITmarkers.png">

单独对clt6做一下leiden
```
sc.tl.leiden(adata_filtered, key_added="leiden_res0_1", resolution=0.1)
sc.tl.leiden(adata_filtered, key_added="leiden_res0_2", resolution=0.2)
sc.tl.leiden(adata_filtered, key_added="leiden_res0_3", resolution=0.3)
sc.tl.leiden(adata_filtered, key_added="leiden_res0_4", resolution=0.4)
sc.tl.leiden(adata_filtered, key_added="leiden_res0_6", resolution=0.6)
sc.tl.leiden(adata_filtered, key_added="leiden_res0_8", resolution=0.8)
sc.tl.leiden(adata_filtered, key_added="leiden_res1", resolution=1)
sc.pl.umap(adata_filtered,color=["leiden_res0_2","leiden_res0_3","leiden_res0_4","leiden_res0_6","leiden_res0_8","leiden_res1"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="CD8T_clt6_resubset.png")
```
```
marker_genes = {
    "Tn": ["CCR7", "LEF1", "SELL", "TCF7", "MAL","IL7R"],
    "Temra": ["CX3CR1","KLRG1","TBX21","FGFBP2","GNLY","GZMH","GZMB","ASCL2","S1PR1","S1PR5"],
    "Tm" :["CXCR4","ZFP36L2"],
    "Tem": ["FOS", "JUN"],
    "Trm": ["ZNF683","CXCR6","CD69"],
    "Tex": ["PDCD1", "HAVCR2", "LAG3", "LAYN", "TIGIT", "ENTPD1", "CXCL13", "CTLA4", "TNFRSF18", "TOX",],
    "MAIT": ["SLC4A10","KLRB1","RORC","ZBTB16","NCR3"],
}
marker_genes_in_data = {}
for ct, markers in marker_genes.items():
    markers_found = []
    for marker in markers:
        if marker in adata.var.index:
            markers_found.append(marker)
    marker_genes_in_data[ct] = markers_found
sc.pl.dotplot(
            adata_filtered,
            groupby="leiden_res0_6",
            var_names=marker_genes_in_data,
            standard_scale="var",
            save="res1_2_epifilter_121_CD8_cl6resubset.png"
        )
```
<img src="..\figures\umapCD8T_clt6_TnMAITmarkers.png">
<img src="..\figures\umapCD8T_clt6_resubset.png">
<img src="..\figures\dotplot_res1_2_epifilter_121_CD8_cl6resubset.png">

0,2应该是Tn-like,1群感觉是Temra-like,
3群应该是MAIT,4群不知道是啥，5群不知道是啥

### 有没有可能是双细胞?
```
marker_genes = ["PTPRC","CD3E","CD79A","COL1A1",
"VWF","PECAM1","EPCAM","SFN","KRT5","LYZ","GCA","TPSAB1","NCAM1","FCGR3A","GNLY","KLRD1"]
sc.pl.umap(adata_filtered,color=marker_genes,vmin=0,vmax="p99",sort_order=False,cmap="Blues",save="121_CD8T_clt6_majormarker.png")
```
<img src=..\figures\umap121_CD8T_clt6_majormarker.png>
不是双细胞

总之这个第6群是有点难搞,就归为Tn-like吧

## 一锤定音
res1.2:
0:Temra,特征基因是FGFBP2,CX3CR1,KLRG1,TBX21,GNLY
1:Tem,特征基因是FOS,JUN,HOOK2,MARCKSL1
2:Tem,特征基因是FOS,JUN,ITM2C,PLEK
3:Tex,特征基因是CXCL13,CTLA4,CXCR6
4:Temra-like,KLRG1,FCRL3
5:Tex,特征基因是CXCL13,CCL3,GNLY,CCL4
6:Tn-like,特征基因是TCF7,LTB,SELL,LEF1
7:Tem:特征基因是FOS,JUN,SNX9,CAPG,HOPX

### 点图
```
adata = sc.read("RData/CD8+T__umapmd0.4_filterepi.h5ad")
cl_annotation = {
    "0": "c10_CD8_Temra_FGFBP2",
    "1": "c08_CD8_Tem_FOS",
    "2": "c07_CD8_Tem_ITM2C",
    "3": "c12_CD8_Tex_CXCL13",
    "4": "c11_CD8_Temra-like_KLRG1",
    "5": "c13_CD8_Tex_CCL3",
    "6": "c06_CD8_Tn-like_TCF7",
    "7": "c09_CD8_Tem_SNX9",
}
adata.obs["minor_celltype"]=adata.obs.leiden_res1_2.map(cl_annotation)
adata.write("RData/CD8T_annotated.h5ad")
marker_genes = ["TCF7","LTB","SELL","LEF1",
                "FOS","JUN","ITM2C",
                "HOOK2","MARCKSL1",
                "HOPX","GZMK","NR4A1",
                "FGFBP2","CX3CR1","TBX21",
                "KLRG1",
                "CXCL13","CTLA4","CXCR6",
                "CCL3","GNLY","CCL4"]
from matplotlib.colors import LinearSegmentedColormap
custom_cmap = LinearSegmentedColormap.from_list("custom", ["#0A3D7C", "#F6F6F6", "#810426"])
dp = sc.pl.dotplot(
        adata,
        groupby="minor_celltype",
        categories_order = ["c06_CD8_Tn-like_TCF7",
                            "c07_CD8_Tem_ITM2C",
                            "c08_CD8_Tem_FOS",
                            "c09_CD8_Tem_SNX9",
                            "c10_CD8_Temra_FGFBP2",
                            "c11_CD8_Temra-like_KLRG1",
                            "c12_CD8_Tex_CXCL13",
                            "c13_CD8_Tex_CCL3",],
        var_names=marker_genes,
        standard_scale="var",
        return_fig=True,
        var_group_labels="        ",
        var_group_positions = [(0,3),(4,7),(8,9),
        (10,12),(13,15),(16,16),(17,19),(20,22)],
        cmap = custom_cmap
    )
dp.style(dot_edge_color='black', dot_edge_lw=1).savefig("plots/dotplot_CD8T.png")
```
