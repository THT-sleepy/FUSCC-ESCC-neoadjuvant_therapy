# T细胞注释过程及结果记录         (What)

* Oct 17, 2025                                 (When)
* liyi /data1/liyi/zhaoyue/FUSCC-ESCC-neoadjuvant_thrapy           (Where)
* 单细胞需要注释后才能进行后续分析 (Why)

## 1 再次检查有无混杂的其它细胞系如上皮的marker表达
* 聚类情况
<img src="..\figures\umapT_umapmd0.4_cluster.png">

```
import scanpy as sc

sc.settings.set_figure_params(
    dpi=300,
    facecolor="white",
    frameon=False
)
sc.settings.figdir = "plots"

adata = sc.read("RData/T__umapmd0.4_addleiden2.h5ad")

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
sc.pl.umap(adata,color=marker_genes_in_data["Fib"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Tsubset_Fib.png")
sc.pl.umap(adata,color=marker_genes_in_data["Endo"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Tsubset_Endo.png")
sc.pl.umap(adata,color=marker_genes_in_data["Epi"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Tsubset_Epi.png")
sc.pl.umap(adata,color=marker_genes_in_data["B"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Tsubset_B.png")
sc.pl.umap(adata,color=marker_genes_in_data["Myeloid"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Tsubset_Myeloid.png")

sc.pl.dotplot(
    adata,
    groupby="leiden_res1_2",
    var_names=marker_genes_in_data,
    standard_scale="var",
    save="Tsubset_大类marker.png"
)
```
<img src="..\figures\umapTsubset_Fib.png">
看起来没有混的纤维细胞<br>

<img src="..\figures\umapTsubset_Endo.png">
好像也没有混的内皮细胞

<img src="..\figures\umapTsubset_Epi.png">
好像也没有混的上皮细胞

<img src="..\figures\umapTsubset_B.png">
B细胞也没有混

<img src="..\figures\umapTsubset_Myeloid.png">
GCA不对劲！这一群(res1.2 cl10)可能是粒细胞混了T细胞

* 点图
<img src="..\figures\dotplot_Tsubset_大类marker.png">
第10群确实是会高表达GCA这个中性粒的marker，去掉吧
```
adata = adata[adata.obs.leiden_res1_2 != '10'].copy()
```

## 2 接下来区分T和ILC-CD4,CD8,DP(去掉),增殖T,NKT和NK及其它可能的ILC
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
sc.pl.umap(adata,color=marker_genes_in_data["T"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Tsubset_T.png")
sc.pl.umap(adata,color=marker_genes_in_data["NK"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Tsubset_NK.png")
sc.pl.dotplot(
    adata,
    groupby="leiden_res1_2",
    var_names=marker_genes_in_data,
    standard_scale="var",
    save="Tsubset_T_ILCmarker.png"
)
```

<img src="..\figures\umapTsubset_T.png">
左上和左下属于ILC
<img src="..\figures\umapTsubset_NK.png">
NK细胞确实分布在左上和左下

<img src="..\figures\dotplot_Tsubset_T_ILCmarker.png">

```
cl_annotation = {
    "0": "CD4+T",
    "1": "NKT/ILC",
    "2": "CD8+T",
    "3": "CD4+T",
    "4": "CD8+T",
    "5": "DP_T?",
    "6": "CD4T",
    "7": "CD8+T",
    "8": "CD8+T/NK",
    "9": "NKT/ILC",
    "11":"DP_T",
}
```
后面结合umap看了下分群有变化
第8群有点麻烦,感觉需要提高一点分辨率(试一下2吧)

```
sc.pl.dotplot(
    adata,
    groupby="leiden_res2",
    var_names=marker_genes_in_data,
    standard_scale="var",
    save="Tsubset_T_ILCmarker_res2.png"
)
```
<img src="..\figures\umapT_umapmd0.4_cluster.png">
<img src="..\figures\umapT_leiden2_cluster.png">
<img src="..\figures\dotplot_Tsubset_T_ILCmarker_res2.png">
其实还是会有T/NK分不开的情况，就把上面的第8群认为是CD8T吧

```
adata = adata[adata.obs.leiden_res1_2 != '11'].copy()
```
把第11群这个明确的DP去掉,记得在method里面补上这个操作
adata.obs["Immune_celltype"]=adata.obs.leiden_res2.map(cl_annotation)

忘记增殖T了,补上
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
sc.pl.umap(adata,color=marker_genes_in_data["Tprf"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Tsubset_Tprf.png")
```
<img src="..\figures\umapTsubset_Tprf.png">
很少的prfT，可以不分这一类


<img src="..\figures\umapT_umapmd0.4_cluster.png">
<img src="..\figures\umapTsubset_T.png">
<img src="..\figures\umapTsubset_NK.png">
<img src="..\figures\dotplot_Tsubset_T_ILCmarker.png">

## 先看看celltypist的结果
<img src="..\figures\umapcelltypist_coarse_T.png">
<img src="..\figures\umapcelltypist_fine_T.png">

```
cl_annotation = {
    "0": "CD4+T",
    "1": "NKT_ILC",
    "2": "CD8+T",
    "3": "CD4+T",
    "4": "CD8+T",
    "5": "CD4+T",
    "6": "CD4+T",
    "7": "CD8+T",
    "8": "CD8+T",
    "9": "NKT_ILC",
    "11":"CD8+T"
}
adata.obs["TILC_celltype"]=adata.obs.leiden_res1_2.map(cl_annotation)
adata.write("RData/T_1st_round_anno.h5ad")
```
### 注意第5群暂时分为CD4+T,实际上有一点CD8 marker表达



### 接下来把NKT/ILC,CD4+T,CD8+T分别从0.6-1再进行亚群聚类
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
celltypes_to_process = ["CD4+T", "CD8+T", "NKT/ILC"]

# 2. 批量运行（可自定义保存目录和图片前缀）
for celltype in celltypes_to_process:
    process_celltype(
        adata=adata,
        celltype_name=celltype,
        save_dir="RData/",  # h5ad文件保存到RData目录
        plot_prefix=""  # 图片名前缀（如scRNA_T_umapmd0.1_cluster.png）
    )
```

# final_anno NK,NKT,ILC
## 先看看celltypist的结果
<img src="..\figures\umapcelltypist_coarse_NKNKT.png">
<img src="..\figures\umapcelltypist_fine_NKNKT.png">

### 利用umap优化marker和确定分辨率
```
import scanpy as sc

sc.settings.set_figure_params(
    dpi=300,
    facecolor="white",
    frameon=False
)
sc.settings.figdir = "plots"
sc.settings.autoshow = False

adata = sc.read("RData/NKT_ILC__umapmd0.4.h5ad") #X是log1p_norm

marker_genes = {
"NK":["NCAM1","NKG7","KLRC1","KLRD1","KLRG1","GNLY", "TYROBP", "FCGR3A"],
"NKT":["KIR","EOMES","TXK","CD3D","CD3E","CD3G","CD4","CD8A","CD8B","NCAM1"],
"ILC":["ID2", "PLCG2", "GNLY", "SYNE1","IL7R"],
}
marker_genes_in_data = {}
for ct, markers in marker_genes.items():
    markers_found = []
    for marker in markers:
        if marker in adata.var.index:
            markers_found.append(marker)
    marker_genes_in_data[ct] = markers_found
sc.pl.umap(adata,color=marker_genes_in_data["NK"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Tresubset_NK.png")
sc.pl.umap(adata,color=marker_genes_in_data["NKT"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Tsubset_NKT.png")
sc.pl.umap(adata,color=marker_genes_in_data["ILC"],vmin=0,vmax="p99",sort_order=False,cmap="Reds",save="Tsubset_ILC.png")
```

* 聚类情况
<img src="..\figures\umapNKT_ILC_umapmd0.4_cluster.png">
* NK
<img src="..\figures\umapTresubset_NK.png">
可以保留KLRC1,FCER1G,NCAM1,KLRD1

* NKT
<img src="..\figures\umapTresubset_NKT.png">
NK的markerEOMES,NCAM1+T细胞marker CD3D,CD3G

* ILC
<img src="..\figures\umapTresubset_ILC.png">
分不出单独的marker

看下来res1的0,2,5,6像是NKT,其它的像是NK，ILC很难分出来

*点图
```
marker_genes = {
"NK":["KLRC1","FCGR3A","NCAM1","KLRD1"],
"NKT":["EOMES","NCAM1","CD3D","CD3G"],
}
sc.pl.dotplot(
        adata,
        groupby="leiden_res1",
        var_names=marker_genes,
        standard_scale="var",
        save="NKNKT亚群.png"
    )
```

<img src="..\figures\dotplot_NKNKT亚群.png">
0,2,5,6归为NK/NKT;1,3,4,7归为NK

* 再根据top10gene选marker
```
##各cluster top10gene
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
  save="NK_NKT_top10gene.png"
)
```
<img src="..\figures\dotplot_NK_NKT_top10gene.png">

### 后面发现0,2,5分不开,所以DEG要重新算一下
```
cl_annotation = {
    "0": "NK/NKT_?",
    "1": "NK_KLRC1",
    "2": "NK/NKT_?",
    "3": "NK_CXXC5",
    "4": "NK_CXXC5",
    "5":"NK/NKT_?",
    "6": "NK/NKT_HSPA1B",
    "7": "NK_NRGN",
}
adata.obs["final_anno_celltype"]=adata.obs.leiden_res1.map(cl_annotation)
sc.tl.rank_genes_groups(
 adata, groupby="final_anno_celltype", method="wilcoxon", key_added="dea_final_anno_celltype"
)
sc.tl.filter_rank_genes_groups(
  adata,
  min_in_group_fraction=0.2,
  max_out_group_fraction=0.2,
  key="dea_final_anno_celltype",
  key_added="dea_final_anno_celltype_filtered",
 )
sc.pl.rank_genes_groups_dotplot(
  adata,
  groupby="final_anno_celltype",
  standard_scale="var",
  n_genes=10,
  key="dea_final_anno_celltype_filtered",
  save="NK_NKT_anno_top10gene.png"
)
```
<img src="..\figures\dotplot_NK_NKT_anno_top10gene.png">




## 最终注释
```
adata = sc.read("RData/NKT_ILC__umapmd0.4.h5ad")
cl_annotation = {
    "0": "NK/NKT_KLRC3",
    "1": "NK_KLRC1",
    "2": "NK/NKT_KLRC3",
    "3": "NK_CXXC5",
    "4": "NK_CXXC5",
    "5":"NK/NKT_KLRC3",
    "6": "NK/NKT_HSPA1B",
    "7": "NK_NRGN",
}
adata.obs["final_anno_celltype"]=adata.obs.leiden_res1.map(cl_annotation)
adata.write("RData/NKNKT_final_anno.h5ad")
```
### final umap和点图
```
sc.pl.umap(adata,color="final_anno_celltype",save="NKNKT_anno.png")

marker_genes = {
"T":["CD3D","CD3G"],
"NK":["KLRC1","NCAM1","FCGR3A","KLRD1"],
"others":["HSPA1B","KLRC3","CXXC5","NRGN"],
}
sc.pl.dotplot(
    adata,
    groupby="final_anno_celltype",
    var_names=marker_genes,
    standard_scale="var",
    save="NKNKT_finalanno.png"
)

```
* final umap
<img src="..\figures\umapNKNKT_anno.png">
* dotplot
<img src="..\figures\dotplot_NKNKT_finalanno.png">







### final anno CD4+T细胞

```
import scanpy as sc

sc.settings.set_figure_params(
    dpi=300,
    facecolor="white",
    frameon=False
)
sc.settings.figdir = "plots"

adata = sc.read("RData/CD4+T__umapmd0.4.h5ad")

```
* 聚类图
<img src="..\figures\umapCD4+T_umapmd0.4_cluster.png">
## 先看看celltypist的结果
<img src="..\figures\umapcelltypist_coarse_CD4T.png">
<img src="..\figures\umapcelltypist_fine_CD4T.png">

## 根据umap优化下marker

```
marker_genes = {
    "Activated Teff":["CD40LG","CD69","CD38"],
    "Treg": ["FOXP3"],
    "Tn_Tcm_Tfh_pre": [
        "TCF7",
        "LEF1",
        "TXK",
        "CCR7",
        "SELL",
        "MAL",
        "CXCR5",
        "ANXA1",
        "ANXA2",
        "GPR183",
        "PRF1",
        "GZMA",
        "ADSL",
        "IL16",
        "IL7R"
    ],
    "Tem": [
        "KLRG1",
        "CX3CR1",
        "TBX21",
        "GZMK",
        "GNLY",
        "IFNG",
        "GZMH",
        "PRF1",
        "NKG7",
        "GZMA",
        "GZMB",
        "FOS",
        "JUN",
        "TNF"
    ],
    "Th17": [
        "RORA",
        "RORC",
        "CCR6",
        "IL23R",
        "IL22",
        "IL17A",
        "IL17F",
        "IL26"
    ],
    "Trm":["ZNF683","CXCR6"],
    "Tex":["PDCD1","TOX","CXCL13","TIGIT",
          "CTLA4","TNFRSF9","HAVCR2","LAG3"],
    "Tfh/Th1": [
        "TOX2",
        "IL21",
        "GNG4",
        "CD200",
        "BCL6",
        "ZBED2",
        "CCL3",
        "CCL4",
        "IFNG",
        "GZMB",
        "LAG3",
    ],
    "Th2": [
        "CCR3",
        "CCR8",
        "CXCR4",
        "ST2",
        "STAT5"
        "STAT6",
        "GATA3",
        "CCR4"
    ],
}
marker_genes_in_data = {}
for ct, markers in marker_genes.items():
    markers_found = []
    for marker in markers:
        if marker in adata.var.index:
            markers_found.append(marker)
    marker_genes_in_data[ct] = markers_found
sc.pl.dotplot(
    adata,
    groupby="leiden_res0_8",
    var_names=marker_genes_in_data,
    standard_scale="var",
    save="Tsubset_CD4Tanno.png"
)
##各cluster top10gene
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

* 聚类图
<img src="..\figures\umapCD4+T_umapmd0.4_cluster.png">
* 点图
<img src="..\figures\dotplot_Tsubset_CD4Tanno.png">
* DEG top10gene
<img src="..\figures\dotplot_CD4T_top10gene.png">
0群:Tn,CD4_Tn_?<br>
1群:Tn+ANXA1,CD4_Tcm_ANXA1<br> 这个ANXA1是参考的张泽民结肠癌的marker
2群:Foxp3+Treg,CD4_Treg_FOXP3<br>
3群：FOS,JUN,IL7R根据胸科marker应该是CD4_Tem_?<br>
4群：FOS,JUN,IL7R根据胸科marker应该是CD4_Tem_?<br>
5群：之前觉得是Tcm，但CCR7是阴性的而且耗竭marker都是阴性的，所以应该不是，按照celltypist注释为Treg吧，CD4_Treg_?<br>

```
cl_annotation = {
    "0": "CD4_Tn_MYC",
    "1": "CD4_Tcm_ANXA1",
    "2": "CD4_Treg_FOXP3",
    "3": "CD4_Tem_KLRB1",
    "4": "CD4_Tem_NR4A1",
    "5": "CD4_Treg_CDC25B",
}
adata.obs["final_anno_celltype"]=adata.obs.leiden_res0_8.map(cl_annotation)
adata.write("RData/CD4T_final_anno.h5ad")
```
后面有修改，3和4分不开

### final umap和点图
```
marker_genes = {
"Activated T" : ["CD69"],
"Treg":["FOXP3"],
"Tn/Tcm":["TCF7","LEF1","TXK","SELL","MAL","CCR7","IL7R"],
"Tcm":["GPR183","ANXA1","ANXA2"],
"Tem":["FOS","JUN"],
"pre Tfh":["CXCR5"],
"Th17":["RORA","RORC","CCR6","IL23R","IL22",
"IL17A","IL17F","IL26"],
"Trm":["ZNF683","CXCR6"],
"Tex":["PDCD1","TOX","CXCL13","TIGIT",
      "CTLA4","TNFRSF9","HAVCR2","LAG3"],
"Tfh/Th1":["TOX2","IL21","GNG4","CD200","BCL6","ZBED2","CCL3","CCL4","IFNG","GZMB","LAG3",],
"Th2":[
    "CCR3",
    "CCR8",
    "CXCR4",
    "ST2",
    "STAT5"
    "STAT6",
    "GATA3",
    "CCR4"
],
}
marker_genes_in_data = {}
for ct, markers in marker_genes.items():
    markers_found = []
    for marker in markers:
        if marker in adata.var.index:
            markers_found.append(marker)
    marker_genes_in_data[ct] = markers_found
sc.pl.dotplot(
        adata,
        groupby="leiden_res0_8",
        var_names=marker_genes_in_data,
        standard_scale="var",
        save="Tsubset_CD4Tfinalanno.png"
    )
```
<img src="..\figures\dotplot_Tsubset_CD4Tfinalanno.png">
3群和4群实际上根据经典marker和topgene都分不开


## 重新选一下特征基因
```
cl_annotation = {
    "0": "CD4_Tn_MYC",
    "1": "CD4_Tcm_ANXA1",
    "2": "CD4_Treg_FOXP3",
    "3": "CD4_Tem_?",
    "4": "CD4_Tem_?",
    "5": "CD4_Treg_CDC25B",
}
adata.obs["final_anno_celltype"]=adata.obs.leiden_res0_8.map(cl_annotation)
sc.tl.rank_genes_groups(
 adata, groupby="final_anno_celltype", method="wilcoxon", key_added="dea_final_anno_celltype"
)
sc.tl.filter_rank_genes_groups(
  adata,
  min_in_group_fraction=0.2,
  max_out_group_fraction=0.2,
  key="dea_final_anno_celltype",
  key_added="dea_final_anno_celltype_filtered",
 )
sc.pl.rank_genes_groups_dotplot(
  adata,
  groupby="final_anno_celltype",
  standard_scale="var",
  n_genes=10,
  key="dea_final_anno_celltype_filtered",
  save="CD4T_2nd_anno_top10gene.png"
)
```

<img src="..\figures\dotplot_CD4T_2nd_anno_top10gene.png">


```
adata = sc.read("RData/CD4+T__umapmd0.4.h5ad")
cl_annotation = {
    "0": "CD4_Tn_TCF7",
    "1": "CD4_Tcm_ANXA1",
    "2": "CD4_Treg_FOXP3",
    "3": "CD4_Tem_NR4A2",
    "4": "CD4_Tem_NR4A2",
    "5": "CD4_Treg_HLA-DPB1",
}
adata.obs["final_anno_celltype"]=adata.obs.leiden_res0_8.map(cl_annotation)
sc.pl.umap(adata,color="final_anno_celltype",save="CD4T_finalanno.png")
marker_genes = {
"Treg":["FOXP3"],
"Tn/Tcm":["TCF7","LEF1","TXK","SELL","MAL","CCR7","IL7R"],
"Tcm":["GPR183","ANXA1","ANXA2"],
"Tem":["NRF4A2","FOS","JUN"],
"Trm":["ZNF683","CXCR6"],
"Tex":["PDCD1","TOX","CXCL13","TIGIT",
      "CTLA4","TNFRSF9","HAVCR2","LAG3"],
"others":["NR4A2","HLA-DPB1"]
}
marker_genes_in_data = {}
for ct, markers in marker_genes.items():
    markers_found = []
    for marker in markers:
        if marker in adata.var.index:
            markers_found.append(marker)
    marker_genes_in_data[ct] = markers_found
sc.pl.dotplot(
        adata,
        groupby="final_anno_celltype",
        var_names=marker_genes_in_data,
        standard_scale="var",
        save="Tresubset_CD4Tfinalanno.png"
    )
```
### final umap
<img src="..\figures\umapCD4T_finalanno.png">

### final点图
<img src="..\figures\dotplot_Tresubset_CD4Tfinalanno.png">







# final anno CD8+T

```
import scanpy as sc

sc.settings.set_figure_params(
    dpi=300,
    facecolor="white",
    frameon=False
)
sc.settings.figdir = "plots"
adata = sc.read("RData/CD8+T__umapmd0.4.h5ad")


```

* 聚类
<img src="..\figures\umapCD8+T_umapmd0.4_cluster.png">
* 先看看celltypist的结果
<img src="..\figures\umapcelltypist_coarse_CD8T.png">
<img src="..\figures\umapcelltypist_fine_CD8T.png">

```
marker_genes = {
    "Tn/Tcm": ["TCF7", "LEF1", "CCR7", "SELL", "MAL"],
    "Tm": ["IL7R",
          "GPR183",
          "ZFP36L2",
          "CXCR4"],
    "Trm": ["ZNF683",
          "CXCR6",
          ],
    "terminal differentiation": ["TBX21", "ASCL2", "CX3CR1", "KLRG1"],
    "Tem": [
         "ID2",
         "GZMK",
         "GZMH",
         "CD27",
         "HLA-DRB1",
         "GNLY",
         "IFNG",
         "PRF1",
         "NKG7",
         "GZMA",
         "GZMB"],
    "Tex": ["PDCD1", "CXCL13", "LAYN","LAG3","CTLA4","TIGIT",
        "HAVCR2"],
}
marker_genes_in_data = {}
for ct, markers in marker_genes.items():
    markers_found = []
    for marker in markers:
        if marker in adata.var.index:
            markers_found.append(marker)
    marker_genes_in_data[ct] = markers_found
sc.pl.dotplot(
    adata,
    groupby="leiden_res0_6",
    var_names=marker_genes_in_data,
    standard_scale="var",
    save="Tsubset_CD8Tanno.png"
)

```
* 聚类
<img src="..\figures\umapCD8+T_umapmd0.4_cluster.png"><br>
<img src="..\figures\dotplot_Tsubset_CD8Tanno.png"><br>


```
sc.pl.dotplot(
    adata,
    groupby="leiden_res0_8",
    var_names=marker_genes_in_data,
    standard_scale="var",
    save="Tsubset_CD8Tanno.png"
)
##各cluster top10gene
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
  save="CD8T_top10gene.png"
)
```

0群:不是Tcm/Tn,也不是Temra,同时有记忆细胞marker的表达,ZNF683/CXCR6也没有高表达所以不能确定为Trm，可以就归到Tem,CD8_Tem_?
1群:效应分子+耗竭分子,CD8_Tex_?,
2群:效应分子+衰老分子+lack of CD27,CD8_Temra_?
3群:记忆细胞marker+GZMK,CD8_Tem_GZMK
4群:记忆细胞marker+GZMK,CD8_Tem_LTB #不太对，Tem应该是激活状态的，IL7R不应该阳性才对吧;感觉这群细胞是很特殊的,又倒回去看来了下文献胸科的Tem也是IL7R高表达的，所以应该没问题
<img src="..\figures\dotplot_CD8T_top10gene.png">

```
adata = sc.read("RData/CD8+T__umapmd0.4.h5ad")
cl_annotation = {
    "0": "CD8_Tem_NR4A1",
    "1": "CD8_Tex_CXCL13",
    "2": "CD8_Temra_CX3CR1",
    "3": "CD8_Tem_GZMK",
    "4": "CD8_Tem_LTB",
}
adata.obs["final_anno_celltype"]=adata.obs.leiden_res0_6.map(cl_annotation)
adata.write("RData/CD8T_final_anno.h5ad")
```
### 第3群有点奇怪的，DEG高表达的基因居然是S100A8,A9和KRT5,14这些..令人担忧,如果自动注释也认为其是其它类型细胞就去除吧，不过celltypist也认为是T

## final umap和点图
sc.pl.umap(adata,color="final_anno_celltype",save="CD8T_finalanno.png")
marker_genes = {
    "Tn/Tcm": ["TCF7", "LEF1", "CCR7", "SELL", "MAL"],
    "Tm": ["IL7R",
          "GPR183",
          "ZFP36L2",
          "CXCR4"],
    "Trm": ["ZNF683",
          "CXCR6",
          ],
    "terminal differentiation": ["TBX21", "ASCL2", "CX3CR1", "KLRG1"],
    "Tem": [
         "ID2",
         "GZMK",
         "GZMH",
         "CD27",
         "HLA-DRB1",
         "GNLY",
         "IFNG",
         "PRF1",
         "NKG7",
         "GZMA",
         "GZMB"],
    "Tex": ["PDCD1", "CXCL13", "LAYN","LAG3","CTLA4","TIGIT",
        "HAVCR2"],
    "others":["NR4A1","LTB"]
}
marker_genes_in_data = {}
for ct, markers in marker_genes.items():
    markers_found = []
    for marker in markers:
        if marker in adata.var.index:
            markers_found.append(marker)
    marker_genes_in_data[ct] = markers_found
sc.pl.dotplot(
        adata,
        groupby="final_anno_celltype",
        var_names=marker_genes_in_data,
        standard_scale="var",
        save="Tresubset_CD8Tfinalanno.png"
    )

## final umap
<img src="..\figures\umapCD8T_finalanno.png">

## final 点图
<img src="..\figures\dotplot_Tresubset_CD8Tfinalanno.png">
