# Figure1和S1的出图代码         (What)

* Nov 28, 2025                                 (When)
* liyi /data1/liyi/zhaoyue/FUSCC-ESCC-neoadjuvant_thrapy           (Where)
* 要有图片才能发文章 (Why)

需要出的图如下:
* 1 umap(color=[大类,大类marker,亚群,(组织来源,treatment stage,患者)])
* 3 每个患者的百分比柱状图+拼在上面的小热图(总的和CD45+各一张)
* 4 Ti vs Pi 散点图(这个目前是缺少每个患者具体的regression rate的值)
* 7 大亚群的基线(pre post都要看)和变化与疗效关系-箱线图
* 8 中亚群的基线(pre post都要看)和变化与疗效关系-箱线图
* 9 小亚群的基线(pre post都要看)和变化与疗效关系-箱线图(总细胞数CD45+和major以及CD4T这种亚major都试一下)
* 10 所有的亚聚类结果和组织来源+markergene点图
* 11 做完把methods给写了


## 使用ggplot2绘图
需要3个东西
1 数据框
cell_id UMAP1 UMAP2 minor_UMAP1 minor_UMAP2 大类注释 中类注释 亚类注释 anatomic_site treatment_stage patient patient_id
2 大类颜色
major_celltype_colors =
{
"T cell":"#E0823A",
"B cell":"#3A6FA4" ,
"Myeloid cell":"#4E9669",
"Fibrocyte":"#81564C",
"Endothelial cell":"#6F5198",
"Epithelial cell":"#BA3931"
}

3 亚类颜色
minor_celltype_colors = {
      'c01_CD4_Tn_TCF7': '#4279AD',
      'c02_CD4_Tcm_S1PR1': '#AFCCD3',
      'c03_CD4_Tcm_ANXA1': '#A25936',
      'c04_CD4_Tem_NR4A1': '#6171A8',
      'c05_CD4_Treg_FOXP3': '#BD5BA5',
      'c06_CD4_TOX': '#D5E1B4',
      'c07_CD8_Tem/Teff_TNFSF9': '#F2D7D0',
      'c08_CD8_Tem/Teff_HSPA1B': '#E8A7BA',
      'c09_CD8_Tem/Teff_GNLY': '#8A7FA0',
      'c10_CD8_Teff_GZMK': '#F4EFA9',
      'c11_CD8_Tex_CXCL13': '#E77772',
      'c12_CD8_Temra_FGFBP2': '#4D4199',
      'c13_CD8_IL7R': '#F1B9AE',
      'c14_Tprf_MKI67': '#C5A4C7',
      'c15_NK_CXXC5': '#DF617B',
      'c16_NK/NKT_LAG3': '#816AA9',
      'c17_NK_KLRC1': '#CBD1E7',
      'c18_NK_NRGN': '#5BA0D4',
      'c19_Bn_TCL1A': '#5AA452',
      'c20_Bm_TNFRSF13B': '#F29696',
      'c21_Bm_TUBA1A': '#D8432E',
      'c22_Bp_MZB1': '#F19B9B',
      'c23_Mono_CD14': '#F4BC82',
      'c24_Mono_HDAC9': '#ABE0ED',
      'c25_Mono_FCGR3A': '#C0FEFC',
      'c26_Mono_CCL5': '#D8C6DF',
      'c27_Macro_C1QC': '#F6C589',
      'c28_Macro_G0S2': '#F2ACAE',
      'c29_cDC1_CADM1': '#EF8E42',
      'c30_cDC2_OTUD1': '#8AA624',
      'c31_cDC2_CD1A': '#D39FC8',
      'c32_cDC2_CD1D': '#E5D9E3',
      'c33_cDC3_LAMP3': '#F2BF7D',
      'c34_pDC_MZB1': '#E9655B',
      'c35_DC_ECE1': '#7DD076',
      'c36_Mast_KIT': '#5D7AB1',
      'c37_VEC_ACKR1': '#D2E9C7',
      'c38_VEC_ITPKC': '#9DD1C8',
      'c39_AEC_GJA4': '#BCDD7A',
      'c40_LEC_PROX1': '#EA8675',
      'c41_capEC_RGCC': '#BCBBD7',
      'c42_EC_CSF2RB': '#9BBBD0',
      'c43_prfEC_MKI67': '#4E8975',
      'c44_Fib_COL11A1': '#E3C567',
      'c45_Fib_MMP1': '#A63D40',
      'c46_Fib_COL6A5': '#7F97A2',
      'c47_proFib_CD34': '#8F913D',
      'c48_apFib_HLA-DRB1': '#884C76',
      'c49_PC_RGS5': '#48C0BF',
      'c50_PC_MYH11': '#D9A779',
      'c51_Epi_Normal': '#7FC696',
      'c52_Epi_Tumor': '#2E4F73'
}

得到绘图df
```
import anndata as ad
import pandas as pd
import scanpy as sc

adata = sc.read("RData/1128_final_escc121.h5ad")

for col in ["major_celltype", "TILC_celltype", "myeloid_celltype"]:
    if col in adata.obs.columns:
        adata.obs[col] = adata.obs[col].astype(str)

# 1. 提取基础列
adata.obs["middle_celltype"] = adata.obs["major_celltype"].copy()
adata.obs.loc[adata.obs["major_celltype"] == "T cell", "middle_celltype"] = adata.obs["TILC_celltype"]
adata.obs.loc[adata.obs["major_celltype"] == "Myeloid cell", "middle_celltype"] = adata.obs["myeloid_celltype"]

df = adata.obs[['major_celltype','middle_celltype', 'minor_celltype', 'sample']].copy()

# 2. 添加细胞ID（行索引即为cell_id）
df['cell_id'] = df.index

# 3. 添加UMAP坐标（假设存储在obsm['X_umap']中）
df[['umap1', 'umap2']] = adata.obsm['X_umap']
df['minor_umap1'] = adata.obs['minor_UMAP1']
df['minor_umap2'] = adata.obs['minor_UMAP2']

# 4. 处理anatomic_site：sample含pbmc则为blood，否则为tumor
df['anatomic_site'] = df['sample'].apply(lambda x: 'blood' if 'pbmc' in x.lower() else 'tumor')

# 5. 处理treatment_stage：根据sample中的关键词判断
def get_treatment_stage(sample):
    sample_lower = sample.lower()  # 统一转为小写，避免大小写问题
    if 'prec1' in sample_lower:
        return 'pre'
    elif 'prec2' in sample_lower:
        return 'on'
    elif 'op' in sample_lower:
        return 'post'
    else:
        return 'unknown'  # 处理未匹配到的情况

df['treatment_stage'] = df['sample'].apply(get_treatment_stage)

# 6. 处理patient：sample中以"_"分隔的第一项加最后一项
df['patient'] = df['sample'].apply(lambda x: '_'.join([x.split('_')[0], x.split('_')[-1]]))
# 有几个患者的住院号是不对的，要调整下
#PDD_2284071应该是PDD_2284017
#ZJM_2306458应该是ZJM_2306548
#YLY_2282682应该是YLY_2286282
replace_map = {
    'PDD_2284071': 'PDD_2284017',
    'ZJM_2306458': 'ZJM_2306548',
    'YLY_2282682': 'YLY_2286282'
}
df['patient'] = df['patient'].replace(replace_map)
df['patient_id'] = df['patient'].str.split('_').str[-1]

# 7 添加Patient_Number
unique_patients = df['patient'].unique()
patient_id_map = {pat: f"P{str(i+1).zfill(2)}" for i, pat in enumerate(unique_patients)}
df['patient_number'] = df['patient'].map(patient_id_map)


# 8. 调整列顺序并命名
df = df[['cell_id', 'umap1', 'umap2','minor_umap1','minor_umap2','major_celltype',
        'middle_celltype','minor_celltype', 'anatomic_site', 'treatment_stage', 'patient','patient_id','patient_number']]

# 9. 保存为TSV文件
df.to_csv('output/df_1128.tsv', sep='\t', index=False)
```

## 除了绘制umap大类的绘图代码在biotrainee上面

## 绘制大类UMAP
```
sc.settings.set_figure_params(
    dpi=300,
    facecolor="white",
    frameon=False
)
sc.settings.figdir = "plots"
marker_genes = {
"Immune":["PTPRC"],
"Fib": ["DCN","COL1A1","COL1A2","FN1","COL3A1","COL6A1"],
"Endo" : ["VWF","PECAM1","ENG","CDH5"],
"Epi" : ["EPCAM","SFN","KRT5","KRT14"],
"T":["CD3D","CD2","CD3E","CD3G"],
"B":["CD79A","JCHAIN","CD19","MZB1"],
"Myeloid":["LYZ","CD14","CD68",
           "HPGDS","TPSAB1",
          "GCA"],}

marker_genes_in_data = {}
for ct, markers in marker_genes.items():
    markers_found = []
    for marker in markers:
        if marker in adata.var.index:
            markers_found.append(marker)
    marker_genes_in_data[ct] = markers_found
sc.pl.umap(adata,color=marker_genes_in_data["Fib"],vmin=0,vmax="p99",sort_order=False,cmap="Blues",save="final_大类marekr_Fib.png")
sc.pl.umap(adata,color=marker_genes_in_data["Endo"],vmin=0,vmax="p99",sort_order=False,cmap="Blues",save="final_大类marekr_Endo.png")
sc.pl.umap(adata,color=marker_genes_in_data["Epi"],vmin=0,vmax="p99",sort_order=False,cmap="Blues",save="final_大类marekr_Epi.png")
sc.pl.umap(adata,color=marker_genes_in_data["T"],vmin=0,vmax="p99",sort_order=False,cmap="Blues",save="final_大类marekr_T.png")
sc.pl.umap(adata,color=marker_genes_in_data["B"],vmin=0,vmax="p99",sort_order=False,cmap="Blues",save="final_大类marekr_B.png")
sc.pl.umap(adata,color=marker_genes_in_data["Myeloid"],vmin=0,vmax="p99",sort_order=False,cmap="Blues",save="final_大类marekr_Myeloid.png")
```

<img src="..\figures\umapfinal_大类marekr_T.png">
CD2+CD3E
<img src="..\figures\umapfinal_大类marekr_B.png">
CD79A
<img src="..\figures\umapfinal_大类marekr_Myeloid.png">
LYZ,TPSAB1
<img src="..\figures\umapfinal_大类marekr_Fib.png">
COL3A1
<img src="..\figures\umapfinal_大类marekr_Endo.png">
VWF,CDH5
<img src="..\figures\umapfinal_大类marekr_Epi.png">
EPCAM,SFN



```
final_major_markers = ["PTPRC","CD2","CD3E",
                      "CD79A","LYZ","TPSAB1",
                      "COL3A1","VWF","CDH5",
                      "EPCAM","SFN"]
```

<img src="..\figures\umapfinal_大类marker.png">
