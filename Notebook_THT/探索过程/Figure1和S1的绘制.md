# Figure1和S1的出图代码         (What)

* Oct 26, 2025                                 (When)
* liyi /data1/liyi/zhaoyue/FUSCC-ESCC-neoadjuvant_thrapy           (Where)
* 要有图片才能发文章 (Why)

需要出的图如下:
* 1 umap(color=[大类,大类marker,亚群,(组织来源,treatment stage,患者)])
* 2 箱线+jitter 图(每个大类MPR vs Non-MPR)-done
* 3 百分比柱状图(pre-on-op)-done
* 3.5换下umap配色-done
* 4 Ti vs Pi 散点图(这个琢磨下再弄)

* 5 每个患者取样的时间(有空再弄)
*



## 我准备使用紫禁城配色
<img src="..\figures\紫禁城配色1.png">
<img src="..\figures\紫禁城配色2.png">
<img src="..\figures\紫禁城配色3.png">


## 使用ggplot2绘图
需要3个东西
1 数据框
cell_id UMAP1 UMAP2 大类注释 亚类注释 anatomic_site treatment_stage patient patient_id
2 大类颜色
<img src="..\figures\大类配色.png">
{
"T":"#ec4c2c", 檎丹
"B":"#fcec6c" , 松花
"Myeloid":"#d4bc44",姚黄
"Fib":"#acd4ec",碧落
"Endo":"#c4dc4c", 清粲
"Epi":"#f4a494" 海天霞
}
3 亚类颜色
cluster_colors = {
    # T细胞相关（基色：T="#DBDBFF"）
    'c01_CD4_Tn_TCF7': '#4279AD',
    'c02_CD4_Tcm_ANXA1': '#AFCCD3',
    'c03_CD4_Treg_FOXP3': '#A25936',
    'c04_CD4_Treg_HLA-DPB1': '#6171A8',
    'c05_CD4_Tem_NR4A2': '#BD5BA5',
    'C08_D8_Tem_LTB': '#D5E1B4',
    'c06_CD8_Tem_NR4A1': '#F2D7D0',
    'c07_CD8_Tem_GZMK': '#E8A7BA',
    'c09_CD8_Temra_CX3CR1': '#8A7FA0',
    'c10_CD8_Tex_CXCL13': '#F4EFA9',
    'c11_NK_KLRC1': '#E77772',
    'c12_NK_CXXC5': '#4D4199',
    'c13_NK_NRGN': '#F1B9AE',
    'c14_NK/NKT_KLRC3': '#C5A4C7',
    'c15_NK/NKT_HSPA1B': '#DF617B',

    # B细胞相关（基色：B="#FECF99"）
    'c16_Bn_FCER2': '#816AA9',
    'c17_Bgc_HMGB2': '#CBD1E7',
    'c18_Bm_TNFRSF13B': '#5BA0D4',
    'c19_Bp_MZB1': '#5AA452',

    # 髓系相关（基色：Myeloid="#FFA430"）
    'c20_Mono_CD14': '#F29696',
    'c21_Mono_FCGR3A': '#D8432E',
    'c22_Macro_C1QC': '#F19B9B',
    'c23_Macro_G0S2': '#F4BC82',
    'c24_Mast_KIT': '#ABE0ED',
    'c25_cDC1/2_CD1C': '#C0FEFC',
    'c26_cDC3_LAMP3': '#D8C6DF',
    'c27_pDC_LILRA4': '#F6C589',

    # 成纤维细胞相关（基色：Fib="#E6CCE7"）
    'c28_Fib_CFD': '#F2ACAE',
    'c29_Fib_MMP1': '#EF8E42',
    'c30_Fib_CST1': '#8AA624',
    'c31_Fib_PTGDS': '#D39FC8',
    'c32_apFib_CD74': '#E5D9E3',

    # 周细胞/平滑肌（归为Fib相关变体）
    'c33_Pericyte_RGS5': '#F2BF7D',
    'c34_Pericyte_MYH11': '#E9655B',

    'c35_SMC_TAGLN': '#7DD076',
    # 内皮细胞相关（基色：Endo="#C999CB"）
    'c36_VEC_IL33': '#5D7AB1',
    'c37_VEC_CSF2RB': '#D2E9C7',
    'c38_AEC_GJA4': '#9DD1C8',
    'c39_LEC_CCL21': '#BCDD7A',
    'c40_capEC_RGCC': '#EA8675',

    # 上皮细胞相关（基色：Epi="#CCCCCC"）
    'c41_Epi_Tumor': '#BCBBD7'
}

得到绘图df
```
import anndata as ad
import pandas as pd
import scanpy as sc

adata = sc.read("RData/final_escc127.h5ad")

# 假设你的AnnData对象名为adata
# 1. 提取基础列
df = adata.obs[['first_round_celltype', 'final_anno_celltype', 'sample']].copy()

# 2. 添加细胞ID（行索引即为cell_id）
df['cell_id'] = df.index

# 3. 添加UMAP坐标（假设存储在obsm['X_umap']中）
df[['umap1', 'umap2']] = adata.obsm['X_umap']

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

# 7. 调整列顺序并命名
df = df[['cell_id', 'UMAP1', 'UMAP2', 'first_round_celltype',
         'final_anno_celltype', 'anatomic_site', 'treatment_stage', 'patient']]

df = df.rename(columns={
    'first_round_celltype': 'MajorCellType',
    'final_anno_celltype': 'SubCellType'
})

# 8 添加Patient_ID
unique_patients = df['patient'].unique()
patient_id_map = {pat: f"P{str(i+1).zfill(2)}" for i, pat in enumerate(unique_patients)}
df['patient_id'] = df['patient'].map(patient_id_map)
df = df[['cell_id', 'UMAP1', 'UMAP2', 'MajorCellType',
         'SubCellType', 'anatomic_site', 'treatment_stage', 'patient', 'patient_id']]


# 9. 保存为TSV文件
df.to_csv('df_1027.tsv', sep='\t', index=False)
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
