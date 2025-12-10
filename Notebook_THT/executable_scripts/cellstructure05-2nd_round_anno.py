import scanpy as sc

## CD4T
adata = sc.read("RData/CD4+T__umapmd0.4.h5ad")
cl_annotation = {
    "0": "c01_CD4_Tn_TCF7",
    "1": "c02_CD4_Tcm_ANXA1",
    "2": "c03_CD4_Treg_FOXP3",
    "3": "c05_CD4_Tem_NR4A2",
    "4": "c05_CD4_Tem_NR4A2",
    "5": "c04_CD4_Treg_HLA-DPB1",
}
adata.obs["final_anno_celltype"]=adata.obs.leiden_res0_8.map(cl_annotation)
adata.write("RData/CD4T_final_anno.h5ad")

## CD8T
adata = sc.read("RData/CD8+T__umapmd0.4.h5ad")
cl_annotation = {
    "0": "c06_CD8_Tem_NR4A1",
    "1": "c10_CD8_Tex_CXCL13",
    "2": "c09_CD8_Temra_CX3CR1",
    "3": "c07_CD8_Tem_GZMK",
    "4": "C08_D8_Tem_LTB",
}
adata.obs["final_anno_celltype"]=adata.obs.leiden_res0_6.map(cl_annotation)
adata.write("RData/CD8T_final_anno.h5ad")


## NKNKT
adata = sc.read("RData/NKT_ILC__umapmd0.4.h5ad")
cl_annotation = {
    "0": "c14_NK/NKT_KLRC3",
    "1": "c11_NK_KLRC1",
    "2": "c14_NK/NKT_KLRC3",
    "3": "c12_NK_CXXC5",
    "4": "c12_NK_CXXC5",
    "5":"c14_NK/NKT_KLRC3",
    "6": "c15_NK/NKT_HSPA1B",
    "7": "c13_NK_NRGN",
}
adata.obs["final_anno_celltype"]=adata.obs.leiden_res1.map(cl_annotation)
adata.write("RData/NKNKT_final_anno.h5ad")

## B
adata = sc.read("RData/B__umapmd0.4.h5ad")
adata = adata[adata.obs.leiden_res0_8 != '4' ].copy()
adata = adata[adata.obs.leiden_res0_8 != '13' ].copy()
cl_annotation = {
    "0": "c18_Bm_TNFRSF13B",
    "1": "c16_Bn_FCER2",
    "2": "c19_Bp_MZB1",
    "3": "c19_Bp_MZB1",
    "5": "c19_Bp_MZB1",
    "6": "c19_Bp_MZB1",
    "7": "c19_Bp_MZB1",
    "8": "c19_Bp_MZB1",
    "9": "c19_Bp_MZB1",
    "10": "c19_Bp_MZB1",
    "11":"c19_Bp_MZB1",
    "12": "c19_Bp_MZB1",
    "14":"c16_Bn_FCER2",
    "15": "c19_Bp_MZB1",
    "16": "c17_Bgc_HMGB2"
}
adata.obs["final_anno_celltype"]=adata.obs.leiden_res0_8.map(cl_annotation)
adata.write("RData/B_final_anno.h5ad")


## Myeloid
adata = sc.read("RData/Myeloid__umapmd0.4.h5ad")
cl_annotation = {
    "0": "c20_Mono_CD14",
    "1": "c22_Macro_C1QC",
    "2": "c20_Mono_CD14",
    "3": "c20_Mono_CD14",
    "4": "c20_Mono_CD14",
    "5":"c23_Macro_G0S2",
    "6": "c20_Mono_CD14",
    "7": "c24_Mast_KIT",
    "8": "c25_cDC1/2_CD1C",
    "9": "c21_Mono_FCGR3A",
    "10": "c27_pDC_LILRA4",
    "11": "c26_cDC3_LAMP3",
}
adata.obs["final_anno_celltype"]=adata.obs.leiden_res1.map(cl_annotation)
adata.write("RData/Myeloid_final_anno.h5ad")

## Fib
adata = sc.read("RData/Fib__umapmd0.4.h5ad")
cl_annotation = {
    "0": "c28_Fib_CFD",
    "1": "c29_Fib_MMP1",
    "2": "c33_Pericyte_RGS5",
    "3": "c30_Fib_CST1",
    "4": "c31_Fib_PTGDS",
    "5":"c35_SMC_TAGLN",
    "6": "c34_Pericyte_MYH11",
    "7": "c32_apFib_CD74",
}
adata.obs["final_anno_celltype"]=adata.obs.leiden_res0_8.map(cl_annotation)
adata.write("RData/Fib_final_anno.h5ad")


## Endo
adata = sc.read("RData/Endo__umapmd0.4.h5ad")
adata = adata[adata.obs.leiden_res0_8 != '5'].copy()
cl_annotation = {
    "0": "c36_VEC_IL33",
    "1": "c40_capEC_RGCC",
    "2": "c36_VEC_IL33",
    "3": "c40_capEC_RGCC",
    "4": "c38_AEC_GJA4",
    "6": "c39_LEC_CCL21",
    "7": "c37_VEC_CSF2RB",
}
adata.obs["final_anno_celltype"]=adata.obs.leiden_res0_8.map(cl_annotation)
adata.write("RData/Endo_final_anno.h5ad")

## Epi
adata = sc.read("RData/Epi__umapmd0.4.h5ad")
adata.obs["final_anno_celltype"]= "Epi"
adata.write("RData/Epi_final_anno.h5ad")


