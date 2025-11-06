## B细胞主要包括Naive B细胞,Memory B细胞,Plasma B细胞以及生发中心B细胞
## Plasmablast,short-lived plasma cells,long-lived plasma cycle_cells
<img src="..\figures\人B细胞发育和细胞标志物.png">

### 张泽明泛癌B
将肿瘤浸润B细胞分为了Naive B(Bn),memory B(Bm),germinal centerB(Bgc)
以及antibody-secreting cells(ASC)
<img src = "..\figures\张泽明泛癌B点图.png">
marker_genes = {
"Bn":["FCER2","TCL1A","IL4R","CD72","BACH2","IGHD","IGHM"],
"Bm":["CD27","TNFRSF13B"],
"B Plasma":["CD38","MZB1","PRDM1","IRF4","XBP1"],
"GCB":["TCL1A","IL4R","BACH2","BCL6","RGS13","AICDA","IL21R","MKI67","STMN1","HMGB2","TOP2A"],
}

### 高强泛癌B
看这篇的情况其实B比较难分得清，搞不清的就直接特征基因了
<img src="..\figures\高强泛癌B-UMAP.png">
<img src="..\figures\高强泛癌B-dotplot.png">
marker_genes={
"Bn":["TCL1A","YBX3","FCER2"],
"Bm":["CRIP2","ITGB1"],
"activated B":["CD69","CD83"], #激活了但还没到GCB阶段
"GCB":["TCL1A","SUGCT","NEIL1","NEF2B","LMO2","GMDS","MARCKSL1","STMN1","HMGB2","MKI67"],
"B Plasma":["XBP1","MZB1","JCHAIN"],
"B Plasmablast":["XBP1","MZB1","JCHAIN","STMN1","HMGB2","MKI67"],
}


### 生信加油站

marker_genes = {
"B Naive":["IGHD","IL4R","FCER2","TCL1A"],
"B Memory":["CD27","IGHG1","AIM2","TNFRSF13B"],
#CD27阳性就不再是初始成熟细胞了
"B Plasma":["PRDM1","MZB1"],
"B Plasmablast":["XBP1","SDC1"],#SDC1 is negative marker
"B Short-lived plasma cells":["SDC1"],#SDC1高表达
"B Long-lived plasma cells":["SDC1","STAT3","IKZF3"],
}

### Liu 2023 CancerCell
marker_genes = {
"B Plasma":["MZB1","XBP1"],
"CD19+ B":["CD19"], #CD19表达于除了浆细胞外的绝大多数B细胞
}

### Single Cell Best Practices
marker_genes = {
"Naive CD20+ B": ["MS4A1", "IL4R", "IGHD", "FCRL1", "IGHM"], #CD20(MS4A1)也是不表达于浆细胞
"B1 B": [
    "MS4A1",
    "SSPN",
    "ITGB1",
    "EPHA4",
    "COL4A4",
    "PRDM1",
    "IRF4",
    "CD38",
    "XBP1",
    "PAX5",
    "BCL11A",
    "BLK",
    "IGHD",
    "IGHM",
    "ZNF215",
],  # Note IGHD and IGHM are negative markers
"Transitional B": ["MME", "CD38", "CD24", "ACSM3", "MSI2"],
#MME是CD10,只在过渡B表达，CD38在初始成熟B和记忆B都低表达，CD24也是只在过渡B高表达，同时
#过渡B还应该是CD20+/MS4A1阳性
"Plasma cells": ["MZB1", "HSP90B1", "FNDC3B", "PRDM1", "IGKC", "JCHAIN"],
"Plasmablast": ["XBP1", "RF4", "PRDM1", "PAX5"],  # Note PAX5 is a negative marker
}

## 各文献分的B细胞类型
林东昕2021 NC ESCC把B细胞分为了GC B,activated B,Resting B,浆细胞
Liu 2024 CancerCell ESCC把B细胞分为了滤泡B细胞,Bm,Bn,浆细胞
Chen 2024 CancerCell 结直肠癌(也有PBMC的样本)把B细胞分为了Bn,Bm,GCB,浆细胞
张泽明泛癌B 把B细胞分为了Bn,Bm,GCB,激活的不典型B细胞，浆细胞

## 总结B细胞注释可能用到的marker
marker_genes = {
"Bn":["TCL1A", #在张泽明泛癌中GCB也会表达这个细胞
      "YBX3", #高强泛癌B中的marker,高表达但不是特别特异
      "FCER2", #preGC和Bm也有可能表达
      "IL4R", #生发中心的也可能表达
      "MS4A1", #CD20,区别浆细胞
      "FCRL1", #SCBP的
      "CD72", #张泽明泛癌,Bm和GCB也可能表达
      "BACH2", #张泽明泛癌,Bm和GCB也可能表达
      "IGHD", #普遍表达但不高
      "IGHM"] ,#普遍表达但不高
"Bm":["CD27", #按理来说是阳性，但在张泽明文章里也没有高表达
      "IGHG1", #IgG,第二次免疫应答才有IgG
      "AIM2", #生信加油站的
      "TNFRSF13B", #这个marker特异性应该可
      "CRIP2", #高强泛癌B中的marker,好像不是很特异
      "ITGB1"], #高强泛癌B中的marker,好像不是很特异
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
                "SDC1"],  # Note PAX5 and SDC1 are negative marker

}
