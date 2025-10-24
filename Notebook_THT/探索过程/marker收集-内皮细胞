# Endo细胞marker收集        (What)

* Oct 23, 2025                                 (When)
* liyi /data1/liyi/zhaoyue/FUSCC-ESCC-neoadjuvant_thrapy           (Where)
* 单细胞需要注释后才能进行后续分析 (Why)

### 先看看各个文献的内皮是怎么分群的
林东昕2021 NC:分为了NEC(normal EC)和TEC(tumor EC)两大类<br> 我们的不是癌和癌旁所以分不了
胸科2023 Cancer Cell,也是分成了NEC和TEC
Chen结直肠癌 Cancer Cell,直接用的特征基因
张泽民泛癌 Arteries,capillaries,hypoxia,lymphatics,tip cell,veins

### markers

### 生信加油站
marker_genes = {
"LEC":["PROX1","PDPN","ALCAM","CD44","VEGFA"],
"VEC":["EPHB4","NR2F2","ACKR1","MMRN1","SELP",
      "VCAM1","POSTN","IGFBP7","CCL14","PRCP"],
"AEC":["GJA4","GJA5","EFNB2","VEGFC","SOX17",
      "DKK2","LTBP4","FBLN5","FN1","MGP",
      "SERPINE2","ENPP2","TFPI","NPR3"],
"TipEC":["CXCL12","CXCR4","ACKR3","LYVE1","DLL4",
        "KCNE3","ESM1","ANGPT2","APLN"],
"CapEC":["CA4","CD36"],
"ISG+EC":["ISG20","IFIT1","IFIT3"],
"prolEC":["MKI67","TOP2A"],
}

### 张泽民泛癌
marker_genes = {
"AEC": ["GJA5", "FBLN5", "GJA4"],
"CapEC": ["CA4", "CD36", "RGCC"],
"Hypoxia": ["MT1X", "MT1E", "MT2A"],
"LEC": ["PROX1", "LYVE1", "CCL21"],
"TipEC": ["COL4A1", "KDR", "ESM1"],
"VEC": ["ACKR1", "SELP", "CLU"]
}

### 合并marker genes
marker_genes = {
    "LEC": ["PROX1", "PDPN", "ALCAM", "CD44", "VEGFA", "LYVE1", "CCL21"],
    "VEC": ["EPHB4", "NR2F2", "ACKR1", "MMRN1", "SELP", "VCAM1", "POSTN", "IGFBP7", "CCL14", "PRCP", "CLU"],
    "AEC": ["GJA4", "GJA5", "EFNB2", "VEGFC", "SOX17", "DKK2", "LTBP4", "FBLN5", "FN1", "MGP", "SERPINE2", "ENPP2", "TFPI", "NPR3"],
    "TipEC": ["CXCL12", "CXCR4", "ACKR3", "LYVE1", "DLL4", "KCNE3", "ESM1", "ANGPT2", "APLN", "COL4A1", "KDR"],
    "CapEC": ["CA4", "CD36", "RGCC"],
    "ISG+EC": ["ISG20", "IFIT1", "IFIT3"],
    "prolEC": ["MKI67", "TOP2A"],
    "Hypoxia": ["MT1X", "MT1E", "MT2A"]
}
