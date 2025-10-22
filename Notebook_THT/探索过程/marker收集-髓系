# Myeloid细胞marker收集        (What)

* Oct 22, 2025                                 (When)
* liyi /data1/liyi/zhaoyue/FUSCC-ESCC-neoadjuvant_thrapy           (Where)
* 单细胞需要注释后才能进行后续分析 (Why)

### 背景知识
髓系主要包括单核细胞/巨噬细胞,肥大细胞,树突状细胞,粒细胞这几个大类

### 各文献注释情况
zheng 2020 NC,把髓系细胞分为了DC,Mono,Macro,MDSC(Myeloid derived suppressor cells)
林东昕 2021 NC,把髓系细胞分为了Mono,cDC(conventional DC),pDC(plasmcytoid DC),
tDC(tolerogenic DC),Mast,TAM(tumor associated macrophage)
胸科 2024 Cancer Cell,把髓系细胞分为了Macro,Mast,Mono,Neutrophils,cDC,pDC,mDC(myeloid DC)
张泽明 泛癌 2021 Cell,Mono/Macro,cDC,pDC,Mast
张泽明 结肠癌,2020 Cell,Mast,pDC,cDC,Mono,Macro,Mono-like,TAM

### marker
### 生信加油站
marker_genes = {
"Mono CD14+":["CD14","S100A9","S100A8"],
"Mono CD16+":["FCGR3A","LST1","LILRB2"],
"pDC":["LILRA4","GZMB","IL3RA"],
"cDC1":["CLEC9A","FLT3","IDO1"],
"cDC2":["CD1C","FCER1A","HLA-DQA1"],
"cDC3":["LAMP3","CCR7","FSCN1"],
"Mast":["KIT","TPSAB1","CPA3"],
"Neutrophil":["CSF3R","S100A9","FGR3B",
              "ALPL","CXCR1"],
"Cycling":["MKI67","STMN1","TOP2A"],
}

### single cell best practices
marker_genes = {
    "Mono CD14+": ["FCN1", "CD14"],
    "Mono CD16+": ["TCF7L2", "FCGR3A", "LYN"],
    "cDC1": ["CLEC9A", "CADM1"],
    "cDC2": [
        "CST3",
        "COTL1",
        "LYZ",
        "DMXL2",
        "CLEC10A",
        "FCER1A",
    ],  # Note: DMXL2 should be negative
    "pDC": ["GZMB", "IL3RA", "COBLL1", "TCF4"],
}

### zheng 2020 NC
marker_genes = {
"mDC":["CD1C","FCER1A"],
"pDC":["CLEC4C"],
"Mast":["TPSB2","CPA3"]
}

### 张泽明泛癌
marker_genes = {
"Mast":["KIT","TPSAB1","CPA3"],
"pDC":["LIRA4","GZMB","IL3RA"],
"cDC1":["CLEC9A","FLT3","IDO1"],
"cDC2":["CD1C","FCER1A","HLA-DQA1"],
"cDC3":["LAMP3","CCR7","FSCN1"],
"Mono CD14+":["FCN1","S100A8","S100A9"],
"Mono CD16+":["FCGR3A","LST1","LILRB2"],
"Macro":["INHBA","IL1RN","CCL4","NLRP3","EREG","IL1B","LYVE1","PLTP","SEPP1","C1QC","C1QA","APOE"]
}

### 胸科
marker_genes = {
"Macro":["CD68"]
}

### 张泽明结肠癌
{"Mast": ["GATA2", "EPAS1", "PBX1","TPSAB1", "CPA3", "TPSD1", "MS4A2", "KIT", "SIGLEC6", "SIGLEC8", "CD22", "CSF1", "TIMP3"],
"pDC": ["IRF4", "IRF7", "IRF8", "SPIB", "SOX4","LILRA4", "SLC15A4", "PLD4", "CCDC50", "IL3RA", "LY9", "SELL", "GAS6"],
"cDC2": ["HMGA1", "PFDN1", "IRF4","CD1C", "CLEC10A", "FCGR2B", "S1RPA", "ADAM8", "FCER1A", "AXL", "ADAM28", "LY86", "TIMM13", "ARAF"],
"cDC1":["BATF3", "ID2", "ETV3","XCR1", "CLEC9A", "PTDSS1", "SCARB1", "IL6ST", "CD40", "TNFRSF10B", "IDO1", "CST7", "CLIC2", "NET1", "ANXA6"],
"Mono CD14+":["CEBPD","FCN1", "CD14", "CD36", "SELL", "S100A8", "S100A12", "CLEC12A", "MS4A6A", "CXCL14"],
"Mono CD16+":["TCF7L2", "KLF3", "IKZF1", "FLI1","FCN1", "FCGR3A", "FCGR3B", "LILRA1B1", "CX3CR1", "IFITM1", "ICAM2", "MTSS1", "CDKN1C", "CDH23", "SLC44A2"],
"Mono-CD14CD16":["CEBPD", "ZFP36L2","FCN1", "CD14", "FCGR3A", "CD300E", "FAM65B"],
}
"Macro": ["NFKB2","NLRP3", "AQP9", "METRNL",
          "GPR132", "ANPEP", "PLAUR", "EREG",
          "VEGFA", "PTGS2","KLF2", "KLF4",
          "EGR1", "CREM", "HIF1A","PLTP",
          "LYVE1", "IL10", "STAB1", "SEPP1",
          "SNX6", "CTSL", "THBS1", "MARCO",
          "CXCL3", "CALR", "CD163", "CD36",
          "HES1", "REL", "BTG2","IL1B",
          "C1QA", "C1QC", "C1QB", "CCL3L3",
          "HLA-DQA1","HLA-DQA2", "HLA-DQB1",
          "GPR183", "CD83", "CLEC10A", "MRC1",
          "MAFB", "ATF3", "MAF", "MEF2A", "PA2G4",
          "C1QA/B/C", "MERTK", "FPR3", "TREM2", "MS4A4A",
          "SLCO2B1", "NRP1", "SLAMF8", "FCGR1A", "ENG", "SIRPA",
          "CEBPB","SPP1", "IL1RN", "OLR1", "CXCL2", "VEGFA",
          "EREG", "C15ORF48", "GPNMB", "PHLDA1", "AQP9", "TNS3", "NDRG1"],
}

### 并集
marker_genes ={
"Mono CD14+": ["CD14", "S100A9", "S100A8",
                "FCN1", "CD36", "SELL", "S100A12",
                "CLEC12A", "MS4A6A", "CXCL14", "CEBPD"],
"Mono CD16+": ["FCGR3A", "LST1", "LILRB2",
              "TCF7L2", "LYN", "FCN1", "FCGR3B",
              "LILRA1B1", "CX3CR1", "IFITM1",
              "ICAM2", "MTSS1", "CDKN1C",
              "CDH23", "SLC44A2", "KLF3",
              "IKZF1", "FLI1"],
"pDC": ["LILRA4", "GZMB", "IL3RA",
        "COBLL1", "TCF4", "CLEC4C",
        "IRF4", "IRF7", "IRF8",
        "SPIB", "SOX4", "SLC15A4",
        "PLD4", "CCDC50", "LY9",
        "SELL", "GAS6", "LIRA4"],
"cDC1": ["CLEC9A", "FLT3", "IDO1",
        "CADM1", "BATF3", "ID2",
        "ETV3", "XCR1", "PTDSS1",
        "SCARB1", "IL6ST", "CD40",
        "TNFRSF10B", "CST7", "CLIC2",
        "NET1", "ANXA6"],
"cDC2": ["CD1C", "FCER1A", "HLA-DQA1",
        "CST3", "COTL1", "LYZ", "DMXL2",
        "CLEC10A", "HMGA1", "PFDN1",
        "IRF4", "FCGR2B", "S1RPA",
        "ADAM8", "AXL", "ADAM28",
        "LY86", "TIMM13", "ARAF"],
"cDC3": ["LAMP3", "CCR7", "FSCN1"],
"Mast": ["KIT", "TPSAB1", "CPA3",
         "TPSB2", "GATA2", "EPAS1",
         "PBX1", "TPSD1", "MS4A2",
         "SIGLEC6", "SIGLEC8", "CD22",
         "CSF1", "TIMP3"],
"Neutrophil": ["CSF3R", "S100A9", "FGR3B",
               "ALPL", "CXCR1"],
"Cycling": ["MKI67", "STMN1", "TOP2A"],
"mDC": ["CD1C", "FCER1A"],
"Macro": ["INHBA", "IL1RN", "CCL4",
          "NLRP3", "EREG", "IL1B",
          "LYVE1", "PLTP", "SEPP1",
          "C1QC", "C1QA", "APOE", "CD68"],
"Mono-CD14CD16": ["CEBPD", "ZFP36L2", "FCN1", "CD14", "FCGR3A", "CD300E", "FAM65B"],
}
