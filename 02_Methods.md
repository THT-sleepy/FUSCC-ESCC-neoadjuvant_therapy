# Methods

## Annotation
去掉了同时表达T和上皮、T和中性粒marker基因的cluster

## scTCR-seq data processiing
The scTCR-seq data was processed using the Cell Ranger toolkit, which aligned reads to the GRCh38 human reference genome and
assembled TCR sequences. We then filtered the preliminary assembled CDR3 nucleotide sequences of TCR α- and β-chain to retain
only those high-confident, full-length, and productive ones, and assigned each with a valid cell barcode and an unambiguous chain
type. For each cell, one pair of α and β chains showing the highest UMI counts was kept for identity assignment. Cells with identical
TCR pairs were regarded to be of the same clonotype.For clonotypes containing both CD4⁺ T cells and CD8⁺ T cells, we eliminated either CD4⁺ or CD8⁺ T cells within the clonotypes following the method described by Guo et al. (2025, Cancer Cell). Specifically, we calculated the ratio of the number of CD4⁺ T cells to that of CD8⁺ T cells in such clonotypes. If the CD4/CD8 ratio was greater than 5, CD8⁺ T cells were excluded; if the ratio was less than 0.2, CD4⁺ T cells were excluded. Otherwise, the entire clonotype was removed from subsequent analyses. After elimination,
78.4% (118203/150719) of T cells with high-quality TCR α-β pair were kept, resulting in 76996 clonotypes. The size of clonotypes ranged from 1 to 23 or 873 for CD4+ and CD8+ T cells,respectively, with 7977 clonotypes comprising more than one cell (Figures ? and ?). We calculated the STARTRAC-expa,STARTRAC-migr and STARTRAC-trans index with   the Startrac.run function from the R package Startrac (version 0.1.0)[citation]

## Differential expression analysis
The differential expression analysis between two cellular populations was performed using the ? function from Scanpy.
Specifically, the wilcox-test was used to evaluate the significance of each gene with multiple hypothesis corrections performed with the Benjamini–Hochberg procedure. Adjusted p values below 0.05 were used to distinguish significantly differentially expressed genes.

## Definition of CD8+ T cell signature scores
We defined the Tex progenitor (Texp) score and the exhaustion score using signatures from ??(liu et al.2022 张泽民 Nature Cancer).Specifically, the exhaustion score was calculated using HAVCR2, ENTPD1, LAYN, and LAG3. The differentially expressed genes of Texp clusters in the abovementioned NSCLC study was used to calculate the Texp score. The cytotoxicity score was calculated according to the expression level of FOS,JUN,PRF1, GZMB,GZMA,GZMH,NKG7, and GNLY according to a published T cell study[Mathewson, N.D.2021 CELL] and 胸科医院Cancer cell. [Wedefined an MHC II related signature score by genes from the GO term ‘‘MHC class II protein complex’’.没用的话就去掉] All scores were calculated by the Scanpy function ??.



## Sample Collection

## Sequence Analysis
