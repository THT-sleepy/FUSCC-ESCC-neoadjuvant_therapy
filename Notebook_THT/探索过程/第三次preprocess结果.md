# 第三次对121个样本做预处理及注释的过程和结果讨论         (What)

* Nov 6, 2025                                 (When)
* liyi /data1/liyi/zhaoyue/FUSCC-ESCC-neoadjuvant_thrapy           (Where)
* 单细胞数据需要preprocess

第二次预处理后，我进行了后续注释。存在两个问题，一个是发现有些样本
实际上是对照组没有使用新辅助治疗(当然这和注释无关)，另一个就是免疫细胞
不太好分得比较细，我觉得有可能的原因是之前mt的标准5%有点严格，所以想着
再看看

```
import scanpy as sc
import matplotlib.pyplot as plt
adata = sc.read("RData/mtremain_merge.h5ad")
plt.figure(figsize=(10, 6))
plt.hist(adata.obs['pct_counts_mt'],bins=40,color='skyblue',edgecolor='black')
plt.xlabel('pct_counts_mt')
plt.ylabel('Frequency')
plt.title('Distribution of pct_counts_mt')
plt.grid(alpha=0.3)
plt.show()
```

<img src="..\figures\mt_counts_直方图.png">
5%没有问题,把多余的样本去了就可以了
