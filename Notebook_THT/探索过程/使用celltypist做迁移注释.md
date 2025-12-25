# CellTypist迁移注释       (What)

* Dec 21, 2025                                 (When)
* liyi /data1/liyi/zhaoyue/FUSCC-ESCC-neoadjuvant_thrapy           (Where)
* 对T/CD8T做一下迁移注释 (Why)

### 导入模块
```python
import scanpy as sc
import celltypist
import time
import numpy as np
```

### 1 ref使用#1,是一个CD8T的单细胞数据

#### 导入我们的CD8T并进行标化
```python
#先把CD8T细胞都提取出来
adata=sc.read("RData/1128_final_escc120.h5ad")
t_cls=[
'c07_CD8_Tem/Teff_TNFSF9',
'c08_CD8_Tem/Teff_HSPA1B',
'c09_CD8_Tem/Teff_GNLY',
'c10_CD8_Tem/Teff_GZMK',
'c11_CD8_Tex_CXCL13',
'c12_CD8_Temra_FGFBP2',
'c13_CD8_IL7R',
'c14_Tprf_MKI67',]
adata_t=adata[adata.obs['minor_celltype'].isin(t_cls)].copy()
adata_t.X = adata_t.layers["soupX_counts"] # set adata.X to raw counts
sc.pp.normalize_total(adata_t, target_sum=10**4 ) # normalize to 10,000 counts per cell
sc.pp.log1p(adata_t) # log-transform
# make .X dense instead of sparse, for compatibility with celltypist:
adata_t.X = adata_t.X.toarray()
```

#### 导入参考数据集并进行标化
```python
adata_andreatta=sc.read("input/ref_TILAtlas_human_cd8t.h5ad") #一共是16k细胞,不用downsample
sc.pp.normalize_total(adata_andreatta, target_sum = 1e4)
sc.pp.log1p(adata_andreatta)
# adata_andreatta.obs['functional.cluster'].unique()
#Categories (9, object): ['CD8_Tex', 'CD8_Tpex', 'CD8_EffectorMemory', 'CD8_EarlyActiv', ...,
#                         'CD4_NaiveLike', 'Tfh', 'Th1', 'Treg']

```
#### 先快速得到一个粗糙的模型
```python
t_start = time.time()
model_fs = celltypist.train(adata_andreatta, 'functional.cluster', n_jobs = 10, max_iter = 5, use_SGD = True)
t_end = time.time()
print(f"Time elapsed: {t_end - t_start} seconds")
```

#### 用这个模型的前300基因作为特征基因(因为一共只有16k细胞，所以设高一点)
```python
gene_index = np.argpartition(np.abs(model_fs.classifier.coef_), -300, axis = 1)[:, -300:]
gene_index = np.unique(gene_index)
```

#### 得到最终的模型并保存
```python
t_start = time.time()
model = celltypist.train(adata_andreatta[:,gene_index], 'functional.cluster', check_expression = False, n_jobs = 10, max_iter = 100)
t_end = time.time()
print(f"Time elapsed: {(t_end - t_start)/60} minutes")
model.write('/home/liyi301/.celltypist/data/models/model_from_Andreatta_2021.pkl')
```

#### 做预测
```python
t_start = time.time()
predictions = celltypist.annotate(adata_t, model = '/home/liyi301/.celltypist/data/models/model_from_Andreatta_2021.pkl', majority_voting = True)
t_end = time.time()
print(f"Time elapsed: {t_end - t_start} seconds")
celltypist.dotplot(predictions, use_as_reference = 'minor_celltype', use_as_prediction = 'majority_voting',save="Andreatta_TA.png")
```
<img src="..\figures\CellTypist_dotplot_Andreatta_TA.png"><br>
结果还挺不错的,但就是说我们是否应该考虑把c10换成Tem/Teff呢？要的(我注释的时候以为KLRG1表达就不能定义为Tem,但实际上不是的)

* #1 Andreatta, M., Corria-Osorio, J., Müller, S. et al. Interpretation of T cell states from single-cell transcriptomics data using reference atlases. Nat Commun 12, 2965 (2021). https://doi.org/10.1038/s41467-021-23324-4

adata = sc.read_10x_mtx(
            "input/MD01-004_MANA0_4/filtered_feature_bc_matrix/",
            var_names='gene_symbols'
        )
