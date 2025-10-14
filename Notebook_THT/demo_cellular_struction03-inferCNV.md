## 修改路径为conda路径


```python
import os
```


```python
print(os.environ['PATH'])
```

    /usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/snap/bin



```python
conda_path = "/home/data/t190513/miniconda3/envs/sc_infercnv/bin:/home/data/t190513/miniconda3/condabin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/snap/bin:/home/data/t190513/miniconda3/bin"
```


```python
os.environ['PATH'] = conda_path + os.environ['PATH']
```


```python
print(os.environ['PATH'])
```

    /home/data/t190513/miniconda3/envs/sc_infercnv/bin:/home/data/t190513/miniconda3/condabin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/snap/bin:/home/data/t190513/miniconda3/bin/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/snap/bin



```python
import scanpy as sc
import infercnvpy as cnv
import matplotlib.pyplot as plt
import warnings
import anndata as ad
import pandas as pd
warnings.simplefilter("ignore")

sc.settings.set_figure_params(figsize=(5, 5))
```


```python
adata = sc.read("RData/mergedemo_annotation.h5ad")
```


```python
adata.var
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>gene_ids</th>
      <th>feature_types</th>
      <th>mt</th>
      <th>ribo</th>
      <th>hb</th>
      <th>n_cells_by_counts</th>
      <th>mean_counts</th>
      <th>log1p_mean_counts</th>
      <th>pct_dropout_by_counts</th>
      <th>total_counts</th>
      <th>log1p_total_counts</th>
      <th>n_cells</th>
      <th>highly_variable</th>
      <th>means</th>
      <th>dispersions</th>
      <th>dispersions_norm</th>
      <th>highly_variable_nbatches</th>
      <th>highly_variable_intersection</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>PEX10</th>
      <td>ENSG00000157911</td>
      <td>Gene Expression</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>2590</td>
      <td>0.314726</td>
      <td>0.273628</td>
      <td>76.940883</td>
      <td>3535.0</td>
      <td>8.170751</td>
      <td>2209</td>
      <td>False</td>
      <td>0.120009</td>
      <td>0.218987</td>
      <td>-0.590220</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>PEX14</th>
      <td>ENSG00000142655</td>
      <td>Gene Expression</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>2347</td>
      <td>0.263711</td>
      <td>0.234052</td>
      <td>79.104345</td>
      <td>2962.0</td>
      <td>7.993958</td>
      <td>1791</td>
      <td>False</td>
      <td>0.145328</td>
      <td>0.346840</td>
      <td>-0.177943</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>PLCH2</th>
      <td>ENSG00000149527</td>
      <td>Gene Expression</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>841</td>
      <td>0.083244</td>
      <td>0.079961</td>
      <td>92.512464</td>
      <td>935.0</td>
      <td>6.841615</td>
      <td>641</td>
      <td>False</td>
      <td>0.043936</td>
      <td>0.466124</td>
      <td>0.109270</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>SPSB1</th>
      <td>ENSG00000171621</td>
      <td>Gene Expression</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>3395</td>
      <td>0.481036</td>
      <td>0.392742</td>
      <td>69.773860</td>
      <td>5403.0</td>
      <td>8.594895</td>
      <td>2570</td>
      <td>False</td>
      <td>0.170279</td>
      <td>0.536394</td>
      <td>0.116207</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>SLC2A5</th>
      <td>ENSG00000142583</td>
      <td>Gene Expression</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>208</td>
      <td>0.018786</td>
      <td>0.018611</td>
      <td>98.148148</td>
      <td>211.0</td>
      <td>5.356586</td>
      <td>154</td>
      <td>False</td>
      <td>0.008491</td>
      <td>0.642713</td>
      <td>0.460389</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>MT-ND4L</th>
      <td>ENSG00000212907</td>
      <td>Gene Expression</td>
      <td>True</td>
      <td>False</td>
      <td>False</td>
      <td>10221</td>
      <td>13.706642</td>
      <td>2.688299</td>
      <td>9.001068</td>
      <td>153953.0</td>
      <td>11.944409</td>
      <td>6843</td>
      <td>False</td>
      <td>1.711212</td>
      <td>1.157528</td>
      <td>-0.112210</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>MT-ND4</th>
      <td>ENSG00000198886</td>
      <td>Gene Expression</td>
      <td>True</td>
      <td>False</td>
      <td>False</td>
      <td>10150</td>
      <td>13.478365</td>
      <td>2.672655</td>
      <td>9.633191</td>
      <td>151389.0</td>
      <td>11.927615</td>
      <td>6839</td>
      <td>False</td>
      <td>1.657537</td>
      <td>1.110053</td>
      <td>-0.246822</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>MT-ND5</th>
      <td>ENSG00000198786</td>
      <td>Gene Expression</td>
      <td>True</td>
      <td>False</td>
      <td>False</td>
      <td>9661</td>
      <td>6.454861</td>
      <td>2.008866</td>
      <td>13.986823</td>
      <td>72501.0</td>
      <td>11.191369</td>
      <td>6414</td>
      <td>False</td>
      <td>1.461008</td>
      <td>1.053399</td>
      <td>-0.191707</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>MT-ND6</th>
      <td>ENSG00000198695</td>
      <td>Gene Expression</td>
      <td>True</td>
      <td>False</td>
      <td>False</td>
      <td>7069</td>
      <td>2.866720</td>
      <td>1.352407</td>
      <td>37.063746</td>
      <td>32199.0</td>
      <td>10.379722</td>
      <td>4793</td>
      <td>False</td>
      <td>1.066389</td>
      <td>1.389778</td>
      <td>0.884135</td>
      <td>1</td>
      <td>False</td>
    </tr>
    <tr>
      <th>MT-CYB</th>
      <td>ENSG00000198727</td>
      <td>Gene Expression</td>
      <td>True</td>
      <td>False</td>
      <td>False</td>
      <td>10456</td>
      <td>26.456998</td>
      <td>3.312621</td>
      <td>6.908832</td>
      <td>297165.0</td>
      <td>12.602046</td>
      <td>7032</td>
      <td>False</td>
      <td>2.396976</td>
      <td>1.615468</td>
      <td>-0.197952</td>
      <td>0</td>
      <td>False</td>
    </tr>
  </tbody>
</table>
<p>15674 rows × 18 columns</p>
</div>



### 需要在var属性有ensg,chromosome,start,end四列


```python
adata.var['ensg'] = adata.var['gene_ids']
```

### 得到基因坐标文件


```python
# 读取 GTF 文件并提取基因坐标
gtf_path = 'tmp/gencode.v47.basic.annotation.gtf.gz'  # 替换为你的文件路径
df = pd.read_csv(gtf_path, sep='\t', comment='#', 
                 names=['chrom', 'source', 'feature', 'start', 'end', 
                        'score', 'strand', 'frame', 'attr'])

# 筛选基因行并提取 ENSG ID 和基因名
df = df[df['feature'] == 'gene'].copy()

# 使用正则表达式提取 gene_id 和 gene_name
df['ensg_with_version'] = df['attr'].str.extract(r'gene_id "([^"]+)"').iloc[:,0]
df['symbol'] = df['attr'].str.extract(r'gene_name "([^"]+)"').iloc[:,0]

# --- 核心步骤：去除版本号 ---
# 使用 str.split() 按 '.' 分割，取第一部分
df['ensg'] = df['ensg_with_version'].str.split('.').str[0]

# 选择需要的列并保存
df = df[['ensg', 'chrom', 'start', 'end', 'symbol']]

# 可选：一个 ENSG (去版本号后) 可能对应多个条目（如不同组装版本），我们取第一个
df = df.drop_duplicates(subset='ensg', keep='first')

df.to_csv('ensembl_gene_coords_no_version.tsv', sep='\t', index=False)
```


```python
df
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>ensg</th>
      <th>chrom</th>
      <th>start</th>
      <th>end</th>
      <th>symbol</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>ENSG00000290825</td>
      <td>chr1</td>
      <td>11121</td>
      <td>24894</td>
      <td>DDX11L16</td>
    </tr>
    <tr>
      <th>6</th>
      <td>ENSG00000223972</td>
      <td>chr1</td>
      <td>12010</td>
      <td>13670</td>
      <td>DDX11L1</td>
    </tr>
    <tr>
      <th>14</th>
      <td>ENSG00000310526</td>
      <td>chr1</td>
      <td>14356</td>
      <td>30744</td>
      <td>WASH7P</td>
    </tr>
    <tr>
      <th>108</th>
      <td>ENSG00000227232</td>
      <td>chr1</td>
      <td>14696</td>
      <td>24886</td>
      <td>WASH7P</td>
    </tr>
    <tr>
      <th>120</th>
      <td>ENSG00000278267</td>
      <td>chr1</td>
      <td>17369</td>
      <td>17436</td>
      <td>MIR6859-1</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>2225559</th>
      <td>ENSG00000198695</td>
      <td>chrM</td>
      <td>14149</td>
      <td>14673</td>
      <td>MT-ND6</td>
    </tr>
    <tr>
      <th>2225564</th>
      <td>ENSG00000210194</td>
      <td>chrM</td>
      <td>14674</td>
      <td>14742</td>
      <td>MT-TE</td>
    </tr>
    <tr>
      <th>2225567</th>
      <td>ENSG00000198727</td>
      <td>chrM</td>
      <td>14747</td>
      <td>15887</td>
      <td>MT-CYB</td>
    </tr>
    <tr>
      <th>2225572</th>
      <td>ENSG00000210195</td>
      <td>chrM</td>
      <td>15888</td>
      <td>15953</td>
      <td>MT-TT</td>
    </tr>
    <tr>
      <th>2225575</th>
      <td>ENSG00000210196</td>
      <td>chrM</td>
      <td>15956</td>
      <td>16023</td>
      <td>MT-TP</td>
    </tr>
  </tbody>
</table>
<p>78724 rows × 5 columns</p>
</div>



### 映射坐标到 AnnData


```python
# 读取处理好的坐标文件
coords = pd.read_csv('ensembl_gene_coords_no_version.tsv', sep='\t')

# 转换坐标列为整数
coords['start'] = coords['start'].astype(int)
coords['end'] = coords['end'].astype(int)

# 映射到 AnnData
# 确保 coords 的索引是 'ensg'
mapping = coords.set_index('ensg')
adata.var = adata.var.join(mapping, on='ensg')

# 重命名列以符合你的要求
adata.var = adata.var.rename(columns={'chrom': 'chromosome'})
```


```python
adata.var.loc[:, ["ensg", "chromosome", "start", "end"]].head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>ensg</th>
      <th>chromosome</th>
      <th>start</th>
      <th>end</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>PEX10</th>
      <td>ENSG00000157911</td>
      <td>chr1</td>
      <td>2403964.0</td>
      <td>2413797.0</td>
    </tr>
    <tr>
      <th>PEX14</th>
      <td>ENSG00000142655</td>
      <td>chr1</td>
      <td>10472288.0</td>
      <td>10630758.0</td>
    </tr>
    <tr>
      <th>PLCH2</th>
      <td>ENSG00000149527</td>
      <td>chr1</td>
      <td>2425980.0</td>
      <td>2505532.0</td>
    </tr>
    <tr>
      <th>SPSB1</th>
      <td>ENSG00000171621</td>
      <td>chr1</td>
      <td>9292894.0</td>
      <td>9369532.0</td>
    </tr>
    <tr>
      <th>SLC2A5</th>
      <td>ENSG00000142583</td>
      <td>chr1</td>
      <td>9035106.0</td>
      <td>9088478.0</td>
    </tr>
  </tbody>
</table>
</div>




```python
sc.pl.umap(adata, color="manual_celltype_level1")
```


    
![png](output_17_0.png)
    



```python
adata
```




    AnnData object with n_obs × n_vars = 28911 × 15674
        obs: 'n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'pct_counts_in_top_20_genes', 'total_counts_mt', 'log1p_total_counts_mt', 'pct_counts_mt', 'total_counts_ribo', 'log1p_total_counts_ribo', 'pct_counts_ribo', 'total_counts_hb', 'log1p_total_counts_hb', 'pct_counts_hb', 'outlier', 'mt_outlier', 'scDblFinder_score', 'scDblFinder_class', 'sample', 'leiden_res0_1', 'leiden_res0_25', 'leiden_res0_5', 'leiden_res1', 'manual_celltype_level1', 'celltypist_cell_label_coarse', 'celltypist_conf_score_coarse', 'celltypist_cell_label_fine', 'celltypist_conf_score_fine'
        var: 'gene_ids', 'feature_types', 'mt', 'ribo', 'hb', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts', 'n_cells', 'highly_variable', 'means', 'dispersions', 'dispersions_norm', 'highly_variable_nbatches', 'highly_variable_intersection', 'ensg', 'chromosome', 'start', 'end', 'symbol'
        uns: 'celltypist_cell_label_coarse_colors', 'celltypist_cell_label_fine_colors', 'dea_leiden_res0_1', 'dea_leiden_res0_1_filtered', 'dendrogram_leiden_res0_1', 'hvg', 'leiden_res0_1', 'leiden_res0_1_colors', 'leiden_res0_1_filtered', 'leiden_res0_25', 'leiden_res0_25_colors', 'leiden_res0_5', 'leiden_res0_5_colors', 'leiden_res1', 'leiden_res1_colors', 'manual_celltype_level1_colors', 'neighbors', 'pca', 'sample_colors', 'scDblFinder_class_colors', 'tsne', 'umap'
        obsm: 'X_pca', 'X_pca_harmony', 'X_tsne', 'X_umap'
        varm: 'PCs'
        layers: 'counts', 'log1p_norm', 'soupX_counts'
        obsp: 'connectivities', 'distances'




```python
# We provide all immune cell types as "normal cells".
cnv.tl.infercnv(
    adata,
    reference_key="manual_celltype_level1",
    reference_cat=[
        "Immune Cell",
        "Fibro Cell",
        "Endo Cell",
    ],
    window_size=250,
)
```

    WARNING: Skipped 146 genes because they don't have a genomic position annotated. 


    100%|██████████| 6/6 [01:13<00:00, 12.28s/it]


### 绘制热图 


```python
cnv.pl.chromosome_heatmap(adata, groupby="manual_celltype_level1")
```


    
![png](output_21_0.png)
    


### Clustering by CNV profiles and identifying tumor cells


```python
cnv.tl.pca(adata)
cnv.pp.neighbors(adata)
cnv.tl.leiden(adata)
```


```python
cnv.pl.chromosome_heatmap(adata, groupby="cnv_leiden", dendrogram=True)
```

    WARNING: dendrogram data not found (using key=dendrogram_cnv_leiden). Running `sc.tl.dendrogram` with default parameters. For fine tuning it is recommended to run `sc.tl.dendrogram` independently.
    WARNING: Groups are not reordered because the `groupby` categories and the `var_group_labels` are different.
    categories: 0, 1, 2, etc.
    var_group_labels: chr1, chr2, chr3, etc.



    
![png](output_24_1.png)
    


### UMAP plot of CNV profiles


```python
cnv.tl.umap(adata)
cnv.tl.cnv_score(adata)
```


```python
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(11, 11))
ax4.axis("off")
cnv.pl.umap(
    adata,
    color="cnv_leiden",
    legend_loc="on data",
    legend_fontoutline=2,
    ax=ax1,
    show=False,
)
cnv.pl.umap(adata, color="cnv_score", ax=ax2, show=False)
cnv.pl.umap(adata, color="manual_celltype_level1", ax=ax3)
```


    
![png](output_27_0.png)
    



```python
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 11), gridspec_kw={"wspace": 0.5})
ax4.axis("off")
sc.pl.umap(adata, color="cnv_leiden", ax=ax1, show=False)
sc.pl.umap(adata, color="cnv_score", ax=ax2, show=False)
sc.pl.umap(adata, color="manual_celltype_level1", ax=ax3)
```


    
![png](output_28_0.png)
    



```python
adata.obs["cnv_status"] = "normal"
adata.obs.loc[adata.obs["cnv_leiden"].isin(["0", "1", "2", "3", "7"]), "cnv_status"] = (
    "tumor"
)
```


```python
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={"wspace": 0.5})
cnv.pl.umap(adata, color="cnv_status", ax=ax1, show=False)
sc.pl.umap(adata, color="cnv_status", ax=ax2)
```


    
![png](output_30_0.png)
    



```python
cnv.pl.chromosome_heatmap(adata[adata.obs["cnv_status"] == "tumor", :])
```


    
![png](output_31_0.png)
    



```python
cnv.pl.chromosome_heatmap(adata[adata.obs["cnv_status"] == "normal", :])
```


    
![png](output_32_0.png)
    

