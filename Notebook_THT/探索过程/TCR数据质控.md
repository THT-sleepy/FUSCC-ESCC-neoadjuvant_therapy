# TCR数据质控和预处理         (What)

* Dec 11, 2025                                 (When)
* liyi /data1/liyi/zhaoyue/FUSCC-ESCC-neoadjuvant_thrapy           (Where)
* 查看TCR数据的大致情况，定义clonotype (Why)

### 导入模块
```python
import warnings
warnings.filterwarnings(    
"ignore",    
".*IProgress not found*",
)
warnings.simplefilter(action="ignore", category=FutureWarning)
import pandas as pd
import scanpy as sc
import scirpy as ir
import os
import muon as mu
import matplotlib.pyplot as plt
warnings.simplefilter(action="ignore", category=pd.errors.DtypeWarning)
```

### 把所有样本整合起来(去除7个不要的)
```python
data_root = "input/ALL_TCR/ALL_TCR"
samples = []
for sample_id in os.listdir(data_root):
    # 样本文件夹的完整路径
    sample_dir = os.path.join(data_root, sample_id)
    # 指定tcr路径
    tcr_fname = "filtered_contig_annotations.csv"
    tcr_path = os.path.join(sample_dir,tcr_fname)
    # 构建样本信息字典
    sample_info = {
                "sample_id": sample_id,
                "tcr_path" : tcr_path,
            }
    # 添加到样本列表
    samples.append(sample_info)
adatas = {}
for sample_dict in samples:
        try:
            sample_id = sample_dict['sample_id']
            tcr_path = sample_dict['tcr_path']
            sample_adata = ir.io.read_10x_vdj(tcr_path)
            sample_adata.obs['sample'] = sample_id
            adatas[sample_id] = sample_adata
            print(f"已读取样本: {sample_id}")
        except Exception as e:
            print(f"读取样本 {sample_id} 失败: {str(e)}")

## 合并
adata_tcr = sc.concat(adatas,label="sample",join="outer")
## 确保细胞名唯一
adata_tcr.obs_names_make_unique()
samples_exclude=[
  'LCR_preC1_2265065',
  'SCH_preC1_2265089',
  'RWX_preC1_2265611',
  'ZJK_preC1_2266252',
  'WFC_preC1_2268478',
  'YDY_preC1_2271160',
  'SZR_op_2279307'
]
adata_tcr=adata_tcr[~adata_tcr.obs['sample'].isin(samples_exclude)]
```

### 计算质控参数
```python
ir.pp.index_chains(adata_tcr)
ir.tl.chain_qc(adata_tcr)
```

### 绘制简单的质控图片
```python
_ = ir.pl.group_abundance(adata_tcr, groupby="sample", target_col="chain_pairing",max_cols=0)
plt.savefig(
    "plots/chain_pairing.png",  # 保存路径+文件名
    dpi=300,  # 分辨率（建议300+）
    bbox_inches="tight"  # 自动裁剪空白边距（避免标签被截断）
)
plt.close()
_ = ir.pl.group_abundance(adata_tcr, groupby="sample", target_col="chain_pairing", normalize=True,max_cols=0)
plt.savefig(
    "plots/chain_pairing_normalize.png",  # 保存路径+文件名
    dpi=300,  # 分辨率（建议300+）
    bbox_inches="tight"  # 自动裁剪空白边距（避免标签被截断）
)
plt.close()

```
<img src="../figures/chain_pairing.png">
<img src="../figures/chain_pairing_normalize.png">
从第二张图可以看出来single pair占比都是超过60%,质量都是不错滴

### 计算clonotype
```python
ir.pp.ir_dist(adata_tcr, sequence="aa")
ir.tl.define_clonotype_clusters(
    adata_tcr, sequence="aa", receptor_arms="all", dual_ir="primary_only"
)
```

### 看下最大的一些clone
```python
ir.tl.clonotype_network(adata_tcr, min_cells=50, sequence="aa")
_ = ir.pl.clonotype_network(
    adata_tcr,
    color="sample",
    base_size=10,
    label_fontsize=9,
    panel_size=(10, 10),
    legend_fontsize=15,
)
plt.savefig(
    "plots/clones.png",  # 保存路径+文件名
    dpi=300,  # 分辨率（建议300+）
    bbox_inches="tight"  # 自动裁剪空白边距（避免标签被截断）
)
plt.close()
```
<img src="..\figures\clones.png">
这个图哈哈哈,看来有好多很大的T细胞亚群


### 保存
```python
adata_tcr.write("RData/TCR_merge_1211.h5ad")
```

### 把转录组数据和TCR数据整合到一起
```python
# Load the TCR data
adata_tcr = sc.read("RData/TCR_merge_1211.h5ad")

# Load the associated transcriptomics data
adata_gex = sc.read("RData/1128_final_escc121.h5ad")

mdata = mu.MuData({"gex": adata_gex, "airr": adata_tcr})
```

### 看下TCR是否和CD3数据重合
```python
fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(10, 4), gridspec_kw={"wspace": 0.5})
mu.pl.embedding(mdata, basis="gex:umap", color=["CD3E"], ax=ax0, show=False)
mu.pl.embedding(mdata, basis="gex:umap", color=["airr:receptor_type"], ax=ax1,show=False)
plt.savefig(
    "plots/TCR_umap.png",  # 保存路径+文件名
    dpi=300,               # 分辨率（越高越清晰）
    bbox_inches="tight",   # 自动裁剪空白边距
    facecolor="white"      # 背景色（默认透明，设为白色更易查看）
)
plt.close()
```
确实是重合的，莫问题
<img src="..\figures\TCR_umap.png">

### 画几张scirpy示例的克隆分析图
```python
t_clusters = [ 'c01_CD4_Tn_TCF7',
  'c02_CD4_Tcm_S1PR1',
  'c03_CD4_Tcm_ANXA1',
  'c04_CD4_Tem_NR4A1',
  'c05_CD4_Treg_FOXP3',
  'c06_CD4_TOX',
  'c07_CD8_Tem/Teff_TNFSF9',
  'c08_CD8_Tem/Teff_HSPA1B',
  'c09_CD8_Tem/Teff_GNLY',
  'c10_CD8_Teff_GZMK',
  'c11_CD8_Tex_CXCL13',
  'c12_CD8_Temra_FGFBP2',
  'c13_CD8_IL7R',
  'c14_Tprf_MKI67',]
mdata = mdata[mdata.obs['gex:minor_celltype'].isin(t_clusters)].copy()
### clonal expansion umap
# using default parameters, `ir_dist` will compute nucleotide sequence identity
ir.pp.ir_dist(mdata,sequence="aa")
ir.tl.define_clonotype_clusters(mdata, receptor_arms="all",sequence="aa", dual_ir="primary_only")

ir.tl.clonal_expansion(mdata,target_col="cc_aa_identity")
mu.pl.embedding(mdata, basis="gex:umap", color=["airr:clonal_expansion","airr:cc_aa_identity_size"],show=False)
plt.savefig(
    "plots/Clonal_expansion_umap.png",  # 保存路径+文件名
    dpi=300,               # 分辨率（越高越清晰）
    bbox_inches="tight",   # 自动裁剪空白边距
    facecolor="white"      # 背景色（默认透明，设为白色更易查看）
)
plt.close()

### clonal expansion barplot
_ = ir.pl.clonal_expansion(mdata, target_col="cc_aa_identity", groupby="gex:minor_celltype", breakpoints=(1, 2, 5), normalize=False)
plt.savefig(
    "plots/Clonal_expansion_barplot.png",  # 保存路径+文件名
    dpi=300,               # 分辨率（越高越清晰）
    bbox_inches="tight",   # 自动裁剪空白边距
    facecolor="white"      # 背景色（默认透明，设为白色更易查看）
)
plt.close()

_ = ir.pl.clonal_expansion(mdata, target_col="cc_aa_identity", groupby="gex:minor_celltype", breakpoints=(1, 2, 5), normalize=True)
plt.savefig(
    "plots/Clonal_expansion_barplot_normalized.png",  # 保存路径+文件名
    dpi=300,               # 分辨率（越高越清晰）
    bbox_inches="tight",   # 自动裁剪空白边距
    facecolor="white"      # 背景色（默认透明，设为白色更易查看）
)
plt.close()

### clonal diversity
_ = ir.pl.alpha_diversity(mdata, target_col="cc_aa_identity",metric="normalized_shannon_entropy", groupby="gex:minor_celltype")
plt.savefig(
    "plots/Clonal_diversity.png",  # 保存路径+文件名
    dpi=300,               # 分辨率（越高越清晰）
    bbox_inches="tight",   # 自动裁剪空白边距
    facecolor="white"      # 背景色（默认透明，设为白色更易查看）
)
plt.close()

### clonal abundance
_ = ir.pl.group_abundance(mdata, groupby="airr:cc_aa_identity", target_col="gex:minor_celltype", max_cols=10)
plt.savefig(
    "plots/clonal_abundance.png",  # 保存路径+文件名
    dpi=300,               # 分辨率（越高越清晰）
    bbox_inches="tight",   # 自动裁剪空白边距
    facecolor="white"      # 背景色（默认透明，设为白色更易查看）
)
plt.close()
_ = ir.pl.group_abundance(
    mdata,
    groupby="airr:cc_aa_identity",
    target_col="gex:minor_celltype",
    max_cols=10,
    normalize="gex:sample",
)
plt.savefig(
    "plots/clonal_abundance_normalize.png",  # 保存路径+文件名
    dpi=300,               # 分辨率（越高越清晰）
    bbox_inches="tight",   # 自动裁剪空白边距
    facecolor="white"      # 背景色（默认透明，设为白色更易查看）
)
plt.close()
mdata['gex'].obs['sample_site'] = mdata['gex'].obs['sample'].apply(lambda x: 'blood' if 'pbmc' in x.lower() else 'tumor')
mdata.update()
_ = ir.pl.group_abundance(mdata, groupby="airr:cc_aa_identity", target_col="gex:sample_site", max_cols=15, figsize=(5, 3))
plt.savefig(
    "plots/clonal_abundance_samplesite.png",  # 保存路径+文件名
    dpi=300,               # 分辨率（越高越清晰）
    bbox_inches="tight",   # 自动裁剪空白边距
    facecolor="white"      # 背景色（默认透明，设为白色更易查看）
)
plt.close()
mdata['gex'].obs['patient'] = mdata['gex'].obs['sample'].apply(lambda x: '_'.join([x.split('_')[0], x.split('_')[-1]]))
mdata.update()
_ = ir.pl.group_abundance(
    mdata,
    groupby="airr:cc_aa_identity",
    target_col="gex:patient",
    max_cols=15,
    figsize=(5, 3),
)
plt.savefig(
    "plots/clonal_abundance_patient.png",  # 保存路径+文件名
    dpi=300,               # 分辨率（越高越清晰）
    bbox_inches="tight",   # 自动裁剪空白边距
    facecolor="white"      # 背景色（默认透明，设为白色更易查看）
)
plt.close()
```
<img src="..\figures\Clonal_expansion_umap.png">
<img src="..\figures\Clonal_expansion_barplot.png">
<img src="..\figures\Clonal_expansion_barplot_normalized.png">
<img src="..\figures\Clonal_diversity.png">
<img src="..\figures\clonal_abundance.png">
<img src="..\figures\clonal_abundance_normalize.png">
<img src="..\figures\clonal_abundance_samplesite.png">
<img src="..\figures\clonal_abundance_patient.png">

### 导出用于startrac分析的数据框
mask = mdata.obs["airr:cc_aa_identity"].notna()
obs_filtered = mdata.obs[mask].copy()

df_extracted = pd.DataFrame({
    "Cell_Name": obs_filtered.index,  # Cell_Name 对应 obs 的索引
    "clone.id": obs_filtered["airr:cc_aa_identity"],  # clone.id 对应 airr:cc_aa_identity
    "patient": obs_filtered["gex:patient"],  # patient 对应 gex:patient
    "loc": obs_filtered["gex:sample_site"],  # loc 对应 gex:sample_site
    "majorCluster": obs_filtered["gex:minor_celltype"],  # majorCluster 对应 gex:minor_celltype
    "sample" : obs_filtered["gex:sample"]
})
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

df_extracted['treatment_stage'] = df_extracted['sample'].apply(get_treatment_stage)

df_extracted = df_extracted.reset_index(drop=True)
df_extracted.to_csv("RData/df_Tclonotype_1211_.csv", index=False, encoding="utf-8")

# 查看结果（前5行）
print("提取的数据框前5行：")
print(df_extracted.head())
print(f"\n符合条件的细胞总数：{len(df_extracted)}")
