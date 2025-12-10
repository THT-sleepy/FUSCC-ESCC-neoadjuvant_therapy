## 导入模块
import scanpy as sc
import infercnvpy as cnv
import matplotlib.pyplot as plt
import warnings
import anndata as ad
import pandas as pd
warnings.simplefilter("ignore")
sc.settings.set_figure_params(figsize=(5, 5))

sc.settings.set_figure_params(
    dpi=300,
    facecolor="white",
    frameon=False
)
sc.settings.figdir = "plots"
sc.settings.autoshow = False

adata = sc.read("RData/2nd_round_anno_umapfiltered.h5ad")

#添加基因坐标
adata.var['ensg'] = adata.var['gene_ids']

## 读取 GTF 文件并提取基因坐标
gtf_path = 'tmp/gencode.v47.basic.annotation.gtf.gz' # 替换为你的文件路径 
df = pd.read_csv(gtf_path, sep='\t', comment='#',names=['chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attr'])
## 筛选基因行并提取 ENSG ID 和基因名
df = df[df['feature'] == 'gene'].copy()
## 使用正则表达式提取 gene_id 和 gene_name
df['ensg_with_version'] = df['attr'].str.extract(r'gene_id "([^"]+)"').iloc[:,0] 
df['symbol'] = df['attr'].str.extract(r'gene_name "([^"]+)"').iloc[:,0]
## --- 核心步骤：去除版本号 --- 使用 str.split() 按 '.' 分割，取第一部分
df['ensg'] = df['ensg_with_version'].str.split('.').str[0]
## 选择需要的列并保存
df = df[['ensg', 'chrom', 'start', 'end', 'symbol']]
## 可选：一个 ENSG (去版本号后) 可能对应多个条目（如不同组装版本），我们取第一个
df = df.drop_duplicates(subset='ensg', keep='first')
df.to_csv('ensembl_gene_coords_no_version.tsv', sep='\t', index=False)
## 读取处理好的坐标文件
coords = pd.read_csv('ensembl_gene_coords_no_version.tsv', sep='\t')
## 转换坐标列为整数
coords['start'] = coords['start'].astype(int) 
coords['end'] = coords['end'].astype(int)
## 映射到 AnnData 确保 coords 的索引是 'ensg'
mapping = coords.set_index('ensg')
adata.var = adata.var.join(mapping, on='ensg')
## 重命名列以符合你的要求
adata.var = adata.var.rename(columns={'chrom': 'chromosome'})

# 把除了上皮外的细胞都设为reference
cnv.tl.infercnv(adata,
		reference_key="major_celltype",
		reference_cat=[ "T&ILC cell","B cell","Myeloid cell","Fibrocyte","Endothelial cell"],
		window_size=250,
)

# 热图
cnv.pl.chromosome_heatmap(adata, groupby="major_celltype",save="1126_cnv_heatmap.pdf")
cnv.pl.chromosome_heatmap(adata, groupby="major_celltype",save="1126_cnv_heatmap.png")

# Clustering by CNV profiles and identifying tumor cells
cnv.tl.pca(adata)
cnv.pp.neighbors(adata)
cnv.tl.leiden(adata)

cnv.pl.chromosome_heatmap(adata, groupby="cnv_leiden", dendrogram=True,save="1126_clustered_cnv_heatmap.pdf")
cnv.pl.chromosome_heatmap(adata, groupby="cnv_leiden", dendrogram=True,save="1126_clustered_cnv_heatmap.png")

# UMAP plot of CNV profiles
cnv.tl.umap(adata)
cnv.tl.cnv_score(adata)

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(11, 11)) 
ax4.axis("off") 
cnv.pl.umap( adata,
	color="cnv_leiden",
	legend_loc="on data",
	legend_fontoutline=2,
	ax=ax1,
	show=False,)
cnv.pl.umap(adata, color="cnv_score", ax=ax2, show=False)
cnv.pl.umap(adata, color="major_celltype", ax=ax3,save="1126_cnv_umap.pdf")

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(11, 11))
ax4.axis("off")
cnv.pl.umap( adata,
        color="cnv_leiden",
        legend_loc="on data",
        legend_fontoutline=2,
        ax=ax1,
        show=False,)
cnv.pl.umap(adata, color="cnv_score", ax=ax2, show=False)
cnv.pl.umap(adata, color="major_celltype", ax=ax3,save="1126_cnv_umap.png")

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 11), gridspec_kw={"wspace": 0.5}) 
ax4.axis("off") 
sc.pl.umap(adata, color="cnv_leiden", ax=ax1, show=False) 
sc.pl.umap(adata, color="cnv_score", ax=ax2, show=False)
sc.pl.umap(adata, color="major_celltype", ax=ax3,save="1126_cnv_scanpyumap.pdf")

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 11), gridspec_kw={"wspace": 0.5})
ax4.axis("off")
sc.pl.umap(adata, color="cnv_leiden", ax=ax1, show=False)
sc.pl.umap(adata, color="cnv_score", ax=ax2, show=False)
sc.pl.umap(adata, color="major_celltype", ax=ax3,save="1126_cnv_scanpyumap.png")

adata.write("RData/infercnv.h5ad")
