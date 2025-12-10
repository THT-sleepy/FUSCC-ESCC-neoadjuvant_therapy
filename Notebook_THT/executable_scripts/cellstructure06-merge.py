import scanpy as sc
import pandas as pd


adata_cd4t = sc.read("RData/CD4T_annotated.h5ad")
adata_cd4t.obs['minor_UMAP1'] = adata_cd4t.obsm['X_umap'][:, 0]
adata_cd4t.obs['minor_UMAP2'] = adata_cd4t.obsm['X_umap'][:, 1]

adata_cd8t = sc.read("RData/CD8T_annotated.h5ad")
adata_cd8t.obs['minor_UMAP1'] = adata_cd8t.obsm['X_umap'][:, 0]
adata_cd8t.obs['minor_UMAP2'] = adata_cd8t.obsm['X_umap'][:, 1]

adata_prfT = sc.read("RData/Tprf_annotated.h5ad")

adata_nknkt = sc.read("RData/ILC_annotated.h5ad")
adata_nknkt.obs['minor_UMAP1'] = adata_nknkt.obsm['X_umap'][:, 0]
adata_nknkt.obs['minor_UMAP2'] = adata_nknkt.obsm['X_umap'][:, 1]

adata_b = sc.read("RData/B_annotated.h5ad")
adata_b.obs['minor_UMAP1'] = adata_b.obsm['X_umap'][:, 0]
adata_b.obs['minor_UMAP2'] = adata_b.obsm['X_umap'][:, 1]

adata_monomacro = sc.read("RData/MonoMacro_annotated.h5ad")
adata_monomacro.obs['minor_UMAP1'] = adata_monomacro.obsm['X_umap'][:, 0]
adata_monomacro.obs['minor_UMAP2'] = adata_monomacro.obsm['X_umap'][:, 1]

adata_dc = sc.read("RData/DC_annotated.h5ad")
adata_dc.obs['minor_UMAP1'] = adata_dc.obsm['X_umap'][:, 0]
adata_dc.obs['minor_UMAP2'] = adata_dc.obsm['X_umap'][:, 1]

adata_mast = sc.read("RData/Mast_annotated.h5ad")

adata_fib = sc.read("RData/Fib_annotated.h5ad")
adata_fib.obs['minor_UMAP1'] = adata_fib.obsm['X_umap'][:, 0]
adata_fib.obs['minor_UMAP2'] = adata_fib.obsm['X_umap'][:, 1]

adata_EC = sc.read("RData/EC_annotated.h5ad")
adata_EC.obs['minor_UMAP1'] = adata_EC.obsm['X_umap'][:, 0]
adata_EC.obs['minor_UMAP2'] = adata_EC.obsm['X_umap'][:, 1]

adata_epi = sc.read("RData/Epithelial cell__umapmd0.4.h5ad")


adatas = {
"CD4T":adata_cd4t,
"CD8T":adata_cd8t,
"Tprf":adata_prfT,
"ILC":adata_nknkt,
"B":adata_b,
"MonoMacro":adata_monomacro,
"DC":adata_dc,
"Mast":adata_mast,
"Fib":adata_fib,
"EC":adata_EC,
"Epi":adata_epi,
}

adata = sc.concat(adatas,join="outer")

## 确保细胞名唯一
adata.obs_names_make_unique()


##进行合并后var属性会丢失我们需要自己添加回来
all_var = [x.var for x in adatas.values()]
all_var = pd.concat(all_var,join="outer")
all_var = all_var[~all_var.index.duplicated()]
adata.var = all_var.loc[adata.var_names]

#adata.write("RData/merge2ndanno.h5ad")

#重新normalize和feature selection
adata.X = adata.layers["soupX_counts"]

##normalize
scales_counts = sc.pp.normalize_total(adata, target_sum=None, inplace=False)
# log1p transform
adata.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)

##feature selection
hvg_n=4000
sc.pp.highly_variable_genes(adata, layer="log1p_norm",n_top_genes=hvg_n,batch_key="sample")

adata.X = adata.layers["log1p_norm"]

## pca
sc.pp.pca(adata, svd_solver="arpack", use_highly_variable=True)


sc.external.pp.bbknn(adata, batch_key="sample",n_pcs=30)  # running bbknn 1.6.0
sc.tl.umap(adata,min_dist=0.4)
output_fpath="RData/merge2ndanno.h5ad"
adata.write(output_fpath)

