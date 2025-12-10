import scanpy as sc
import celltypist
from celltypist import models


sc.settings.set_figure_params(
    dpi=300,
    facecolor="white",
    frameon=False
)
sc.settings.figdir = "plots"
sc.settings.autoshow = False

##
adata = sc.read("RData/NKT_ILC__umapmd0.4_regressoutmtribo.h5ad")

##需要标化到10000
adata_celltypist = adata.copy() # make a copy of our adata 
adata_celltypist.X = adata.layers["soupX_counts"] # set adata.X to raw counts 
sc.pp.normalize_total(adata_celltypist, target_sum=10**4 ) # normalize to 10,000 counts per cell 
sc.pp.log1p(adata_celltypist) # log-transform
# make .X dense instead of sparse, for compatibility with celltypist:
adata_celltypist.X = adata_celltypist.X.toarray()

model_low = models.Model.load(model="Immune_All_Low.pkl") 
model_high = models.Model.load(model="Immune_All_High.pkl") 

predictions_high = celltypist.annotate(
    adata_celltypist, model=model_high, majority_voting=True ) 

predictions_high_adata = predictions_high.to_adata() 

adata.obs["celltypist_cell_label_coarse"] = predictions_high_adata.obs.loc[adata.obs.index, "majority_voting" ] 

adata.obs["celltypist_conf_score_coarse"] = predictions_high_adata.obs.loc[adata.obs.index, "conf_score"] 

predictions_low = celltypist.annotate(
    adata_celltypist, model=model_low, majority_voting=True
)

predictions_low_adata = predictions_low.to_adata()

adata.obs["celltypist_cell_label_fine"] = predictions_low_adata.obs.loc[ adata.obs.index, "majority_voting"]

adata.obs["celltypist_conf_score_fine"] = predictions_low_adata.obs.loc[adata.obs.index, "conf_score"]

sc.pl.umap(adata, color=["celltypist_cell_label_coarse", "celltypist_conf_score_coarse"], frameon=False, sort_order=False, wspace=1,save="121celltypist_coarse_ILC.png") 

sc.pl.umap(adata, color=["celltypist_cell_label_fine", "celltypist_conf_score_fine"], frameon=False, sort_order=False, wspace=1,save="121celltypist_fine_ILC.png")
