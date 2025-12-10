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
adata = sc.read("RData/Endothelial cell__umapmd0.4.h5ad")

##需要标化到10000
adata_celltypist = adata.copy() # make a copy of our adata 
adata_celltypist.X = adata.layers["soupX_counts"] # set adata.X to raw counts 
sc.pp.normalize_total(adata_celltypist, target_sum=10**4 ) # normalize to 10,000 counts per cell 
sc.pp.log1p(adata_celltypist) # log-transform
# make .X dense instead of sparse, for compatibility with celltypist:
adata_celltypist.X = adata_celltypist.X.toarray()

model = models.Model.load(model="Cells_Intestinal_Tract.pkl")

predictions = celltypist.annotate(
    adata_celltypist, model=model, majority_voting=True ) 

predictions_adata = predictions.to_adata() 

adata.obs["celltypist_cell_label"] = predictions_adata.obs.loc[adata.obs.index, "majority_voting" ] 

adata.obs["celltypist_conf_score"] = predictions_adata.obs.loc[adata.obs.index, "conf_score"]

sc.pl.umap(adata, color=["celltypist_cell_label", "celltypist_conf_score"], frameon=False, sort_order=False, wspace=1,save="1121_celltypist__Endo.png")

