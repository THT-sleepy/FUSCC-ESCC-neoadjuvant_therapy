##导入模块，设置参数
import sys
import scanpy as sc
sc.settings.verbosity = 0
sc.settings.set_figure_params(dpi=80, facecolor="white", frameon=False)
sc.settings.figdir = "plots"
sc.settings.autoshow = False

##读取1
input_path="RData/hvg4000_PC30_batchdoubletremoved_umapmd0.4.h5ad"
adata = sc.read(input_path)

##聚类1
sc.tl.leiden(adata, key_added="leiden_res0_2", resolution=0.2)
sc.tl.leiden(adata, key_added="leiden_res0_6", resolution=0.6)
sc.tl.leiden(adata, key_added="leiden_res1", resolution=1.0)
sc.tl.leiden(adata, key_added="leiden_res1_2", resolution=1.2)
sc.tl.leiden(adata, key_added="leiden_res2", resolution=2)



##绘图1
plot_path="hvg"+sys.argv[1]+"_cluster.png"
sc.pl.umap( adata, color=["leiden_res0_2", "leiden_res0_6", "leiden_res1","leiden_res1_2","leiden_res2"], legend_loc="on data",save=plot_path)

output_path="RData/hvg"+ sys.argv[1] +"_cluster.h5ad"
adata.write(output_path)
