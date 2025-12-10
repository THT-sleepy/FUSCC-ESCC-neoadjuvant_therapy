##导入模块，设置参数
import sys
import scanpy as sc
import re

sc.settings.verbosity = 0
sc.settings.set_figure_params(dpi=80, facecolor="white", frameon=False)
sc.settings.figdir = "plots"
sc.settings.autoshow = False

##读取1
input_path="RData/CD8+T__umapmd0.4.h5ad"
adata = sc.read(input_path)



##聚类1
sc.tl.leiden(adata, key_added="leiden_res0_8", resolution=0.8)
sc.tl.leiden(adata, key_added="leiden_res1", resolution=1)
sc.tl.leiden(adata, key_added="leiden_res1_2", resolution=1.2)
sc.tl.leiden(adata, key_added="leiden_res1_4", resolution=1.4)
sc.tl.leiden(adata, key_added="leiden_res1_6", resolution=1.6)
sc.tl.leiden(adata, key_added="leiden_res1_8", resolution=1.8)


##绘图1
plot_path="CD8_cluster_res1_2to1_8.png"
sc.pl.umap( adata, color=["leiden_res1_2","leiden_res1_4","leiden_res1_6","leiden_res1_8"], legend_loc="on data",save=plot_path)

output_path="RData/CD8+T__umapmd0.4.h5ad"
adata.write(output_path)
