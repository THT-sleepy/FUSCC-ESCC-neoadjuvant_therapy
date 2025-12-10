import scanpy as sc
import pandas as pd

sc.settings.verbosity = 0
sc.settings.set_figure_params(
    dpi=300,
    facecolor="white",
    frameon=False
)
sc.settings.figdir = "plots"
sc.settings.autoshow = False

adata = sc.read("RData/Myeloid_1stround_annotated.h5ad") #X是log1p_norm

def process_celltype(adata, celltype_name, save_dir="RData/", plot_prefix=""):
    """
    批量处理特定细胞类型的单细胞数据（归一化→PCA→BBKNN→UMAP→Leiden聚类）

    参数：
    adata: 原始AnnData对象（需包含first_round_celltype列和soupX_counts层）
    celltype_name: 目标细胞类型（如"T"、"B"、"Epi"）
    save_dir: h5ad文件保存目录（默认RData/）
    plot_prefix: 图片文件名前缀（默认空，可加项目标识）
    """
    # 1. 提取目标细胞类型子集
    adata_sub = adata[adata.obs.myeloid_celltype == celltype_name].copy()
    print(f"开始处理 {celltype_name} 细胞，共 {adata_sub.n_obs} 个细胞")

    # 1.5 去掉细胞数小于10个的样本
    sample_counts = adata_sub.obs['sample'].value_counts()
    keep_samples = sample_counts[sample_counts >= 10].index
    adata_sub = adata_sub[adata_sub.obs['sample'].isin(keep_samples)].copy()

    # 2. 重新normalize和feature selection
    adata_sub.X = adata_sub.layers["soupX_counts"]  # 用soupX校正后的counts
    scales_counts = sc.pp.normalize_total(adata_sub, target_sum=None, inplace=False)
    adata_sub.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)  # log1p转换
    sc.pp.highly_variable_genes(  # 筛选高变基因
        adata_sub, layer="log1p_norm",batch_key="sample",n_top_genes=2000
    )
    adata_sub.X = adata_sub.layers["log1p_norm"]  # 将log1p后的数据设为默认X
    adata_sub.layers["regress"] = adata_sub.layers["log1p_norm"]
    adata_sub.X = adata_sub.layers["regress"]  # 将log1p后的数据设为默认X


    #2.5 regressout mt ribo
    sc.pp.regress_out(adata_sub, ['pct_counts_mt','pct_counts_ribo'])

    # 3. PCA降维
    sc.pp.pca(adata_sub, svd_solver="arpack", use_highly_variable=True)

    # 4. 循环运行不同min_dist的UMAP+Leiden聚类（0.4）
    for min_dist in [0.4]:
        # 复制数据避免相互干扰
        adata_md = adata_sub.copy()

        # BBKNN批次校正 + UMAP降维
        sc.external.pp.bbknn(adata_md, batch_key="sample", n_pcs=30)
        sc.tl.umap(adata_md, min_dist=min_dist)

        # 多分辨率Leiden聚类
        sc.tl.leiden(adata_md, key_added="leiden_res0_6", resolution=0.6)
        sc.tl.leiden(adata_md, key_added="leiden_res0_8", resolution=0.8)
        sc.tl.leiden(adata_md, key_added="leiden_res1", resolution=1.0)

        # 保存UMAP聚类图
        plot_path = f"{plot_prefix}{celltype_name}_umapmd{min_dist}_cluster.png"
        sc.pl.umap(
            adata_md,
            color=["leiden_res0_6", "leiden_res0_8", "leiden_res1"],
            legend_loc="on data",
            save=plot_path,
            show=False  # 避免弹出图片窗口，仅保存
        )

        # 保存AnnData对象
        output_fpath = f"{save_dir}{celltype_name}__umapmd{min_dist}.h5ad"
        adata_md.X = adata_md.layers["log1p_norm"]
        adata_md.write(output_fpath)
        print(f"  {celltype_name} (min_dist={min_dist}) 已保存：{output_fpath}")

    print(f"{celltype_name} 细胞处理完成！\n")
# 假设你的原始AnnData对象已加载为 'adata'
# 1. 定义要处理的细胞类型列表（从T开始，按需求调整顺序）
celltypes_to_process = ["Mono_Macro","DC","Mast"]

# 2. 批量运行（可自定义保存目录和图片前缀）
for celltype in celltypes_to_process:
    process_celltype(
        adata=adata,
        celltype_name=celltype,
        save_dir="RData/",  # h5ad文件保存到RData目录
        plot_prefix=""  # 图片名前缀（如scRNA_T_umapmd0.1_cluster.png）
    )
