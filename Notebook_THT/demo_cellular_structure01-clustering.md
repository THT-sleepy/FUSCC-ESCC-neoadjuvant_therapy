## 修改路径为conda路径


```python
import os
```


```python
print(os.environ['PATH'])
```

    /usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/snap/bin



```python
conda_path = "/home/data/t190513/miniconda3/envs/sc_clustering/bin:/home/data/t190513/miniconda3/condabin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/snap/bin:/home/data/t190513/miniconda3/bin"
```


```python
os.environ['PATH'] = conda_path + os.environ['PATH']
```


```python
print(os.environ['PATH'])
```

    /home/data/t190513/miniconda3/envs/sc_clustering/bin:/home/data/t190513/miniconda3/condabin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/snap/bin:/home/data/t190513/miniconda3/bin/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/snap/bin


## 导入包


```python
import scanpy as sc

sc.settings.verbosity = 0
sc.settings.set_figure_params(dpi=80, facecolor="white", frameon=False)
```


```python
adata = sc.read("RData/mergedemo_BatchHarmony.h5ad")
```


```python
#sc.pp.neighbors(adata)
#sc.tl.umap(adata)
```

    /home/data/t190513/miniconda3/envs/sc_clustering/lib/python3.10/site-packages/tqdm/auto.py:21: TqdmWarning: IProgress not found. Please update jupyter and ipywidgets. See https://ipywidgets.readthedocs.io/en/stable/user_install.html
      from .autonotebook import tqdm as notebook_tqdm



```python
sc.tl.leiden(adata, key_added="leiden_res0_1", resolution=0.1)
sc.tl.leiden(adata, key_added="leiden_res0_25", resolution=0.25)
sc.tl.leiden(adata, key_added="leiden_res0_5", resolution=0.5)
sc.tl.leiden(adata, key_added="leiden_res1", resolution=1.0)
```

    /tmp/ipykernel_1162421/242196680.py:1: FutureWarning: In the future, the default backend for leiden will be igraph instead of leidenalg.
    
     To achieve the future defaults please pass: flavor="igraph" and n_iterations=2.  directed must also be False to work with igraph's implementation.
      sc.tl.leiden(adata, key_added="leiden_res0_1", resolution=0.1)



```python
sc.pl.umap(
    adata,
    color=["leiden_res0_1","leiden_res0_25", "leiden_res0_5", "leiden_res1"],
    legend_loc="on data",
)
```


    
![png](output_11_0.png)
    



```python
sc.pl.umap(
    adata,
    color=['scDblFinder_class'],
    legend_loc="on data",
)
```


    
![png](output_12_0.png)
    



```python
adata.write("RData/mergedemo_clustered.h5ad")
```


```python
adata
```




    AnnData object with n_obs × n_vars = 28911 × 15674
        obs: 'n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'pct_counts_in_top_20_genes', 'total_counts_mt', 'log1p_total_counts_mt', 'pct_counts_mt', 'total_counts_ribo', 'log1p_total_counts_ribo', 'pct_counts_ribo', 'total_counts_hb', 'log1p_total_counts_hb', 'pct_counts_hb', 'outlier', 'mt_outlier', 'scDblFinder_score', 'scDblFinder_class', 'sample', 'leiden_res0_1', 'leiden_res0_25', 'leiden_res0_5', 'leiden_res1'
        var: 'gene_ids', 'feature_types', 'mt', 'ribo', 'hb', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts', 'n_cells', 'highly_variable', 'means', 'dispersions', 'dispersions_norm', 'highly_variable_nbatches', 'highly_variable_intersection'
        uns: 'hvg', 'neighbors', 'pca', 'sample_colors', 'scDblFinder_class_colors', 'tsne', 'umap', 'leiden_res0_1', 'leiden_res0_25', 'leiden_res0_5', 'leiden_res1', 'leiden_res0_1_colors', 'leiden_res0_25_colors', 'leiden_res0_5_colors', 'leiden_res1_colors'
        obsm: 'X_pca', 'X_pca_harmony', 'X_tsne', 'X_umap'
        varm: 'PCs'
        layers: 'counts', 'log1p_norm', 'soupX_counts'
        obsp: 'connectivities', 'distances'




```python
?sc.pp.neighbors
```


    [0;31mSignature:[0m
    [0msc[0m[0;34m.[0m[0mpp[0m[0;34m.[0m[0mneighbors[0m[0;34m([0m[0;34m[0m
    [0;34m[0m    [0madata[0m[0;34m:[0m [0;34m'AnnData'[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0mn_neighbors[0m[0;34m:[0m [0;34m'int'[0m [0;34m=[0m [0;36m15[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0mn_pcs[0m[0;34m:[0m [0;34m'int | None'[0m [0;34m=[0m [0;32mNone[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0;34m*[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0muse_rep[0m[0;34m:[0m [0;34m'str | None'[0m [0;34m=[0m [0;32mNone[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0mknn[0m[0;34m:[0m [0;34m'bool'[0m [0;34m=[0m [0;32mTrue[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0mmethod[0m[0;34m:[0m [0;34m'_Method'[0m [0;34m=[0m [0;34m'umap'[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0mtransformer[0m[0;34m:[0m [0;34m'KnnTransformerLike | _KnownTransformer | None'[0m [0;34m=[0m [0;32mNone[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0mmetric[0m[0;34m:[0m [0;34m'_Metric | _MetricFn'[0m [0;34m=[0m [0;34m'euclidean'[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0mmetric_kwds[0m[0;34m:[0m [0;34m'Mapping[str, Any]'[0m [0;34m=[0m [0mmappingproxy[0m[0;34m([0m[0;34m{[0m[0;34m}[0m[0;34m)[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0mrandom_state[0m[0;34m:[0m [0;34m'_LegacyRandom'[0m [0;34m=[0m [0;36m0[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0mkey_added[0m[0;34m:[0m [0;34m'str | None'[0m [0;34m=[0m [0;32mNone[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0mcopy[0m[0;34m:[0m [0;34m'bool'[0m [0;34m=[0m [0;32mFalse[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m[0;34m)[0m [0;34m->[0m [0;34m'AnnData | None'[0m[0;34m[0m[0;34m[0m[0m
    [0;31mDocstring:[0m
    Compute the nearest neighbors distance matrix and a neighborhood graph of observations :cite:p:`McInnes2018`.
    
    The neighbor search efficiency of this heavily relies on UMAP :cite:p:`McInnes2018`,
    which also provides a method for estimating connectivities of data points -
    the connectivity of the manifold (`method=='umap'`). If `method=='gauss'`,
    connectivities are computed according to :cite:t:`Coifman2005`, in the adaption of
    :cite:t:`Haghverdi2016`.
    
    Parameters
    ----------
    adata : 'AnnData'
        Annotated data matrix.
    n_neighbors : 'int', optional (default: 15)
        The size of local neighborhood (in terms of number of neighboring data
        points) used for manifold approximation. Larger values result in more
        global views of the manifold, while smaller values result in more local
        data being preserved. In general values should be in the range 2 to 100.
        If `knn` is `True`, number of nearest neighbors to be searched. If `knn`
        is `False`, a Gaussian kernel width is set to the distance of the
        `n_neighbors` neighbor.
    
        *ignored if ``transformer`` is an instance.*
    n_pcs : 'int | None', optional (default: None)
        Use this many PCs. If `n_pcs==0` use `.X` if `use_rep is None`.
    use_rep : 'str | None', optional (default: None)
        Use the indicated representation. `'X'` or any key for `.obsm` is valid.
        If `None`, the representation is chosen automatically:
        For `.n_vars` < :attr:`~scanpy._settings.ScanpyConfig.N_PCS` (default: 50), `.X` is used, otherwise 'X_pca' is used.
        If 'X_pca' is not present, it’s computed with default parameters or `n_pcs` if present.
    knn : 'bool', optional (default: True)
        If `True`, use a hard threshold to restrict the number of neighbors to
        `n_neighbors`, that is, consider a knn graph. Otherwise, use a Gaussian
        Kernel to assign low weights to neighbors more distant than the
        `n_neighbors` nearest neighbor.
    method : '_Method', optional (default: 'umap')
        Use 'umap' :cite:p:`McInnes2018` or 'gauss' (Gauss kernel following :cite:t:`Coifman2005`
        with adaptive width :cite:t:`Haghverdi2016`) for computing connectivities.
    transformer : 'KnnTransformerLike | _KnownTransformer | None', optional (default: None)
        Approximate kNN search implementation following the API of
        :class:`~sklearn.neighbors.KNeighborsTransformer`.
        See :doc:`/how-to/knn-transformers` for more details.
        Also accepts the following known options:
    
        `None` (the default)
            Behavior depends on data size.
            For small data, we will calculate exact kNN, otherwise we use
            :class:`~pynndescent.pynndescent_.PyNNDescentTransformer`
        `'pynndescent'`
            :class:`~pynndescent.pynndescent_.PyNNDescentTransformer`
        `'rapids'`
            A transformer based on :class:`cuml.neighbors.NearestNeighbors`.
    
            .. deprecated:: 1.10.0
               Use :func:`rapids_singlecell.pp.neighbors` instead.
    metric : '_Metric | _MetricFn', optional (default: 'euclidean')
        A known metric’s name or a callable that returns a distance.
    
        *ignored if ``transformer`` is an instance.*
    metric_kwds : 'Mapping[str, Any]', optional (default: mappingproxy({}))
        Options for the metric.
    
        *ignored if ``transformer`` is an instance.*
    random_state : '_LegacyRandom', optional (default: 0)
        A numpy random seed.
    
        *ignored if ``transformer`` is an instance.*
    key_added : 'str | None', optional (default: None)
        If not specified, the neighbors data is stored in `.uns['neighbors']`,
        distances and connectivities are stored in `.obsp['distances']` and
        `.obsp['connectivities']` respectively.
        If specified, the neighbors data is added to .uns[key_added],
        distances are stored in `.obsp[key_added+'_distances']` and
        connectivities in `.obsp[key_added+'_connectivities']`.
    copy : 'bool', optional (default: False)
        Return a copy instead of writing to adata.
    
    Returns
    -------
    Returns `None` if `copy=False`, else returns an `AnnData` object. Sets the following fields:
    
    `adata.obsp['distances' | key_added+'_distances']` : :class:`scipy.sparse.csr_matrix` (dtype `float`)
        Distance matrix of the nearest neighbors search. Each row (cell) has `n_neighbors`-1 non-zero entries. These are the distances to their `n_neighbors`-1 nearest neighbors (excluding the cell itself).
    `adata.obsp['connectivities' | key_added+'_connectivities']` : :class:`scipy.sparse._csr.csr_matrix` (dtype `float`)
        Weighted adjacency matrix of the neighborhood graph of data
        points. Weights should be interpreted as connectivities.
    `adata.uns['neighbors' | key_added]` : :class:`dict`
        neighbors parameters.
    
    Examples
    --------
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> # Basic usage
    >>> sc.pp.neighbors(adata, 20, metric="cosine")
    >>> # Provide your own transformer for more control and flexibility
    >>> from sklearn.neighbors import KNeighborsTransformer
    >>> transformer = KNeighborsTransformer(
    ...     n_neighbors=10, metric="manhattan", algorithm="kd_tree"
    ... )
    >>> sc.pp.neighbors(adata, transformer=transformer)
    >>> # now you can e.g. access the index: `transformer._tree`
    
    See Also
    --------
    :doc:`/how-to/knn-transformers`
    [0;31mFile:[0m      ~/miniconda3/envs/sc_clustering/lib/python3.10/site-packages/scanpy/neighbors/__init__.py
    [0;31mType:[0m      function

