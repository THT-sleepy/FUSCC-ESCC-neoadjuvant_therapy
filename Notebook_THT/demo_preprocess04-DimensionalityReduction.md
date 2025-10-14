## 修改路径为conda路径


```python
import os
```


```python
print(os.environ['PATH'])
```

    /usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/snap/bin



```python
conda_path = "/home/data/t190513/miniconda3/envs/sc_preprocess/bin:/home/data/t190513/miniconda3/condabin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/snap/bin:/home/data/t190513/miniconda3/bin"
```


```python
os.environ['PATH'] = conda_path + os.environ['PATH']
```


```python
print(os.environ['PATH'])
```

    /home/data/t190513/miniconda3/envs/sc_preprocess/bin:/home/data/t190513/miniconda3/condabin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/snap/bin:/home/data/t190513/miniconda3/bin/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/snap/bin


## 导入包及设置参数


```python
import scanpy as sc

sc.settings.verbosity = 0
sc.settings.set_figure_params(
    dpi=80,
    facecolor="white",
    frameon=False
)
sc.settings.figdir = "plots"
sc.settings.autoshow = False 
```


```python
?sc.pp.neighbors
```


    [31mSignature:[39m
    sc.pp.neighbors(
        adata: [33m'AnnData'[39m,
        n_neighbors: [33m'int'[39m = [32m15[39m,
        n_pcs: [33m'int | None'[39m = [38;5;28;01mNone[39;00m,
        *,
        use_rep: [33m'str | None'[39m = [38;5;28;01mNone[39;00m,
        knn: [33m'bool'[39m = [38;5;28;01mTrue[39;00m,
        method: [33m'_Method'[39m = [33m'umap'[39m,
        transformer: [33m'KnnTransformerLike | _KnownTransformer | None'[39m = [38;5;28;01mNone[39;00m,
        metric: [33m'_Metric | _MetricFn'[39m = [33m'euclidean'[39m,
        metric_kwds: [33m'Mapping[str, Any]'[39m = mappingproxy({}),
        random_state: [33m'_LegacyRandom'[39m = [32m0[39m,
        key_added: [33m'str | None'[39m = [38;5;28;01mNone[39;00m,
        copy: [33m'bool'[39m = [38;5;28;01mFalse[39;00m,
    ) -> [33m'AnnData | None'[39m
    [31mDocstring:[39m
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
    [31mFile:[39m      ~/miniconda3/envs/sc_preprocess/lib/python3.12/site-packages/scanpy/neighbors/__init__.py
    [31mType:[39m      function



```python
adata = sc.read(
    filename="RData/mergedemo_feature_selection.h5ad"
)
```


```python
adata.X = adata.layers["log1p_norm"]
```


```python
# setting highly variable as highly deviant to use scanpy 'use_highly_variable' argument in sc.pp.pca
#adata.var["highly_variable"] = adata.var["highly_deviant"]
sc.pp.pca(adata, svd_solver="arpack", use_highly_variable=True) 
```

    /home/data/t190513/miniconda3/envs/sc_preprocess/lib/python3.12/site-packages/scanpy/preprocessing/_pca/__init__.py:438: FutureWarning: Argument `use_highly_variable` is deprecated, consider using the mask argument. Use_highly_variable=True can be called through mask_var="highly_variable". Use_highly_variable=False can be called through mask_var=None
      warn(msg, FutureWarning)



```python
sc.pl.pca_scatter(adata, color="total_counts",save="mergedemo_pca.pdf")
```




    <Axes: title={'center': 'total_counts'}, xlabel='PC1', ylabel='PC2'>




```python
sc.tl.tsne(adata, use_rep="X_pca") 
```


```python
sc.pl.tsne(adata, color="total_counts")
```




    <Axes: title={'center': 'total_counts'}, xlabel='tSNE1', ylabel='tSNE2'>




    
![png](output_14_1.png)
    



```python
sc.pp.neighbors(adata) ##默认使用X_pca,PC默认是50，neighbors默认是15，按需调整
sc.tl.umap(adata) 
```


```python
sc.pl.umap(adata, color="total_counts",save="demo_umap.pdf")
```




    <Axes: title={'center': 'total_counts'}, xlabel='UMAP1', ylabel='UMAP2'>




```python
sc.pl.umap(
    adata,
    color=["total_counts", "pct_counts_mt", "scDblFinder_score", "scDblFinder_class"],
)
```




    [<Axes: title={'center': 'total_counts'}, xlabel='UMAP1', ylabel='UMAP2'>,
     <Axes: title={'center': 'pct_counts_mt'}, xlabel='UMAP1', ylabel='UMAP2'>,
     <Axes: title={'center': 'scDblFinder_score'}, xlabel='UMAP1', ylabel='UMAP2'>,
     <Axes: title={'center': 'scDblFinder_class'}, xlabel='UMAP1', ylabel='UMAP2'>]




    
![png](output_17_1.png)
    



```python
adata.write("RData/mergedemo_dimensionality_reduction.h5ad")
```

### 这里代码有一个可以优化的地方是不切换adata.X而是直接在函数中指明layer
