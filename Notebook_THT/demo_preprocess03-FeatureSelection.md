## 修改conda路径


```python
import os
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


## 导入包


```python
import logging

import matplotlib.pyplot as plt
import numpy as np
np.float_ = np.float64
import anndata2ri
import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro
import scanpy as sc
import seaborn as sns

sc.settings.verbosity = 0
sc.settings.set_figure_params(
    dpi=80,
    facecolor="white",
    frameon=False,
)

rcb.logger.setLevel(logging.ERROR)
ro.pandas2ri.activate()
anndata2ri.activate()

%load_ext rpy2.ipython
```

## 挑选高变基因


```python
adata = sc.read("RData/merge_demo.normalization.h5ad")
```


```python
sc.pp.highly_variable_genes(adata, layer="log1p_norm",n_top_genes=1000,batch_key="sample")
```


```python
adata.write("RData/mergedemo_feature_selection.h5ad")
```
