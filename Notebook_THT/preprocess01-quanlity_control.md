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


## 导入库


```python
import numpy as np
np.float_ = np.float64
import scanpy as sc
import seaborn as sns
import pandas as pd
from scipy.stats import median_abs_deviation

sc.settings.verbosity = 0
sc.settings.set_figure_params(
    dpi=80,
    facecolor="white",
    frameon=False,
)
```


```python
data_root = "input/final_cellranger_data"
```


```python
samples = []
```


```python
for sample_id in os.listdir(data_root):
    # 样本文件夹的完整路径
    sample_dir = os.path.join(data_root, sample_id)
    
    # 构建filtered特征矩阵文件夹路径
    filtered_matrix_dir = os.path.join(sample_dir, "sample_filtered_feature_bc_matrix")
    # 构建raw 特征矩阵文件路径
    raw_matrix_path = os.path.join(sample_dir, "raw_feature_bc_matrix.h5")
    # 构建cluster信息文件路径
    cluster_path = os.path.join(sample_dir,"clusters.csv")
    # 指定质控后h5ad文件输出路径
    quanlity_control_h5_path = os.path.join(sample_dir,"quanlity_control.h5ad")
    # 指定质控信息输出路径
    quanlity_control_info_path = os.path.join(sample_dir,"quanlity_control_info.txt")
    # 构建样本信息字典
    sample_info = {
                "sample_id": sample_id,
                "filtered_matrix_dir_path": filtered_matrix_dir,
                "raw_matrix_h5_path": raw_matrix_path,
                "cluster_path": cluster_path,
                "quanlity_control_h5_path" : quanlity_control_h5_path,
                "quanlity_control_info_path": quanlity_control_info_path,
            }
    # 添加到样本列表
    samples.append(sample_info)
```


```python
sample = samples[1]
```


```python
sample
```




    {'sample_id': 'YXY_preC2_2282374',
     'filtered_matrix_dir_path': 'input/final_cellranger_data/YXY_preC2_2282374/sample_filtered_feature_bc_matrix',
     'raw_matrix_h5_path': 'input/final_cellranger_data/YXY_preC2_2282374/raw_feature_bc_matrix.h5',
     'cluster_path': 'input/final_cellranger_data/YXY_preC2_2282374/clusters.csv',
     'quanlity_control_h5_path': 'input/final_cellranger_data/YXY_preC2_2282374/quanlity_control.h5ad',
     'quanlity_control_info_path': 'input/final_cellranger_data/YXY_preC2_2282374/quanlity_control_info.txt'}




```python
samples
```




    [{'sample_id': 'XZM_preC1_2285332',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/XZM_preC1_2285332/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/XZM_preC1_2285332/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/XZM_preC1_2285332/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/XZM_preC1_2285332/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/XZM_preC1_2285332/quanlity_control_info.txt'},
     {'sample_id': 'YXY_preC2_2282374',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/YXY_preC2_2282374/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/YXY_preC2_2282374/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/YXY_preC2_2282374/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/YXY_preC2_2282374/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/YXY_preC2_2282374/quanlity_control_info.txt'},
     {'sample_id': 'GZP_preC1_2271757',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/GZP_preC1_2271757/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/GZP_preC1_2271757/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/GZP_preC1_2271757/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/GZP_preC1_2271757/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/GZP_preC1_2271757/quanlity_control_info.txt'},
     {'sample_id': 'THX_op_2293786',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/THX_op_2293786/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/THX_op_2293786/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/THX_op_2293786/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/THX_op_2293786/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/THX_op_2293786/quanlity_control_info.txt'},
     {'sample_id': 'XSY_pbmc_op_2308306',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/XSY_pbmc_op_2308306/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/XSY_pbmc_op_2308306/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/XSY_pbmc_op_2308306/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/XSY_pbmc_op_2308306/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/XSY_pbmc_op_2308306/quanlity_control_info.txt'},
     {'sample_id': 'XDQ_preC1_2286869',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/XDQ_preC1_2286869/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/XDQ_preC1_2286869/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/XDQ_preC1_2286869/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/XDQ_preC1_2286869/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/XDQ_preC1_2286869/quanlity_control_info.txt'},
     {'sample_id': 'WDY_op_2302792',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/WDY_op_2302792/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/WDY_op_2302792/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/WDY_op_2302792/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/WDY_op_2302792/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/WDY_op_2302792/quanlity_control_info.txt'},
     {'sample_id': 'DJZ_pbmc_op_2301100',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/DJZ_pbmc_op_2301100/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/DJZ_pbmc_op_2301100/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/DJZ_pbmc_op_2301100/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/DJZ_pbmc_op_2301100/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/DJZ_pbmc_op_2301100/quanlity_control_info.txt'},
     {'sample_id': 'SML_op_2290669',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/SML_op_2290669/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/SML_op_2290669/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/SML_op_2290669/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/SML_op_2290669/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/SML_op_2290669/quanlity_control_info.txt'},
     {'sample_id': 'WJJ_preC1_2285056',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/WJJ_preC1_2285056/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/WJJ_preC1_2285056/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/WJJ_preC1_2285056/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/WJJ_preC1_2285056/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/WJJ_preC1_2285056/quanlity_control_info.txt'},
     {'sample_id': 'CDZ_preC2_2282421',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/CDZ_preC2_2282421/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/CDZ_preC2_2282421/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/CDZ_preC2_2282421/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/CDZ_preC2_2282421/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/CDZ_preC2_2282421/quanlity_control_info.txt'},
     {'sample_id': 'OBC_pbmc_op_2292682',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/OBC_pbmc_op_2292682/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/OBC_pbmc_op_2292682/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/OBC_pbmc_op_2292682/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/OBC_pbmc_op_2292682/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/OBC_pbmc_op_2292682/quanlity_control_info.txt'},
     {'sample_id': 'LEG_pbmc_op_2282331',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/LEG_pbmc_op_2282331/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/LEG_pbmc_op_2282331/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/LEG_pbmc_op_2282331/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/LEG_pbmc_op_2282331/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/LEG_pbmc_op_2282331/quanlity_control_info.txt'},
     {'sample_id': 'GZP_pbmc_op_2271757',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/GZP_pbmc_op_2271757/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/GZP_pbmc_op_2271757/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/GZP_pbmc_op_2271757/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/GZP_pbmc_op_2271757/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/GZP_pbmc_op_2271757/quanlity_control_info.txt'},
     {'sample_id': 'MQW_preC1_2268090',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/MQW_preC1_2268090/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/MQW_preC1_2268090/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/MQW_preC1_2268090/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/MQW_preC1_2268090/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/MQW_preC1_2268090/quanlity_control_info.txt'},
     {'sample_id': 'LEG_preC2_2282331',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/LEG_preC2_2282331/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/LEG_preC2_2282331/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/LEG_preC2_2282331/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/LEG_preC2_2282331/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/LEG_preC2_2282331/quanlity_control_info.txt'},
     {'sample_id': 'LCR_preC1_2265065',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/LCR_preC1_2265065/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/LCR_preC1_2265065/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/LCR_preC1_2265065/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/LCR_preC1_2265065/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/LCR_preC1_2265065/quanlity_control_info.txt'},
     {'sample_id': 'THX_preC1_2293786',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/THX_preC1_2293786/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/THX_preC1_2293786/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/THX_preC1_2293786/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/THX_preC1_2293786/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/THX_preC1_2293786/quanlity_control_info.txt'},
     {'sample_id': 'YLY_pbmc_op_2282682',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/YLY_pbmc_op_2282682/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/YLY_pbmc_op_2282682/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/YLY_pbmc_op_2282682/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/YLY_pbmc_op_2282682/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/YLY_pbmc_op_2282682/quanlity_control_info.txt'},
     {'sample_id': 'YLB_pbmc_op_2277686',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/YLB_pbmc_op_2277686/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/YLB_pbmc_op_2277686/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/YLB_pbmc_op_2277686/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/YLB_pbmc_op_2277686/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/YLB_pbmc_op_2277686/quanlity_control_info.txt'},
     {'sample_id': 'XDQ_preC2_2286869',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/XDQ_preC2_2286869/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/XDQ_preC2_2286869/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/XDQ_preC2_2286869/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/XDQ_preC2_2286869/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/XDQ_preC2_2286869/quanlity_control_info.txt'},
     {'sample_id': 'WKD_op_2274363',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/WKD_op_2274363/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/WKD_op_2274363/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/WKD_op_2274363/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/WKD_op_2274363/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/WKD_op_2274363/quanlity_control_info.txt'},
     {'sample_id': 'ZMZ_preC2_2267536',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/ZMZ_preC2_2267536/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/ZMZ_preC2_2267536/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/ZMZ_preC2_2267536/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/ZMZ_preC2_2267536/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/ZMZ_preC2_2267536/quanlity_control_info.txt'},
     {'sample_id': 'WJJ_preC2_2285056',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/WJJ_preC2_2285056/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/WJJ_preC2_2285056/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/WJJ_preC2_2285056/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/WJJ_preC2_2285056/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/WJJ_preC2_2285056/quanlity_control_info.txt'},
     {'sample_id': 'XDQ_op_2286869',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/XDQ_op_2286869/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/XDQ_op_2286869/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/XDQ_op_2286869/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/XDQ_op_2286869/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/XDQ_op_2286869/quanlity_control_info.txt'},
     {'sample_id': 'RWX_preC1_2265611',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/RWX_preC1_2265611/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/RWX_preC1_2265611/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/RWX_preC1_2265611/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/RWX_preC1_2265611/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/RWX_preC1_2265611/quanlity_control_info.txt'},
     {'sample_id': 'ZHM_op_2301356',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/ZHM_op_2301356/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/ZHM_op_2301356/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/ZHM_op_2301356/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/ZHM_op_2301356/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/ZHM_op_2301356/quanlity_control_info.txt'},
     {'sample_id': 'SML_pbmc_op_2290669',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/SML_pbmc_op_2290669/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/SML_pbmc_op_2290669/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/SML_pbmc_op_2290669/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/SML_pbmc_op_2290669/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/SML_pbmc_op_2290669/quanlity_control_info.txt'},
     {'sample_id': 'GZP_op_2271757',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/GZP_op_2271757/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/GZP_op_2271757/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/GZP_op_2271757/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/GZP_op_2271757/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/GZP_op_2271757/quanlity_control_info.txt'},
     {'sample_id': 'GHC_pbmc_op_2297515',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/GHC_pbmc_op_2297515/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/GHC_pbmc_op_2297515/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/GHC_pbmc_op_2297515/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/GHC_pbmc_op_2297515/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/GHC_pbmc_op_2297515/quanlity_control_info.txt'},
     {'sample_id': 'PDD_pbmc_op_2284071',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/PDD_pbmc_op_2284071/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/PDD_pbmc_op_2284071/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/PDD_pbmc_op_2284071/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/PDD_pbmc_op_2284071/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/PDD_pbmc_op_2284071/quanlity_control_info.txt'},
     {'sample_id': 'PDD_op_2284071',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/PDD_op_2284071/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/PDD_op_2284071/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/PDD_op_2284071/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/PDD_op_2284071/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/PDD_op_2284071/quanlity_control_info.txt'},
     {'sample_id': 'XJZ_preC1_2296878',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/XJZ_preC1_2296878/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/XJZ_preC1_2296878/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/XJZ_preC1_2296878/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/XJZ_preC1_2296878/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/XJZ_preC1_2296878/quanlity_control_info.txt'},
     {'sample_id': 'WDY_preC2_2302792',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/WDY_preC2_2302792/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/WDY_preC2_2302792/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/WDY_preC2_2302792/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/WDY_preC2_2302792/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/WDY_preC2_2302792/quanlity_control_info.txt'},
     {'sample_id': 'MHY_pbmc_op_2268689',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/MHY_pbmc_op_2268689/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/MHY_pbmc_op_2268689/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/MHY_pbmc_op_2268689/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/MHY_pbmc_op_2268689/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/MHY_pbmc_op_2268689/quanlity_control_info.txt'},
     {'sample_id': 'THX_preC2_2293786',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/THX_preC2_2293786/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/THX_preC2_2293786/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/THX_preC2_2293786/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/THX_preC2_2293786/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/THX_preC2_2293786/quanlity_control_info.txt'},
     {'sample_id': 'OBC_preC2_2292682',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/OBC_preC2_2292682/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/OBC_preC2_2292682/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/OBC_preC2_2292682/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/OBC_preC2_2292682/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/OBC_preC2_2292682/quanlity_control_info.txt'},
     {'sample_id': 'LYG_op_2290266',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/LYG_op_2290266/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/LYG_op_2290266/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/LYG_op_2290266/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/LYG_op_2290266/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/LYG_op_2290266/quanlity_control_info.txt'},
     {'sample_id': 'LGZ_op_2300309',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/LGZ_op_2300309/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/LGZ_op_2300309/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/LGZ_op_2300309/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/LGZ_op_2300309/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/LGZ_op_2300309/quanlity_control_info.txt'},
     {'sample_id': 'YXY_preC1_2282374',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/YXY_preC1_2282374/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/YXY_preC1_2282374/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/YXY_preC1_2282374/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/YXY_preC1_2282374/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/YXY_preC1_2282374/quanlity_control_info.txt'},
     {'sample_id': 'YLB_op_2277686',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/YLB_op_2277686/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/YLB_op_2277686/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/YLB_op_2277686/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/YLB_op_2277686/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/YLB_op_2277686/quanlity_control_info.txt'},
     {'sample_id': 'ZJM_op_2306458',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/ZJM_op_2306458/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/ZJM_op_2306458/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/ZJM_op_2306458/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/ZJM_op_2306458/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/ZJM_op_2306458/quanlity_control_info.txt'},
     {'sample_id': 'CDZ_op_2282421',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/CDZ_op_2282421/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/CDZ_op_2282421/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/CDZ_op_2282421/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/CDZ_op_2282421/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/CDZ_op_2282421/quanlity_control_info.txt'},
     {'sample_id': 'ZHM_preC1_2301356',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/ZHM_preC1_2301356/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/ZHM_preC1_2301356/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/ZHM_preC1_2301356/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/ZHM_preC1_2301356/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/ZHM_preC1_2301356/quanlity_control_info.txt'},
     {'sample_id': 'ZJM_preC2_2306458',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/ZJM_preC2_2306458/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/ZJM_preC2_2306458/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/ZJM_preC2_2306458/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/ZJM_preC2_2306458/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/ZJM_preC2_2306458/quanlity_control_info.txt'},
     {'sample_id': 'LDS_preC2_2304172',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/LDS_preC2_2304172/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/LDS_preC2_2304172/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/LDS_preC2_2304172/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/LDS_preC2_2304172/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/LDS_preC2_2304172/quanlity_control_info.txt'},
     {'sample_id': 'LGZ_pbmc_op_2300309',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/LGZ_pbmc_op_2300309/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/LGZ_pbmc_op_2300309/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/LGZ_pbmc_op_2300309/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/LGZ_pbmc_op_2300309/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/LGZ_pbmc_op_2300309/quanlity_control_info.txt'},
     {'sample_id': 'CJC_pbmc_op_2274041',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/CJC_pbmc_op_2274041/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/CJC_pbmc_op_2274041/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/CJC_pbmc_op_2274041/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/CJC_pbmc_op_2274041/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/CJC_pbmc_op_2274041/quanlity_control_info.txt'},
     {'sample_id': 'LDS_op_2304172',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/LDS_op_2304172/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/LDS_op_2304172/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/LDS_op_2304172/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/LDS_op_2304172/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/LDS_op_2304172/quanlity_control_info.txt'},
     {'sample_id': 'ZMZ_op_2267536',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/ZMZ_op_2267536/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/ZMZ_op_2267536/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/ZMZ_op_2267536/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/ZMZ_op_2267536/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/ZMZ_op_2267536/quanlity_control_info.txt'},
     {'sample_id': 'MHY_preC2_2268689',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/MHY_preC2_2268689/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/MHY_preC2_2268689/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/MHY_preC2_2268689/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/MHY_preC2_2268689/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/MHY_preC2_2268689/quanlity_control_info.txt'},
     {'sample_id': 'ZJM_preC1_2306458',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/ZJM_preC1_2306458/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/ZJM_preC1_2306458/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/ZJM_preC1_2306458/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/ZJM_preC1_2306458/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/ZJM_preC1_2306458/quanlity_control_info.txt'},
     {'sample_id': 'LEG_op_2282331',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/LEG_op_2282331/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/LEG_op_2282331/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/LEG_op_2282331/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/LEG_op_2282331/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/LEG_op_2282331/quanlity_control_info.txt'},
     {'sample_id': 'GHC_op_2297515',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/GHC_op_2297515/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/GHC_op_2297515/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/GHC_op_2297515/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/GHC_op_2297515/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/GHC_op_2297515/quanlity_control_info.txt'},
     {'sample_id': 'LYG_preC1_2290266',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/LYG_preC1_2290266/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/LYG_preC1_2290266/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/LYG_preC1_2290266/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/LYG_preC1_2290266/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/LYG_preC1_2290266/quanlity_control_info.txt'},
     {'sample_id': 'LDS_pbmc_op_2304172',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/LDS_pbmc_op_2304172/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/LDS_pbmc_op_2304172/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/LDS_pbmc_op_2304172/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/LDS_pbmc_op_2304172/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/LDS_pbmc_op_2304172/quanlity_control_info.txt'},
     {'sample_id': 'XZM_pbmc_op_2285332',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/XZM_pbmc_op_2285332/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/XZM_pbmc_op_2285332/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/XZM_pbmc_op_2285332/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/XZM_pbmc_op_2285332/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/XZM_pbmc_op_2285332/quanlity_control_info.txt'},
     {'sample_id': 'PDD_preC2_2284071',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/PDD_preC2_2284071/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/PDD_preC2_2284071/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/PDD_preC2_2284071/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/PDD_preC2_2284071/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/PDD_preC2_2284071/quanlity_control_info.txt'},
     {'sample_id': 'DJZ_preC2_2301100',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/DJZ_preC2_2301100/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/DJZ_preC2_2301100/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/DJZ_preC2_2301100/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/DJZ_preC2_2301100/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/DJZ_preC2_2301100/quanlity_control_info.txt'},
     {'sample_id': 'WDY_pbmc_op_2302792',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/WDY_pbmc_op_2302792/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/WDY_pbmc_op_2302792/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/WDY_pbmc_op_2302792/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/WDY_pbmc_op_2302792/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/WDY_pbmc_op_2302792/quanlity_control_info.txt'},
     {'sample_id': 'QBG_preC1_2300416',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/QBG_preC1_2300416/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/QBG_preC1_2300416/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/QBG_preC1_2300416/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/QBG_preC1_2300416/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/QBG_preC1_2300416/quanlity_control_info.txt'},
     {'sample_id': 'GZP_preC2_2271757',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/GZP_preC2_2271757/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/GZP_preC2_2271757/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/GZP_preC2_2271757/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/GZP_preC2_2271757/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/GZP_preC2_2271757/quanlity_control_info.txt'},
     {'sample_id': 'WDY_preC1_2302792',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/WDY_preC1_2302792/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/WDY_preC1_2302792/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/WDY_preC1_2302792/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/WDY_preC1_2302792/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/WDY_preC1_2302792/quanlity_control_info.txt'},
     {'sample_id': 'SML_preC2_2290669',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/SML_preC2_2290669/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/SML_preC2_2290669/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/SML_preC2_2290669/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/SML_preC2_2290669/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/SML_preC2_2290669/quanlity_control_info.txt'},
     {'sample_id': 'LYG_preC2_2290266',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/LYG_preC2_2290266/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/LYG_preC2_2290266/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/LYG_preC2_2290266/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/LYG_preC2_2290266/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/LYG_preC2_2290266/quanlity_control_info.txt'},
     {'sample_id': 'XZM_preC2_2285332',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/XZM_preC2_2285332/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/XZM_preC2_2285332/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/XZM_preC2_2285332/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/XZM_preC2_2285332/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/XZM_preC2_2285332/quanlity_control_info.txt'},
     {'sample_id': 'WYH_preC1_2296283',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/WYH_preC1_2296283/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/WYH_preC1_2296283/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/WYH_preC1_2296283/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/WYH_preC1_2296283/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/WYH_preC1_2296283/quanlity_control_info.txt'},
     {'sample_id': 'YXY_op_2282374',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/YXY_op_2282374/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/YXY_op_2282374/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/YXY_op_2282374/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/YXY_op_2282374/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/YXY_op_2282374/quanlity_control_info.txt'},
     {'sample_id': 'WYH_preC2_2296283',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/WYH_preC2_2296283/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/WYH_preC2_2296283/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/WYH_preC2_2296283/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/WYH_preC2_2296283/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/WYH_preC2_2296283/quanlity_control_info.txt'},
     {'sample_id': 'SML_preC1_2290669',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/SML_preC1_2290669/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/SML_preC1_2290669/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/SML_preC1_2290669/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/SML_preC1_2290669/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/SML_preC1_2290669/quanlity_control_info.txt'},
     {'sample_id': 'ZJM_pbmc_op_2306458',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/ZJM_pbmc_op_2306458/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/ZJM_pbmc_op_2306458/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/ZJM_pbmc_op_2306458/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/ZJM_pbmc_op_2306458/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/ZJM_pbmc_op_2306458/quanlity_control_info.txt'},
     {'sample_id': 'SZR_preC1_2279307',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/SZR_preC1_2279307/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/SZR_preC1_2279307/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/SZR_preC1_2279307/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/SZR_preC1_2279307/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/SZR_preC1_2279307/quanlity_control_info.txt'},
     {'sample_id': 'XSY_op_2308306',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/XSY_op_2308306/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/XSY_op_2308306/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/XSY_op_2308306/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/XSY_op_2308306/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/XSY_op_2308306/quanlity_control_info.txt'},
     {'sample_id': 'MHY_op_2268689',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/MHY_op_2268689/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/MHY_op_2268689/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/MHY_op_2268689/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/MHY_op_2268689/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/MHY_op_2268689/quanlity_control_info.txt'},
     {'sample_id': 'HLN_op_2295618',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/HLN_op_2295618/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/HLN_op_2295618/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/HLN_op_2295618/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/HLN_op_2295618/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/HLN_op_2295618/quanlity_control_info.txt'},
     {'sample_id': 'CJC_op_2274041',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/CJC_op_2274041/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/CJC_op_2274041/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/CJC_op_2274041/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/CJC_op_2274041/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/CJC_op_2274041/quanlity_control_info.txt'},
     {'sample_id': 'OBC_op_2292682',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/OBC_op_2292682/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/OBC_op_2292682/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/OBC_op_2292682/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/OBC_op_2292682/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/OBC_op_2292682/quanlity_control_info.txt'},
     {'sample_id': 'SCH_preC1_2265089',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/SCH_preC1_2265089/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/SCH_preC1_2265089/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/SCH_preC1_2265089/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/SCH_preC1_2265089/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/SCH_preC1_2265089/quanlity_control_info.txt'},
     {'sample_id': 'PLR_pbmc_op_2269964',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/PLR_pbmc_op_2269964/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/PLR_pbmc_op_2269964/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/PLR_pbmc_op_2269964/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/PLR_pbmc_op_2269964/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/PLR_pbmc_op_2269964/quanlity_control_info.txt'},
     {'sample_id': 'YLB_preC2_2277686',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/YLB_preC2_2277686/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/YLB_preC2_2277686/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/YLB_preC2_2277686/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/YLB_preC2_2277686/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/YLB_preC2_2277686/quanlity_control_info.txt'},
     {'sample_id': 'preC1_2291490',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/preC1_2291490/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/preC1_2291490/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/preC1_2291490/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/preC1_2291490/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/preC1_2291490/quanlity_control_info.txt'},
     {'sample_id': 'LCL_preC2_2284745',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/LCL_preC2_2284745/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/LCL_preC2_2284745/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/LCL_preC2_2284745/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/LCL_preC2_2284745/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/LCL_preC2_2284745/quanlity_control_info.txt'},
     {'sample_id': 'HJJ_preC1_2297052',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/HJJ_preC1_2297052/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/HJJ_preC1_2297052/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/HJJ_preC1_2297052/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/HJJ_preC1_2297052/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/HJJ_preC1_2297052/quanlity_control_info.txt'},
     {'sample_id': 'YXY_pbmc_op_2282374',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/YXY_pbmc_op_2282374/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/YXY_pbmc_op_2282374/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/YXY_pbmc_op_2282374/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/YXY_pbmc_op_2282374/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/YXY_pbmc_op_2282374/quanlity_control_info.txt'},
     {'sample_id': 'LGZ_preC1_2300309',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/LGZ_preC1_2300309/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/LGZ_preC1_2300309/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/LGZ_preC1_2300309/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/LGZ_preC1_2300309/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/LGZ_preC1_2300309/quanlity_control_info.txt'},
     {'sample_id': 'XSY_preC1_2308306',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/XSY_preC1_2308306/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/XSY_preC1_2308306/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/XSY_preC1_2308306/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/XSY_preC1_2308306/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/XSY_preC1_2308306/quanlity_control_info.txt'},
     {'sample_id': 'WXY_preC2_2295902',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/WXY_preC2_2295902/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/WXY_preC2_2295902/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/WXY_preC2_2295902/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/WXY_preC2_2295902/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/WXY_preC2_2295902/quanlity_control_info.txt'},
     {'sample_id': 'DJZ_preC1_2301100',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/DJZ_preC1_2301100/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/DJZ_preC1_2301100/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/DJZ_preC1_2301100/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/DJZ_preC1_2301100/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/DJZ_preC1_2301100/quanlity_control_info.txt'},
     {'sample_id': 'PLR_preC1_2269964',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/PLR_preC1_2269964/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/PLR_preC1_2269964/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/PLR_preC1_2269964/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/PLR_preC1_2269964/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/PLR_preC1_2269964/quanlity_control_info.txt'},
     {'sample_id': 'XZM_op_2285332',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/XZM_op_2285332/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/XZM_op_2285332/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/XZM_op_2285332/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/XZM_op_2285332/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/XZM_op_2285332/quanlity_control_info.txt'},
     {'sample_id': 'ZMZ_pbmc_op_2267536',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/ZMZ_pbmc_op_2267536/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/ZMZ_pbmc_op_2267536/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/ZMZ_pbmc_op_2267536/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/ZMZ_pbmc_op_2267536/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/ZMZ_pbmc_op_2267536/quanlity_control_info.txt'},
     {'sample_id': 'PDD_preC1_2284071',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/PDD_preC1_2284071/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/PDD_preC1_2284071/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/PDD_preC1_2284071/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/PDD_preC1_2284071/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/PDD_preC1_2284071/quanlity_control_info.txt'},
     {'sample_id': 'WXY_op_2295902',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/WXY_op_2295902/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/WXY_op_2295902/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/WXY_op_2295902/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/WXY_op_2295902/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/WXY_op_2295902/quanlity_control_info.txt'},
     {'sample_id': 'HLN_preC1_2295618',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/HLN_preC1_2295618/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/HLN_preC1_2295618/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/HLN_preC1_2295618/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/HLN_preC1_2295618/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/HLN_preC1_2295618/quanlity_control_info.txt'},
     {'sample_id': 'GHC_preC1_2297515',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/GHC_preC1_2297515/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/GHC_preC1_2297515/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/GHC_preC1_2297515/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/GHC_preC1_2297515/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/GHC_preC1_2297515/quanlity_control_info.txt'},
     {'sample_id': 'WKD_preC1_2274363',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/WKD_preC1_2274363/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/WKD_preC1_2274363/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/WKD_preC1_2274363/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/WKD_preC1_2274363/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/WKD_preC1_2274363/quanlity_control_info.txt'},
     {'sample_id': 'TWH_op_2280553',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/TWH_op_2280553/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/TWH_op_2280553/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/TWH_op_2280553/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/TWH_op_2280553/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/TWH_op_2280553/quanlity_control_info.txt'},
     {'sample_id': 'WXY_pbmc_op_2295902',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/WXY_pbmc_op_2295902/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/WXY_pbmc_op_2295902/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/WXY_pbmc_op_2295902/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/WXY_pbmc_op_2295902/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/WXY_pbmc_op_2295902/quanlity_control_info.txt'},
     {'sample_id': 'WXY_preC1_2295902',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/WXY_preC1_2295902/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/WXY_preC1_2295902/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/WXY_preC1_2295902/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/WXY_preC1_2295902/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/WXY_preC1_2295902/quanlity_control_info.txt'},
     {'sample_id': 'SZR_pbmc_op_2279307',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/SZR_pbmc_op_2279307/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/SZR_pbmc_op_2279307/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/SZR_pbmc_op_2279307/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/SZR_pbmc_op_2279307/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/SZR_pbmc_op_2279307/quanlity_control_info.txt'},
     {'sample_id': 'THX_pbmc_op_2293786',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/THX_pbmc_op_2293786/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/THX_pbmc_op_2293786/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/THX_pbmc_op_2293786/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/THX_pbmc_op_2293786/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/THX_pbmc_op_2293786/quanlity_control_info.txt'},
     {'sample_id': 'WJJ_pbmc_op_2285056',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/WJJ_pbmc_op_2285056/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/WJJ_pbmc_op_2285056/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/WJJ_pbmc_op_2285056/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/WJJ_pbmc_op_2285056/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/WJJ_pbmc_op_2285056/quanlity_control_info.txt'},
     {'sample_id': 'ZMZ_preC1_2267536',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/ZMZ_preC1_2267536/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/ZMZ_preC1_2267536/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/ZMZ_preC1_2267536/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/ZMZ_preC1_2267536/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/ZMZ_preC1_2267536/quanlity_control_info.txt'},
     {'sample_id': 'CJC_preC1_2274041',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/CJC_preC1_2274041/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/CJC_preC1_2274041/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/CJC_preC1_2274041/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/CJC_preC1_2274041/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/CJC_preC1_2274041/quanlity_control_info.txt'},
     {'sample_id': 'HJJ_op_2297052',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/HJJ_op_2297052/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/HJJ_op_2297052/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/HJJ_op_2297052/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/HJJ_op_2297052/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/HJJ_op_2297052/quanlity_control_info.txt'},
     {'sample_id': 'YLB_preC1_2277686',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/YLB_preC1_2277686/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/YLB_preC1_2277686/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/YLB_preC1_2277686/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/YLB_preC1_2277686/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/YLB_preC1_2277686/quanlity_control_info.txt'},
     {'sample_id': 'WJJ_op_2285056',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/WJJ_op_2285056/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/WJJ_op_2285056/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/WJJ_op_2285056/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/WJJ_op_2285056/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/WJJ_op_2285056/quanlity_control_info.txt'},
     {'sample_id': 'DJZ_op_2301100',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/DJZ_op_2301100/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/DJZ_op_2301100/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/DJZ_op_2301100/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/DJZ_op_2301100/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/DJZ_op_2301100/quanlity_control_info.txt'},
     {'sample_id': 'QBG_preC2_2300416',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/QBG_preC2_2300416/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/QBG_preC2_2300416/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/QBG_preC2_2300416/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/QBG_preC2_2300416/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/QBG_preC2_2300416/quanlity_control_info.txt'},
     {'sample_id': 'YDY_preC1_2271160',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/YDY_preC1_2271160/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/YDY_preC1_2271160/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/YDY_preC1_2271160/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/YDY_preC1_2271160/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/YDY_preC1_2271160/quanlity_control_info.txt'},
     {'sample_id': 'WFC_preC1_2268478',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/WFC_preC1_2268478/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/WFC_preC1_2268478/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/WFC_preC1_2268478/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/WFC_preC1_2268478/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/WFC_preC1_2268478/quanlity_control_info.txt'},
     {'sample_id': 'GHC_preC2_2297515',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/GHC_preC2_2297515/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/GHC_preC2_2297515/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/GHC_preC2_2297515/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/GHC_preC2_2297515/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/GHC_preC2_2297515/quanlity_control_info.txt'},
     {'sample_id': 'LDS_preC1_2304172',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/LDS_preC1_2304172/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/LDS_preC1_2304172/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/LDS_preC1_2304172/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/LDS_preC1_2304172/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/LDS_preC1_2304172/quanlity_control_info.txt'},
     {'sample_id': 'ZJK_preC1_2266252',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/ZJK_preC1_2266252/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/ZJK_preC1_2266252/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/ZJK_preC1_2266252/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/ZJK_preC1_2266252/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/ZJK_preC1_2266252/quanlity_control_info.txt'},
     {'sample_id': 'MHY_preC1_2268689',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/MHY_preC1_2268689/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/MHY_preC1_2268689/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/MHY_preC1_2268689/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/MHY_preC1_2268689/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/MHY_preC1_2268689/quanlity_control_info.txt'},
     {'sample_id': 'LCL_preC1_2284745',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/LCL_preC1_2284745/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/LCL_preC1_2284745/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/LCL_preC1_2284745/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/LCL_preC1_2284745/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/LCL_preC1_2284745/quanlity_control_info.txt'},
     {'sample_id': 'HLN_preC2_2295618',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/HLN_preC2_2295618/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/HLN_preC2_2295618/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/HLN_preC2_2295618/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/HLN_preC2_2295618/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/HLN_preC2_2295618/quanlity_control_info.txt'},
     {'sample_id': 'XSY_preC2_2308306',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/XSY_preC2_2308306/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/XSY_preC2_2308306/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/XSY_preC2_2308306/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/XSY_preC2_2308306/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/XSY_preC2_2308306/quanlity_control_info.txt'},
     {'sample_id': 'ZHM_preC2_2301356',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/ZHM_preC2_2301356/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/ZHM_preC2_2301356/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/ZHM_preC2_2301356/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/ZHM_preC2_2301356/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/ZHM_preC2_2301356/quanlity_control_info.txt'},
     {'sample_id': 'LDL_preC1_2266790',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/LDL_preC1_2266790/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/LDL_preC1_2266790/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/LDL_preC1_2266790/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/LDL_preC1_2266790/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/LDL_preC1_2266790/quanlity_control_info.txt'},
     {'sample_id': 'PLR_preC2_2269964',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/PLR_preC2_2269964/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/PLR_preC2_2269964/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/PLR_preC2_2269964/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/PLR_preC2_2269964/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/PLR_preC2_2269964/quanlity_control_info.txt'},
     {'sample_id': 'LGZ_preC2_2300309',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/LGZ_preC2_2300309/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/LGZ_preC2_2300309/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/LGZ_preC2_2300309/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/LGZ_preC2_2300309/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/LGZ_preC2_2300309/quanlity_control_info.txt'},
     {'sample_id': 'YLY_op_2282682',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/YLY_op_2282682/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/YLY_op_2282682/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/YLY_op_2282682/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/YLY_op_2282682/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/YLY_op_2282682/quanlity_control_info.txt'},
     {'sample_id': 'CDZ_pbmc_op_2282421',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/CDZ_pbmc_op_2282421/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/CDZ_pbmc_op_2282421/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/CDZ_pbmc_op_2282421/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/CDZ_pbmc_op_2282421/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/CDZ_pbmc_op_2282421/quanlity_control_info.txt'},
     {'sample_id': 'XDQ_pbmc_op_2286869',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/XDQ_pbmc_op_2286869/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/XDQ_pbmc_op_2286869/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/XDQ_pbmc_op_2286869/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/XDQ_pbmc_op_2286869/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/XDQ_pbmc_op_2286869/quanlity_control_info.txt'},
     {'sample_id': 'HLC_op_2306665',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/HLC_op_2306665/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/HLC_op_2306665/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/HLC_op_2306665/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/HLC_op_2306665/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/HLC_op_2306665/quanlity_control_info.txt'},
     {'sample_id': 'SZR_op_2279307',
      'filtered_matrix_dir_path': 'input/final_cellranger_data/SZR_op_2279307/sample_filtered_feature_bc_matrix',
      'raw_matrix_h5_path': 'input/final_cellranger_data/SZR_op_2279307/raw_feature_bc_matrix.h5',
      'cluster_path': 'input/final_cellranger_data/SZR_op_2279307/clusters.csv',
      'quanlity_control_h5_path': 'input/final_cellranger_data/SZR_op_2279307/quanlity_control.h5ad',
      'quanlity_control_info_path': 'input/final_cellranger_data/SZR_op_2279307/quanlity_control_info.txt'}]




```python
print(f"找到 {len(samples)} 个有效样本")
```

    找到 127 个有效样本


### 统计下每个样本的细胞数,污染比例,soupx过滤后细胞数,去除低质量细胞后细胞数，双细胞数目过滤双细胞后剩余细胞数


```python
sample_qc_info = {
                "sample_id" : sample["sample_id"],
                "cellranger_filtered_cell_count": "",
                "contamination_fraction_bysoupx": "",
                "cell_count_after_soupx": "",
                "cell_count_after_lowquanlity_filter": "",
                "double_cell_numbers":"",
            }
```

### 读取filtered count matrix


```python
adata = sc.read_10x_mtx(
            sample["filtered_matrix_dir_path"],
            var_names='gene_symbols'
        )
# 确保基因名唯一
adata.var_names_make_unique()
sample_qc_info["cellranger_filtered_cell_count"] = adata.n_obs
```

    /home/data/t190513/miniconda3/envs/sc_preprocess/lib/python3.12/site-packages/scanpy/readwrite.py:585: UserWarning: Suffix used (-[0-9]+) to deduplicate index values may make index values difficult to interpret. There values with a similar suffixes in the index. Consider using a different delimiter by passing `join={delimiter}`. Example key collisions generated by the make_index_unique algorithm: ['SNORD116-1', 'SNORD116-2', 'SNORD116-3', 'SNORD116-4', 'SNORD116-5']
      adata = _read_10x_mtx(


## 质控步骤1 去掉cell-free RNA-soupX


```python
import logging
```


```python
import anndata2ri
```


```python
import rpy2.rinterface_lib.callbacks as rcb
```


```python
import rpy2.robjects as ro
```


```python
rcb.logger.setLevel(logging.ERROR)
```


```python
ro.pandas2ri.activate()
```


```python
anndata2ri.activate()
```


```python
%load_ext rpy2.ipython
```


```r
%%R
library(SoupX)
```

    
        WARNING: The R package "reticulate" only fixed recently
        an issue that caused a segfault when used with rpy2:
        https://github.com/rstudio/reticulate/pull/1188
        Make sure that you use a version of that package that includes
        the fix.
        


```python
### SoupX需要使用的信息包括简单的cluster信息,cellranger的raw count以及filtered count
# 来自cellranger 的couster信息
soupx_groups = pd.read_csv(sample["cluster_path"], encoding='utf-8')
soupx_groups = soupx_groups.set_index("Barcode")["Cluster"]


#来自cellranger 的filtered count matrix
cells = adata.obs_names
genes = adata.var_names
data = adata.X.T ##因为SoupX要的矩阵是features x barcodes所以要转置一下

#来自cellranger 的raw count matrix
adata_raw = sc.read_10x_h5(
    sample["raw_matrix_h5_path"]
)
adata_raw.var_names_make_unique()
data_tod = adata_raw.X.T ###tod指的是table of droplets也就是cellranger的原始矩阵

```

    /home/data/t190513/miniconda3/envs/sc_preprocess/lib/python3.12/site-packages/anndata/_core/anndata.py:1793: UserWarning: Variable names are not unique. To make them unique, call `.var_names_make_unique`.
      utils.warn_names_duplicates("var")
    /home/data/t190513/miniconda3/envs/sc_preprocess/lib/python3.12/site-packages/anndata/_core/anndata.py:1793: UserWarning: Variable names are not unique. To make them unique, call `.var_names_make_unique`.
      utils.warn_names_duplicates("var")
    /tmp/ipykernel_3184233/2525819892.py:16: UserWarning: Suffix used (-[0-9]+) to deduplicate index values may make index values difficult to interpret. There values with a similar suffixes in the index. Consider using a different delimiter by passing `join={delimiter}`. Example key collisions generated by the make_index_unique algorithm: ['SNORD116-1', 'SNORD116-2', 'SNORD116-3', 'SNORD116-4', 'SNORD116-5']
      adata_raw.var_names_make_unique()



```r
%%R -i data -i data_tod -i genes -i cells -i soupx_groups -o out -o contamination_fraction

# specify row and column names of data
rownames(data) = genes
colnames(data) = cells
# ensure correct sparse format for table of counts and table of droplets
data <- as(data, "sparseMatrix")
data_tod <- as(data_tod, "sparseMatrix")

# Generate SoupChannel Object for SoupX 
sc = SoupChannel(data_tod, data, calcSoupProfile = FALSE)

# Add extra meta data to the SoupChannel object
soupProf = data.frame(row.names = rownames(data), est = rowSums(data)/sum(data), counts = rowSums(data))
sc = setSoupProfile(sc, soupProf)

# Set cluster information in SoupChannel
sc = setClusters(sc, soupx_groups)

# Estimate contamination fraction
sc  = autoEstCont(sc, doPlot=FALSE)
contamination_fraction = sc$fit$rhoEst
# Infer corrected table of counts and rount to integer
out = adjustCounts(sc, roundToInt = TRUE)
```


```python
adata.layers["counts"] = adata.X
adata.layers["soupX_counts"] = out.T
adata.X = adata.layers["soupX_counts"]
```


```python
sample_qc_info["contamination_fraction_bysoupx"] = contamination_fraction.item()
sample_qc_info["cell_count_after_soupx"] = adata.n_obs
```


```python
sample_qc_info
```




    {'sample_id': 'YXY_preC2_2282374',
     'cellranger_filtered_cell_count': 12216,
     'contamination_fraction_bysoupx': 0.019,
     'cell_count_after_soupx': 12216,
     'cell_count_after_lowquanlity_filter': '',
     'double_cell_numbers': ''}



## 质控步骤2 去掉低质量细胞


```python
# mitochondrial genes
adata.var["mt"] = adata.var_names.str.startswith("MT-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes.
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")
```


```python
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True
)
adata
```




    AnnData object with n_obs × n_vars = 12216 × 62700
        obs: 'n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'pct_counts_in_top_20_genes', 'total_counts_mt', 'log1p_total_counts_mt', 'pct_counts_mt', 'total_counts_ribo', 'log1p_total_counts_ribo', 'pct_counts_ribo', 'total_counts_hb', 'log1p_total_counts_hb', 'pct_counts_hb'
        var: 'gene_ids', 'feature_types', 'mt', 'ribo', 'hb', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts'
        layers: 'counts', 'soupX_counts'




```python
def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier
```


```python
adata.obs["outlier"] = (
    is_outlier(adata, "log1p_total_counts", 5)
    | is_outlier(adata, "log1p_n_genes_by_counts", 5)
    | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
)
adata.obs.outlier.value_counts()
```




    outlier
    False    11094
    True      1122
    Name: count, dtype: int64




```python
adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (
    adata.obs["pct_counts_mt"] > 8
)
adata.obs.mt_outlier.value_counts()
```




    mt_outlier
    False    8445
    True     3771
    Name: count, dtype: int64




```python
adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()
```


```python
sample_qc_info["cell_count_after_lowquanlity_filter"] = adata.n_obs
```


```python
print(f"Total number of genes: {adata.n_vars}")

# Min 20 cells - filters out 0 count genes
sc.pp.filter_genes(adata, min_cells=20)
print(f"Number of genes after cell filter: {adata.n_vars}")
```

    Total number of genes: 62700
    Number of genes after cell filter: 19772


## 质控步骤3 标记双细胞


```r
%%R
library(Seurat)
library(scater)
library(scDblFinder)
library(BiocParallel)
```


```python
data_mat = adata.X.T
```


```r
%%R -i data_mat -o doublet_score -o doublet_class

set.seed(123)
sce = scDblFinder(
SingleCellExperiment(
    list(counts=data_mat),
) 
)
doublet_score = sce$scDblFinder.score
doublet_class = sce$scDblFinder.class
```


```python
adata.obs["scDblFinder_score"] = doublet_score
adata.obs["scDblFinder_class"] = doublet_class
adata.obs.scDblFinder_class.value_counts()
```




    scDblFinder_class
    singlet    7314
    doublet     747
    Name: count, dtype: int64




```python
sample_qc_info["double_cell_numbers"] = int((adata.obs['scDblFinder_class'] == 'doublet').sum())
sample_qc_info
```




    {'sample_id': 'YXY_preC2_2282374',
     'cellranger_filtered_cell_count': 12216,
     'contamination_fraction_bysoupx': 0.019,
     'cell_count_after_soupx': 12216,
     'cell_count_after_lowquanlity_filter': 8061,
     'double_cell_numbers': 747}




```python
adata.write(sample["quanlity_control_h5_path"])
# 保存到文本文件
with open(sample["quanlity_control_info_path"], "w") as f:
    f.write(str(sample_qc_info))  # 将元组转为字符串写入
```
