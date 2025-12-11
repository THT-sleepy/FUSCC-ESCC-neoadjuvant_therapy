# TCR数据质控和预处理         (What)

* Dec 11, 2025                                 (When)
* liyi /data1/liyi/zhaoyue/FUSCC-ESCC-neoadjuvant_thrapy           (Where)
* 查看TCR数据的大致情况，定义clonotype (Why)

### 把所有TCR数据整合成一个文件

```linux
# 第一步：提取表头（不变）
find . -type f -name "filtered_contig_annotations.csv" -print -quit | xargs head -n 1 > merged_127samples_tcr.csv

# 第二步：加 -q 禁用分隔符，跳过表头追加数据
find . -type f -name "filtered_contig_annotations.csv" | xargs tail -q -n +2 >> merged_127samples_tcr.csv

```

### 导入模块和文件路径
```
import warnings
warnings.filterwarnings(    
"ignore",    
".*IProgress not found*",
)
warnings.simplefilter(action="ignore", category=FutureWarning)
import pandas as pd
import scanpy as sc
import scirpy as ir
warnings.simplefilter(action="ignore", category=pd.errors.DtypeWarning)
path_tcr_input = "input/merged_127samples_tcr.csv"
```

### 需要先把rawdata导入进来，后面才能够添加sample等信息(因为这些信息不会被scirpy自动导入)
```
df_tcr_raw = pd.read_csv(path_tcr_input)
# The column 'productive' contains mixed data types which are not compatible with downstream tools.
# We correct them by casting them to strings.
df_tcr_raw["productive"] = df_tcr_raw["productive"].astype(str)
print(f"Total measurements: {len(df_tcr_raw)}")
df_tcr_raw.head(5)
df_tcr_raw["full_length"].value_counts()
df_tcr_raw["productive"].value_counts()

patient_information = [
    "barcode",
    "sample"
]
df_patient = df_tcr_raw[patient_information].copy()

df_patient = df_patient.drop_duplicates(subset=["barcode"], keep="first").reset_index(drop=True)

# Assigning barcode as index（此时索引已无重复）
df_patient.index = df_patient.pop("barcode")
df_patient.index.name = None
adata_tcr.obs[df_patient.columns] = df_patient

```


### 导入数据并去掉不要的7个样本(6个非化免,1个质控不合格)
```
adata_tcr = ir.io.read_10x_vdj(path_tcr_input)
ir.pp.index_chains(mdata)
ir.tl.chain_qc(mdata)
```
