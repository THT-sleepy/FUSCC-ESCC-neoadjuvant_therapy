#!/usr/bin/env python
# coding: utf-8

import sys
import os
import subprocess
import tempfile
import json

import matplotlib.pyplot as plt
import numpy as np
np.float_ = np.float64
import scanpy as sc
import seaborn as sns

sc.settings.verbosity = 0
sc.settings.set_figure_params(
    dpi=80,
    facecolor="white",
    frameon=False,
)

# ## 挑选高变基因

input_fpath = "RData/"+"mad"+sys.argv[1]+"_"+"normalization.h5ad"
adata = sc.read(input_fpath)

# 使用subprocess调用R脚本
r_script = '''
library(scry)
library(zellkonverter)
library(SingleCellExperiment)

args <- commandArgs(trailingOnly=TRUE)
input_fpath <- args[1]
output_fpath <- args[2]

adata <- readH5AD(input_fpath, use_hdf5=TRUE)
sce <- devianceFeatureSelection(adata, assay="X")

binomial_deviance <- rowData(sce)$binomial_deviance
writeLines(as.character(binomial_deviance), con=output_fpath)
'''

# 创建临时文件
with tempfile.NamedTemporaryFile(mode='w', suffix='.R', delete=False) as r_file:
    r_file.write(r_script)
    r_script_path = r_file.name

with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as output_file:
    output_path = output_file.name

try:
    # 执行R脚本
    subprocess.run([
        'Rscript', r_script_path, 
        input_fpath, 
        output_path
    ], check=True)
    
    # 读取结果
    with open(output_path, 'r') as f:
        binomial_deviance = np.array([float(x) for x in f.read().splitlines()])
        
finally:
    # 清理临时文件
    os.unlink(r_script_path)
    os.unlink(output_path)

hvg_n = int(sys.argv[2])
# 其余处理代码...
idx = binomial_deviance.argsort()[-hvg_n:]
mask = np.zeros(adata.var_names.shape, dtype=bool)
mask[idx] = True

adata.var["highly_deviant"] = mask
adata.var["binomial_deviance"] = binomial_deviance

sc.pp.highly_variable_genes(adata, layer="log1p_norm")

output_fpath = "RData/"+"mad"+sys.argv[1]+"_"+"hvgn"+sys.argv[2]+"_feature_selection.h5ad"
adata.write(output_fpath)

print("脚本执行完成")
