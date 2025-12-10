#!/bin/bash
##mad=5 hvg=4000
seq 0 120 | xargs -I {} -P 10 sh -c 'ipython preprocess01-quanlity_control.py {}  1>tmp/{}_qc.log 2>&1'
ipython preprocess02-1-merge.py  1>tmp/merge.log 2>&1
ipython preprocess02-normalization.py  1>tmp/normalization.log 2>&1
ipython preprocess03-FeatureSelection.py  4000 1>tmp/feature_selection.log 2>&1
ipython preprocess04-DimensionalityReduction.py  4000 1>tmp/DimensionalityReduction.log 2>&1
python preprocess05-Batch_Doublet_Removing.py  4000 1>tmp/Batch_Doublet_Removing.log 2>&1
python cellstructure01-clustering.py  4000 1>tmp/clustering.log 2>&1

##mad=5 hvg=3000
#seq 0 126 | xargs -I {} -P 10 sh -c 'ipython preprocess01-quanlity_control.py  {}  1>tmp/{}_qc.log 2>&1'
#ipython preprocess02-normalization.py  1>tmp/normalization.log 2>&1
#ipython preprocess03-FeatureSelection.py  3000 1>tmp/feature_selection.log 2>&1
#ipython preprocess04-DimensionalityReduction.py  3000 1>tmp/DimensionalityReduction.log 2>&1
#python preprocess05-Batch_Doublet_Removing.py  3000 1>tmp/Batch_Doublet_Removing.log 2>&1
#python cellstructure01-clustering.py 3000 1>tmp/clustering.log 2>&1

##mad=5 hvg=2000
#seq 0 126 | xargs -I {} -P 10 sh -c 'ipython preprocess01-quanlity_control.py  {}  1>tmp/{}_qc.log 2>&1'
#ipython preprocess02-normalization.py  1>tmp/normalization.log 2>&1
#ipython preprocess03-FeatureSelection.py  2000 1>tmp/feature_selection.log 2>&1
#ipython preprocess04-DimensionalityReduction.py  2000 1>tmp/DimensionalityReduction.log 2>&1
#python preprocess05-Batch_Doublet_Removing.py  2000 1>tmp/Batch_Doublet_Removing.log 2>&1
#python cellstructure01-clustering.py 2000 1>tmp/clustering.log 2>&1










##mad=2
#seq 0 126 | xargs -I {} -P 10 sh -c 'ipython preprocess01-quanlity_control.py 2 {}  1>tmp/{}_qc.log 2>&1'
#ipython preprocess02-normalization.py 2 1>tmp/normalization.log 2>&1
#ipython preprocess03-FeatureSelection.py 2 1>tmp/feature_selection.log 2>&1
#ipython preprocess04-DimensionalityReduction.py 2 1>tmp/DimensionalityReduction.log 2>&1
