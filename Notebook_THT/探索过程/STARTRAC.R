### loading packages
library("Startrac")
library("tictoc")
library("ggpubr")
library("ggplot2")
library("ComplexHeatmap")
library("RColorBrewer")
library("circlize")
library("data.table")

### loading example files
in.dat <- read.delim("input/1211_filtered_Tclonotype.tsv")
in.dat$patient <- paste(in.dat$patient,in.dat$treatment_stage,sep="_")
in.dat <- as.data.table(in.dat)
head(in.dat)

### run pipeline
tic("Startrac.run")
out <- Startrac.run(in.dat, proj="ESCC",verbose=F)
toc()

### visualization
#cluster level index of all data
Startrac::plot(out,index.type="cluster.all",byPatient=F)

#cluster level index by patients.
Startrac::plot(out,index.type="cluster.all",byPatient=T)

#pairwise index of all data
Startrac::plot(out,index.type="pairwise.migr",byPatient=F)

Startrac::plot(out,index.type="pairwise.tran",byPatient=T)

### 保存结果
write.table(out@cluster.data,"output/1211_startrac_clusterindex.tsv",sep = "\t",quote = F,row.names = F)
write.table(as.data.table(out@cluster.sig.data)[aid!=out@proj,],"output/1211_startrac_clusterindex_bypatient.tsv",sep = "\t",quote = F,row.names = F)
write.table(out@pIndex.migr,"output/1211_startrac_migr.tsv",sep = "\t",quote = F,row.names = F)
write.table(out@pIndex.tran,"output/1211_startrac_tran.tsv",sep = "\t",quote = F,row.names = F)

### 仅计算肿瘤样本
in.dat <- read.delim("input/1211_filtered_Tclonotype.tsv")
in.dat$patient <- paste(in.dat$patient,in.dat$treatment_stage,sep="_")
in.dat <- as.data.table(in.dat)
tic("Startrac.run")
out <- Startrac.run(in.dat[in.dat$loc=="tumor",], proj="ESCC",verbose=F)
toc()
write.table(out@cluster.data,"output/1211_startrac_clusterindex_onlytumor.tsv",sep = "\t",quote = F,row.names = F)
write.table(as.data.table(out@cluster.sig.data)[aid!=out@proj,],"output/1211_startrac_clusterindex_bypatient_onlytumor.tsv",sep = "\t",quote = F,row.names = F)



