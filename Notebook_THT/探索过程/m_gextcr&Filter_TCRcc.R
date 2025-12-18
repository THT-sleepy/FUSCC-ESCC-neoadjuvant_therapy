### 这个代码首先把转录组和TCR数据整合起来


### 这个代码其次是处理一下CD4和CD8共享的clone
### 按照董晨2025肝癌cancer cell的做法，首先
### 计算克隆里CD4细胞数/CD8细胞数的值,>5则去
### 除克隆里的CD8细胞,<0.2去除CD4细胞，否则
### 去除该clone

### 去除后做了一些基本的统计

### loading packages
library(dplyr)
library(ggpubr)
library(ggplot2)
library(scales)

### loading files
df_gex <- read.delim("input/df_1216.tsv")
df_tcr <- read.delim("input/df_cc_1216.tsv")
df_gex$raw_cellid <- sapply(
  strsplit(as.character(df_gex$cell_id), "-"), # 按"-"分割（需转字符型避免因子）
  function(x) paste(head(x, 2), collapse = "-") # 取前2项+重新连接
)
df_tcr$raw_cellid <- sapply(
  strsplit(as.character(df_tcr$cell_id), "-"), # 按"-"分割（需转字符型避免因子）
  function(x) paste(head(x, 2), collapse = "-") # 取前2项+重新连接
)
df_tcr$cell_id <- NULL
df_merge <- df_gex %>%
  left_join(df_tcr,by=c("sample","raw_cellid"))


#sum(!is.na(df_merge$cc_aa_identity)) 141733个细胞
df_cc <- df_merge[!is.na(df_merge$cc_aa_identity),]
#table(df_cc$minor_celltype)
#将注释非T细胞对应的cc设为NA(可能是双细胞什么的)
t_clusters = c( 'c01_CD4_Tn_TCF7',
               'c02_CD4_Tcm_S1PR1',
               'c03_CD4_Tcm_ANXA1',
               'c04_CD4_Tem_NR4A1',
               'c05_CD4_Treg_FOXP3',
               'c06_CD4_TOX',
               'c07_CD8_Tem/Teff_TNFSF9',
               'c08_CD8_Tem/Teff_HSPA1B',
               'c09_CD8_Tem/Teff_GNLY',
               'c10_CD8_Teff_GZMK',
               'c11_CD8_Tex_CXCL13',
               'c12_CD8_Temra_FGFBP2',
               'c13_CD8_IL7R',
               'c14_Tprf_MKI67')
df_cc <- df_cc[df_cc$minor_celltype %in% t_clusters,] #136305个T细胞

### filter
df <- df_cc
df$type <- sapply(strsplit(as.character(df$minor_celltype), "_"), `[`, 2)
df_filtered <- df %>%
  # 按clone.id分组，统计CD4/CD8的数量（按行数统计：n()；按计数和：sum(count)）
  group_by(cc_aa_identity) %>%
  mutate(
    cd4_n = sum(type == "CD4"),  # 该clone的CD4数量（行数）
    cd8_n = sum(type == "CD8"),  # 该clone的CD8数量（行数）
    ratio = ifelse(cd8_n == 0, Inf, cd4_n / cd8_n),  # 计算CD4/CD8比例，处理分母为0
    # 标记每个clone的处理规则
    rule = case_when(
      ratio > 5 ~ "remove_CD8",    # 比例>5：删CD8
      ratio < 0.2 ~ "remove_CD4",  # 比例<0.2：删CD4
      cd4_n > 0 & cd8_n > 0 ~ "remove_clone",  # 0.2~5：删整个clone
      TRUE ~ "keep"  # 仅CD4/仅CD8：保留
    )
  ) %>%
  ungroup() %>%
  # 按规则筛选行：保留符合条件的obs
  filter(
    rule == "keep" |
      (rule == "remove_CD8" & type != "CD8") |
      (rule == "remove_CD4" & type != "CD4")
  ) %>%
  # 清理临时列（可选，保留也不影响分析）
  select(-cd4_n, -cd8_n, -ratio, -rule)
#计算filter后的clone大小
df_filtered1 <- df_filtered %>%
  mutate(clone.id=cc_aa_identity) %>%
  group_by(clone.id) %>%
  mutate(clone_size = n()) %>%
  ungroup() %>%
  select(clone.id,clone_size,sample,raw_cellid)

# 将数据合并到df_gex中
df_merge2 <- df_gex %>%
  left_join(df_filtered1,by=c('sample','raw_cellid'))

write.table(df_merge2,"input/gex_tcr_combined_1216.tsv",sep = "\t",quote = F,row.names = F)

### 计算cd4和cd8各有多少克隆，克隆大小的范围
df_filtered$clone.id <- df_filtered$cc_aa_identity
clone_stats <- df_filtered %>%
  
  # 步骤1：按【clone.id + type】分组，计算每个克隆在对应类型下的大小（每行=1个细胞）
  group_by(clone.id, type) %>%
  summarise(
    clone_size = n(),  # 单个克隆的大小（细胞数）
    .groups = "drop"   # 完全取消分组，转为普通数据框
  ) %>%
  # 步骤2：仅按【type】分组，统计该类型下所有克隆的整体信息
  group_by(type) %>%
  summarise(
    克隆数量 = n_distinct(clone.id),  # CD4/CD8各自的克隆总数（唯一clone.id的数量）
    最小克隆大小 = min(clone_size),   # 该类型下所有克隆的最小大小
    最大克隆大小 = max(clone_size),   # 该类型下所有克隆的最大大小
    平均克隆大小 = round(mean(clone_size), 2),  # 平均值（保留2位小数）
    中位克隆大小 = median(clone_size) # 中位数
  ) %>%
  ungroup()
cat("===== CD4/CD8 克隆统计结果 =====\n")
print(clone_stats)  # 长格式（简洁）

### 绘图
clone_unique <- df_filtered %>%
  group_by(clone.id) %>%
  mutate(clone_size = n()) %>%
  ungroup() %>%
  distinct(clone.id, .keep_all = TRUE) 

# 步骤2：统计每个clone_size的克隆数量（即count）
clone_size_count <- clone_unique %>%
  count(clone_size, name = "count",type) %>%  # 按clone_size分组，统计数量→列名count
  arrange(clone_size)  # 按克隆大小升序排列
  
  
##cd8
clone_size_count %>%
  filter(type=="CD8") %>%
  ggplot(aes(x = clone_size,y=count)) + 
  geom_bar(stat = "identity",width=3,fill = "#8A281E",color = "#8A281E")+
  scale_y_log10()+
  theme_pubr()+
  labs(title="CD8+ T cell clone size distribution",y="Number of clonotypes",x="Clone size")+
  theme(legend.position = "none",
        axis.text.x = element_text(size=15), 
        axis.text.y = element_text(size=15),
        axis.title = element_text(size = 20),)
##cd4
clone_size_count %>%
  filter(type=="CD4") %>%
  ggplot(aes(x = clone_size,y=count)) + 
  geom_bar(stat = "identity",width=0.5,fill = "#5778A2",color = "#5778A2")+
  scale_y_log10()+
  theme_pubr()+
  labs(title="CD4+ T cell clone size distribution",y="Number of clonotypes",x="Clone size")+
  theme(legend.position = "none",
        axis.text.x = element_text(size=15), 
        axis.text.y = element_text(size=15),
        axis.title = element_text(size = 20),)

if(0){
## 
target_celltype <- "c11_CD8_Tex_CXCL13"

# 提取该细胞类型对应的所有唯一clone.id（过滤NA/空克隆）
tex_clone_ids <- df_filtered %>%
  filter(
    minor_celltype == target_celltype,  # 匹配目标细胞类型
  ) %>%
  distinct(clone.id) %>%  # 去重，得到该细胞类型关联的所有克隆ID
  pull(clone.id)           # 转为向量，方便后续筛选

# 验证：打印提取的克隆ID
cat("c11_CD8_Tex_CXCL13关联的克隆数：", length(tex_clone_ids), "\n")

tex_clone_all_cells <- df_filtered %>%
  filter(clone.id %in% target_clone_ids) 

table(tex_clone_all_cells[tex_clone_all_cells$minor_celltype!="c11_CD8_Tex_CXCL13" & 
                          tex_clone_all_cells$response_mpr==1,]$treatment_stage)
#pre on post
#2337 1540 851
table(df_filtered[df_filtered$response_mpr==1 &
                  df_filtered$minor_celltype!="c11_CD8_Tex_CXCL13",]$treatment_stage)
#pre on post
#14186 10029 40989

#MPR者Texp变化是pre0.165-on0.085-post0.038

table(tex_clone_all_cells[tex_clone_all_cells$minor_celltype!="c11_CD8_Tex_CXCL13" & 
                            tex_clone_all_cells$response_mpr==0,]$treatment_stage)
#pre on post
#1068 1089 1643
table(df_filtered[df_filtered$response_mpr==0 &
                    df_filtered$minor_celltype!="c11_CD8_Tex_CXCL13",]$treatment_stage)
#pre on post
#9313 8664 36622

#Non-MPR者Texp变化是pre0.114-on0.126-post0.045

table(tex_clone_all_cells[tex_clone_all_cells$minor_celltype!="c11_CD8_Tex_CXCL13" & 
                            tex_clone_all_cells$response_pcr==1,]$treatment_stage)
#pre on post
#1865 229 721
table(df_filtered[df_filtered$response_pcr==1 &
                    df_filtered$minor_celltype!="c11_CD8_Tex_CXCL13",]$treatment_stage)
#pre on post
#6647 1174 11377
#pCR者Texp变化是pre0.28-on0.195-post0.144
}



