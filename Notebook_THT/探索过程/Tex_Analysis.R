### 此脚本主要是模仿张泽民CRC的文章对Tex进行详细分析

# Library and data loading
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(tidydr)
library(ggrepel)
require(Seurat)
require(ggpubr)


df = read.delim('input/df_1128.tsv',header = T)
#去掉SZR_2279307 post tumor
df$sample <- paste(df$patient_number,df$patient_id,df$treatment_stage,df$anatomic_site,sep="_")
df <- df[df$sample != "P29_2279307_post_tumor",]

df_clin <- read.csv("input/clin_metadata.csv")
load('RData/major_cell_colors.rda')
load('RData/subtype_col.rda')

df_startrac <- read.delim("output/1211_startrac_clusterindex_bypatient_onlytumor.tsv")
df_startrac_expa <- df_startrac %>%
  mutate(minor_celltype=majorCluster,patient_id=as.integer(str_split_fixed(aid, "_", n = Inf)[, 2])) %>%
  filter(index=="expa") %>%
  select(-patient_id)

df_merge <- df %>%
  left_join(df_clin,by=c("patient_id")) %>%
  mutate(Response=ifelse(MPR==0,"Non-MPR",
                  ifelse(MPR==1,"MPR",NA))) %>%
  mutate(aid=paste(patient_id,treatment_stage,sep="_")) %>%
  left_join(df_startrac_expa,by=c("aid","minor_celltype")) 

## UMAP
df_merge %>%
  mutate(UMAP1=minor_umap1,UMAP2=minor_umap2) %>%
  filter(anatomic_site=="tumor",minor_celltype %in% c("c07_CD8_Tem/Teff_TNFSF9","c08_CD8_Tem/Teff_HSPA1B", 
                                                      "c09_CD8_Tem/Teff_GNLY","c10_CD8_Teff_GZMK",
                                                      "c11_CD8_Tex_CXCL13","c12_CD8_Temra_FGFBP2", 
                                                      "c13_CD8_IL7R")) %>%
  ggplot() + 
  geom_point(aes(x=UMAP1, y=UMAP2, color=minor_celltype), size=.01, alpha=1) + 
  scale_color_manual(values = subtype_col, name='T sub-types in tumor') + 
  #     geom_text(aes(x=m1, y=m2, label=anno_brief), data = text_anno, size=5) +
  theme_dr() +
  theme(plot.title = element_text(size=22),
        strip.text = element_text(size = 18),
        strip.background =element_rect(fill='white', color = 'black', size = .1), 
        panel.grid = element_blank(),
        axis.ticks = element_blank(),       # 去掉刻度线
        axis.text = element_blank(),
        axis.title = element_text(size = 18),
        axis.line = element_line(size = .8),
        legend.position = 'right',
        legend.key.size = unit(1.4,'line'),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15, face = 'bold'),
        plot.margin = margin(b=0, l=10, t=10, r=10)) +
  guides(colour = guide_legend(override.aes = list(size=5, alpha = 1)))


df_merge %>%
  mutate(UMAP1=minor_umap1,UMAP2=minor_umap2) %>%
  filter(anatomic_site=="tumor",
         minor_celltype %in% c("c07_CD8_Tem/Teff_TNFSF9","c08_CD8_Tem/Teff_HSPA1B", 
                                                      "c09_CD8_Tem/Teff_GNLY","c10_CD8_Teff_GZMK",
                                                      "c11_CD8_Tex_CXCL13","c12_CD8_Temra_FGFBP2", 
                                                      "c13_CD8_IL7R"),
         !is.na(Response)) %>%
  ggplot(aes(x=UMAP1, y=UMAP2)) + 
  facet_grid(cols=vars(factor(treatment_stage, levels=c('pre','on','post'))),
             rows=vars(Response)) +
  geom_point(aes(x=UMAP1, y=UMAP2, color=minor_celltype), size=.01, alpha=1) + 
  scale_color_manual(values = subtype_col)+
  theme_pubr() +
  theme(plot.title = element_text(size=22),
        strip.text = element_text(size = 18),
        strip.background =element_rect(fill=NA, color = NA, size = .1), 
        axis.text = element_blank(),
        axis.title = element_text(size = 18),
        axis.line = element_line(size = .8),
        legend.position = 'none',
        legend.key.size = unit(1.4,'line'),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15, face = 'bold'),
        plot.margin = margin(b=0, l=10, t=10, r=10))


## Proportion

global_df = df_merge %>%
  filter(middle_celltype %in% c("CD4+T","CD8+T","Tprf")) %>%
  mutate(Tissue=anatomic_site,Patient=patient_id,Stage=treatment_stage) %>%
  group_by(Tissue, Patient, Stage) %>% summarise(n1 = n()) 

sub_df = df_merge %>%
  filter(middle_celltype %in% c("CD4+T","CD8+T","Tprf")) %>%
  mutate(Tissue=anatomic_site,Patient=patient_id,Stage=treatment_stage,annotation=minor_celltype) %>%
  group_by(Tissue, annotation, Patient, Stage) %>% 
  summarise(n = n())


sub_freq_df =
  df_merge %>%
  filter(middle_celltype %in% c("CD4+T","CD8+T","Tprf")) %>%
  mutate(Tissue=anatomic_site,Patient=patient_id,Stage=treatment_stage,annotation=minor_celltype) %>%
  mutate(sample_stage = paste(Patient, Stage, Tissue, sep = '-')) %>% 
  tidyr::expand(sample_stage,annotation) %>%
  separate(sample_stage, into = c('Patient','Stage','Tissue'), sep = '-') %>%
  mutate(Patient=as.integer(Patient)) %>%
  left_join(sub_df, by = c('Tissue', 'Patient', 'Stage','annotation')) %>% 
  left_join(global_df, by = c('Tissue', 'Patient', 'Stage')) %>%
  replace_na(list(n=0, n1=0)) %>%
  mutate(freq = n/n1) %>%
  mutate(patient_id=Patient) %>%
  left_join(df_clin,by = c("patient_id")) %>%
  mutate(Response=ifelse(MPR==0,"Non-MPR",
                                ifelse(MPR==1,"MPR",NA)))


sub_freq_df$Stage = factor(sub_freq_df$Stage, levels = c('pre','on','post'))

# head(sub_freq_df)
# dim(sub_freq_df)

### Line chart

cls = 'c11_CD8_Tex_CXCL13'

input = sub_freq_df %>% 
  filter(Tissue == 'tumor') %>%
  filter(annotation==cls,!is.na(Response))
input %>%
  mutate(group_label = paste(Response, Stage, sep = "_")) %>%
  ggplot(aes(x=Stage, y=freq, color=Response, group=Patient)) + 
  geom_point( size=2, stroke=2) +
  geom_line(alpha=.5) + 
  scale_color_manual(values = get_palette('nejm', 3)[c(1,2)], name='Response') +
  xlab('Sampling stage') + ylab('Frequency(%T cells)') +
  labs(title="Cellular abundance of c11_CD8_Tex_CXCL13")+
  facet_wrap(~Response) +
  theme_pubr() +
  theme(strip.text.x = element_text(size = 18),
        strip.background =element_rect(fill=NA, color = NA, size = 1), 
        axis.text.x = element_text(size=15, angle = 0, hjust = .5),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size = 18),
        axis.line = element_line(size = .8),
        legend.position = 'none',
        legend.key.size = unit(1.4,'line'),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15, face = 'bold'),
        plot.margin = margin(b = 10, t=10, l=10, r=10),
        plot.title = element_text(size = 20))


###计算p值
groups <- c("MPR_on", "MPR_post", "MPR_pre", "Non-MPR_on", "Non-MPR_post", "Non-MPR_pre")
comparison_pairs <- combn(groups, 2, simplify = FALSE)

# 创建结果列表
results <- list()

for (pair in comparison_pairs) {
  group1_data <- input %>%
    filter(paste(Response, Stage, sep = "_") == pair[1]) %>%
    pull(freq)
  
  group2_data <- input %>%
    filter(paste(Response, Stage, sep = "_") == pair[2]) %>%
    pull(freq)
  
  if (length(group1_data) >= 2 && length(group2_data) >= 2) {
    test_res <- wilcox.test(group1_data, group2_data, exact = FALSE)
    p_val <- test_res$p.value
  } else {
    p_val <- NA
  }
  
  results[[paste(pair[1], "vs", pair[2])]] <- data.frame(
    group1 = pair[1],
    group2 = pair[2],
    p_value = p_val
  )
}
p_values_simple <- do.call(rbind, results)
rownames(p_values_simple) <- NULL

print(p_values_simple)





## Expansion 


# head(expa_df)


### line chart

### 前面的df_startrac_expa已经是仅肿瘤样本了，这里不用再筛选
input=df_startrac_expa %>%
  mutate(Patient=str_split_fixed(aid, "_", n = Inf)[, 2],
         patient_id=as.integer(str_split_fixed(aid, "_", n = Inf)[, 2]),
         Stage=str_split_fixed(aid, "_", n = Inf)[, 3]) %>%
  rename(annotation=majorCluster,exp=value) %>%
  filter(annotation=='c11_CD8_Tex_CXCL13') %>%
  left_join(df_clin,by=c("patient_id")) %>%
  mutate(Response=ifelse(MPR==0,"Non-MPR",
                         ifelse(MPR==1,"MPR",NA))) %>%
  filter(!is.na(Response))

input %>%
  ggplot(aes(x=factor(Stage, levels=c('pre','on','post')), y=exp,
             group=Patient, col=Response)) + 
  geom_point(size=2, stroke=2) +
  geom_line(alpha=.5) + 
  scale_color_manual(values = get_palette('nejm', 3)[c(1,2)], name='Response')+
  xlab('') + ylab(paste0('STARTRAC-expa')) +
  facet_wrap(~Response) +
  theme_pubr() +
  theme(strip.text.x = element_text(size = 18),
        strip.background =element_rect(fill=NA, color = NA, size = 1), 
        axis.text.x = element_text(size=15, angle = 0, hjust = .5),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size = 18),
        axis.line = element_line(size = .8),
        legend.position = 'right',
        legend.key.size = unit(1.4,'line'),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15, face = 'bold'),
        plot.margin = margin(b = 10, t=10, l=10, r=10)) +
  guides(colour = guide_legend(override.aes = list(size=3, alpha = 1),ncol=1))

###计算p值
groups <- c("MPR_on", "MPR_post", "MPR_pre", "Non-MPR_on", "Non-MPR_post", "Non-MPR_pre")
comparison_pairs <- combn(groups, 2, simplify = FALSE)

# 创建结果列表
results <- list()

for (pair in comparison_pairs) {
  group1_data <- input %>%
    filter(paste(Response, Stage, sep = "_") == pair[1]) %>%
    pull(exp)
  
  group2_data <- input %>%
    filter(paste(Response, Stage, sep = "_") == pair[2]) %>%
    pull(exp)
  
  if (length(group1_data) >= 2 && length(group2_data) >= 2) {
    test_res <- wilcox.test(group1_data, group2_data, exact = FALSE)
    p_val <- test_res$p.value
  } else {
    p_val <- NA
  }
  
  results[[paste(pair[1], "vs", pair[2])]] <- data.frame(
    group1 = pair[1],
    group2 = pair[2],
    p_value = p_val
  )
}
p_values_simple <- do.call(rbind, results)
rownames(p_values_simple) <- NULL

print(p_values_simple)

## boxplot
## pre-expa MPR vs Non-MPR
## boxplot startrac-expa(MPR vs Non-MPR)
p <- input %>%
  select(Response,Stage,exp,patient_id) %>%
  filter(!is.na(Response)) %>%
  filter(Stage=="pre") %>%
  ggplot(aes(x=Response, y=exp)) +
  geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
  geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
  theme_pubr() + 
  labs(x="",y="STARTRAC-expa") +
  scale_color_manual(values = c("Non-MPR"='#5778A2',"MPR"='#8A281E'), name='Response') +
  scale_x_discrete(labels = c("Non-MPR" = "Non-MPR\n(n=16)", "MPR" = "MPR\n(n=19)")) +
  theme(strip.text.x = element_text(size = 20),
        strip.background =element_rect(fill=NA, color = 'white', size = 1), 
        axis.text.x = element_text(size=15), 
        axis.text.y = element_text(size=15),
        axis.title = element_text(size = 20),
        axis.line = element_line(size = .5),
        #plot.title = element_text(hjust = 0.5,size = 15),
        legend.position = 'none',
        #legend.key.size = unit(1.4,'line'),
        #legend.text = element_text(size = 13),
        #legend.title = element_text(size = 15, face = 'bold'),
        plot.margin = margin(t=30, l=10)
  )
my_comparisons <- list(c("Non-MPR","MPR"))
p1 <- p+stat_compare_means(comparisons = my_comparisons)+
  theme() 
p1
## on-expa MPR vs Non-MPR
p <- input %>%
  select(Response,Stage,exp,patient_id) %>%
  filter(!is.na(Response)) %>%
  filter(Stage=="on") %>%
  ggplot(aes(x=Response, y=exp)) +
  geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
  geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
  theme_pubr() + 
  labs(x="",y="STARTRAC-expa") +
  scale_color_manual(values = c("Non-MPR"='#5778A2',"MPR"='#8A281E'), name='Response') +
  scale_x_discrete(labels = c("Non-MPR" = "Non-MPR\n(n=16)", "MPR" = "MPR\n(n=19)")) +
  theme(strip.text.x = element_text(size = 20),
        strip.background =element_rect(fill=NA, color = 'white', size = 1), 
        axis.text.x = element_text(size=15), 
        axis.text.y = element_text(size=15),
        axis.title = element_text(size = 20),
        axis.line = element_line(size = .5),
        #plot.title = element_text(hjust = 0.5,size = 15),
        legend.position = 'none',
        #legend.key.size = unit(1.4,'line'),
        #legend.text = element_text(size = 13),
        #legend.title = element_text(size = 15, face = 'bold'),
        plot.margin = margin(t=30, l=10)
  )
my_comparisons <- list(c("Non-MPR","MPR"))
p1 <- p+stat_compare_means(comparisons = my_comparisons)+
  theme() 
p1
## post-expa MPR vs Non-MPR
p <- input %>%
  select(Response,Stage,exp,patient_id) %>%
  filter(!is.na(Response)) %>%
  filter(Stage=="post") %>%
  ggplot(aes(x=Response, y=exp)) +
  geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
  geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
  theme_pubr() + 
  labs(x="",y="STARTRAC-expa") +
  scale_color_manual(values = c("Non-MPR"='#5778A2',"MPR"='#8A281E'), name='Response') +
  scale_x_discrete(labels = c("Non-MPR" = "Non-MPR\n(n=16)", "MPR" = "MPR\n(n=19)")) +
  theme(strip.text.x = element_text(size = 20),
        strip.background =element_rect(fill=NA, color = 'white', size = 1), 
        axis.text.x = element_text(size=15), 
        axis.text.y = element_text(size=15),
        axis.title = element_text(size = 20),
        axis.line = element_line(size = .5),
        #plot.title = element_text(hjust = 0.5,size = 15),
        legend.position = 'none',
        #legend.key.size = unit(1.4,'line'),
        #legend.text = element_text(size = 13),
        #legend.title = element_text(size = 15, face = 'bold'),
        plot.margin = margin(t=30, l=10)
  )
my_comparisons <- list(c("Non-MPR","MPR"))
p1 <- p+stat_compare_means(comparisons = my_comparisons)+
  theme() 
p1

## Differential expression analysis
### Plot
### MPR post vs pre
df_deg <- read.csv("input/MPR_post_vs_pre_deg_results.csv")

df_deg %>% 
  mutate(rank=rank(padj, ties.method = "min")) %>%
  mutate(label=ifelse((abs(log2FC)>1) & (rank<2000) &(padj <=0.05) & (!str_detect(gene, "ENSG|LINC|MIR|-AS|IGL")),gene,NA)) %>%
  filter(padj!=1,pct>0.25,abs(log2FC)<5) %>%
  ggplot() + 
  theme_pubr() +
  xlab('Log2(fold change)') + 
  ylab('-Log10 (adj. P value)') +
  geom_point(aes(x=log2FC, y = -log10(padj), 
                 color = abs(log2FC) > 0 & padj <=0.01), size=1) + 
  scale_color_manual(values = c('grey', '#BD4739'), name='') + 
  geom_text_repel(aes(x=log2FC, y = -log10(padj), label=label), 
                  size=4, max.overlaps=100, fontface = "italic") +
  theme(strip.text.x = element_text(size = 18),
        strip.background =element_rect(fill='white', color = 'black', size = 1), 
        axis.text.x = element_text(size=20, angle = 0, hjust = .5),
        axis.text.y = element_text(size=20),
        axis.title = element_text(size = 20),
        axis.line = element_line(size = .8),
        legend.position = 'none',
        legend.key.size = unit(1.4,'line'),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15, face = 'bold'),
        plot.margin = margin(t=10, l=10,r=10)) +
  guides(colour = guide_legend(override.aes = list(size=3, alpha = 1)))


