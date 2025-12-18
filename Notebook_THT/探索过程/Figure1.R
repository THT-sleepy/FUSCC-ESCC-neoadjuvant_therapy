#---- Figure 1 ----#
setwd("~/FUSCC-ESCC-neoadjuvant_thrapy/")

## Library and data loading

library(data.table)
library(dplyr)
library(tidyr)
library(Seurat)
library(ggpubr)
library(ggplot2)
library(tidydr)

df = read.delim('input/df_1128.tsv',header = T)
#去掉SZR_2279307 post tumor
df$sample <- paste(df$patient_number,df$patient_id,df$treatment_stage,df$anatomic_site,sep="_")
df <- df[df$sample != "P29_2279307_post_tumor",]

df_clin <- read.csv("input/clin_metadata.csv")
load('RData/major_cell_colors.rda')
load('RData/subtype_col.rda')
# dim(df)
# 储存umap图片的时候，png是1600x800(大图例)，1000x800(小图例)，800x800(无图例)，pdf是10x8或8x8或16x8



if(0){
major_cell_colors <- c("B cell"="#3A6FA4",
                       "T&ILC cell"="#E0823A",
                       "Myeloid cell"="#4E9669",
                       "Epithelial cell"= "#BA3931",
                       "Endothelial cell" = "#6F5198",
                       "Fibrocyte" = "#81564C"
)
save(major_cell_colors,file="RData/major_cell_colors.rda")
subtype_col <- c(
  'c01_CD4_Tn_TCF7'= '#4279AD',
  'c02_CD4_Tcm_S1PR1'= '#AFCCD3',
  'c03_CD4_Tcm_ANXA1'= '#A25936',
  'c04_CD4_Tem_NR4A1'= '#6171A8',
  'c05_CD4_Treg_FOXP3'= '#BD5BA5',
  'c06_CD4_TOX'= '#D5E1B4',
  'c07_CD8_Tem/Teff_TNFSF9'= '#F2D7D0',
  'c08_CD8_Tem/Teff_HSPA1B'= '#E8A7BA',
  'c09_CD8_Tem/Teff_GNLY'= '#8A7FA0',
  'c10_CD8_Teff_GZMK'= '#F4EFA9',
  'c11_CD8_Tex_CXCL13'= '#E77772',
  'c12_CD8_Temra_FGFBP2'= '#4D4199',
  'c13_CD8_IL7R'= '#F1B9AE',
  'c14_Tprf_MKI67'= '#C5A4C7',
  'c15_NK_CXXC5'= '#DF617B',
  'c16_NK/NKT_LAG3'= '#816AA9',
  'c17_NK_KLRC1'= '#CBD1E7',
  'c18_NK_NRGN'= '#5BA0D4',
  'c19_Bn_TCL1A'= '#5AA452',
  'c20_Bm_TNFRSF13B'= '#F29696',
  'c21_Bm_TUBA1A'= '#D8432E',
  'c22_Bp_MZB1'= '#F19B9B',
  'c23_Mono_CD14'= '#F4BC82',
  'c24_Mono_HDAC9'= '#ABE0ED',
  'c25_Mono_FCGR3A'= '#C0FEFC',
  'c26_Mono_CCL5'= '#D8C6DF',
  'c27_Macro_C1QC'= '#F6C589',
  'c28_Macro_G0S2'= '#F2ACAE',
  'c29_cDC1_CADM1'= '#EF8E42',
  'c30_cDC2_OTUD1'= '#8AA624',
  'c31_cDC2_CD1A'= '#D39FC8',
  'c32_cDC2_CD1D'= '#E5D9E3',
  'c33_cDC3_LAMP3'= '#F2BF7D',
  'c34_pDC_MZB1'= '#E9655B',
  'c35_DC_ECE1'= '#7DD076',
  'c36_Mast_KIT'= '#5D7AB1',
  'c37_VEC_ACKR1'= '#D2E9C7',
  'c38_VEC_ITPKC'= '#9DD1C8',
  'c39_AEC_GJA4'= '#BCDD7A',
  'c40_LEC_PROX1'= '#EA8675',
  'c41_capEC_RGCC'= '#BCBBD7',
  'c42_EC_CSF2RB'= '#9BBBD0',
  'c43_prfEC_MKI67'= '#4E8975',
  'c44_Fib_COL11A1'= '#E3C567',
  'c45_Fib_MMP1'= '#A63D40',
  'c46_Fib_COL6A5'= '#7F97A2',
  'c47_proFib_CD34'= '#8F913D',
  'c48_apFib_HLA-DRB1'= '#884C76',
  'c49_PC_RGS5'= '#48C0BF',
  'c50_PC_MYH11'= '#D9A779',
  'c51_Epi_Normal'= '#7FC696',
  'c52_Epi_Tumor'= '#2E4F73'
)
save(subtype_col,file ="RData/subtype_col.rda")
}
## umap
### Cell subtypes 带图例
df %>% 
  ggplot(aes(x=umap1, y=umap2, color=minor_celltype)) + 
  geom_point(size=.1, shape=16, stroke=0) + 
  theme_dr() + #小箭头坐标轴 
  scale_color_manual(values = subtype_col, name='') + 
  labs(x="UMAP1",y="UMAP2")+
  theme( panel.grid = element_blank(),
         axis.line = element_blank(),        # 去掉坐标轴线条
         axis.ticks = element_blank(),       # 去掉刻度线
         axis.text = element_blank(),        # 去掉刻度标签（如数值）
         axis.title = element_text(size = 15, face = 'bold'),
        legend.key.size = unit(1.4,'line'),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15, face = 'bold'),
        legend.position = 'right')+
  guides(colour = guide_legend(override.aes = list(size=5, alpha = 1)))
### Major cell types

df %>% 
  ggplot(aes(x=umap1, y=umap2, color=major_celltype)) + 
  geom_point(size=.1, shape=16, stroke=0) + 
  theme_dr() + 
  scale_color_manual(values = major_cell_colors, name='') + 
  labs(x="UMAP1",y="UMAP2")+
  theme(panel.grid = element_blank(),
        axis.line = element_blank(),        # 去掉坐标轴线条
        axis.ticks = element_blank(),       # 去掉刻度线
        axis.text = element_blank(),        # 去掉刻度标签（如数值）
        axis.title =element_text(size = 15, face = 'bold'),         # 去掉轴标题（如“umap1”“umap2”）
        legend.key.size = unit(1.4,'line'),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15, face = 'bold'),
        legend.position = 'none')+
  guides(colour = guide_legend(override.aes = list(size=5, alpha = 1)))

### 组织来源
df %>% 
  ggplot(aes(x=umap1, y=umap2, color=anatomic_site)) + 
  geom_point(size=.1, shape=16, stroke=0) + 
  theme_dr() +
  labs(x="UMAP1",y="UMAP2")+
  scale_color_manual(values = c("tumor"="#8A281E","blood"="#5778A2"))+
  theme(panel.grid = element_blank(),
        axis.title =element_text(size = 15, face = 'bold'), 
        legend.key.size = unit(1.4,'line'),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15, face = 'bold'),
        legend.position = 'right')+
  guides(colour = guide_legend(override.aes = list(size=5, alpha = 1)))


### 治疗阶段
df %>% 
  ggplot(aes(x=umap1, y=umap2, color=treatment_stage)) + 
  geom_point(size=.1, shape=16, stroke=0) + 
  scale_color_manual(values=c("pre"="#E7872B","on"="#925A44","post"="#C03830"))+
  theme_dr() + 
  labs(x="UMAP1",y="UMAP2")+
  theme(panel.grid = element_blank(),
        axis.title =element_text(size = 15, face = 'bold'), 
        legend.key.size = unit(1.4,'line'),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15, face = 'bold'),
        legend.position = 'right')+
  guides(colour = guide_legend(override.aes = list(size=5, alpha = 1)))

### patient
df %>% 
  ggplot(aes(x=umap1, y=umap2, color=patient_number)) + 
  geom_point(size=.1, shape=16, stroke=0) + 
  scale_color_manual(values= c(
    "P01" = "#4279AD",
    "P02" = "#AFCCD3",
    "P03" = "#A25936",
    "P04" = "#6171A8",
    "P05" = "#BD5BA5",
    "P06" = "#D5E1B4",
    "P07" = "#F2D7D0",
    "P08" = "#E8A7BA",
    "P09" = "#8A7FA0",
    "P10" = "#F4EFA9",
    "P11" = "#E77772",
    "P12" = "#4D4199",
    "P13" = "#F1B9AE",
    "P14" = "#C5A4C7",
    "P15" = "#DF617B",
    "P16" = "#816AA9",
    "P17" = "#CBD1E7",
    "P18" = "#5BA0D4",
    "P19" = "#5AA452",
    "P20" = "#F29696",
    "P21" = "#D8432E",
    "P22" = "#F19B9B",
    "P23" = "#F4BC82",
    "P24" = "#ABE0ED",
    "P25" = "#C0FEFC",
    "P26" = "#D8C6DF",
    "P27" = "#F6C589",
    "P28" = "#F2ACAE",
    "P29" = "#EF8E42",
    "P30" = "#8AA624",
    "P31" = "#D39FC8",
    "P32" = "#E5D9E3",
    "P33" = "#F2BF7D",
    "P34" = "#E9655B",
    "P35" = "#7DD076",
    "P36" = "#5D7AB1",
    "P37" = "#D2E9C7",
    "P38" = "#9DD1C8",
    "P39" = "#BCDD7A",
    "P40" = "#EA8675"
  ))+
  theme_dr() + 
  labs(x="UMAP1",y="UMAP2")+
  theme(panel.grid = element_blank(),
        axis.title =element_text(size = 15, face = 'bold'), 
        legend.key.size = unit(1.4,'line'),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15, face = 'bold'),
        legend.position = 'right')+
  guides(colour = guide_legend(override.aes = list(size=5, alpha = 1)))







## Proportion of major celltypes over all immune cells in tumor
### Calculating abundance matrix
df$major_celltype <- ifelse((df$major_celltype == "T&ILC cell") & (df$middle_celltype %in% c("CD4+T","CD8+T","Tprf")),"T cell",
                            ifelse(df$minor_celltype %in% c("c15_NK_CXXC5","c17_NK_KLRC1","c18_NK_NRGN"),"NK cell",
                            ifelse(df$minor_celltype == "c16_NK/NKT_LAG3","NK/NKT cell",df$major_celltype)))

imm_cls = c('T cell','B cell','Myeloid cell',"NK cell","NK/NKT cell")
meta = df %>% filter(anatomic_site == 'tumor')
anno_df = data.frame()
global_df = meta %>% 
  filter(major_celltype %in% imm_cls) %>%
  group_by(anatomic_site,patient, treatment_stage) %>% 
  summarise(n = n())

sub_df = meta %>% 
  filter(major_celltype %in% imm_cls) %>%
  group_by(major_celltype, anatomic_site, patient, treatment_stage) %>% 
  summarise(n1 = n())

anno_df= meta %>% 
  filter(major_celltype %in% imm_cls) %>%
  tidyr::expand(major_celltype, patient,treatment_stage, anatomic_site) %>%
  select(major_celltype, patient, treatment_stage, anatomic_site) %>% 
  left_join(sub_df, by=c('major_celltype', 'patient', 'treatment_stage','anatomic_site')) %>% 
  left_join(global_df, by = c('patient', 'treatment_stage','anatomic_site')) %>% 
  replace_na(list(n1=0,n=0)) %>%
  mutate(freq = n1/n) 

tmp = anno_df %>% 
  group_by(major_celltype, patient) %>% 
  summarise(freq_dnmc = freq[which(treatment_stage=='post')] - freq[which(treatment_stage=='pre')])
anno_df = anno_df %>% left_join(tmp, by = c('major_celltype','patient'))
# head(anno_df)

### Plotting
#### pre-treatment proportion
df_clin$patient_id <- as.character(df_clin$patient_id)
input = anno_df %>% 
  select(-c(n1,n)) %>% 
  distinct() %>% 
  filter(treatment_stage == 'pre') %>% 
  filter(major_celltype != "NK/NKT cell") %>%
  mutate(patient_id=sapply(stringr::str_split(patient, "_"), function(x) if (length(x)>=2) x[2] else NA)) %>%
  left_join(df_clin,by=c("patient_id")) %>%
  select(-patient_id) %>%
  filter(!is.na(MPR)) %>%#去掉未手术者
  mutate(Response=ifelse(MPR==1,"Well responders","Poor responders"))
input$major_celltype <- factor(input$major_celltype,levels=c("T cell","B cell","Myeloid cell","NK cell","NK/NKT cell"))
input$Response <- factor(input$Response,levels=c("Poor responders","Well responders"))
input %>%
  ggplot(aes(x=Response, y=freq)) +
  geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
  geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
  theme_pubr() + 
  facet_wrap(~major_celltype, ncol=2, scales = 'free_y') +
  labs(x="",y="Frequency at baseline (%CD45)") +
  scale_color_manual(values = c("Poor responders"='#5778A2',"Well responders"='#8A281E'), name='Response') +
  scale_x_discrete(labels = c("Poor responders" = "Poor\n responders\n(n=13)", "Well responders" = "Well\n responders\n(n=16)")) +
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

if(0){
##### T
p <- input %>%
  filter(major_celltype == "T cell") %>%
  ggplot(aes(x=Response, y=freq)) +
  geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
  geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
  theme_pubr() + 
  labs(x="",y="Frequency at baseline (%CD45)",title="T cell") +
  scale_y_continuous(limits=c(0,1)) +
  scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
  theme(strip.text.x = element_text(size = 15),
        strip.background =element_rect(fill=NA, color = 'black', size = 1), 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size=13),
        axis.title = element_text(size = 15),
        axis.line = element_line(size = .5),
        plot.title = element_text(hjust = 0.5,size = 15),
        legend.position = 'right',
        legend.key.size = unit(1.4,'line'),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15, face = 'bold'),
        plot.margin = margin(t=30, l=10)
        )
my_comparisons <- list(c("Poor responders","Well responders"))
p1 <- p+stat_compare_means(comparisons = my_comparisons)+
  theme() 

##### B
p <- input %>%
  filter(major_celltype == "B cell") %>%
  ggplot(aes(x=Response, y=freq)) +
  geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
  geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
  theme_pubr() + 
  labs(x="",y="Frequency at baseline (%CD45)",title="B cell") +
  scale_y_continuous(limits=c(0,1)) +
  scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
  theme(strip.text.x = element_text(size = 15),
        strip.background =element_rect(fill=NA, color = 'black', size = 1), 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size=13),
        axis.title = element_text(size = 15),
        axis.line = element_line(size = .5),
        plot.title = element_text(hjust = 0.5,size = 15),
        legend.position = 'right',
        legend.key.size = unit(1.4,'line'),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15, face = 'bold'),
        plot.margin = margin(t=30, l=10)
  )
my_comparisons <- list(c("Poor responders","Well responders"))
p2 <- p+stat_compare_means(comparisons = my_comparisons)+
  theme() 

##### Myeloid
p <- input %>%
  filter(major_celltype == "Myeloid cell") %>%
  ggplot(aes(x=Response, y=freq)) +
  geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
  geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
  theme_pubr() + 
  labs(x="",y="Frequency at baseline (%CD45)",title="Myeloid cell") +
  scale_y_continuous(limits=c(0,1)) +
  scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
  theme(strip.text.x = element_text(size = 15),
        strip.background =element_rect(fill=NA, color = 'black', size = 1), 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size=13),
        axis.title = element_text(size = 15),
        axis.line = element_line(size = .5),
        plot.title = element_text(hjust = 0.5,size = 15),
        legend.position = 'right',
        legend.key.size = unit(1.4,'line'),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15, face = 'bold'),
        plot.margin = margin(t=30, l=10)
  )
my_comparisons <- list(c("Poor responders","Well responders"))
p3 <- p+stat_compare_means(comparisons = my_comparisons)+
  theme() 


##### NK
p <- input %>%
  filter(major_celltype == "NK cell") %>%
  ggplot(aes(x=Response, y=freq)) +
  geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
  geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
  theme_pubr() + 
  labs(x="",y="Frequency at baseline (%CD45)",title="NK cell") +
  scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
  theme(strip.text.x = element_text(size = 15),
        strip.background =element_rect(fill=NA, color = 'black', size = 1), 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size=13),
        axis.title = element_text(size = 15),
        axis.line = element_line(size = .5),
        plot.title = element_text(hjust = 0.5,size = 15),
        legend.position = 'right',
        legend.key.size = unit(1.4,'line'),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15, face = 'bold'),
        plot.margin = margin(t=30, l=10)
  )
my_comparisons <- list(c("Poor responders","Well responders"))
p4 <- p+stat_compare_means(comparisons = my_comparisons)+
  theme() 
p_baseline_imm <- p1+p2+p3+p4

} #带p值画的代码

#### post-treatment proportion
df_clin$patient_id <- as.character(df_clin$patient_id)
input = anno_df %>% 
  select(-c(n1,n)) %>% 
  distinct() %>% 
  filter(treatment_stage == 'post') %>% 
  filter(major_celltype != "NK/NKT cell") %>%
  mutate(patient_id=sapply(stringr::str_split(patient, "_"), function(x) if (length(x)>=2) x[2] else NA)) %>%
  left_join(df_clin,by=c("patient_id")) %>%
  select(-patient_id) %>%
  filter(!is.na(MPR)) %>%#去掉未手术者
  mutate(Response=ifelse(MPR==1,"Well responders","Poor responders"))
input$major_celltype <- factor(input$major_celltype,levels=c("T cell","B cell","Myeloid cell","NK cell","NK/NKT cell"))
input$Response <- factor(input$Response,levels=c("Poor responders","Well responders"))
input %>%
  ggplot(aes(x=Response, y=freq)) +
  geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
  geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
  theme_pubr() + 
  facet_wrap(~major_celltype, ncol=2, scales = 'free_y') +
  labs(x="",y="Frequency at post-treatment (%CD45)") +
  scale_color_manual(values = c("Poor responders"='#5778A2',"Well responders"='#8A281E'), name='Response') +
  scale_x_discrete(labels = c("Poor responders" = "Poor\n responders\n(n=15)", "Well responders" = "Well\n responders\n(n=16)")) +
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

if(0){
  ##### T
  p <- input %>%
    filter(major_celltype == "T cell") %>%
    ggplot(aes(x=Response, y=freq)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency at post-treatment (%CD45)",title="T cell") +
    scale_y_continuous(limits=c(0,1)) +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p1 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  ##### B
  p <- input %>%
    filter(major_celltype == "B cell") %>%
    ggplot(aes(x=Response, y=freq)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency at post-treatment (%CD45)",title="B cell") +
    scale_y_continuous(limits=c(0,1)) +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p2 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  ##### Myeloid
  p <- input %>%
    filter(major_celltype == "Myeloid cell") %>%
    ggplot(aes(x=Response, y=freq)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency at post-treatment (%CD45)",title="Myeloid cell") +
    scale_y_continuous(limits=c(0,1)) +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p3 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  
  ##### NK
  p <- input %>%
    filter(major_celltype == "NK cell") %>%
    ggplot(aes(x=Response, y=freq)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency at post-treatment (%CD45)",title="NK cell") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p4 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  p_posttreatment_imm <- p1+p2+p3+p4
  p_posttreatment_imm
} #带p值画的代码


#### Proportion dynamics
df_clin$patient_id <- as.character(df_clin$patient_id)
input = anno_df %>% 
  select(-c(n1,n)) %>% 
  distinct() %>% 
  filter(treatment_stage == 'pre') %>% 
  filter(major_celltype != "NK/NKT cell") %>%
  mutate(patient_id=sapply(stringr::str_split(patient, "_"), function(x) if (length(x)>=2) x[2] else NA)) %>%
  left_join(df_clin,by=c("patient_id")) %>%
  select(-patient_id) %>%
  filter(!is.na(MPR)) %>%#去掉未手术者
  mutate(Response=ifelse(MPR==1,"Well responders","Poor responders"))
input$major_celltype <- factor(input$major_celltype,levels=c("T cell","B cell","Myeloid cell","NK cell","NK/NKT cell"))
input$Response <- factor(input$Response,levels=c("Poor responders","Well responders"))

input %>%
  ggplot(aes(x=Response, y=freq_dnmc)) +
  geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
  geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
  theme_pubr() + 
  facet_wrap(~major_celltype, ncol=2, scales = 'free_y') +
  labs(x="",y="Frequency change (%CD45)") +
  scale_color_manual(values = c("Poor responders"='#5778A2',"Well responders"='#8A281E'), name='Response') +
  scale_x_discrete(labels = c("Poor responders" = "Poor\n responders\n(n=12)", "Well responders" = "Well\n responders\n(n=13)")) +
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
if(0){
  ##### T
  p <- input %>%
    filter(major_celltype == "T cell") %>%
    ggplot(aes(x=Response, y=freq_dnmc)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency change (%CD45)",title="T cell") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p1 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  ##### B
  p <- input %>%
    filter(major_celltype == "B cell") %>%
    ggplot(aes(x=Response, y=freq_dnmc)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency change (%CD45)",title="B cell") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p2 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  ##### Myeloid
  p <- input %>%
    filter(major_celltype == "Myeloid cell") %>%
    ggplot(aes(x=Response, y=freq_dnmc)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency change (%CD45)",title="Myeloid cell") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p3 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  
  ##### NK
  p <- input %>%
    filter(major_celltype == "NK cell") %>%
    ggplot(aes(x=Response, y=freq_dnmc)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency change (%CD45)",title="NK cell") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p4 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  p_change_imm <- p1+p2+p3+p4
  p_change_imm
} #带p值画的代码

## Proportion of middle celltypes over all immune cells in tumor
imm_cls = c('CD4+T','CD8+T',"Tprf",'B cell','DC','Mono_Macro','Mast',"NK_ILC")
meta = df %>% filter(anatomic_site == 'tumor')
anno_df = data.frame()
global_df = meta %>% 
  filter(middle_celltype %in% imm_cls) %>%
  group_by(anatomic_site,patient, treatment_stage) %>% 
  summarise(n = n())

sub_df = meta %>% 
  filter(middle_celltype %in% imm_cls) %>%
  group_by(middle_celltype, anatomic_site, patient, treatment_stage) %>% 
  summarise(n1 = n())

anno_df= meta %>% 
  filter(middle_celltype %in% imm_cls) %>%
  tidyr::expand(middle_celltype, patient,treatment_stage, anatomic_site) %>%
  select(middle_celltype, patient, treatment_stage, anatomic_site) %>% 
  left_join(sub_df, by=c('middle_celltype', 'patient', 'treatment_stage','anatomic_site')) %>% 
  left_join(global_df, by = c('patient', 'treatment_stage','anatomic_site')) %>% 
  replace_na(list(n1=0,n=0)) %>%
  mutate(freq = n1/n) 

tmp = anno_df %>% 
  group_by(middle_celltype, patient) %>% 
  summarise(freq_dnmc = freq[which(treatment_stage=='post')] - freq[which(treatment_stage=='pre')])
anno_df = anno_df %>% left_join(tmp, by = c('middle_celltype','patient'))
# head(anno_df)

### Plotting
#### pre-treatment proportion
df_clin$patient_id <- as.character(df_clin$patient_id)
input = anno_df %>% 
  select(-c(n1,n)) %>% 
  distinct() %>% 
  filter(treatment_stage == 'pre') %>% 
  filter(middle_celltype %in% c("CD4+T","CD8+T","DC","Mono_Macro")) %>%
  mutate(patient_id=sapply(stringr::str_split(patient, "_"), function(x) if (length(x)>=2) x[2] else NA)) %>%
  left_join(df_clin,by=c("patient_id")) %>%
  select(-patient_id) %>%
  filter(!is.na(MPR)) %>%#去掉未手术者
  mutate(Response=ifelse(MPR==1,"Well responders","Poor responders"))
input$Response <- factor(input$Response,levels=c("Poor responders","Well responders"))
input %>%
  ggplot(aes(x=Response, y=freq)) +
  geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
  geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
  theme_pubr() + 
  facet_wrap(~middle_celltype, ncol=2, scales = 'free_y') +
  labs(x="",y="Frequency at baseline (%CD45)") +
  scale_color_manual(values = c("Poor responders"='#5778A2',"Well responders"='#8A281E'), name='Response') +
  scale_x_discrete(labels = c("Poor responders" = "Poor\n responders\n(n=13)", "Well responders" = "Well\n responders\n(n=16)")) +
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
if(0){
  ##### CD4+T
  p <- input %>%
    filter(middle_celltype == "CD4+T") %>%
    ggplot(aes(x=Response, y=freq)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency at baseline (%CD45)",title="CD4+T") +
    scale_y_continuous(limits=c(0,1)) +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p1 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  ##### CD8+T
  p <- input %>%
    filter(middle_celltype == "CD8+T") %>%
    ggplot(aes(x=Response, y=freq)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency at baseline (%CD45)",title="CD8+T") +
    scale_y_continuous(limits=c(0,1)) +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p2 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  ##### DC
  p <- input %>%
    filter(middle_celltype == "DC") %>%
    ggplot(aes(x=Response, y=freq)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency at baseline (%CD45)",title="DC") +
    scale_y_continuous(limits=c(0,1)) +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p3 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  
  ##### Mono/Macro
  p <- input %>%
    filter(middle_celltype == "Mono_Macro") %>%
    ggplot(aes(x=Response, y=freq)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency at baseline (%CD45)",title="Mono/Macro") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p4 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  p_baseline_imm <- p1+p2+p3+p4
  
} #带p值画的代码

#### post-treatment proportion
df_clin$patient_id <- as.character(df_clin$patient_id)
input = anno_df %>% 
  select(-c(n1,n)) %>% 
  distinct() %>% 
  filter(treatment_stage == 'post') %>% 
  filter(middle_celltype %in% c("CD4+T","CD8+T","DC","Mono_Macro")) %>%
  mutate(patient_id=sapply(stringr::str_split(patient, "_"), function(x) if (length(x)>=2) x[2] else NA)) %>%
  left_join(df_clin,by=c("patient_id")) %>%
  select(-patient_id) %>%
  filter(!is.na(MPR)) %>%#去掉未手术者
  mutate(Response=ifelse(MPR==1,"Well responders","Poor responders"))
input$Response <- factor(input$Response,levels=c("Poor responders","Well responders"))
input %>%
  ggplot(aes(x=Response, y=freq)) +
  geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
  geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
  theme_pubr() + 
  facet_wrap(~middle_celltype, ncol=2, scales = 'free_y') +
  labs(x="",y="Frequency at post-treatment (%CD45)") +
  scale_color_manual(values = c("Poor responders"='#5778A2',"Well responders"='#8A281E'), name='Response') +
  scale_x_discrete(labels = c("Poor responders" = "Poor\n responders\n(n=15)", "Well responders" = "Well\n responders\n(n=16)")) +
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
if(0){
  ##### CD4+T
  p <- input %>%
    filter(middle_celltype == "CD4+T") %>%
    ggplot(aes(x=Response, y=freq)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency at post-treatment (%CD45)",title="CD4+T") +
    scale_y_continuous(limits=c(0,1)) +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p1 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  ##### CD8+T
  p <- input %>%
    filter(middle_celltype == "CD8+T") %>%
    ggplot(aes(x=Response, y=freq)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency at post-treatment (%CD45)",title="CD8+T") +
    scale_y_continuous(limits=c(0,1)) +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p2 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  ##### DC
  p <- input %>%
    filter(middle_celltype == "DC") %>%
    ggplot(aes(x=Response, y=freq)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency at post-treatment (%CD45)",title="DC") +
    scale_y_continuous(limits=c(0,1)) +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p3 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  
  ##### Mono/Macro
  p <- input %>%
    filter(middle_celltype == "Mono_Macro") %>%
    ggplot(aes(x=Response, y=freq)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency at post-treatment (%CD45)",title="Mono/Macro") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p4 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  p_posttreatment_imm <- p1+p2+p3+p4
  
} #带p值画的代码

#### Proportion dynamics
df_clin$patient_id <- as.character(df_clin$patient_id)
input = anno_df %>% 
  select(-c(n1,n)) %>% 
  distinct() %>% 
  filter(treatment_stage == 'pre') %>% 
  filter(middle_celltype %in% c("CD4+T","CD8+T","DC","Mono_Macro")) %>%
  mutate(patient_id=sapply(stringr::str_split(patient, "_"), function(x) if (length(x)>=2) x[2] else NA)) %>%
  left_join(df_clin,by=c("patient_id")) %>%
  select(-patient_id) %>%
  filter(!is.na(MPR)) %>%#去掉未手术者
  mutate(Response=ifelse(MPR==1,"Well responders","Poor responders"))
input$Response <- factor(input$Response,levels=c("Poor responders","Well responders"))

input %>%
  ggplot(aes(x=Response, y=freq_dnmc)) +
  geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
  geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
  theme_pubr() + 
  facet_wrap(~middle_celltype, ncol=2, scales = 'free_y') +
  labs(x="",y="Frequency change (%CD45)") +
  scale_color_manual(values = c("Poor responders"='#5778A2',"Well responders"='#8A281E'), name='Response') +
  scale_x_discrete(labels = c("Poor responders" = "Poor\n responders\n(n=12)", "Well responders" = "Well\n responders\n(n=13)")) +
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
if(0){
  ##### CD4+T
  p <- input %>%
    filter(middle_celltype == "CD4+T") %>%
    ggplot(aes(x=Response, y=freq_dnmc)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency change (%CD45)",title="CD4+T") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p1 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  ##### CD8+T
  p <- input %>%
    filter(middle_celltype == "CD8+T") %>%
    ggplot(aes(x=Response, y=freq_dnmc)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency change (%CD45)",title="CD8+T") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p2 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  ##### DC
  p <- input %>%
    filter(middle_celltype == "DC") %>%
    ggplot(aes(x=Response, y=freq_dnmc)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency change (%CD45)",title="DC") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p3 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  
  ##### Mono/Macro
  p <- input %>%
    filter(middle_celltype == "Mono_Macro") %>%
    ggplot(aes(x=Response, y=freq_dnmc)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency change (%CD45)",title="Mono/Macro") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p4 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  p_change_imm <- p1+p2+p3+p4
  p_change_imm
} #带p值画的代码

## Proportion of minor celltypes over all immune cells in tumor
imm_cls = c('T cell','B cell','Myeloid cell',"NK cell",'NK/NKT cell')
meta = df %>% filter(anatomic_site == 'tumor')
anno_df = data.frame()
global_df = meta %>% 
  filter(major_celltype %in% imm_cls) %>%
  group_by(anatomic_site,patient, treatment_stage) %>% 
  summarise(n = n())
sub_df = meta %>% 
  filter(major_celltype %in% imm_cls) %>%
  group_by(minor_celltype, anatomic_site, patient, treatment_stage) %>% 
  summarise(n1 = n())
anno_df= meta %>% 
  filter(major_celltype %in% imm_cls) %>%
  tidyr::expand(minor_celltype,patient,treatment_stage, anatomic_site) %>%
  select(minor_celltype, patient, treatment_stage, anatomic_site) %>% 
  left_join(sub_df, by=c('minor_celltype', 'patient', 'treatment_stage','anatomic_site')) %>% 
  left_join(global_df, by = c('patient', 'treatment_stage','anatomic_site')) %>% 
  replace_na(list(n1=0,n=0,n2=0)) %>%
  mutate(freq = n1/n) 
tmp = anno_df %>% 
  group_by(minor_celltype, patient) %>% 
  summarise(freq_dnmc = freq[which(treatment_stage=='post')] - freq[which(treatment_stage=='pre')])
anno_df = anno_df %>% left_join(tmp, by = c('minor_celltype','patient'))

#### pre-treatment proportion
df_clin$patient_id <- as.character(df_clin$patient_id)
input = anno_df %>% 
  select(-c(n1,n)) %>% 
  distinct() %>% 
  filter(treatment_stage == 'pre') %>% 
  mutate(patient_id=sapply(stringr::str_split(patient, "_"), function(x) if (length(x)>=2) x[2] else NA)) %>%
  left_join(df_clin,by=c("patient_id")) %>%
  select(-patient_id) %>%
  filter(!is.na(MPR)) %>%#去掉未手术者
  mutate(Response=ifelse(MPR==1,"Well responders","Poor responders"))
input$Response <- factor(input$Response,levels=c("Poor responders","Well responders"))
input %>%
  ggplot(aes(x=Response, y=freq)) +
  geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
  geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
  theme_pubr() + 
  facet_wrap(~minor_celltype, ncol=7, scales = 'free_y') +
  labs(x="",y="Frequency at baseline (%CD45)") +
  scale_color_manual(values = c("Poor responders"='#5778A2',"Well responders"='#8A281E'), name='Response') +
  scale_x_discrete(labels = c("Poor responders" = "Poor\n responders\n(n=13)", "Well responders" = "Well\n responders\n(n=16)")) +
  theme(strip.text.x = element_text(size = 15),
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
#### post-treatment proportion
df_clin$patient_id <- as.character(df_clin$patient_id)
input = anno_df %>% 
  select(-c(n1,n)) %>% 
  distinct() %>% 
  filter(treatment_stage == 'post') %>% 
  mutate(patient_id=sapply(stringr::str_split(patient, "_"), function(x) if (length(x)>=2) x[2] else NA)) %>%
  left_join(df_clin,by=c("patient_id")) %>%
  select(-patient_id) %>%
  filter(!is.na(MPR)) %>%#去掉未手术者
  mutate(Response=ifelse(MPR==1,"Well responders","Poor responders"))
input$Response <- factor(input$Response,levels=c("Poor responders","Well responders"))
input %>%
  ggplot(aes(x=Response, y=freq)) +
  geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
  geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
  theme_pubr() + 
  facet_wrap(~minor_celltype, ncol=7, scales = 'free_y') +
  labs(x="",y="Frequency at post-treatment (%CD45)") +
  scale_color_manual(values = c("Poor responders"='#5778A2',"Well responders"='#8A281E'), name='Response') +
  scale_x_discrete(labels = c("Poor responders" = "Poor\n responders\n(n=15)", "Well responders" = "Well\n responders\n(n=16)")) +
  theme(strip.text.x = element_text(size = 15),
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

#### Proportion dynamics
df_clin$patient_id <- as.character(df_clin$patient_id)
input = anno_df %>% 
  select(-c(n1,n)) %>% 
  distinct() %>% 
  filter(treatment_stage == 'pre') %>% 
  mutate(patient_id=sapply(stringr::str_split(patient, "_"), function(x) if (length(x)>=2) x[2] else NA)) %>%
  left_join(df_clin,by=c("patient_id")) %>%
  select(-patient_id) %>%
  filter(!is.na(MPR)) %>%#去掉未手术者
  mutate(Response=ifelse(MPR==1,"Well responders","Poor responders"))
input$Response <- factor(input$Response,levels=c("Poor responders","Well responders"))

input %>%
  ggplot(aes(x=Response, y=freq_dnmc)) +
  geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
  geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
  theme_pubr() + 
  facet_wrap(~minor_celltype, ncol=7, scales = 'free_y') +
  labs(x="",y="Frequency change (%CD45)") +
  scale_color_manual(values = c("Poor responders"='#5778A2',"Well responders"='#8A281E'), name='Response') +
  scale_x_discrete(labels = c("Poor responders" = "Poor\n responders\n(n=12)", "Well responders" = "Well\n responders\n(n=13)")) +
  theme(strip.text.x = element_text(size = 15),
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

## Proportion of T minor celltypes over CD4/CD8 cells in tumor
##### CD4T
imm_cls = c('CD4+T')
meta = df %>% filter(anatomic_site == 'tumor')
anno_df = data.frame()
global_df = meta %>% 
  filter(middle_celltype %in% imm_cls) %>%
  group_by(anatomic_site,patient, treatment_stage) %>% 
  summarise(n = n())
sub_df = meta %>% 
  filter(middle_celltype %in% imm_cls) %>%
  group_by(minor_celltype, anatomic_site, patient, treatment_stage) %>% 
  summarise(n1 = n())
anno_df= meta %>% 
  filter(middle_celltype %in% imm_cls) %>%
  tidyr::expand(minor_celltype,patient,treatment_stage, anatomic_site) %>%
  select(minor_celltype, patient, treatment_stage, anatomic_site) %>% 
  left_join(sub_df, by=c('minor_celltype', 'patient', 'treatment_stage','anatomic_site')) %>% 
  left_join(global_df, by = c('patient', 'treatment_stage','anatomic_site')) %>% 
  replace_na(list(n1=0,n=0)) %>%
  mutate(freq = n1/n) 
tmp = anno_df %>% 
  group_by(minor_celltype, patient) %>% 
  summarise(freq_dnmc = freq[which(treatment_stage=='post')] - freq[which(treatment_stage=='pre')])
anno_df = anno_df %>% left_join(tmp, by = c('minor_celltype','patient'))

#### pre-treatment proportion
df_clin$patient_id <- as.character(df_clin$patient_id)
input = anno_df %>% 
  select(-c(n1,n)) %>% 
  distinct() %>% 
  filter(treatment_stage == 'pre') %>% 
  mutate(patient_id=sapply(stringr::str_split(patient, "_"), function(x) if (length(x)>=2) x[2] else NA)) %>%
  left_join(df_clin,by=c("patient_id")) %>%
  select(-patient_id) %>%
  filter(!is.na(MPR)) %>%#去掉未手术者
  mutate(Response=ifelse(MPR==1,"Well responders","Poor responders"))
input$Response <- factor(input$Response,levels=c("Poor responders","Well responders"))
input %>%
  ggplot(aes(x=Response, y=freq)) +
  geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
  geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
  theme_pubr() + 
  facet_wrap(~minor_celltype, ncol=3, scales = 'free_y') +
  labs(x="",y="Frequency at baseline (%CD4+T)") +
  scale_color_manual(values = c("Poor responders"='#5778A2',"Well responders"='#8A281E'), name='Response') +
  scale_x_discrete(labels = c("Poor responders" = "Poor\n responders\n(n=13)", "Well responders" = "Well\n responders\n(n=16)")) +
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
if(0){
  ##### c01_CD4_Tn_TCF7
  p <- input %>%
    filter(minor_celltype == "c01_CD4_Tn_TCF7") %>%
    ggplot(aes(x=Response, y=freq)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency at baseline (%CD4+T)",title="c01_CD4_Tn_TCF7") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p1 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  ##### c02_CD4_Tcm_S1PR1
  p <- input %>%
    filter(minor_celltype == "c02_CD4_Tcm_S1PR1") %>%
    ggplot(aes(x=Response, y=freq)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency at baseline (%CD4+T)",title="c02_CD4_Tcm_S1PR1") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p2 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  ##### c03_CD4_Tcm_ANXA1
  p <- input %>%
    filter(minor_celltype == "c03_CD4_Tcm_ANXA1") %>%
    ggplot(aes(x=Response, y=freq)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency at baseline (%CD4+T)",title="c03_CD4_Tcm_ANXA1") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p3 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  ##### c04_CD4_Tem_NR4A1
  p <- input %>%
    filter(minor_celltype == "c04_CD4_Tem_NR4A1") %>%
    ggplot(aes(x=Response, y=freq)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency at baseline (%CD4+T)",title="c04_CD4_Tem_NR4A1") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p4 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  ##### c05_CD4_Treg_FOXP3
  p <- input %>%
    filter(minor_celltype == "c05_CD4_Treg_FOXP3") %>%
    ggplot(aes(x=Response, y=freq)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency at baseline (%CD4+T)",title="c05_CD4_Treg_FOXP3") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p5 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  ##### c06_CD4_TOX
  p <- input %>%
    filter(minor_celltype == "c06_CD4_TOX") %>%
    ggplot(aes(x=Response, y=freq)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency at baseline (%CD4+T)",title="c06_CD4_TOX") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p6 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  p_pre_imm <- (p1+p2+p3)/(p4+p5+p6)
  p_pre_imm
} #带p值画的代码


#### post-treatment proportion
df_clin$patient_id <- as.character(df_clin$patient_id)
input = anno_df %>% 
  select(-c(n1,n)) %>% 
  distinct() %>% 
  filter(treatment_stage == 'post') %>% 
  mutate(patient_id=sapply(stringr::str_split(patient, "_"), function(x) if (length(x)>=2) x[2] else NA)) %>%
  left_join(df_clin,by=c("patient_id")) %>%
  select(-patient_id) %>%
  filter(!is.na(MPR)) %>%#去掉未手术者
  mutate(Response=ifelse(MPR==1,"Well responders","Poor responders"))
input$Response <- factor(input$Response,levels=c("Poor responders","Well responders"))
input %>%
  ggplot(aes(x=Response, y=freq)) +
  geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
  geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
  theme_pubr() + 
  facet_wrap(~minor_celltype, ncol=3, scales = 'free_y') +
  labs(x="",y="Frequency at post-treatment (%CD4+T)") +
  scale_color_manual(values = c("Poor responders"='#5778A2',"Well responders"='#8A281E'), name='Response') +
  scale_x_discrete(labels = c("Poor responders" = "Poor\n responders\n(n=15)", "Well responders" = "Well\n responders\n(n=16)")) +
  theme(strip.text.x = element_text(size = 15),
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
if(0){
  ##### c01_CD4_Tn_TCF7
  p <- input %>%
    filter(minor_celltype == "c01_CD4_Tn_TCF7") %>%
    ggplot(aes(x=Response, y=freq)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency change (%CD4+T)",title="c01_CD4_Tn_TCF7") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p1 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  ##### c02_CD4_Tcm_S1PR1
  p <- input %>%
    filter(minor_celltype == "c02_CD4_Tcm_S1PR1") %>%
    ggplot(aes(x=Response, y=freq)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency change (%CD4+T)",title="c02_CD4_Tcm_S1PR1") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p2 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  ##### c03_CD4_Tcm_ANXA1
  p <- input %>%
    filter(minor_celltype == "c03_CD4_Tcm_ANXA1") %>%
    ggplot(aes(x=Response, y=freq)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency change (%CD4+T)",title="c03_CD4_Tcm_ANXA1") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p3 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  ##### c04_CD4_Tem_NR4A1
  p <- input %>%
    filter(minor_celltype == "c04_CD4_Tem_NR4A1") %>%
    ggplot(aes(x=Response, y=freq)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency change (%CD4+T)",title="c04_CD4_Tem_NR4A1") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p4 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  ##### c05_CD4_Treg_FOXP3
  p <- input %>%
    filter(minor_celltype == "c05_CD4_Treg_FOXP3") %>%
    ggplot(aes(x=Response, y=freq)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency change (%CD4+T)",title="c05_CD4_Treg_FOXP3") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p5 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  ##### c06_CD4_TOX
  p <- input %>%
    filter(minor_celltype == "c06_CD4_TOX") %>%
    ggplot(aes(x=Response, y=freq)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency change (%CD4+T)",title="c06_CD4_TOX") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p6 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  p_post_imm <- (p1+p2+p3)/(p4+p5+p6)
  p_post_imm
} #带p值画的代码


#### Proportion dynamics
df_clin$patient_id <- as.character(df_clin$patient_id)
input = anno_df %>% 
  select(-c(n1,n)) %>% 
  distinct() %>% 
  filter(treatment_stage == 'pre') %>% 
  mutate(patient_id=sapply(stringr::str_split(patient, "_"), function(x) if (length(x)>=2) x[2] else NA)) %>%
  left_join(df_clin,by=c("patient_id")) %>%
  select(-patient_id) %>%
  filter(!is.na(MPR)) %>%#去掉未手术者
  mutate(Response=ifelse(MPR==1,"Well responders","Poor responders"))
input$Response <- factor(input$Response,levels=c("Poor responders","Well responders"))

input %>%
  ggplot(aes(x=Response, y=freq_dnmc)) +
  geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
  geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
  theme_pubr() + 
  facet_wrap(~minor_celltype, ncol=3, scales = 'free_y') +
  labs(x="",y="Frequency change (%CD4+T)") +
  scale_color_manual(values = c("Poor responders"='#5778A2',"Well responders"='#8A281E'), name='Response') +
  scale_x_discrete(labels = c("Poor responders" = "Poor\n responders\n(n=12)", "Well responders" = "Well\n responders\n(n=13)")) +
  theme(strip.text.x = element_text(size = 15),
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

if(0){
  ##### c01_CD4_Tn_TCF7
  p <- input %>%
    filter(minor_celltype == "c01_CD4_Tn_TCF7") %>%
    ggplot(aes(x=Response, y=freq_dnmc)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="freq_dnmcuency change (%CD4+T)",title="c01_CD4_Tn_TCF7") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p1 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  ##### c02_CD4_Tcm_S1PR1
  p <- input %>%
    filter(minor_celltype == "c02_CD4_Tcm_S1PR1") %>%
    ggplot(aes(x=Response, y=freq_dnmc)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Freqcuency change (%CD4+T)",title="c02_CD4_Tcm_S1PR1") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p2 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  ##### c03_CD4_Tcm_ANXA1
  p <- input %>%
    filter(minor_celltype == "c03_CD4_Tcm_ANXA1") %>%
    ggplot(aes(x=Response, y=freq_dnmc)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency change (%CD4+T)",title="c03_CD4_Tcm_ANXA1") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p3 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  ##### c04_CD4_Tem_NR4A1
  p <- input %>%
    filter(minor_celltype == "c04_CD4_Tem_NR4A1") %>%
    ggplot(aes(x=Response, y=freq_dnmc)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency change (%CD4+T)",title="c04_CD4_Tem_NR4A1") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p4 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  ##### c05_CD4_Treg_FOXP3
  p <- input %>%
    filter(minor_celltype == "c05_CD4_Treg_FOXP3") %>%
    ggplot(aes(x=Response, y=freq_dnmc)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency change (%CD4+T)",title="c05_CD4_Treg_FOXP3") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p5 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  ##### c06_CD4_TOX
  p <- input %>%
    filter(minor_celltype == "c06_CD4_TOX") %>%
    ggplot(aes(x=Response, y=freq_dnmc)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency change (%CD4+T)",title="c06_CD4_TOX") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p6 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  p_change_imm <- (p1+p2+p3)/(p4+p5+p6)
  p_change_imm
} #带p值画的代码

##### CD8T
imm_cls = c('CD8+T')
meta = df %>% filter(anatomic_site == 'tumor')
anno_df = data.frame()
global_df = meta %>% 
  filter(middle_celltype %in% imm_cls) %>%
  group_by(anatomic_site,patient, treatment_stage) %>% 
  summarise(n = n())
sub_df = meta %>% 
  filter(middle_celltype %in% imm_cls) %>%
  group_by(minor_celltype, anatomic_site, patient, treatment_stage) %>% 
  summarise(n1 = n())
anno_df= meta %>% 
  filter(middle_celltype %in% imm_cls) %>%
  tidyr::expand(minor_celltype,patient,treatment_stage, anatomic_site) %>%
  select(minor_celltype, patient, treatment_stage, anatomic_site) %>% 
  left_join(sub_df, by=c('minor_celltype', 'patient', 'treatment_stage','anatomic_site')) %>% 
  left_join(global_df, by = c('patient', 'treatment_stage','anatomic_site')) %>% 
  replace_na(list(n1=0,n=0)) %>%
  mutate(freq = n1/n) 
tmp = anno_df %>% 
  group_by(minor_celltype, patient) %>% 
  summarise(freq_dnmc = freq[which(treatment_stage=='post')] - freq[which(treatment_stage=='pre')])
anno_df = anno_df %>% left_join(tmp, by = c('minor_celltype','patient'))

#### pre-treatment proportion
df_clin$patient_id <- as.character(df_clin$patient_id)
input = anno_df %>% 
  select(-c(n1,n)) %>% 
  distinct() %>% 
  filter(treatment_stage == 'pre') %>% 
  mutate(patient_id=sapply(stringr::str_split(patient, "_"), function(x) if (length(x)>=2) x[2] else NA)) %>%
  left_join(df_clin,by=c("patient_id")) %>%
  select(-patient_id) %>%
  filter(!is.na(MPR)) %>%#去掉未手术者
  mutate(Response=ifelse(MPR==1,"Well responders","Poor responders"))
input$Response <- factor(input$Response,levels=c("Poor responders","Well responders"))
input %>%
  ggplot(aes(x=Response, y=freq)) +
  geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
  geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
  theme_pubr() + 
  facet_wrap(~minor_celltype, ncol=3, scales = 'free_y') +
  labs(x="",y="Frequency at baseline (%CD8+T)") +
  scale_color_manual(values = c("Poor responders"='#5778A2',"Well responders"='#8A281E'), name='Response') +
  scale_x_discrete(labels = c("Poor responders" = "Poor\n responders\n(n=13)", "Well responders" = "Well\n responders\n(n=16)")) +
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
if(0){
  ##### c07_CD8_Tem/Teff_TNFSF9
  p <- input %>%
    filter(minor_celltype == "c07_CD8_Tem/Teff_TNFSF9") %>%
    ggplot(aes(x=Response, y=freq)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency at baseline (%CD8+T)",title="c07_CD8_Tem/Teff_TNFSF9") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p1 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  ##### c08_CD8_Tem/Teff_HSPA1B
  p <- input %>%
    filter(minor_celltype == "c08_CD8_Tem/Teff_HSPA1B") %>%
    ggplot(aes(x=Response, y=freq)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency at baseline (%CD8+T)",title="c08_CD8_Tem/Teff_HSPA1B") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p2 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  ##### c09_CD8_Tem/Teff_GNLY
  p <- input %>%
    filter(minor_celltype == "c09_CD8_Tem/Teff_GNLY") %>%
    ggplot(aes(x=Response, y=freq)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency at baseline (%CD8+T)",title="c09_CD8_Tem/Teff_GNLY") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p3 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  ##### c10_CD8_Teff_GZMK
  p <- input %>%
    filter(minor_celltype == "c10_CD8_Teff_GZMK") %>%
    ggplot(aes(x=Response, y=freq)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency at baseline (%CD8+T)",title="c10_CD8_Teff_GZMK") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p4 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  ##### c11_CD8_Tex_CXCL13
  p <- input %>%
    filter(minor_celltype == "c11_CD8_Tex_CXCL13") %>%
    ggplot(aes(x=Response, y=freq)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency at baseline (%CD8+T)",title="c11_CD8_Tex_CXCL13") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p5 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  ##### c12_CD8_Temra_FGFBP2
  p <- input %>%
    filter(minor_celltype == "c12_CD8_Temra_FGFBP2") %>%
    ggplot(aes(x=Response, y=freq)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency at baseline (%CD8+T)",title="c12_CD8_Temra_FGFBP2") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p6 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  ##### c13_CD8_IL7R
  p <- input %>%
    filter(minor_celltype == "c13_CD8_IL7R") %>%
    ggplot(aes(x=Response, y=freq)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency at baseline (%CD8+T)",title="c13_CD8_IL7R") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p7 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  p_pre_imm <- (p1+p2+p3)/(p4+p5+p6)/(p7+p7+p7)
  p_pre_imm
} #带p值画的代码


#### post-treatment proportion
df_clin$patient_id <- as.character(df_clin$patient_id)
input = anno_df %>% 
  select(-c(n1,n)) %>% 
  distinct() %>% 
  filter(treatment_stage == 'post') %>% 
  mutate(patient_id=sapply(stringr::str_split(patient, "_"), function(x) if (length(x)>=2) x[2] else NA)) %>%
  left_join(df_clin,by=c("patient_id")) %>%
  select(-patient_id) %>%
  filter(!is.na(MPR)) %>%#去掉未手术者
  mutate(Response=ifelse(MPR==1,"Well responders","Poor responders"))
input$Response <- factor(input$Response,levels=c("Poor responders","Well responders"))
input %>%
  ggplot(aes(x=Response, y=freq)) +
  geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
  geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
  theme_pubr() + 
  facet_wrap(~minor_celltype, ncol=3, scales = 'free_y') +
  labs(x="",y="Frequency at post-treatment (%CD8+T)") +
  scale_color_manual(values = c("Poor responders"='#5778A2',"Well responders"='#8A281E'), name='Response') +
  scale_x_discrete(labels = c("Poor responders" = "Poor\n responders\n(n=15)", "Well responders" = "Well\n responders\n(n=16)")) +
  theme(strip.text.x = element_text(size = 15),
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
if(0){
  ##### c07_CD8_Tem/Teff_TNFSF9
  p <- input %>%
    filter(minor_celltype == "c07_CD8_Tem/Teff_TNFSF9") %>%
    ggplot(aes(x=Response, y=freq)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency at post-treatment (%CD8+T)",title="c07_CD8_Tem/Teff_TNFSF9") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p1 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  ##### c08_CD8_Tem/Teff_HSPA1B
  p <- input %>%
    filter(minor_celltype == "c08_CD8_Tem/Teff_HSPA1B") %>%
    ggplot(aes(x=Response, y=freq)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency at post-treatment (%CD8+T)",title="c08_CD8_Tem/Teff_HSPA1B") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p2 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  ##### c09_CD8_Tem/Teff_GNLY
  p <- input %>%
    filter(minor_celltype == "c09_CD8_Tem/Teff_GNLY") %>%
    ggplot(aes(x=Response, y=freq)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency at post-treatment (%CD8+T)",title="c09_CD8_Tem/Teff_GNLY") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p3 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  ##### c10_CD8_Teff_GZMK
  p <- input %>%
    filter(minor_celltype == "c10_CD8_Teff_GZMK") %>%
    ggplot(aes(x=Response, y=freq)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency at post-treatment (%CD8+T)",title="c10_CD8_Teff_GZMK") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p4 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  ##### c11_CD8_Tex_CXCL13
  p <- input %>%
    filter(minor_celltype == "c11_CD8_Tex_CXCL13") %>%
    ggplot(aes(x=Response, y=freq)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency at post-treatment (%CD8+T)",title="c11_CD8_Tex_CXCL13") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p5 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  ##### c12_CD8_Temra_FGFBP2
  p <- input %>%
    filter(minor_celltype == "c12_CD8_Temra_FGFBP2") %>%
    ggplot(aes(x=Response, y=freq)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency at post-treatment (%CD8+T)",title="c12_CD8_Temra_FGFBP2") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p6 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  ##### c13_CD8_IL7R
  p <- input %>%
    filter(minor_celltype == "c13_CD8_IL7R") %>%
    ggplot(aes(x=Response, y=freq)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency at post-treatment (%CD8+T)",title="c13_CD8_IL7R") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p7 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  p_post_imm <- (p1+p2+p3)/(p4+p5+p6)/(p7+p7+p7)
  p_post_imm
} #带p值画的代码


#### Proportion dynamics
df_clin$patient_id <- as.character(df_clin$patient_id)
input = anno_df %>% 
  select(-c(n1,n)) %>% 
  distinct() %>% 
  filter(treatment_stage == 'pre') %>% 
  mutate(patient_id=sapply(stringr::str_split(patient, "_"), function(x) if (length(x)>=2) x[2] else NA)) %>%
  left_join(df_clin,by=c("patient_id")) %>%
  select(-patient_id) %>%
  filter(!is.na(MPR)) %>%#去掉未手术者
  mutate(Response=ifelse(MPR==1,"Well responders","Poor responders"))
input$Response <- factor(input$Response,levels=c("Poor responders","Well responders"))

input %>%
  ggplot(aes(x=Response, y=freq_dnmc)) +
  geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
  geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
  theme_pubr() + 
  facet_wrap(~minor_celltype, ncol=3, scales = 'free_y') +
  labs(x="",y="Frequency change (%CD8+T)") +
  scale_color_manual(values = c("Poor responders"='#5778A2',"Well responders"='#8A281E'), name='Response') +
  scale_x_discrete(labels = c("Poor responders" = "Poor\n responders\n(n=12)", "Well responders" = "Well\n responders\n(n=13)")) +
  theme(strip.text.x = element_text(size = 15),
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

if(0){
  ##### c07_CD8_Tem/Teff_TNFSF9
  p <- input %>%
    filter(minor_celltype == "c07_CD8_Tem/Teff_TNFSF9") %>%
    ggplot(aes(x=Response, y=freq_dnmc)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency change (%CD8+T)",title="c07_CD8_Tem/Teff_TNFSF9") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p1 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  ##### c08_CD8_Tem/Teff_HSPA1B
  p <- input %>%
    filter(minor_celltype == "c08_CD8_Tem/Teff_HSPA1B") %>%
    ggplot(aes(x=Response, y=freq_dnmc)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency change (%CD8+T)",title="c08_CD8_Tem/Teff_HSPA1B") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p2 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  ##### c09_CD8_Tem/Teff_GNLY
  p <- input %>%
    filter(minor_celltype == "c09_CD8_Tem/Teff_GNLY") %>%
    ggplot(aes(x=Response, y=freq_dnmc)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency change (%CD8+T)",title="c09_CD8_Tem/Teff_GNLY") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p3 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  ##### c10_CD8_Teff_GZMK
  p <- input %>%
    filter(minor_celltype == "c10_CD8_Teff_GZMK") %>%
    ggplot(aes(x=Response, y=freq_dnmc)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency change (%CD8+T)",title="c10_CD8_Teff_GZMK") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p4 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  ##### c11_CD8_Tex_CXCL13
  p <- input %>%
    filter(minor_celltype == "c11_CD8_Tex_CXCL13") %>%
    ggplot(aes(x=Response, y=freq_dnmc)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency change (%CD8+T)",title="c11_CD8_Tex_CXCL13") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p5 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  ##### c12_CD8_Temra_FGFBP2
  p <- input %>%
    filter(minor_celltype == "c12_CD8_Temra_FGFBP2") %>%
    ggplot(aes(x=Response, y=freq_dnmc)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency change (%CD8+T)",title="c12_CD8_Temra_FGFBP2") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p6 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  ##### c13_CD8_IL7R
  p <- input %>%
    filter(minor_celltype == "c13_CD8_IL7R") %>%
    ggplot(aes(x=Response, y=freq_dnmc)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    labs(x="",y="Frequency change (%CD8+T)",title="c13_CD8_IL7R") +
    scale_color_manual(values = c('#5778A2','#8A281E'), name='Response') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)
    )
  my_comparisons <- list(c("Poor responders","Well responders"))
  p7 <- p+stat_compare_means(comparisons = my_comparisons)+
    theme() 
  
  p_change_imm <- (p1+p2+p3)/(p4+p5+p6)/(p7+p7+p7)
  p_change_imm
} #带p值画的代码


## Proportion of major cell types in tumor samples
if(0){
library(ComplexHeatmap)

##### matrix
mat1 <- df_clin %>%
  left_join(df,by=c("patient_id")) %>%
  filter(anatomic_site=="tumor") %>%
  select(pCR,MPR,patient_id,patient_number,treatment_stage,anatomic_site) %>%
  distinct() 
mat1$treatment_stage <- factor(mat1$treatment_stage,levels=c("pre","on","post"))
mat1 <- mat1 %>%
  arrange(desc(MPR),patient_number,desc(anatomic_site),treatment_stage,) %>%
  mutate(pCR=ifelse(pCR==0,"Non-pCR","pCR"),MPR=ifelse(MPR==0,"Non-MPR","MPR"))
rownames(mat1) <- paste(mat1$patient_number,mat1$patient_id,mat1$treatment_stage,mat1$anatomic_site,sep="_")
#split level determines show order
col_split <- factor(mat1$patient_number,levels=unique(mat1$patient_number))
row_split <- factor(c("pCR","MPR","anatomic_site","treatment_stage"),levels=c("MPR","pCR","anatomic_site","treatment_stage"))
mat1 <- mat1[,c(1,2,6,5)]
mat1 <- t(mat1)
col <- c(
  "Non-pCR"="#FFEBD9",
  "pCR"="#FCB07E",
  "Non-MPR"="#5778A2",
  "MPR"="#8A281E",
  "pre"="#E7872B",
  "on"="#925A44",
  "post"="#C03830",
  "tumor"="#65B7B3",
  "blood"="#EEB4D6"
)
##### legend
lgd_MPR = Legend(
                 nrow=1,
              title = "Resonpse-MPR", 
              title_position = "topleft",
              labels = c("MPR",
                         "Non-MPR",
                         "NA"),
              grid_height = unit(0.4, "cm"), # 区别于前面的legend_height
              grid_width = unit(4, "mm"),
              border = "white",
              legend_gp = gpar(fill = c("#8A281E",
                                        "#5778A2",
                                        "grey"
                                        )))
lgd_pCR = Legend(nrow=1,
                 title = "Response-pCR", 
                 title_position = "topleft",
                 labels = c("pCR",
                            "Non-pCR",
                            "NA"),
                 grid_height = unit(0.4, "cm"), # 区别于前面的legend_height
                 grid_width = unit(4, "mm"),
                 border = "white",
                 legend_gp = gpar(fill = c("#FCB07E",
                                           "#FFEBD9",
                                           "grey"
                 )))
lgd_treatmentstage = Legend(nrow = 1,
                 title = "Sample Stage", 
                 title_position = "topleft",
                 labels = c("Pre",
                            "On",
                            "Post"),
                 grid_height = unit(0.4, "cm"), # 区别于前面的legend_height
                 grid_width = unit(4, "mm"),
                 border = "white",
                 legend_gp = gpar(fill = c("#E7872B",
                                           "#925A44",
                                           "#C03830"
                 )))
lgd_tissue = Legend(nrow = 1,
                            title = "Tissue", 
                            title_position = "topleft",
                            labels = c("Tumor",
                                       "Blood"),
                            grid_height = unit(0.4, "cm"), # 区别于前面的legend_height
                            grid_width = unit(4, "mm"),
                            border = "white",
                            legend_gp = gpar(fill = c("#65B7B3",
                                                      "#EEB4D6")))
lgd_majorcelltype = Legend(
  by_row = T,
  title = "", 
  title_position = "topleft",
  labels = c("B cell","Myeloid cell","T&NK cell","Epithelial cell","Endothelial cell","Fibrocyte"),
  grid_height = unit(0.4, "cm"), # 区别于前面的legend_height
  grid_width = unit(4, "mm"),
  border = "white",
  legend_gp = gpar(fill = c("#3A6FA4",
                            "#4E9669",
                            "#E0823A",
                            "#BA3931",
                            "#6F5198",
                            "#81564C")))
##### column
m <- df %>%
  mutate(sample=paste(patient_number,patient_id,treatment_stage,anatomic_site,sep="_")) %>%
  group_by(sample, major_celltype) %>%
  summarise(count = n(), .groups = "drop_last") %>%  
  mutate(ratio = count / sum(count)) %>%
  select(-count) %>%
  pivot_wider(
    id_cols = sample,
    names_from = major_celltype,
    values_from = ratio,
    values_fill = 0 
  )
m <- as.data.frame(m)
rownames(m) <- m$sample
m$sample <- NULL
m <- m[colnames(mat1),]
m <- m[,c(1,5,6,3,2,4)]
bottom_ha <- HeatmapAnnotation(foo = anno_barplot(m, gp = gpar(fill =c("#3A6FA4",
                                                                       "#4E9669",
                                                                       "#E0823A",
                                                                       "#BA3931",
                                                                       "#6F5198",
                                                                       "#81564C"
                                                                       )),
                                                  bar_width = 1,
                                                  height = unit(3, "cm")),
                               show_annotation_name = FALSE)
ht1 <- Heatmap(mat1,
               row_split=row_split,
               row_labels = c("Response-pCR","Response-MPR","Tissue","Sample Stage"),
               column_split = col_split,
               col=col,
               row_title = NULL,
               row_names_gp = gpar(fontsize =12),
               show_column_names = FALSE,
               column_title_gp = gpar(fontsize = 12),
               column_title_side = "bottom",
               column_title_rot = 90,
               width = ncol(mat1)*unit(3, "mm"),
               height = nrow(mat1)*unit(4, "mm"),
               row_gap = unit(2, "mm"),
               column_gap = unit(1.5, "mm"),
               show_heatmap_legend = FALSE,
               bottom_annotation = bottom_ha,
               na_col = "grey",
               )
pd = packLegend(lgd_MPR,lgd_pCR,lgd_treatmentstage,lgd_tissue,lgd_majorcelltype,
                direction = "vertical", # "vertical", "horizontal"
                max_height = unit(20, "cm"), 
                row_gap = unit(2, "mm") 
) 
draw(ht1,annotation_legend_list=pd)  
}

## Proportion of immune major cell types in each sample
library(ComplexHeatmap)

##### matrix
mat1 <- df_clin %>%
  left_join(df,by=c("patient_id")) %>%
  select(pCR,MPR,patient_id,patient_number,treatment_stage,anatomic_site) %>%
  distinct() 
mat1$treatment_stage <- factor(mat1$treatment_stage,levels=c("pre","on","post"))
mat1 <- mat1 %>%
  arrange(desc(MPR),patient_number,desc(anatomic_site),treatment_stage,) %>%
  mutate(pCR=ifelse(pCR==0,"Non-pCR","pCR"),MPR=ifelse(MPR==0,"Non-MPR","MPR"))
rownames(mat1) <- paste(mat1$patient_number,mat1$patient_id,mat1$treatment_stage,mat1$anatomic_site,sep="_")
#split level determines show order
col_split <- factor(mat1$patient_number,levels=unique(mat1$patient_number))
row_split <- factor(c("pCR","MPR","anatomic_site","treatment_stage"),levels=c("MPR","pCR","anatomic_site","treatment_stage"))
mat1 <- mat1[,c(1,2,6,5)]
mat1 <- t(mat1)
col <- c(
  "Non-pCR"="#FFEBD9",
  "pCR"="#FCB07E",
  "Non-MPR"="#5778A2",
  "MPR"="#8A281E",
  "pre"="#E7872B",
  "on"="#925A44",
  "post"="#C03830",
  "tumor"="#65B7B3",
  "blood"="#EEB4D6"
)
##### legend
lgd_MPR = Legend(
  nrow=1,
  title = "Resonpse-MPR", 
  title_position = "topleft",
  labels = c("MPR",
             "Non-MPR",
             "NA"),
  grid_height = unit(0.4, "cm"), # 区别于前面的legend_height
  grid_width = unit(4, "mm"),
  border = "white",
  legend_gp = gpar(fill = c("#8A281E",
                            "#5778A2",
                            "grey"
  )))
lgd_pCR = Legend(nrow=1,
                 title = "Response-pCR", 
                 title_position = "topleft",
                 labels = c("pCR",
                            "Non-pCR",
                            "NA"),
                 grid_height = unit(0.4, "cm"), # 区别于前面的legend_height
                 grid_width = unit(4, "mm"),
                 border = "white",
                 legend_gp = gpar(fill = c("#FCB07E",
                                           "#FFEBD9",
                                           "grey"
                 )))
lgd_treatmentstage = Legend(nrow = 1,
                            title = "Sample Stage", 
                            title_position = "topleft",
                            labels = c("Pre",
                                       "On",
                                       "Post"),
                            grid_height = unit(0.4, "cm"), # 区别于前面的legend_height
                            grid_width = unit(4, "mm"),
                            border = "white",
                            legend_gp = gpar(fill = c("#E7872B",
                                                      "#925A44",
                                                      "#C03830"
                            )))
lgd_tissue = Legend(nrow = 1,
                    title = "Tissue", 
                    title_position = "topleft",
                    labels = c("Tumor",
                               "Blood"),
                    grid_height = unit(0.4, "cm"), # 区别于前面的legend_height
                    grid_width = unit(4, "mm"),
                    border = "white",
                    legend_gp = gpar(fill = c("#65B7B3",
                                              "#EEB4D6")))
lgd_majorcelltype = Legend(
  by_row = T,
  title = "", 
  title_position = "topleft",
  labels = c("B cell","Myeloid cell","T&NK cell"),
  grid_height = unit(0.4, "cm"), # 区别于前面的legend_height
  grid_width = unit(4, "mm"),
  border = "white",
  legend_gp = gpar(fill = c("#3A6FA4",
                            "#4E9669",
                            "#E0823A")))
##### column
m <- df %>%
  filter(major_celltype %in% c("B cell","Myeloid cell","T&ILC cell")) %>%
  mutate(sample=paste(patient_number,patient_id,treatment_stage,anatomic_site,sep="_")) %>%
  group_by(sample, major_celltype) %>%
  summarise(count = n(), .groups = "drop_last") %>%  
  mutate(ratio = count / sum(count)) %>%
  select(-count) %>%
  pivot_wider(
    id_cols = sample,
    names_from = major_celltype,
    values_from = ratio,
    values_fill = 0 
  )
m <- as.data.frame(m)
rownames(m) <- m$sample
m$sample <- NULL
m <- m[colnames(mat1),]
m <- m[,c(3,2,1)]
bottom_ha <- HeatmapAnnotation(foo = anno_barplot(m, gp = gpar(fill =c("#E0823A",
                                                                       "#4E9669",
                                                                       "#3A6FA4"
                                                                       
)),
bar_width = 1,
height = unit(3, "cm")),
show_annotation_name = FALSE)
ht1 <- Heatmap(mat1,
               row_split=row_split,
               row_labels = c("Response-pCR","Response-MPR","Tissue","Sample Stage"),
               column_split = col_split,
               col=col,
               row_title = NULL,
               row_names_gp = gpar(fontsize =12),
               show_column_names = FALSE,
               column_title_gp = gpar(fontsize = 12),
               column_title_side = "bottom",
               column_title_rot = 90,
               width = ncol(mat1)*unit(3, "mm"),
               height = nrow(mat1)*unit(4, "mm"),
               row_gap = unit(2, "mm"),
               column_gap = unit(1.5, "mm"),
               show_heatmap_legend = FALSE,
               bottom_annotation = bottom_ha,
               na_col = "grey",
)
pd = packLegend(lgd_MPR,lgd_pCR,lgd_treatmentstage,lgd_tissue,lgd_majorcelltype,
                direction = "vertical", # "vertical", "horizontal"
                max_height = unit(20, "cm"), 
                row_gap = unit(2, "mm") 
) 
draw(ht1,annotation_legend_list=pd)  















