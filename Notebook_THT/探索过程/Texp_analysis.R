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




