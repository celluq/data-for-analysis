##### 下载整合临床数据 -----
rm(list = ls())

# 安装加载需要的R包
library(pacman)
p_load(data.table, tidyverse, magrittr)

# 设置工作目录
projectPath = "/Users/prf/SCIrepro/PARPBP"
dataPath = paste(projectPath, "Data", sep = "/")
setwd(dataPath)

#### 下载临床数据【UCSC Xena】----
# 临床信息
# https://tcga.xenahubs.net/download/TCGA.LIHC.sampleMap/LIHC_clinicalMatrix
# 生存信息
# https://tcga.xenahubs.net/download/survival/LIHC_survival.txt.gz

#### 读取临床数据 ----
# 临床信息
clinic_data = read.table("TCGA.UCEC.sampleMap_UCEC_clinicalMatrix", header = T, sep = "\t")
# 生存信息
surv_data = read.table("survival_UCEC_survival.txt", header = T, sep = "\t")
# 合并
clinic_data = merge(clinic_data, surv_data, by= "sampleID")

#### 读取表达数据 ----
# FPKM-UQ用于table和预后
exp_data_fpkm = fread("LIHC_Portal_RNA_FPKM-UQ.txt", header = T, sep = "\t") %>% 
    column_to_rownames("ensembl_gene_id") %>% 
    dplyr::select(ends_with("-01A")) %>%    #筛选出肿瘤样本中的表达水平
    set_colnames(gsub("A$", "", colnames(.))) # 要把末尾的字符去掉才能与临床样本向吻合

# Counts用于临床因素分组间差异表达分析
exp_data_counts = fread("LIHC_Portal_RNA_Counts.txt", h = T, sep = "\t", check.names = F) %>% 
    column_to_rownames("ensembl_gene_id") %>% 
    dplyr::select(ends_with("-01A")) %>%    #筛选出肿瘤样本中的表达水平
    set_colnames(gsub("A$", "", colnames(.)))

#### 样本过滤 ----
# 样本过滤一般只要求生存时间完整，也可适当添加其他限制因素

# 过滤，根据要求的临床条件筛选出相应的肿瘤样本
clinic_data %>%   #生存时间必须大于0且生存状态必须不为空
    filter(OS.time > 0 & !is.na(OS) & str_detect(clinical_stage, "^Stage")) %>% 
    droplevels    #取出多余的因子
# 筛选出临床样本和表达样本中共同的样本数
Common = intersect(as.character(as.matrix(clinic_data$sampleID)), colnames(exp_data_fpkm)) # exp_data_counts同样
length(Common) # [1] 339

# 过滤后的临床数据，筛选出共同样本的临床数据
clinic_data %<>% filter(sampleID %in% Common)
write.table(clinic_data, "clinic_data_clean.txt", row.names = F, sep = "\t", quote = F)

# 过滤后的表达数据（只有肿瘤样本数据），筛选出共同样本的表达数据
exp_data_fpkm = exp_data_fpkm["ENSG00000185480",Common]   #筛选出目标基因在共有样本的表达值
write.table(data.frame(Symbol = "PARPBP", exp_data_fpkm, check.names = F), "Target_fpkm-uq_tumor.txt", row.names = F, sep = "\t", quote = F)
exp_data_counts = exp_data_counts[,Common]    #筛选出共同样本中的表达值
fwrite(data.frame(Symbol = rownames(exp_data_counts), exp_data_counts, check.names = F), "All_counts_tumor.txt", row.names = F, sep = "\t", quote = F)
