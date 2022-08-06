##### ���������ٴ����� -----
rm(list = ls())

# ��װ������Ҫ��R��
library(pacman)
p_load(data.table, tidyverse, magrittr)

# ���ù���Ŀ¼
projectPath = "/Users/prf/SCIrepro/PARPBP"
dataPath = paste(projectPath, "Data", sep = "/")
setwd(dataPath)

#### �����ٴ����ݡ�UCSC Xena��----
# �ٴ���Ϣ
# https://tcga.xenahubs.net/download/TCGA.LIHC.sampleMap/LIHC_clinicalMatrix
# ������Ϣ
# https://tcga.xenahubs.net/download/survival/LIHC_survival.txt.gz

#### ��ȡ�ٴ����� ----
# �ٴ���Ϣ
clinic_data = read.table("TCGA.UCEC.sampleMap_UCEC_clinicalMatrix", header = T, sep = "\t")
# ������Ϣ
surv_data = read.table("survival_UCEC_survival.txt", header = T, sep = "\t")
# �ϲ�
clinic_data = merge(clinic_data, surv_data, by= "sampleID")

#### ��ȡ������� ----
# FPKM-UQ����table��Ԥ��
exp_data_fpkm = fread("LIHC_Portal_RNA_FPKM-UQ.txt", header = T, sep = "\t") %>% 
    column_to_rownames("ensembl_gene_id") %>% 
    dplyr::select(ends_with("-01A")) %>%    #ɸѡ�����������еı��ˮƽ
    set_colnames(gsub("A$", "", colnames(.))) # Ҫ��ĩβ���ַ�ȥ���������ٴ��������Ǻ�

# Counts�����ٴ����ط������������
exp_data_counts = fread("LIHC_Portal_RNA_Counts.txt", h = T, sep = "\t", check.names = F) %>% 
    column_to_rownames("ensembl_gene_id") %>% 
    dplyr::select(ends_with("-01A")) %>%    #ɸѡ�����������еı��ˮƽ
    set_colnames(gsub("A$", "", colnames(.)))

#### �������� ----
# ��������һ��ֻҪ������ʱ��������Ҳ���ʵ����������������

# ���ˣ�����Ҫ����ٴ�����ɸѡ����Ӧ����������
clinic_data %>%   #����ʱ��������0������״̬���벻Ϊ��
    filter(OS.time > 0 & !is.na(OS) & str_detect(clinical_stage, "^Stage")) %>% 
    droplevels    #ȡ�����������
# ɸѡ���ٴ������ͱ�������й�ͬ��������
Common = intersect(as.character(as.matrix(clinic_data$sampleID)), colnames(exp_data_fpkm)) # exp_data_countsͬ��
length(Common) # [1] 339

# ���˺���ٴ����ݣ�ɸѡ����ͬ�������ٴ�����
clinic_data %<>% filter(sampleID %in% Common)
write.table(clinic_data, "clinic_data_clean.txt", row.names = F, sep = "\t", quote = F)

# ���˺�ı�����ݣ�ֻ�������������ݣ���ɸѡ����ͬ�����ı������
exp_data_fpkm = exp_data_fpkm["ENSG00000185480",Common]   #ɸѡ��Ŀ������ڹ��������ı��ֵ
write.table(data.frame(Symbol = "PARPBP", exp_data_fpkm, check.names = F), "Target_fpkm-uq_tumor.txt", row.names = F, sep = "\t", quote = F)
exp_data_counts = exp_data_counts[,Common]    #ɸѡ����ͬ�����еı��ֵ
fwrite(data.frame(Symbol = rownames(exp_data_counts), exp_data_counts, check.names = F), "All_counts_tumor.txt", row.names = F, sep = "\t", quote = F)
