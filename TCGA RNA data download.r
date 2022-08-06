##### RNA��������������-----
rm(list = ls())

# ��װ������Ҫ��R��
install.packages("pacman", "TCGAbiolinks","tidyverse","magrittr","data.table","biomaRt")
library(pacman)
p_load(TCGAbiolinks, tidyverse, magrittr, data.table, biomaRt)

# ���ù���Ŀ¼
setwd("D:\\R\\��ϰ")

#### ʹ��TCGAbiolinks����Counts��PKM-UQ����----
#https://www.bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html
#��һ�� ����count����
query_count = GDCquery(project = "TCGA-UCEC", legacy = FALSE, experimental.strategy = "RNA-Seq", data.category = "Transcriptome Profiling", data.type = "Gene Expression Quantification", workflow.type = "HTSeq - Counts")
GDCdownload(query_count, files.per.chunk = 100)   # ��������
# ����Ԥ����
dataAssy_count = GDCprepare(query_count, summarizedExperiment = F)
# ����������Ա�����
filename_count = paste0("TCGA-UCEC.RNA.query.","Counts", ".RData")
save(query_count, dataAssy_count, file = filename_count)
load(filename_count)

#�ڶ��� ����FKPM����
query_fpkm = GDCquery(project = "TCGA-UCEC", legacy = FALSE, experimental.strategy = "RNA-Seq", data.category = "Transcriptome Profiling", data.type = "Gene Expression Quantification", workflow.type = "HTSeq - FPKM")
GDCdownload(query_fpkm, files.per.chunk = 100)   # ��������
# 
#����Ԥ����
dataAssy_fpkm = GDCprepare(query_fpkm, summarizedExperiment = F)
#     
# # # ����������Ա�����
filename_fpkm = paste0("TCGA-UCEC.RNA.query.","FPKM", ".RData")
save(query_fpkm, dataAssy_fpkm, file = filename_fpkm)
load(filename_fpkm)
#     
# ��������
exp_data = dataAssy_fpkm %>% 
    filter(str_detect(X1, "^ENSG")) %>% # ��ȡENSG��ͷ�Ļ���
    mutate(X1 = gsub("\\..*", "", X1)) %>% # ɾ���汾��
    rename(ensembl_gene_id = X1) %>% # "X1"������λ"ensembl_gene_id"
    set_colnames(str_replace(colnames(.), "(-\\w*){3}$", "")) # �޸�TCGA����������3����-�ַ���ɾ��
       
exp_data[1:3,1:3]
    
    # if(data_type == "HTSeq - FPKM-UQ"){
    #     
    #     # log2ת��
    #     exp_data[,-1] = log2(exp_data[,-1] + 1)
        
# ���Ŀ����������������еı��ֵ
# DGAT2->ENSG00000062282: https://asia.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000185480;r=12:102120185-102197520
DGAT2_exp_fpkm = exp_data %>% 
    filter(ensembl_gene_id == "ENSG00000062282") %>% 
    mutate(Symbol = "DGAT2") %>%                 #����һ��symbol
    dplyr::select(Symbol, starts_with("TCGA-"))
    write.table(DGAT2_exp_data, "DGAT2_fpkm_allsamples.txt", row.names = F, sep = "\t", quote = F)
        
    
# ���
outfile = paste0("LIHC_Portal_RNA_","count", ".txt")
fwrite(exp_data, outfile, row.names = F, sep = "\t", quote = F)
    
# #### ����ע�� ----
# # biomaRt������nsembl���ݿ�
ensembl = biomaRt::useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
genes_info = biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name"), filters = "ensembl_gene_id", values = exp_data$ensembl_gene_id, mart = ensembl) # �޶�Ӧsymbol������ɾ��
# # symbolע��
exp_data_final = inner_join(exp_data, genes_info, by = "ensembl_gene_id", keep = F)
outfile = paste0("UCEC_RNA_","count", "_addSymbol.txt")
write.table(exp_data_final, outfile, row.names = F, sep = "\t", quote = F)
