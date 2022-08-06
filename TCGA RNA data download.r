##### RNA数据下载与整理-----
rm(list = ls())

# 安装加载需要的R包
install.packages("pacman", "TCGAbiolinks","tidyverse","magrittr","data.table","biomaRt")
library(pacman)
p_load(TCGAbiolinks, tidyverse, magrittr, data.table, biomaRt)

# 设置工作目录
setwd("D:\\R\\练习")

#### 使用TCGAbiolinks下载Counts和PKM-UQ数据----
#https://www.bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html
#第一步 下载count数据
query_count = GDCquery(project = "TCGA-UCEC", legacy = FALSE, experimental.strategy = "RNA-Seq", data.category = "Transcriptome Profiling", data.type = "Gene Expression Quantification", workflow.type = "HTSeq - Counts")
GDCdownload(query_count, files.per.chunk = 100)   # 分批下载
# 数据预处理
dataAssy_count = GDCprepare(query_count, summarizedExperiment = F)
# 保存变量，以备后用
filename_count = paste0("TCGA-UCEC.RNA.query.","Counts", ".RData")
save(query_count, dataAssy_count, file = filename_count)
load(filename_count)

#第二步 下载FKPM数据
query_fpkm = GDCquery(project = "TCGA-UCEC", legacy = FALSE, experimental.strategy = "RNA-Seq", data.category = "Transcriptome Profiling", data.type = "Gene Expression Quantification", workflow.type = "HTSeq - FPKM")
GDCdownload(query_fpkm, files.per.chunk = 100)   # 分批下载
# 
#数据预处理
dataAssy_fpkm = GDCprepare(query_fpkm, summarizedExperiment = F)
#     
# # # 保存变量，以备后用
filename_fpkm = paste0("TCGA-UCEC.RNA.query.","FPKM", ".RData")
save(query_fpkm, dataAssy_fpkm, file = filename_fpkm)
load(filename_fpkm)
#     
# 数据整理
exp_data = dataAssy_fpkm %>% 
    filter(str_detect(X1, "^ENSG")) %>% # 提取ENSG开头的基因
    mutate(X1 = gsub("\\..*", "", X1)) %>% # 删除版本号
    rename(ensembl_gene_id = X1) %>% # "X1"重命名位"ensembl_gene_id"
    set_colnames(str_replace(colnames(.), "(-\\w*){3}$", "")) # 修改TCGA名，将后面3个“-字符”删除
       
exp_data[1:3,1:3]
    
    # if(data_type == "HTSeq - FPKM-UQ"){
    #     
    #     # log2转化
    #     exp_data[,-1] = log2(exp_data[,-1] + 1)
        
# 输出目标基因在所在样本中的表达值
# DGAT2->ENSG00000062282: https://asia.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000185480;r=12:102120185-102197520
DGAT2_exp_fpkm = exp_data %>% 
    filter(ensembl_gene_id == "ENSG00000062282") %>% 
    mutate(Symbol = "DGAT2") %>%                 #增加一列symbol
    dplyr::select(Symbol, starts_with("TCGA-"))
    write.table(DGAT2_exp_data, "DGAT2_fpkm_allsamples.txt", row.names = F, sep = "\t", quote = F)
        
    
# 输出
outfile = paste0("LIHC_Portal_RNA_","count", ".txt")
fwrite(exp_data, outfile, row.names = F, sep = "\t", quote = F)
    
# #### 基因注释 ----
# # biomaRt包访问nsembl数据库
ensembl = biomaRt::useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
genes_info = biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name"), filters = "ensembl_gene_id", values = exp_data$ensembl_gene_id, mart = ensembl) # 无对应symbol的自行删除
# # symbol注释
exp_data_final = inner_join(exp_data, genes_info, by = "ensembl_gene_id", keep = F)
outfile = paste0("UCEC_RNA_","count", "_addSymbol.txt")
write.table(exp_data_final, outfile, row.names = F, sep = "\t", quote = F)
