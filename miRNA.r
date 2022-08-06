##### 成熟体miRNA的下载-----
rm(list = ls())
library(pacman)
p_load(TCGAbiolinks, tidyverse, magrittr, data.table, miRBaseVersions.db, miRNAmeConverter)

#### TCGAbiolinks ----
repeat {
    query =  try(GDCquery(project = "TCGA-UCEC", legacy = FALSE, data.category = "Transcriptome Profiling", data.type = "Isoform Expression Quantification"), silent = TRUE)	
    if(class(query) != "try-error") {
        break
    }else{
        print("Try to connect to GDC...")
    }
}

GDCdownload(query, files.per.chunk = 100)
dataAssy = GDCprepare(query)
# save(query, dataAssy, file = "TCGA-LIHC.miRNA.Isoform.query.raw.RData")
# load("TCGA-LIHC.miRNA.Isoform.query.raw.RData")

# 分别整理count和RPM值
for(item in c("read_count", "reads_per_million_miRNA_mapped")){
    
    # 配置输出文件名
    name = ifelse(item == "read_count", "Counts", "RPM")
    
    # 输出数据, miRNA_region为所有的miRNA名
    data = dataAssy %>% filter(str_detect(miRNA_region, "mature,")) %>% # 提取成熟体的行
        dplyr::select(!!item, "miRNA_region", "barcode") %>% # 提取3列
        # 构建表达谱矩阵
        group_by(barcode, miRNA_region) %>% # 以样本名和miRNA为分组
        summarize(!!paste0(item, "_sum") := sum(get(item), na.rm = TRUE)) %>% #同一样本同一miRNA值叠加
        spread(key = barcode, !!paste0(item, "_sum"), fill = 0) %>% # spread构建矩阵，空值设为0
        mutate(miRNA_region = gsub("mature,", "", miRNA_region)) %>% # 重命名，去除mature字符
        column_to_rownames(var = "miRNA_region") %>% # miRNA_region列转为行名
        set_colnames(str_match(colnames(.), "(TCGA-[^-]*-[^-]*-[^-]*)")[,2])   # 设置列名
    if(name == "RPM") data = log2(data+1)
    
    #### miRNA名字转换
    items = miRBaseVersions.db::select(miRBaseVersions.db, keys = rownames(data), keytype = "MIMAT", columns = "*")
    names = items[items$VERSION == 22, c("ACCESSION", "NAME")] # miRBase最新为v22
    data_final = data[names$ACCESSION,]
    rownames(data_final) = names$NAME
    
    outfile = paste0("LIHC_Portal_miRNA.mature_", name, ".txt")
    fwrite(data.frame(Symbol = rownames(data_final), data_final, check.names = F), outfile, sep = "\t", row.names = F, quote = F)

}
