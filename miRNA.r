##### ������miRNA������-----
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

# �ֱ�����count��RPMֵ
for(item in c("read_count", "reads_per_million_miRNA_mapped")){
    
    # ��������ļ���
    name = ifelse(item == "read_count", "Counts", "RPM")
    
    # �������, miRNA_regionΪ���е�miRNA��
    data = dataAssy %>% filter(str_detect(miRNA_region, "mature,")) %>% # ��ȡ���������
        dplyr::select(!!item, "miRNA_region", "barcode") %>% # ��ȡ3��
        # ��������׾���
        group_by(barcode, miRNA_region) %>% # ����������miRNAΪ����
        summarize(!!paste0(item, "_sum") := sum(get(item), na.rm = TRUE)) %>% #ͬһ����ͬһmiRNAֵ����
        spread(key = barcode, !!paste0(item, "_sum"), fill = 0) %>% # spread�������󣬿�ֵ��Ϊ0
        mutate(miRNA_region = gsub("mature,", "", miRNA_region)) %>% # ��������ȥ��mature�ַ�
        column_to_rownames(var = "miRNA_region") %>% # miRNA_region��תΪ����
        set_colnames(str_match(colnames(.), "(TCGA-[^-]*-[^-]*-[^-]*)")[,2])   # ��������
    if(name == "RPM") data = log2(data+1)
    
    #### miRNA����ת��
    items = miRBaseVersions.db::select(miRBaseVersions.db, keys = rownames(data), keytype = "MIMAT", columns = "*")
    names = items[items$VERSION == 22, c("ACCESSION", "NAME")] # miRBase����Ϊv22
    data_final = data[names$ACCESSION,]
    rownames(data_final) = names$NAME
    
    outfile = paste0("LIHC_Portal_miRNA.mature_", name, ".txt")
    fwrite(data.frame(Symbol = rownames(data_final), data_final, check.names = F), outfile, sep = "\t", row.names = F, quote = F)

}
