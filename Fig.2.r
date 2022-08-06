##### 临床分组表达分析 -----
rm(list = ls())

# 安装需要加载的包
library(pacman)
p_load(tidyverse, magrittr, data.table, edgeR, ggbeeswarm, ggsignif)

# 设置工作目录
projectPath = "/Users/prf/SCIrepro/PARPBP"
deaPath = paste(projectPath, "DEA_Clinic", sep = "/")
if(!dir.exists(deaPath)) dir.create(deaPath)  # 
dataPath = paste(projectPath, "Data", sep = "/")
setwd(dataPath)

#### 读取数据 ----
exp_data = fread("clinic_data_clean_counts.txt", h = T, sep = "\t", check.names = F) %>% column_to_rownames("gene_symbol")
clinic_data = read.table("clinic_data_clean.txt", header = T, sep = "\t")

# 目标基因
target = "ENSG00000141753"
symbol = "PHGDH"

#### 定义差异分析和绘图函数 ----
# 蜂群图绘图函数
beeswarm_plot = function(target_gene_exp, group, pairs_list, anno, map_signif_level = FALSE, title = NULL, parse = FALSE, symbol){
    
    data_for_plot = data.frame(exp = as.numeric(target_gene_exp), group)
    nG = length(levels(data_for_plot$group))
    # y轴范围
    ymax = max(data_for_plot$exp)
    ymin = min(data_for_plot$exp)
    p = ggplot(data_for_plot, aes(x = group, y = exp)) +
        geom_beeswarm(size = 1.2, cex = 1.3, alpha = 0.8, aes(col = group, shape = group),groupOnX = T) +  
        labs(title = title, x = "", y = bquote(italic(.(symbol))~expression~(log[2]))) +
        scale_x_discrete(labels = paste(levels(data_for_plot$group), "\n", "(N=", table(data_for_plot$group), ")", sep = "")) +
        scale_colour_manual(values = c("#125488", "#E96463", "#74af96", "#e39c5c", "#E8487F")[seq_len(nG)]) + 
        scale_shape_manual(values = c(16,17,15,18,25,8)[seq_len(nG)]) + 
        stat_summary(fun = median, geom = "point", size = 1, col = "black") +  # median处两个点
        stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", size = 0.3, width = 0.4, col = "black") +  
        stat_summary(fun.data = function(x) median_hilow(x, 0.5), geom = "errorbar", width = 0.25, col = "black") +  
        theme_classic() + theme(legend.position = 0)
    
    # 显著性标记
    if(str_sub(map_signif_level, 1, 1) == "T"){
        anno = case_when(anno < 0.001 ~ "***",
                         anno < 0.01 ~ "**",
                         anno < 0.05 ~ "*",
                         TRUE ~ "NS")
    }else{
        anno = paste0("list(~italic(p)==", signif(anno, 3), ")")
        parse = TRUE
    }

    # 绘图
    step = 0.2
    # https://github.com/cran/ggsignif/blob/master/R/significance_annotation.R
    p + ggsignif::geom_signif(comparisons = pairs_list, # 比较组，列表
                              test = NULL, # 不做显著性检验，用自己统计学方法得到的annotations!
                              annotations = anno, # 注释内容，用*表示,配合test = NULL使用
                              step_increase = step, # 每级跨度, step倍的y轴范围
                              margin_top = 0.1, # 最低标签到y轴最大值的距离, 0.1倍的y
                              vjust = 0.05, # 标签到横坐标的距离
                              size = 0.5, # 横线段的粗细
                              textsize = 4, # 标签字体大小
                              parse = parse,
                              tip_length = 0.02) + # 小竖线的长度
        ylim(ymin-0.1, ymax+0.1*(ymax-ymin)+step*(ymax-ymin)*length(pairs_list)) +
        theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14)) + 
        theme(plot.title = element_text(size = 15, hjust = 0.5), plot.margin = unit(c(0.3, 0.2, -0.3, 0.2), "cm"))
    ggsave(paste0(title, ".pdf"), width = 4, height = 3)
    
}

# edgeR差异基因表达函数
edgeR_dea = function(exp, group, title, symbol){
    
    group = factor(group)
    dge = DGEList(counts = exp, group = group)
    
    # filter
    keep = filterByExpr(dge)
    dge_filter = dge[keep, keep.lib.sizes = FALSE]
    
    # TMM normalization
    dge_norm = calcNormFactors(dge_filter)  # 默认method = "TMM"
    logcpm = cpm(dge_norm, log = TRUE, prior.count = 1)
    exp_only_target = t(as.matrix(logcpm[symbol,]))

    # Estimating the dispersion
    design = model.matrix(~0+group)
    colnames(design) = levels(group)
    y = estimateDisp(dge_norm, design, robust = TRUE)
    fit = glmQLFit(y, design)
    
    # 所有的比较
    pairs_tmp = t(combn(levels(group), 2))[,2:1]
    if(is.matrix(pairs_tmp)){
        pairs_list = lapply(seq_len(nrow(pairs_tmp)), function(p) pairs_tmp[p,])
        pairs = tidyr::unite(data.frame(pairs_tmp), col = pairs, sep = "-")$pairs
    }else{
        pairs_list = list(pairs_tmp)
        pairs = paste(pairs_tmp, collapse = "-")
    }

    # 目标基因的差异分析
    result = sapply(seq_len(length((pairs))), function(i) {
        res = glmQLFTest(fit, contrast = makeContrasts(contrasts = pairs[i], levels = design))
        topTags(res, n = Inf)$table[symbol,]
    })
    
    dea_final = data.frame(t(result), row.names = pairs)
    
    # 绘图
    beeswarm_plot(target_gene_exp = exp_only_target, group = group, pairs_list = pairs_list, anno = as.numeric(dea_final$PValue), map_signif_level = T, title = title, symbol = symbol)
    
}

### FIGO stage ----
clinic_data_stage = clinic_data %>% 
    dplyr::select(c("sampleID", "clinical_stage")) %>% 
    mutate(clinical_stage = str_match(clinical_stage, "Stage ([^ABC]*)")[,2]) %>% 
    filter(!is.na(clinical_stage))

edgeR_dea(exp = exp_data[,clinic_data_stage$sampleID], group = clinic_data_stage$clinical_stage, title = "FIGO stage", symbol = symbol) 

### Histologic grade ----
clinic_data_grade = clinic_data %>% 
    dplyr::select(c("sampleID", "neoplasm_histologic_grade")) %>%
    mutate(neoplasm_histologic_grade = 
               case_when(neoplasm_histologic_grade %in% c("G1","G2") ~ "G1G2",
                         neoplasm_histologic_grade %in% c("G3","High Grade") ~ "G3")) %>% 
    filter(!is.na(neoplasm_histologic_grade) & neoplasm_histologic_grade != "") %>% 
    droplevels
edgeR_dea(exp = exp_data[,clinic_data_grade$sampleID], group = clinic_data_grade$neoplasm_histologic_grade, title = "Histologic grade", symbol = symbol)

### Diabetes ----
clinic_diabetes = clinic_data %>% 
    dplyr::select(c("sampleID", "diabetes")) %>% 
    filter(!is.na(diabetes) & diabetes != "") %>% 
    droplevels
edgeR_dea(exp = exp_data[,clinic_diabetes$sampleID], group = clinic_diabetes$diabetes, title = "Diabetes",symbol = symbol)

###historical type ----
clinic_data_historical = clinic_data %>% 
    dplyr::select(c("sampleID", "histological_type")) %>% 
     mutate(histological_type = 
               case_when(histological_type == "Endometrioid endometrial adenocarcinoma" ~ "EEC",
                         histological_type== "Serous endometrial adenocarcinoma"~ "SC",
                         histological_type== "Mixed serous and endometrioid"~ "Mix")) %>% 
    #mutate(histological_type = 
               #factor(histological_type, levels = c("Endometrioid endometrial adenocarcinoma","Mixed serous and endometrioid",
                #                                    "Serous endometrial adenocarcinoma"),labels = c("EEC", "Mix", "SC"))) %>% 
    filter(!is.na(histological_type))
edgeR_dea(exp = exp_data[,clinic_data_historical$sampleID], group = clinic_data_historical$histological_type, title = "Histological Type new", symbol = symbol)

### BMI ----
clinic_data_bmi = clinic_data %>% 
    dplyr::select(c("sampleID", "BMI")) %>% 
    mutate(BMI = 
               case_when(BMI < 30 ~ "None_Obesity",
                         BMI >= 30 ~ "Obsity")) %>% 
    filter(!is.na(BMI))
edgeR_dea(exp = exp_data[,clinic_data_bmi$sampleID], group = clinic_data_bmi$BMI, title = "BMI", symbol = symbol)

### Disease status ----
clinic_data_dfi = clinic_data %>% 
    dplyr::select(c("sampleID", "DFI")) %>% 
    mutate(DFI = ifelse(DFI == 1, "YES", "NO")) %>% 
    filter(!is.na(DFI))
edgeR_dea(exp = exp_data[,clinic_data_dfi$sampleID], group = clinic_data_dfi$DFI, title = "Disease status(Recurrence)", symbol = symbol)

### Living status ----
clinic_data_os = clinic_data %>% 
    dplyr::select(c("sampleID", "OS")) %>% 
    mutate(OS = ifelse(OS == 1, "Dead", "Alive")) %>% 
    filter(!is.na(OS))
edgeR_dea(exp = exp_data[,clinic_data_os$sampleID], group = clinic_data_os$OS, title = "Living status", symbol = symbol)

###TCGA ----
clinic_data_TCGA = clinic_data %>% 
    dplyr::select(c("sampleID", "CDE_ID_3226963")) %>% 
    mutate(CDE_ID_3226963 = 
               factor(CDE_ID_3226963, levels = c("MSI-H","MSI-L","MSS"),
                      labels = c("MSI_H","MSI_L","MSS"))) %>% 
    filter(!is.na(CDE_ID_3226963))
edgeR_dea(exp = exp_data[,clinic_data_TCGA$sampleID], group = clinic_data_TCGA$CDE_ID_3226963, title = "TCGA", symbol = symbol)

###TCGA MSS/MSI----
clinic_data_TCGA_MSI_MSS = clinic_data %>% 
    dplyr::select(c("sampleID", "CDE_ID_3226963")) %>% 
    mutate(CDE_ID_3226963 = 
             case_when(CDE_ID_3226963 == "MSI-H" ~ "MSI",
                       CDE_ID_3226963 %in% c("MSS","MSI-L")~ "MSS")) %>% 
    filter(!is.na(CDE_ID_3226963))
edgeR_dea(exp = exp_data[,clinic_data_TCGA_MSI_MSS$sampleID], group = clinic_data_TCGA_MSI_MSS$CDE_ID_3226963, title = "MSI_MSS", symbol = symbol)

#Age
clinic_data_age = clinic_data %>% 
    dplyr::select(c("sampleID", "age_at_initial_pathologic_diagnosis")) %>% 
    mutate(age_at_initial_pathologic_diagnosis = 
               case_when(age_at_initial_pathologic_diagnosis < 60 ~ "younger",
                         age_at_initial_pathologic_diagnosis >= 60 ~ "elder")) %>% 
    filter(!is.na(age_at_initial_pathologic_diagnosis))
edgeR_dea(exp = exp_data[,clinic_data_age$sampleID], group = clinic_data_age$age_at_initial_pathologic_diagnosis, title = "Age", symbol = symbol)

#menopause
clinic_data_menopause = clinic_data %>% 
    dplyr::select(c("sampleID", "menopause_status")) %>% 
    mutate(menopause_status=
               factor(menopause_status,levels = c("Post (prior bilateral ovariectomy OR >12 mo since LMP with no prior hysterectomy)","Pre (<6 months since LMP AND no prior bilateral ovariectomy AND not on estrogen replacement)"),labels = c("Yes","No"))) %>% 
    filter(!is.na(menopause_status))
edgeR_dea(exp = exp_data[,clinic_data_menopause$sampleID], group = clinic_data_menopause$menopause_status, title = "Menopause", symbol = symbol)

