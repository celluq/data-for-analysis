
rm(list = ls())
library(pacman)
p_load(tidyverse, survival, survminer)

# 设置路径
projectPath = "/Users/prf/SCIrepro/PARPBP"
survPath = paste(projectPath, "Gene_Surv", sep = "/")
if(!dir.exists(survPath)) dir.create(survPath)  #
dataPath = paste(projectPath, "Data", sep = "/")
setwd(dataPath)

#### 读入数据 ----
exp_data = read.table("clinic_data_clean_fpkm.txt", header = T, sep = "\t", check.names = F)
clinic_data = read.table("clinic_data_clean.txt", header = T, sep = "\t", row.names = 1)

#### 生存分析----
# 定义生存函数,分别用于OS和DFI的计算
do_surv = function(surv_type, subset = NULL){ # subset用于设置Fig.4样本子集: subset = c("I", "II")
    
    # surv_type = "OS"
    # subset = c("II","III", "IV")
    if(!is.null(subset)){ 
        # 基于分期的样本子集Fig.4
        stage = str_match(clinic_data$clinical_stage, "Stage ([^ABC]*)")[,2]
        index = which(stage %in% subset)
        clinic_data = clinic_data[index,]
        exp_data = exp_data[,index]
        
        #基于分级的样本子集
        grade<-clinic_data$neoplasm_histologic_grade
        subset = c("G3","High Grade")
        index = which(grade %in% subset)
        clinic_data = clinic_data[index,]
        exp_data = exp_data[,index]
        
        #基于病理类型的样本子集
        historical_type<-clinic_data$histological_type
        subset = c("Endometrioid endometrial adenocarcinoma")
        index = which(historical_type %in% subset)
        clinic_data = clinic_data[index,]
        exp_data = exp_data[,index]
        
        subsurvPath = paste(projectPath, "Clinic_Surv", sep = "/")
        if(!dir.exists(subsurvPath)) dir.create(subsurvPath) 
        setwd(subsurvPath)
        name = paste(paste0(subset, collapse = "-"), surv_type, sep = "_")
    }else{
        setwd(survPath)
        name = surv_type
    }
    
    # 整合
    data_for_surv = data.frame(exp = as.numeric(exp_data_mul[1,]), surv_time = clinic_data[,paste0(surv_type, ".time")]/365, surv_status = clinic_data[,surv_type])
    # 过滤
    index = which(!is.na(data_for_surv$surv_status) & data_for_surv$surv_time > 0)
    data_for_surv = data_for_surv[index,]
    
    # 计算分类点
    surv.cut = surv_cutpoint(data_for_surv, time = "surv_time", event = "surv_status", variables = "exp", minprop = 0.3)
    # 最优分类点
    opticut = surv.cut$cutpoint$cutpoint 
    # 样本分组
    data_for_surv_final = surv_categorize(surv.cut, labels = c(0, 1)) # 0为"low", 1为"high"
    
    # 生存分析
    # survfit object
    surv.obj = as.formula("Surv(surv_time, surv_status)~exp")
    # K-M survival curves
    fit = surv_fit(surv.obj, data = data_for_surv_final) # 鐢╯urvival::survfit鎶ラ敊: https://github.com/kassambara/survminer/issues/403
    # log-rank p http://rpkgs.datanovia.com/survminer/reference/surv_pvalue.html
    pval = surv_pvalue(fit)$pval
    # Cox proportional hazards regression model
    cox = coxph(surv.obj, data = data_for_surv_final)
    cox_summary = summary(cox)
    cox_res = c(cox_summary$coefficients[,c(1:2,5)], cox_summary$conf.int[,c(3,4)])
    # 提取结果，显示表格
    #cox %>% gtsummary::tbl_regression(exp = TRUE, label = list(exp ~ paste0(symbol, " (high vs. low)")))  # label瑕佹眰list鏍煎紡,~鍚庝负鏂板懡鍚?
    
    # 提取结果，输出
    surv_res = matrix(c(as.numeric(table(data_for_surv_final$exp)), pval, cox_res), nrow = 1)
    final_data<-data.frame(surv_res)
    
    #循环
    for (i in c(2:nrow(exp_data_mul))) {
        data_for_surv = data.frame(exp = as.numeric(exp_data_mul[i,]), surv_time = clinic_data[,paste0(surv_type, ".time")]/365, surv_status = clinic_data[,surv_type])
        # 过滤
        index = which(!is.na(data_for_surv$surv_status) & data_for_surv$surv_time > 0)
        data_for_surv = data_for_surv[index,]
        
        # 计算分类点
        surv.cut = surv_cutpoint(data_for_surv, time = "surv_time", event = "surv_status", variables = "exp", minprop = 0.3)
        # 最优分类点
        opticut = surv.cut$cutpoint$cutpoint 
        # 样本分组
        data_for_surv_final = surv_categorize(surv.cut, labels = c(0, 1)) # 0为"low", 1为"high"
        
        # 生存分析
        # survfit object
        surv.obj = as.formula("Surv(surv_time, surv_status)~exp")
        # K-M survival curves
        fit = surv_fit(surv.obj, data = data_for_surv_final) # 鐢╯urvival::survfit鎶ラ敊: https://github.com/kassambara/survminer/issues/403
        # log-rank p http://rpkgs.datanovia.com/survminer/reference/surv_pvalue.html
        pval = surv_pvalue(fit)$pval
        # Cox proportional hazards regression model
        cox = coxph(surv.obj, data = data_for_surv_final)
        cox_summary = summary(cox)
        cox_res = c(cox_summary$coefficients[,c(1:2,5)], cox_summary$conf.int[,c(3,4)])
        # 提取结果，显示表格
        #cox %>% gtsummary::tbl_regression(exp = TRUE, label = list(exp ~ paste0(symbol, " (high vs. low)")))  # label瑕佹眰list鏍煎紡,~鍚庝负鏂板懡鍚?
        
        # 提取结果，输出
        surv_res = data.frame(matrix(c(as.numeric(table(data_for_surv_final$exp)), pval, cox_res), nrow = 1))
        final_data<-final_data %>% add_row(surv_res)
        
    }
    colnames(final_data) = c("Num_Low", "Num_High", "Log-rank p", "coef", "HR", "Cox p", "Lower95%", "Upper95%")
    colnames(surv_res) = c("Num_Low", "Num_High", "Log-rank p", "coef", "HR", "Cox p", "Lower95%", "Upper95%")
    write.csv(data.frame(Symbol = symbol, surv_res), paste0(name, "_", symbol, "_surv_res.csv"), row.names = F, fileEncoding = "GBK")
    
    #### 绘图 ----
    # x轴范围
    maxX = ceiling(max(data_for_surv_final$surv_time))*1.05
    # log-rank p值格式化
    pText = ifelse(pval < 0.01, formatC(pval, digits = 2, format = "E"), round(pval, digits = 3))
    # https://github.com/kassambara/survminer
    # http://www.sthda.com/english/wiki/survminer-r-package-survival-data-analysis-and-visualization
    # ggsurvplot
    p = ggsurvplot(fit, data = data_for_surv_final, xlab = "Time(Years)", ylab = "Survival Probability (%)",  # ylab = "Survival probability (%)"
                   size = 1.2, palette = c("navyblue", "firebrick1"), # 曲线设置
                   xlim = c(0, maxX), ylim = c(0, 101),
                   legend.title = 'Category', legend.labs = c("Low","High"), font.legend = 14, legend = c(0.8,0.85), # 图例设置
                   #pval = paste0("HR = ", round(surv_res[5], 3), "\nlog-rank p = ", pText), pval.size = 5, pval.coord = c(0.6, 18), # p值设置
                   pval = F, # 用annotate设置文本
                   parse = T,
                   break.time.by = 2, # x轴刻度间隔
                   risk.table = T, risk.table.height = 0.2, risk.table.fontsize = 5, risk.table.col = "strata", tables.y.text = F, # 风险表设置
                   font.x = 16, font.y = 16, font.xtickslab = 13, font.ytickslab = 13,
                   fun = function(y) y*100)  # 纵坐标刻度显示为百分数
    # 设置"p"值斜体! 
    p$plot = p$plot + 
        annotate("text", x = 0.6, y = c(20, 10), label = c(paste0("HR = ", round(surv_res[5], 3)), "log-rank"), size = 5, hjust = 0, vjust = 1) + # hjust设置左右:0-左对齐, 0.5-居中, 1-右对齐；vjust设置上下对齐0-1
        # 科学计数法显示
        # https://stackoverflow.com/questions/25108062/how-to-write-after-decimal-zero-in-ggplot-geom-text
        # https://rstudio-pubs-static.s3.amazonaws.com/136237_170402e5f0b54561bf7605bdea98267a.html
        annotate("text", x = 2.5, y = 10, label = as.character(as.expression(substitute(italic(p)~"="~pText, list(pText = formatC(pText, format = "e", digits = 2))))), size = 5, hjust = 0, vjust = 1, parse = T)
    p$table = p$table + theme(axis.ticks = element_blank(), axis.line = element_blank(), # 去除轴和刻度
                              axis.title.x = element_blank(), axis.text.x = element_blank(), # 去除标题和刻度
                              axis.title.y = element_blank()) #, axis.text.y = element_blank() 与 tables.y.text = T 搭配使用
    ggsave(plot = print(p, newpage = FALSE), paste0(name, "_", symbol, "_Pvalue_", pval, ".pdf"), height = 5, width = 6)

}
ggsave(plot = print(p, newpage = FALSE), filename = "OS.StromalScore.pdf",height = 5, width = 10)
# 调用
do_surv("OS")
do_surv("DFI")

# 分期子集
do_surv("OS", subset = c("I", "II"))
do_surv("OS", subset = c("III", "IV"))
do_surv("DFI", subset = c("I", "II"))
do_surv("DFI", subset = c("III", "IV"))
