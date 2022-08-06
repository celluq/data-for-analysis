
rm(list = ls())
library(pacman)
p_load(tidyverse, survival, survminer)

# ����·��
projectPath = "/Users/prf/SCIrepro/PARPBP"
survPath = paste(projectPath, "Gene_Surv", sep = "/")
if(!dir.exists(survPath)) dir.create(survPath)  #
dataPath = paste(projectPath, "Data", sep = "/")
setwd(dataPath)

#### �������� ----
exp_data = read.table("clinic_data_clean_fpkm.txt", header = T, sep = "\t", check.names = F)
clinic_data = read.table("clinic_data_clean.txt", header = T, sep = "\t", row.names = 1)

#### �������----
# �������溯��,�ֱ�����OS��DFI�ļ���
do_surv = function(surv_type, subset = NULL){ # subset��������Fig.4�����Ӽ�: subset = c("I", "II")
    
    # surv_type = "OS"
    # subset = c("II","III", "IV")
    if(!is.null(subset)){ 
        # ���ڷ��ڵ������Ӽ�Fig.4
        stage = str_match(clinic_data$clinical_stage, "Stage ([^ABC]*)")[,2]
        index = which(stage %in% subset)
        clinic_data = clinic_data[index,]
        exp_data = exp_data[,index]
        
        #���ڷּ��������Ӽ�
        grade<-clinic_data$neoplasm_histologic_grade
        subset = c("G3","High Grade")
        index = which(grade %in% subset)
        clinic_data = clinic_data[index,]
        exp_data = exp_data[,index]
        
        #���ڲ������͵������Ӽ�
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
    
    # ����
    data_for_surv = data.frame(exp = as.numeric(exp_data_mul[1,]), surv_time = clinic_data[,paste0(surv_type, ".time")]/365, surv_status = clinic_data[,surv_type])
    # ����
    index = which(!is.na(data_for_surv$surv_status) & data_for_surv$surv_time > 0)
    data_for_surv = data_for_surv[index,]
    
    # ��������
    surv.cut = surv_cutpoint(data_for_surv, time = "surv_time", event = "surv_status", variables = "exp", minprop = 0.3)
    # ���ŷ����
    opticut = surv.cut$cutpoint$cutpoint 
    # ��������
    data_for_surv_final = surv_categorize(surv.cut, labels = c(0, 1)) # 0Ϊ"low", 1Ϊ"high"
    
    # �������
    # survfit object
    surv.obj = as.formula("Surv(surv_time, surv_status)~exp")
    # K-M survival curves
    fit = surv_fit(surv.obj, data = data_for_surv_final) # 用survival::survfit报错: https://github.com/kassambara/survminer/issues/403
    # log-rank p http://rpkgs.datanovia.com/survminer/reference/surv_pvalue.html
    pval = surv_pvalue(fit)$pval
    # Cox proportional hazards regression model
    cox = coxph(surv.obj, data = data_for_surv_final)
    cox_summary = summary(cox)
    cox_res = c(cox_summary$coefficients[,c(1:2,5)], cox_summary$conf.int[,c(3,4)])
    # ��ȡ�������ʾ���
    #cox %>% gtsummary::tbl_regression(exp = TRUE, label = list(exp ~ paste0(symbol, " (high vs. low)")))  # label要求list格式,~后为新命�?
    
    # ��ȡ��������
    surv_res = matrix(c(as.numeric(table(data_for_surv_final$exp)), pval, cox_res), nrow = 1)
    final_data<-data.frame(surv_res)
    
    #ѭ��
    for (i in c(2:nrow(exp_data_mul))) {
        data_for_surv = data.frame(exp = as.numeric(exp_data_mul[i,]), surv_time = clinic_data[,paste0(surv_type, ".time")]/365, surv_status = clinic_data[,surv_type])
        # ����
        index = which(!is.na(data_for_surv$surv_status) & data_for_surv$surv_time > 0)
        data_for_surv = data_for_surv[index,]
        
        # ��������
        surv.cut = surv_cutpoint(data_for_surv, time = "surv_time", event = "surv_status", variables = "exp", minprop = 0.3)
        # ���ŷ����
        opticut = surv.cut$cutpoint$cutpoint 
        # ��������
        data_for_surv_final = surv_categorize(surv.cut, labels = c(0, 1)) # 0Ϊ"low", 1Ϊ"high"
        
        # �������
        # survfit object
        surv.obj = as.formula("Surv(surv_time, surv_status)~exp")
        # K-M survival curves
        fit = surv_fit(surv.obj, data = data_for_surv_final) # 用survival::survfit报错: https://github.com/kassambara/survminer/issues/403
        # log-rank p http://rpkgs.datanovia.com/survminer/reference/surv_pvalue.html
        pval = surv_pvalue(fit)$pval
        # Cox proportional hazards regression model
        cox = coxph(surv.obj, data = data_for_surv_final)
        cox_summary = summary(cox)
        cox_res = c(cox_summary$coefficients[,c(1:2,5)], cox_summary$conf.int[,c(3,4)])
        # ��ȡ�������ʾ���
        #cox %>% gtsummary::tbl_regression(exp = TRUE, label = list(exp ~ paste0(symbol, " (high vs. low)")))  # label要求list格式,~后为新命�?
        
        # ��ȡ��������
        surv_res = data.frame(matrix(c(as.numeric(table(data_for_surv_final$exp)), pval, cox_res), nrow = 1))
        final_data<-final_data %>% add_row(surv_res)
        
    }
    colnames(final_data) = c("Num_Low", "Num_High", "Log-rank p", "coef", "HR", "Cox p", "Lower95%", "Upper95%")
    colnames(surv_res) = c("Num_Low", "Num_High", "Log-rank p", "coef", "HR", "Cox p", "Lower95%", "Upper95%")
    write.csv(data.frame(Symbol = symbol, surv_res), paste0(name, "_", symbol, "_surv_res.csv"), row.names = F, fileEncoding = "GBK")
    
    #### ��ͼ ----
    # x�᷶Χ
    maxX = ceiling(max(data_for_surv_final$surv_time))*1.05
    # log-rank pֵ��ʽ��
    pText = ifelse(pval < 0.01, formatC(pval, digits = 2, format = "E"), round(pval, digits = 3))
    # https://github.com/kassambara/survminer
    # http://www.sthda.com/english/wiki/survminer-r-package-survival-data-analysis-and-visualization
    # ggsurvplot
    p = ggsurvplot(fit, data = data_for_surv_final, xlab = "Time(Years)", ylab = "Survival Probability (%)",  # ylab = "Survival probability (%)"
                   size = 1.2, palette = c("navyblue", "firebrick1"), # ��������
                   xlim = c(0, maxX), ylim = c(0, 101),
                   legend.title = 'Category', legend.labs = c("Low","High"), font.legend = 14, legend = c(0.8,0.85), # ͼ������
                   #pval = paste0("HR = ", round(surv_res[5], 3), "\nlog-rank p = ", pText), pval.size = 5, pval.coord = c(0.6, 18), # pֵ����
                   pval = F, # ��annotate�����ı�
                   parse = T,
                   break.time.by = 2, # x��̶ȼ��
                   risk.table = T, risk.table.height = 0.2, risk.table.fontsize = 5, risk.table.col = "strata", tables.y.text = F, # ���ձ�����
                   font.x = 16, font.y = 16, font.xtickslab = 13, font.ytickslab = 13,
                   fun = function(y) y*100)  # ������̶���ʾΪ�ٷ���
    # ����"p"ֵб��! 
    p$plot = p$plot + 
        annotate("text", x = 0.6, y = c(20, 10), label = c(paste0("HR = ", round(surv_res[5], 3)), "log-rank"), size = 5, hjust = 0, vjust = 1) + # hjust��������:0-�����, 0.5-����, 1-�Ҷ��룻vjust�������¶���0-1
        # ��ѧ��������ʾ
        # https://stackoverflow.com/questions/25108062/how-to-write-after-decimal-zero-in-ggplot-geom-text
        # https://rstudio-pubs-static.s3.amazonaws.com/136237_170402e5f0b54561bf7605bdea98267a.html
        annotate("text", x = 2.5, y = 10, label = as.character(as.expression(substitute(italic(p)~"="~pText, list(pText = formatC(pText, format = "e", digits = 2))))), size = 5, hjust = 0, vjust = 1, parse = T)
    p$table = p$table + theme(axis.ticks = element_blank(), axis.line = element_blank(), # ȥ����Ϳ̶�
                              axis.title.x = element_blank(), axis.text.x = element_blank(), # ȥ������Ϳ̶�
                              axis.title.y = element_blank()) #, axis.text.y = element_blank() �� tables.y.text = T ����ʹ��
    ggsave(plot = print(p, newpage = FALSE), paste0(name, "_", symbol, "_Pvalue_", pval, ".pdf"), height = 5, width = 6)

}
ggsave(plot = print(p, newpage = FALSE), filename = "OS.StromalScore.pdf",height = 5, width = 10)
# ����
do_surv("OS")
do_surv("DFI")

# �����Ӽ�
do_surv("OS", subset = c("I", "II"))
do_surv("OS", subset = c("III", "IV"))
do_surv("DFI", subset = c("I", "II"))
do_surv("DFI", subset = c("III", "IV"))
