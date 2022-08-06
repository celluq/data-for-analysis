rm(list = ls())
library(pacman)
p_load(survival, rms, tidyverse)

#### 读入数据 ----
exp_data = read.table("PHGDH_fpkm_tumor.txt", header = T, sep = "\t", check.names = F)
symbol = rownames(exp_data)
clinic_data = read.table("clinic_data_clean.txt", header = T, sep = "\t")

#### 绘图 ----
setwd(nomoPath)
for(surv_type in c("OS", "DFI")){
  
  # surv_type = "OS"
  data_for_nomo = clinic_data %>% 
    mutate(!!symbol := as.numeric(exp_data), # 添加表达量为新的一列
           Stage = str_match(clinical_stage, "Stage ([^ABC]*)")[,2],
           Histologic_Grade = case_when(neoplasm_histologic_grade %in% c("G1") ~ "G1",
                                         neoplasm_histologic_grade  %in% c("G2") ~ "G2",
                                         neoplasm_histologic_grade  %in% c("G3","High Grade") ~ "G3" ),
           Histological_type = case_when(histological_type %in% c("Endometrioid endometrial adenocarcinoma") ~ "EEC",
                       histological_type  %in% c("Serous endometrial adenocarcinoma") ~ "SC",
                       histological_type  %in% c("Mixed serous and endometrioid" ) ~ "Mixed" ),
           Age = ifelse(age_at_initial_pathologic_diagnosis < 60, "<60", "≥60")) %>% # 
    rename(surv_time := !!paste0(surv_type, ".time"), surv_status := !!surv_type) %>% 
    dplyr::select(sampleID, surv_time, surv_status, all_of(symbol),Stage,Histologic_Grade,Histological_type,Age) %>%
    filter(!is.na(surv_type) & surv_time > 0)
  
  # 预测变量的分布汇总
  dfn = datadist(data_for_nomo)
  options(datadist = "dfn") 
  
  # 构建cox比例风险回归模型
  surv.obj = as.formula(paste0("Surv(surv_time, surv_status)~", symbol, "+Stage+Histologic_Grade+Age"))
  # psm和cph函数的区别: https://stackoverflow.com/questions/56161660/what-is-the-difference-between-psm-and-cph-in-rms-package-in-r
  fit = cph(surv.obj, data = data_for_nomo, x = T, y = T, surv = T)
  # 生存函数
  surv = Survival(fit) # rms包中的函数
  times = c(3,5,8) # 时间点，生存时间数据以年为单位
  # 图中显示的标签
  labes = lapply(times, function(time) paste0(time, "-year ", surv_type))
  funcs = lapply(times, function(time) eval(parse(text = paste0("function(x) surv(", time, "*365, x)"))))
  
  # 绘制列线图
  nom = nomogram(
    fit, # 生存模型
    lp = FALSE, # 不显示线性预测指标linear predictor
    fun = funcs, # an optional function to transform the linear predictors, and to plot on another axis.
    funlabel = unlist(labes), # 标签，如1-year OS"
    fun.at = seq(0, 1, 0.1) # 计算生存率的店
  )
  pdf(paste0(surv_type, "_Nomo.pdf"), width = 6, height = 4)
  par(mar = c(1, 0, 1, 0)) 
  plot(nom, xfrac = 0.4) # xfrac调整变量名到主图的距离
  dev.off()
  
  # The predictive accuracy and discriminative ability of the nomogram were determined by concordance index (C-index) and calibration curve
  sum.surv = summary(coxph(surv.obj, data = data_for_nomo))
  c_index_se = sum.surv$concordance
  # C-index
  c_index = c_index_se[1] # 0.6694121
  # C-index 95%置信区间
  c_index.ci_low = c_index - c_index_se[2] 
  c_index.ci_hig = c_index + c_index_se[2]
  
  # calibration curve
  # 配置绘图函数
  # 基于bootstrap训练集重抽样的方法验证列线图的预测准确性
  data_for_plot = lapply(times, function(time) {
    fit_for_plot = cph(surv.obj, x = T, y = T, data = data_for_nomo, surv = T, time.inc = time*365)
    # time.inc要与模型中设置的time.inc一致
    cal_for_plot = calibrate(fit_for_plot, cmethod = 'KM',
                             method = 'boot', # bootstrap抽样
                             u = time*365, # 时间要与模型中的time.inc一致
                             m = ceiling(nrow(data_for_nomo)/4), # 每组的样本量数
                             B = nrow(data_for_nomo)) # 最大再抽样的样本量
    return(cal_for_plot)
  })
  
  ## 绘图
  # Model performance is shown by the plot, relative to the 45-degree line, which represents perfect prediction.
  # 细节设置: https://stackoverflow.com/questions/57834244/changing-the-colour-of-a-calibration-plot
  colors = c("#69AC40", "#0095B0", "#9D1535")
  pdf(paste0(surv_type, "_Calib_new_1.pdf"), width = 5, height = 3)
  par(mar = c(3.2, 3.2, 1, 0.4), mgp = c(1.8, 0.5, 0)) 
  plot(data_for_plot[[1]], bty = "l", lty = 1, lwd = 2, xlim = c(0,0.95), ylim= c(0,1), 
       #errbar.col = c("#00468BFF"),
       errbar = F, # 如果绘制多条曲线，尽量不绘制errbar
       xlab = "Nomogram predicted survival probability",
       ylab = "Actual survival probability",
       riskdist = F, # 去除顶部的竖线
       subtitles = F, # 去除左下以及右下的注释
       par.corrected = list(col = "white")) # 不画X(overfitting-corrected estimates)
  lines(data_for_plot[[1]][,c('mean.predicted',"KM")], 
        type = 'o', # (overplotted: lines over points); "b"(both: points and lines)
        lwd = 2, col = colors[1], pch = 16)
  abline(0, 1, lty = 3, lwd = 2, col = c("#224444"))
  # 加曲线
  lines(data_for_plot[[2]][,c('mean.predicted',"KM")], type = 'o', lwd = 2, col = colors[2], pch = 16)
  lines(data_for_plot[[3]][,c('mean.predicted',"KM")], type = 'o', lwd = 2, col = colors[3], pch = 16)
  legend("bottomright", paste0(times, "-Year Survival"), col = colors[1:3], adj = 0, cex = 0.8, lwd = 2, bty = "n", pch = 16)
  dev.off()
  
}
