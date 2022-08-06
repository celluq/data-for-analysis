rm(list = ls())
library(pacman)
p_load(survival, rms, tidyverse)

#### �������� ----
exp_data = read.table("PHGDH_fpkm_tumor.txt", header = T, sep = "\t", check.names = F)
symbol = rownames(exp_data)
clinic_data = read.table("clinic_data_clean.txt", header = T, sep = "\t")

#### ��ͼ ----
setwd(nomoPath)
for(surv_type in c("OS", "DFI")){
  
  # surv_type = "OS"
  data_for_nomo = clinic_data %>% 
    mutate(!!symbol := as.numeric(exp_data), # ��ӱ����Ϊ�µ�һ��
           Stage = str_match(clinical_stage, "Stage ([^ABC]*)")[,2],
           Histologic_Grade = case_when(neoplasm_histologic_grade %in% c("G1") ~ "G1",
                                         neoplasm_histologic_grade  %in% c("G2") ~ "G2",
                                         neoplasm_histologic_grade  %in% c("G3","High Grade") ~ "G3" ),
           Histological_type = case_when(histological_type %in% c("Endometrioid endometrial adenocarcinoma") ~ "EEC",
                       histological_type  %in% c("Serous endometrial adenocarcinoma") ~ "SC",
                       histological_type  %in% c("Mixed serous and endometrioid" ) ~ "Mixed" ),
           Age = ifelse(age_at_initial_pathologic_diagnosis < 60, "<60", "��60")) %>% # 
    rename(surv_time := !!paste0(surv_type, ".time"), surv_status := !!surv_type) %>% 
    dplyr::select(sampleID, surv_time, surv_status, all_of(symbol),Stage,Histologic_Grade,Histological_type,Age) %>%
    filter(!is.na(surv_type) & surv_time > 0)
  
  # Ԥ������ķֲ�����
  dfn = datadist(data_for_nomo)
  options(datadist = "dfn") 
  
  # ����cox�������ջع�ģ��
  surv.obj = as.formula(paste0("Surv(surv_time, surv_status)~", symbol, "+Stage+Histologic_Grade+Age"))
  # psm��cph����������: https://stackoverflow.com/questions/56161660/what-is-the-difference-between-psm-and-cph-in-rms-package-in-r
  fit = cph(surv.obj, data = data_for_nomo, x = T, y = T, surv = T)
  # ���溯��
  surv = Survival(fit) # rms���еĺ���
  times = c(3,5,8) # ʱ��㣬����ʱ����������Ϊ��λ
  # ͼ����ʾ�ı�ǩ
  labes = lapply(times, function(time) paste0(time, "-year ", surv_type))
  funcs = lapply(times, function(time) eval(parse(text = paste0("function(x) surv(", time, "*365, x)"))))
  
  # ��������ͼ
  nom = nomogram(
    fit, # ����ģ��
    lp = FALSE, # ����ʾ����Ԥ��ָ��linear predictor
    fun = funcs, # an optional function to transform the linear predictors, and to plot on another axis.
    funlabel = unlist(labes), # ��ǩ����1-year OS"
    fun.at = seq(0, 1, 0.1) # ���������ʵĵ�
  )
  pdf(paste0(surv_type, "_Nomo.pdf"), width = 6, height = 4)
  par(mar = c(1, 0, 1, 0)) 
  plot(nom, xfrac = 0.4) # xfrac��������������ͼ�ľ���
  dev.off()
  
  # The predictive accuracy and discriminative ability of the nomogram were determined by concordance index (C-index) and calibration curve
  sum.surv = summary(coxph(surv.obj, data = data_for_nomo))
  c_index_se = sum.surv$concordance
  # C-index
  c_index = c_index_se[1] # 0.6694121
  # C-index 95%��������
  c_index.ci_low = c_index - c_index_se[2] 
  c_index.ci_hig = c_index + c_index_se[2]
  
  # calibration curve
  # ���û�ͼ����
  # ����bootstrapѵ�����س����ķ�����֤����ͼ��Ԥ��׼ȷ��
  data_for_plot = lapply(times, function(time) {
    fit_for_plot = cph(surv.obj, x = T, y = T, data = data_for_nomo, surv = T, time.inc = time*365)
    # time.incҪ��ģ�������õ�time.incһ��
    cal_for_plot = calibrate(fit_for_plot, cmethod = 'KM',
                             method = 'boot', # bootstrap����
                             u = time*365, # ʱ��Ҫ��ģ���е�time.incһ��
                             m = ceiling(nrow(data_for_nomo)/4), # ÿ�����������
                             B = nrow(data_for_nomo)) # ����ٳ�����������
    return(cal_for_plot)
  })
  
  ## ��ͼ
  # Model performance is shown by the plot, relative to the 45-degree line, which represents perfect prediction.
  # ϸ������: https://stackoverflow.com/questions/57834244/changing-the-colour-of-a-calibration-plot
  colors = c("#69AC40", "#0095B0", "#9D1535")
  pdf(paste0(surv_type, "_Calib_new_1.pdf"), width = 5, height = 3)
  par(mar = c(3.2, 3.2, 1, 0.4), mgp = c(1.8, 0.5, 0)) 
  plot(data_for_plot[[1]], bty = "l", lty = 1, lwd = 2, xlim = c(0,0.95), ylim= c(0,1), 
       #errbar.col = c("#00468BFF"),
       errbar = F, # ������ƶ������ߣ�����������errbar
       xlab = "Nomogram predicted survival probability",
       ylab = "Actual survival probability",
       riskdist = F, # ȥ������������
       subtitles = F, # ȥ�������Լ����µ�ע��
       par.corrected = list(col = "white")) # ����X(overfitting-corrected estimates)
  lines(data_for_plot[[1]][,c('mean.predicted',"KM")], 
        type = 'o', # (overplotted: lines over points); "b"(both: points and lines)
        lwd = 2, col = colors[1], pch = 16)
  abline(0, 1, lty = 3, lwd = 2, col = c("#224444"))
  # ������
  lines(data_for_plot[[2]][,c('mean.predicted',"KM")], type = 'o', lwd = 2, col = colors[2], pch = 16)
  lines(data_for_plot[[3]][,c('mean.predicted',"KM")], type = 'o', lwd = 2, col = colors[3], pch = 16)
  legend("bottomright", paste0(times, "-Year Survival"), col = colors[1:3], adj = 0, cex = 0.8, lwd = 2, bty = "n", pch = 16)
  dev.off()
  
}
