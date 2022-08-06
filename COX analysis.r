##### COX������-�����ط��� -----
rm(list = ls())
library(pacman)
p_load(gtsummary, gt, tidyverse, finalfit, survival)
setwd("D:\\R\\��ϰ\\file1")
#### ��ȡ���� ----
exp_data = read.table("clinic_data_clean_fpkm.txt", header = T, sep = "\t", check.names = F)
symbol = rownames(exp_data)
# �Ѿ���table1�������������
clinic_data_tabl = read.table("clinic_data_tab.txt", header = T, sep = "\t") %>% 
  dplyr::select(-c("OS", "DFI"))
# ��ȡ��������
clinic_data_final = read.table("clinic_data_clean.txt", header = T, sep = "\t")
#### COX----
# ���庯��
do_cox = function(surv_type){
  # surv_type = "OS"
  clinic_data = clinic_data_final %>% dplyr::select(sampleID, all_of(surv_type), paste0(surv_type, ".time")) %>% 
    left_join(clinic_data_tabl, ., by = "sampleID") %>% 
    rename("surv_status" := !!surv_type, "surv_time" := !!paste0(surv_type, ".time")) # ������
  
  # ����
  data_for_cox = clinic_data %>% mutate(exp = as.numeric(exp_data_hubgene[1,])) %>% column_to_rownames("sampleID")
  
  # ������������ĺ���
  data_clean = function(data){
    data %>% 
      filter(str_detect(explanatory, "Unknown", negate = T)) %>% 
      mutate_at(c("HR", "L95", "U95"), ~round(.,2)) %>% # ������λ��Ч����
      mutate(HR = ifelse(is.na(HR), "-", paste0(HR, "(", L95, ",", U95, ")")), # ����HR 95CI��һ��
             p = ifelse(is.na(p), "-", p_tidy(p, 3, prefix = ""))) %>% # pֵ��ʽ��
      dplyr::select(explanatory, HR, p) %>% # ��ȡ�ٴ����ء�HR��P
      rename("HR (95% CI)" = "HR")
  }
  
  #### COX�����ط��� ----
  # https://finalfit.org/reference/coxphuni.html
  explanatory_uni = setdiff(colnames(data_for_cox), c("surv_status", "surv_time"))  #����Ҫ���е����ط����ı���
  dependent = "Surv(surv_time, surv_status)"  #����y
  coxphuni_res = data_for_cox %>%
    coxphuni(dependent, explanatory_uni) %>%
    fit2df(condense = F) %>% # ������ϲ�
    data_clean() 
  final<-coxphuni_res[12,]

#ѭ��
for (i in c(2:nrow(exp_data))) 
    {data_for_cox = clinic_data %>% mutate(exp = as.numeric(exp_data[i,])) %>% column_to_rownames("sampleID")
    explanatory_uni = setdiff(colnames(data_for_cox), c("surv_status", "surv_time"))  #����Ҫ���е����ط����ı���
    dependent = "Surv(surv_time, surv_status)"  #����y
    coxphuni_res = data_for_cox %>%
      coxphuni(dependent, explanatory_uni) %>%
      fit2df(condense = F) %>% # ������ϲ�
      data_clean() 
    final<-final%>% add_row(coxphuni_res[12,])
  }
  # ѡ��Ҫ��������ط���������
  explanatory_multi = coxphuni_res %>% filter(p<0.1) %>% pull(explanatory) %>% str_extract(., paste0(explanatory_uni, collapse = "|")) %>% unique

  
  #### COX�����ط��� ----
  # https://finalfit.org/reference/coxphmulti.html
  data_for_cox = clinic_data %>% mutate(exp = as.numeric(exp_sig[1,])) %>% column_to_rownames("sampleID")
  explanatory_multi=c("clinical_stage","neoplasm_histologic_grade","histological_type","exp")
  coxphmulti_res_old = clinic_data_hub %>%
    coxphmulti(dependent, explanatory_multi) %>%
    fit2df(condense = F) %>% # ������ϲ�
    data_clean()
  final_mul<-coxphmulti_res[8,]
  
  #ѭ��
  for (i in c(2:nrow(exp_sig))) 
  {data_for_cox = clinic_data %>% mutate(exp = as.numeric(exp_sig[i,])) %>% column_to_rownames("sampleID")
  explanatory_multi=c("clinical_stage","neoplasm_histologic_grade","histological_type","exp")
  coxphmulti_res = data_for_cox %>%
    coxphmulti(dependent, explanatory_multi) %>%
    fit2df(condense = F) %>% # ������ϲ�
    data_clean()
  final_mul<-final_mul%>% add_row(coxphmulti_res[8,])
  } 
  #����ROC
  multi_COX<-coxph(Surv(surv_time, surv_status) ~ GAL + CXCR3 + clinical_stage + neoplasm_histologic_grade, data=clinic_data_hub)
  RiskScore<-predict(multi_COX,type = "risk",newdata =clinic_data_hub)
  risk_group<-ifelse(RiskScore>=median(RiskScore),'high','low')
  #����һ���µ����ݿ�
  KM.input<-data.frame(clinic_data_hub,RiskScore=RiskScore,risk_group=risk_group)
  library(survminer)
  Survival_ROC_input<-KM.input
  survival_ROC<-survivalROC(Stime=Survival_ROC_input$surv_time, #����ʱ�䣬Event time or censoring time for subjects
                            status=Survival_ROC_input$surv_status, #����״̬,dead or alive
                            marker=Survival_ROC_input$RiskScore, #���յ÷֣�Predictor or marker value
                            predict.time=5*365, #Ԥ��5�������ʱ��
                            method="KM") #ʹ��KM��������ϣ�Ĭ�ϵķ�����method="NNE")
  survival_ROC_AUC<-round(survival_ROC$AUC,3)
  pdf("ROC.pdf", width = 5, height = 4)
  par(mar = c(3.5, 3.5, 1, 1), mgp = c(2, 0.6, 0), las = 1)
  plot(survival_ROC$FP,survival_ROC$TP,type="l",xlim=c(0,1),ylim=c(0,1),
       xlab="False positive rate",  
       ylab="True positive rate",
       main=paste0("ROC Curve", " (", "AUC = ",survival_ROC_AUC," )"),  #����
       cex.main=1.5,#�Ŵ����
       cex.lab=1.3,#�������ǩ�����ƣ������ű���
       cex.axis=1.3, font=1.2, #�Ŵ����ϵ�����
       lwd=1.5, #ָ���������
       col="red")  #��ɫ)
       abline(0, 1, lty = 2, lwd = 2, col = "#0026FA")
       dev.off()
  
  #��һ�ֻ�ͼ����
       time_ROC_input<-KM.input
       time_ROC<-timeROC(T=time_ROC_input$surv_time, #����ʱ��(dead��alive������ʱ��).
                         delta=time_ROC_input$surv_status, #�����֣�Censored������������0��ʾ
                         marker=time_ROC_input$RiskScore, #Ԥ��ı����������Ƿ������֣���û����������£�ֵԽ��Խ���׷���Σ���¼�
                         cause=1, #���Խ�ֵĸ�ֵ��������1��2����Ҳ����dead�ĸ�ֵ������dead��1��ʾ��
                         weighting = "marginal", #Ȩ�ؼ��㷽��������Ĭ�Ϸ�����ѡ��KM����ɾʧ�ֲ���weighting="aalen" [ѡ��COX]��weighting="aalen" [ѡ��Aalen]
                         times = c(3*365,5*365,8*365), #����3��5��10���ROC����
                         ROC=TRUE,
                         iid=TRUE #����AUC
       )
       time_ROC #�鿴��������Կ�������������SE
       #����ROC������
       summary(time_ROC) #����12������
       time_ROC$AUC
       library(ggplot2)
       time_ROC$TP
       summary(time_ROC)
       #��������ݿ�gpplot��Ҫ���ݿ�
       time_ROC.res<-data.frame(TP_3year=time_ROC$TP[,1], #��ȡ3���ROC��TP
                                FP_3year=time_ROC$FP[,1],  #��ȡ3���ROC��FP
                                TP_5year=time_ROC$TP[,2],  #��ȡ5���ROC��TP
                                FP_5year=time_ROC$FP[,2], #��ȡ5���ROC��FP
                                TP_8year=time_ROC$TP[,3], #��ȡ10���ROC��TP
                                FP_8year=time_ROC$FP[,3]) #��ȡ10���ROC��FP
       time_ROC$AUC #�������3��5��10���AUC�������ֱ�ȡ�Ӽ�time_ROC$AUC[[1]]��time_ROC$AUC[[2]]��time_ROC$AUC[[3]]
       TimeROC_plot<-ggplot()+
         geom_line(data=time_ROC.res,aes(x=FP_3year,y=TP_3year),size=1,color="red")+
         geom_line(data=time_ROC.res,aes(x=FP_5year,y=TP_5year),size=1,color="blue")+
         geom_line(data=time_ROC.res,aes(x=FP_8year,y=TP_8year),size=1,color="black")+
         geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey",size=1, linetype = 2 #�������
         )+
         theme_bw()+
         annotate("text",x=0.75,y=0.25,size=4.5,
                  label=paste0("AUC at 3 years = ",round(time_ROC$AUC[[1]],3)),color="red")+
         annotate("text",x=0.75,y=0.15,size=4.5,
                  label=paste0("AUC at 5 years = ",round(time_ROC$AUC[[2]],3)),color="blue")+
         annotate("text",x=0.75,y=0.05,size=4.5,
                  label=paste0("AUC at 10 years = ",round(time_ROC$AUC[[3]],3)),color="black")+
         labs(x="False positive rate",y="True positive rate")+
         theme(axis.text=element_text(face="bold", size=11,  color="black"),#�Ӵ̶ֿȱ�ǩ
               axis.title=element_text(face="bold", size=14, color="black")) #�Ӵ�xy���ǩ����
       TimeROC_plot  #�ǳ�����
  # finalfit�?
  # fmla = as.formula(paste0("Surv(surv_time, surv_status) ~",paste0(explanatory_uni, collapse = '+')))
  # fit = coxph(fmla, data = data_for_cox)
  # fit %>% tbl_regression(exp = T)
  
  #### ���չʾ ----
  # �����Ҫ��ʾ���ٴ���������labs��Ҫ��coxphuni_res��һ�£�
  labs = c("Age (>=60 vs. < 60)", "Menopause (Yes vs. No)", 
           "FIGO stage (II vs. I)", "FIGO stage (III vs. I)", "FIGO stage (IV vs. I)", "Histologic grade (G2 vs. G1)", "Histologic grade (G3 vs. G1)",
           "Histological Type (Mixed vs. EEC)","Histological Type (SC vs. EEC)", "Diabetes (Yes vs. No)", "BMI (>=30 vs. < 30)", 
           "PHGDH (high vs. low)")
  
  cox_res = full_join(coxphuni_res, coxphmulti_res, by = "explanatory") %>% 
    mutate(explanatory = labs) %>% 
    rename("Variables" = "explanatory") %>% 
    expss::if_na(.,"-")
  
  # �������
  # https://gt.rstudio.com/index.html
  cox_res %>% gt() %>% 
    tab_spanner(label = md("**Univariate analysis**"), columns = ends_with(".x")) %>%
    tab_spanner(label = md("**Multivariate analysis**"), columns = ends_with(".y")) %>%
    # ������ʾ������
    cols_label("HR (95% CI).x" = "HR (95% CI)", "p.x" = md("***p***"), 
               "HR (95% CI).y" = "HR (95% CI)", "p.y" = md("***p***")) %>%
    tab_style(style = list(cell_fill(color = "white"), # �ױ���
                           cell_borders(color = "white")), # ������ޱ߿�
              locations = cells_body()) %>%
    tab_style(style = cell_text(weight = "bold"), # ����pֵ�Ӵ�
              locations = list(cells_body(columns = gt::vars(p.x), rows = p.x < 0.05 & p.x != "-"), # 如果"Error: The object to `data` is not a `gt_tbl` object", 则在vars函数前加包名为gt::vars
                               cells_body(columns = gt::vars(p.y), rows = p.y < 0.05 & p.y != "-"))) %>%
    gtsave(paste0(surv_type, "_cox.html")) # ���
  
}

# ����
do_cox("OS") 
do_cox("DFI")
