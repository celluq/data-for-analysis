##### COX单因素-多因素分析 -----
rm(list = ls())
library(pacman)
p_load(gtsummary, gt, tidyverse, finalfit, survival)
setwd("D:\\R\\练习\\file1")
#### 读取数据 ----
exp_data = read.table("clinic_data_clean_fpkm.txt", header = T, sep = "\t", check.names = F)
symbol = rownames(exp_data)
# 已经在table1中整理过的数据
clinic_data_tabl = read.table("clinic_data_tab.txt", header = T, sep = "\t") %>% 
  dplyr::select(-c("OS", "DFI"))
# 提取生存数据
clinic_data_final = read.table("clinic_data_clean.txt", header = T, sep = "\t")
#### COX----
# 定义函数
do_cox = function(surv_type){
  # surv_type = "OS"
  clinic_data = clinic_data_final %>% dplyr::select(sampleID, all_of(surv_type), paste0(surv_type, ".time")) %>% 
    left_join(clinic_data_tabl, ., by = "sampleID") %>% 
    rename("surv_status" := !!surv_type, "surv_time" := !!paste0(surv_type, ".time")) # 重命名
  
  # 整合
  data_for_cox = clinic_data %>% mutate(exp = as.numeric(exp_data_hubgene[1,])) %>% column_to_rownames("sampleID")
  
  # 定义数据整理的函数
  data_clean = function(data){
    data %>% 
      filter(str_detect(explanatory, "Unknown", negate = T)) %>% 
      mutate_at(c("HR", "L95", "U95"), ~round(.,2)) %>% # 保留两位有效数字
      mutate(HR = ifelse(is.na(HR), "-", paste0(HR, "(", L95, ",", U95, ")")), # 整合HR 95CI于一列
             p = ifelse(is.na(p), "-", p_tidy(p, 3, prefix = ""))) %>% # p值格式化
      dplyr::select(explanatory, HR, p) %>% # 提取临床因素、HR、P
      rename("HR (95% CI)" = "HR")
  }
  
  #### COX单因素分析 ----
  # https://finalfit.org/reference/coxphuni.html
  explanatory_uni = setdiff(colnames(data_for_cox), c("surv_status", "surv_time"))  #设置要进行单因素分析的变量
  dependent = "Surv(surv_time, surv_status)"  #设置y
  coxphuni_res = data_for_cox %>%
    coxphuni(dependent, explanatory_uni) %>%
    fit2df(condense = F) %>% # 结果不合并
    data_clean() 
  final<-coxphuni_res[12,]

#循环
for (i in c(2:nrow(exp_data))) 
    {data_for_cox = clinic_data %>% mutate(exp = as.numeric(exp_data[i,])) %>% column_to_rownames("sampleID")
    explanatory_uni = setdiff(colnames(data_for_cox), c("surv_status", "surv_time"))  #设置要进行单因素分析的变量
    dependent = "Surv(surv_time, surv_status)"  #设置y
    coxphuni_res = data_for_cox %>%
      coxphuni(dependent, explanatory_uni) %>%
      fit2df(condense = F) %>% # 结果不合并
      data_clean() 
    final<-final%>% add_row(coxphuni_res[12,])
  }
  # 选择要进入多因素分析的因素
  explanatory_multi = coxphuni_res %>% filter(p<0.1) %>% pull(explanatory) %>% str_extract(., paste0(explanatory_uni, collapse = "|")) %>% unique

  
  #### COX多因素分析 ----
  # https://finalfit.org/reference/coxphmulti.html
  data_for_cox = clinic_data %>% mutate(exp = as.numeric(exp_sig[1,])) %>% column_to_rownames("sampleID")
  explanatory_multi=c("clinical_stage","neoplasm_histologic_grade","histological_type","exp")
  coxphmulti_res_old = clinic_data_hub %>%
    coxphmulti(dependent, explanatory_multi) %>%
    fit2df(condense = F) %>% # 结果不合并
    data_clean()
  final_mul<-coxphmulti_res[8,]
  
  #循环
  for (i in c(2:nrow(exp_sig))) 
  {data_for_cox = clinic_data %>% mutate(exp = as.numeric(exp_sig[i,])) %>% column_to_rownames("sampleID")
  explanatory_multi=c("clinical_stage","neoplasm_histologic_grade","histological_type","exp")
  coxphmulti_res = data_for_cox %>%
    coxphmulti(dependent, explanatory_multi) %>%
    fit2df(condense = F) %>% # 结果不合并
    data_clean()
  final_mul<-final_mul%>% add_row(coxphmulti_res[8,])
  } 
  #绘制ROC
  multi_COX<-coxph(Surv(surv_time, surv_status) ~ GAL + CXCR3 + clinical_stage + neoplasm_histologic_grade, data=clinic_data_hub)
  RiskScore<-predict(multi_COX,type = "risk",newdata =clinic_data_hub)
  risk_group<-ifelse(RiskScore>=median(RiskScore),'high','low')
  #构建一个新的数据框
  KM.input<-data.frame(clinic_data_hub,RiskScore=RiskScore,risk_group=risk_group)
  library(survminer)
  Survival_ROC_input<-KM.input
  survival_ROC<-survivalROC(Stime=Survival_ROC_input$surv_time, #生存时间，Event time or censoring time for subjects
                            status=Survival_ROC_input$surv_status, #生存状态,dead or alive
                            marker=Survival_ROC_input$RiskScore, #风险得分，Predictor or marker value
                            predict.time=5*365, #预测5年的生存时间
                            method="KM") #使用KM法进行拟合，默认的方法是method="NNE")
  survival_ROC_AUC<-round(survival_ROC$AUC,3)
  pdf("ROC.pdf", width = 5, height = 4)
  par(mar = c(3.5, 3.5, 1, 1), mgp = c(2, 0.6, 0), las = 1)
  plot(survival_ROC$FP,survival_ROC$TP,type="l",xlim=c(0,1),ylim=c(0,1),
       xlab="False positive rate",  
       ylab="True positive rate",
       main=paste0("ROC Curve", " (", "AUC = ",survival_ROC_AUC," )"),  #标题
       cex.main=1.5,#放大标题
       cex.lab=1.3,#坐标轴标签（名称）的缩放倍数
       cex.axis=1.3, font=1.2, #放大轴上的数字
       lwd=1.5, #指定线条宽度
       col="red")  #红色)
       abline(0, 1, lty = 2, lwd = 2, col = "#0026FA")
       dev.off()
  
  #另一种绘图方法
       time_ROC_input<-KM.input
       time_ROC<-timeROC(T=time_ROC_input$surv_time, #生存时间(dead和alive的生存时间).
                         delta=time_ROC_input$surv_status, #生存结局，Censored的样本必须用0表示
                         marker=time_ROC_input$RiskScore, #预测的变量，这里是风险评分，在没有特殊情况下，值越大，越容易发生危险事件
                         cause=1, #阳性结局的赋值（必须是1或2），也就是dead的赋值，这里dead是1表示的
                         weighting = "marginal", #权重计算方法，这是默认方法，选择KM计算删失分布，weighting="aalen" [选用COX]，weighting="aalen" [选用Aalen]
                         times = c(3*365,5*365,8*365), #计算3、5、10年的ROC曲线
                         ROC=TRUE,
                         iid=TRUE #计算AUC
       )
       time_ROC #查看结果，可以看到，还包括了SE
       #绘制ROC曲线啦
       summary(time_ROC) #返回12个参数
       time_ROC$AUC
       library(ggplot2)
       time_ROC$TP
       summary(time_ROC)
       #整理成数据框，gpplot需要数据框
       time_ROC.res<-data.frame(TP_3year=time_ROC$TP[,1], #获取3年的ROC的TP
                                FP_3year=time_ROC$FP[,1],  #获取3年的ROC的FP
                                TP_5year=time_ROC$TP[,2],  #获取5年的ROC的TP
                                FP_5year=time_ROC$FP[,2], #获取5年的ROC的FP
                                TP_8year=time_ROC$TP[,3], #获取10年的ROC的TP
                                FP_8year=time_ROC$FP[,3]) #获取10年的ROC的FP
       time_ROC$AUC #这里放了3，5，10年的AUC，后续分别取子集time_ROC$AUC[[1]]，time_ROC$AUC[[2]]，time_ROC$AUC[[3]]
       TimeROC_plot<-ggplot()+
         geom_line(data=time_ROC.res,aes(x=FP_3year,y=TP_3year),size=1,color="red")+
         geom_line(data=time_ROC.res,aes(x=FP_5year,y=TP_5year),size=1,color="blue")+
         geom_line(data=time_ROC.res,aes(x=FP_8year,y=TP_8year),size=1,color="black")+
         geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey",size=1, linetype = 2 #添加虚线
         )+
         theme_bw()+
         annotate("text",x=0.75,y=0.25,size=4.5,
                  label=paste0("AUC at 3 years = ",round(time_ROC$AUC[[1]],3)),color="red")+
         annotate("text",x=0.75,y=0.15,size=4.5,
                  label=paste0("AUC at 5 years = ",round(time_ROC$AUC[[2]],3)),color="blue")+
         annotate("text",x=0.75,y=0.05,size=4.5,
                  label=paste0("AUC at 10 years = ",round(time_ROC$AUC[[3]],3)),color="black")+
         labs(x="False positive rate",y="True positive rate")+
         theme(axis.text=element_text(face="bold", size=11,  color="black"),#加粗刻度标签
               axis.title=element_text(face="bold", size=14, color="black")) #加粗xy轴标签名字
       TimeROC_plot  #非常美丽
  # finalfit鍖?
  # fmla = as.formula(paste0("Surv(surv_time, surv_status) ~",paste0(explanatory_uni, collapse = '+')))
  # fit = coxph(fmla, data = data_for_cox)
  # fit %>% tbl_regression(exp = T)
  
  #### 结果展示 ----
  # 表格中要显示的临床变量名（labs中要和coxphuni_res中一致）
  labs = c("Age (>=60 vs. < 60)", "Menopause (Yes vs. No)", 
           "FIGO stage (II vs. I)", "FIGO stage (III vs. I)", "FIGO stage (IV vs. I)", "Histologic grade (G2 vs. G1)", "Histologic grade (G3 vs. G1)",
           "Histological Type (Mixed vs. EEC)","Histological Type (SC vs. EEC)", "Diabetes (Yes vs. No)", "BMI (>=30 vs. < 30)", 
           "PHGDH (high vs. low)")
  
  cox_res = full_join(coxphuni_res, coxphmulti_res, by = "explanatory") %>% 
    mutate(explanatory = labs) %>% 
    rename("Variables" = "explanatory") %>% 
    expss::if_na(.,"-")
  
  # 制作表格
  # https://gt.rstudio.com/index.html
  cox_res %>% gt() %>% 
    tab_spanner(label = md("**Univariate analysis**"), columns = ends_with(".x")) %>%
    tab_spanner(label = md("**Multivariate analysis**"), columns = ends_with(".y")) %>%
    # 更改显示的列名
    cols_label("HR (95% CI).x" = "HR (95% CI)", "p.x" = md("***p***"), 
               "HR (95% CI).y" = "HR (95% CI)", "p.y" = md("***p***")) %>%
    tab_style(style = list(cell_fill(color = "white"), # 白背景
                           cell_borders(color = "white")), # 表格区无边框
              locations = cells_body()) %>%
    tab_style(style = cell_text(weight = "bold"), # 显著p值加粗
              locations = list(cells_body(columns = gt::vars(p.x), rows = p.x < 0.05 & p.x != "-"), # 濡傛灉"Error: The object to `data` is not a `gt_tbl` object", 鍒欏湪vars鍑芥暟鍓嶅姞鍖呭悕涓篻t::vars
                               cells_body(columns = gt::vars(p.y), rows = p.y < 0.05 & p.y != "-"))) %>%
    gtsave(paste0(surv_type, "_cox.html")) # 输出
  
}

# 调用
do_cox("OS") 
do_cox("DFI")
