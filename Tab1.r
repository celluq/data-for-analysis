##### 基因表达与临床因素相关性 -----
rm(list = ls())

# 安装加载需要的R包
library(pacman)
p_load(data.table, tidyverse, expss, table1, survminer)

# 设置工作目录
setwd("D:\\R\\paper")

#### 读取数据----
exp_data<-read.table("clinic_data_clean_fpkm.txt", header = T, sep = "\t", row.names = 1, check.names = F)
rownames(exp_data)<-"symbol"    #将行名定义为symbol
clinic_data<-read.table("clinic_data_clean.txt", header = T, sep = "\t", row.names = 1)
# 整合
data_for_surv = data.frame(exp = as.numeric(exp_data), surv_time = clinic_data$OS.time, surv_status = clinic_data$OS)   #目标基因在各样本中的表达数据，各样本的生存时间，生存状态整合为数据框
# 计算最优分类点
os.cut = surv_cutpoint(data_for_surv, time = "surv_time", event = "surv_status", variables = "exp", minprop = 0.3)
opticut = os.cut$cutpoint$cutpoint # 表达量分界点
os_cat_data = surv_categorize(os.cut, labels = c("Low", "High")) # 高低表达分类

#### 制表准备 ----
# 将要分析的临床变量
vars = c("age_at_initial_pathologic_diagnosis", "menopause_status", 
         "clinical_stage", "neoplasm_histologic_grade","histological_type",
         "diabetes", "BMI", "OS", "DFI")

# 与vars变量中要对应显示的临床变量名
labs = c("Age (year)", "Menopause","FIGO stage", "Histologic grade",
         "Histological type", "Diabetes", "BMI","Living status", "Disease status")

# 构建对应的关系列表〃(格式要求)
labs_list = lapply(labs, function(lab) lab)  #将labs按顺序排列成列表
names(labs_list) = vars   #将labs_list行名与vars对应

# 先看各个变量的取值情况，确定是否需要转化以及如何转化
#table(clinic_data$pathologic_stage) # 查看分期

#  对于分类变量，将Unknown在表中最后一行显示
orderfac = function(x){
    if(is.character(x) | is.factor(x)){
        x = factor(x, exclude = c("",NA))
        if("Unknown" %in% levels(x)) x = fct_relevel(x, "Unknown", after = Inf)
    }
    return(x)
}

# 数据转换
clinic_data_for_tab = clinic_data %>% dplyr::select(vars) %>% 
    mutate(age_at_initial_pathologic_diagnosis = 
               ifelse(age_at_initial_pathologic_diagnosis < 60, "<60", "≥60"),
           menopause_status=
               factor(menopause_status,levels = c("Post (prior bilateral ovariectomy OR >12 mo since LMP with no prior hysterectomy)","Pre (<6 months since LMP AND no prior bilateral ovariectomy AND not on estrogen replacement)"),labels = c("Yes","No")),
           BMI = 
               ifelse(BMI>=30, "Obesity", "Non-Obesity"),
           diabetes =
               factor(diabetes,levels = c("YES","NO"),labels = c("Yes","No")),
           clinical_stage = 
               str_match(clinical_stage, "Stage ([^ABC]*)")[,2],
           neoplasm_histologic_grade = 
               case_when(neoplasm_histologic_grade %in% c("G1") ~ "G1",
                         neoplasm_histologic_grade %in% c("G2") ~ "G2",
                         neoplasm_histologic_grade %in% c("G3","High Grade") ~ "G3"),
           histological_type = 
               case_when(histological_type %in% c("Endometrioid endometrial adenocarcinoma") ~ "EEC",
                         histological_type  %in% c("Serous endometrial adenocarcinoma") ~ "SC",
                         histological_type  %in% c("Mixed serous and endometrioid" ) ~ "Mixed" ),
           OS = 
               ifelse(OS == 1, "Dead", "Alive"),
           DFI = 
               ifelse(DFI == 1, "YES", "NO")) %>% 
    na_if(c("", "[Discrepancy]")) %>% # 空值转化为NA
    if_na("Unknown") %>% # NA值替换为Unknown
    mutate_all(orderfac) #对所有新形成的行进行排序
# 输出
write.table(data.frame(sampleID = clinic_data$sampleID, clinic_data_for_tab), "clinic_data_tab.new.txt", row.names = F, sep = "\t", quote = F)

# 构建数据(临床数据和分类因子)
data_for_table = data.frame(clinic_data_for_tab, os_cat_data$exp)
colnames(data_for_table)[ncol(data_for_table)] = "PHGDH"  #将最后一列的列名改为symbol

# https://cran.r-project.org/web/packages/table1/vignettes/table1-examples.html
# 变量名设置
labels = list(
    # variables = list(age_at_initial_pathologic_diagnosis = "Age (year)",
    #           gender = "Gender", 
    #           relative_family_cancer_history = "Family history of cancer",
    #           pathologic_stage = "TNM stage",...),
    variables = labs_list, # 
    groups = list("", "PHGDH Expression", "")) # 表头大组

# 设置分组
data_for_table[,"PHGDH"]   #发现未分组，还不是factor
data_for_table[,"PHGDH"]<-as.factor(data_for_table[,"PHGDH"]) #转化为factor
data_for_table[,"PHGDH"]

levels(data_for_table[,"PHGDH"]) = c("High","Low","P-value")
strata = c(list(Total = data_for_table), split(data_for_table, data_for_table[,"PHGDH"]))

# 设置函数用于添加p值
rndr = function(x, name, ...) {
    # 将data_for_table1改为自己需要的文件
    if (length(x) == 0) {
        y <- data_for_table[[name]]
        s <- rep("", length(render.default(x=y, name=name, ...)))
        if (is.numeric(y)) {
            p <- t.test(y ~ data_for_table[,"PHGDH"])$p.value  # symbol为用于分组的基因名
        } else {
            p <- chisq.test(table(y, droplevels(data_for_table[,"PHGDH"])))$p.value
        }
        s[2] <- sub("<", "&lt;", format.pval(p, digits=3, eps=0.001))
        s
    } else {
        render.default(x=x, name=name, ...)
    }
}
rndr.strat <- function(label, n, ...) {
    ifelse(n==0, label, render.strat.default(label, n, ...))
}

#显示表格
tab1 = table1(strata, labels, groupspan = c(1, 2, 1), droplevels = F, render = rndr, render.strat = rndr.strat)
