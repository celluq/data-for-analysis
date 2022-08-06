##### ���������ٴ���������� -----
rm(list = ls())

# ��װ������Ҫ��R��
library(pacman)
p_load(data.table, tidyverse, expss, table1, survminer)

# ���ù���Ŀ¼
setwd("D:\\R\\paper")

#### ��ȡ����----
exp_data<-read.table("clinic_data_clean_fpkm.txt", header = T, sep = "\t", row.names = 1, check.names = F)
rownames(exp_data)<-"symbol"    #����������Ϊsymbol
clinic_data<-read.table("clinic_data_clean.txt", header = T, sep = "\t", row.names = 1)
# ����
data_for_surv = data.frame(exp = as.numeric(exp_data), surv_time = clinic_data$OS.time, surv_status = clinic_data$OS)   #Ŀ������ڸ������еı�����ݣ�������������ʱ�䣬����״̬����Ϊ���ݿ�
# �������ŷ����
os.cut = surv_cutpoint(data_for_surv, time = "surv_time", event = "surv_status", variables = "exp", minprop = 0.3)
opticut = os.cut$cutpoint$cutpoint # ������ֽ��
os_cat_data = surv_categorize(os.cut, labels = c("Low", "High")) # �ߵͱ�����

#### �Ʊ�׼�� ----
# ��Ҫ�������ٴ�����
vars = c("age_at_initial_pathologic_diagnosis", "menopause_status", 
         "clinical_stage", "neoplasm_histologic_grade","histological_type",
         "diabetes", "BMI", "OS", "DFI")

# ��vars������Ҫ��Ӧ��ʾ���ٴ�������
labs = c("Age (year)", "Menopause","FIGO stage", "Histologic grade",
         "Histological type", "Diabetes", "BMI","Living status", "Disease status")

# ������Ӧ�Ĺ�ϵ�б�(��ʽҪ��)
labs_list = lapply(labs, function(lab) lab)  #��labs��˳�����г��б�
names(labs_list) = vars   #��labs_list������vars��Ӧ

# �ȿ�����������ȡֵ�����ȷ���Ƿ���Ҫת���Լ����ת��
#table(clinic_data$pathologic_stage) # �鿴����

#  ���ڷ����������Unknown�ڱ������һ����ʾ
orderfac = function(x){
    if(is.character(x) | is.factor(x)){
        x = factor(x, exclude = c("",NA))
        if("Unknown" %in% levels(x)) x = fct_relevel(x, "Unknown", after = Inf)
    }
    return(x)
}

# ����ת��
clinic_data_for_tab = clinic_data %>% dplyr::select(vars) %>% 
    mutate(age_at_initial_pathologic_diagnosis = 
               ifelse(age_at_initial_pathologic_diagnosis < 60, "<60", "��60"),
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
    na_if(c("", "[Discrepancy]")) %>% # ��ֵת��ΪNA
    if_na("Unknown") %>% # NAֵ�滻ΪUnknown
    mutate_all(orderfac) #���������γɵ��н�������
# ���
write.table(data.frame(sampleID = clinic_data$sampleID, clinic_data_for_tab), "clinic_data_tab.new.txt", row.names = F, sep = "\t", quote = F)

# ��������(�ٴ����ݺͷ�������)
data_for_table = data.frame(clinic_data_for_tab, os_cat_data$exp)
colnames(data_for_table)[ncol(data_for_table)] = "PHGDH"  #�����һ�е�������Ϊsymbol

# https://cran.r-project.org/web/packages/table1/vignettes/table1-examples.html
# ����������
labels = list(
    # variables = list(age_at_initial_pathologic_diagnosis = "Age (year)",
    #           gender = "Gender", 
    #           relative_family_cancer_history = "Family history of cancer",
    #           pathologic_stage = "TNM stage",...),
    variables = labs_list, # 
    groups = list("", "PHGDH Expression", "")) # ��ͷ����

# ���÷���
data_for_table[,"PHGDH"]   #����δ���飬������factor
data_for_table[,"PHGDH"]<-as.factor(data_for_table[,"PHGDH"]) #ת��Ϊfactor
data_for_table[,"PHGDH"]

levels(data_for_table[,"PHGDH"]) = c("High","Low","P-value")
strata = c(list(Total = data_for_table), split(data_for_table, data_for_table[,"PHGDH"]))

# ���ú����������pֵ
rndr = function(x, name, ...) {
    # ��data_for_table1��Ϊ�Լ���Ҫ���ļ�
    if (length(x) == 0) {
        y <- data_for_table[[name]]
        s <- rep("", length(render.default(x=y, name=name, ...)))
        if (is.numeric(y)) {
            p <- t.test(y ~ data_for_table[,"PHGDH"])$p.value  # symbolΪ���ڷ���Ļ�����
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

#��ʾ���
tab1 = table1(strata, labels, groupspan = c(1, 2, 1), droplevels = F, render = rndr, render.strat = rndr.strat)
