
rm(list = ls())
library(pacman)
p_load(data.table, tidyverse, magrittr, ggpubr)

#### 读取数据 ----
# 目标基因表达数据
rna_exp_data = read.table("PHGDH_fpkm_allsamples.txt", h = T, sep = "\t", check.names = F)
# miRNA
mir_exp_data = fread("UCEC_Portal_miRNA.mature_RPM_jiaozheng.txt", h = T, sep = "\t", check.names = F) %>% 
    filter(str_detect(Symbol, "hsa-miR-27b-3p")) %>%
    column_to_rownames("Symbol")

####数据整理 ----
mir_exp_data_final1<-mir_exp_data_final[10,]
#去除表达量为0的样本
samples_mir = colnames(mir_exp_data_final1)[mir_exp_data_final1 != 0]
samples_rna = colnames(rna_exp_data)[rna_exp_data != 0]

samples = intersect(samples_rna, samples_mir)
data4plot = as.data.frame(t(rbind(rna_exp_data[,samples], mir_exp_data_final1[,samples])))
colnames(data4plot) = c("mrna", "mirna") # 

#### 绘图 ----
setwd(corPath)
p = ggscatter(data4plot, x = "mrna", y = "mirna", 
              add = "reg.line",  # Add regressin line
              add.params = list(color = "blue", fill = "gray60"), # Customize reg. line
              conf.int = TRUE, # Add confidence interval
              ggtheme = theme_bw())
# Add correlation coefficient
p + labs(x = bquote(italic("PHGDH")~expression~(log[2])), y = bquote(italic("hsa-miR-27b-3p")~expression~(log[2]))) +
    stat_cor(method = "pearson", cor.coef.name = "cor", label.x.npc = 0.5, label.y.npc = 1)
ggsave("corplot.pdf", width = 4.5, height = 3.5)
