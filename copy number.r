##### ��������CNV�ͼ׻���֮��Ĺ�ϵ -----
rm(list = ls())
library(pacman)
p_load_gh("cBioPortal/cgdsr")
p_load(data.table, tidyverse, ggsci, ggpubr)
setwd("D:\\R\\��ʿ����\\DGAT2\\DGAT2��׻���")

#### ʹ��cgdsr������BioPortal���� ----
# Create CGDS object
mycgds = CGDS("http://www.cbioportal.org/")

# Test the CGDS endpoint URL using a few simple API tests
test(mycgds)
# Retrieve case lists for a specific cancer study
sample_list = getCaseLists(mycgds, "ucec_tcga") 
sample_list[,1:3] # ��ȡ����

#### �������뿽�������� ----
data = getProfileData(mycgds, "PHGDH", c("ucec_tcga_gistic", "ucec_tcga_rna_seq_v2_mrna"), "ucec_tcga_all")
data = data[complete.cases(data),]
data_for_plot = data.frame(exp = log2(data$ucec_tcga_rna_seq_v2_mrna + 1), cnv = factor(data$ucec_tcga_gistic))
levels(data_for_plot$cnv) = c("Shallow Deletion", "Diploid", "Gain", "Amplification")
group = data_for_plot$cnv
ymax = max(data_for_plot$exp)
ymin = min(data_for_plot$exp)
#С����ͼ
p1 = ggplot(data_for_plot, aes(x = group, y = exp, col = group)) +
    geom_violin(trim = FALSE) + geom_boxplot(width = 0.3) + 
    labs(x = "", y = "PHGDH expression(log2)") +
    scale_x_discrete(labels = paste(levels(group), "\n", "(N=", table(group), ")", sep = "")) +
    scale_colour_lancet() +   # ���ñ߿���ɫ
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
    theme(axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 12)) +
    theme(plot.margin = unit(c(0.3, 0.2, -0.3, 0.2), "cm"))

# ��Ⱥͼ
p2 = ggplot(data_for_plot, aes(x = group, y = exp, col = group)) +
    geom_beeswarm(size = 1, cex = 1.5, alpha = 0.5) +  # cex���÷�ɢ�Ŀ�ȣ�alpha���õ��͸����
    ylim(1.1*ymin, 1.3*ymax) + 
    labs(x = "", y = "PHGDH expression(log2)") +    #.(symbol) ������ʾ����
    scale_x_discrete(labels = paste(levels(group), "\n", "(N=", table(group), ")", sep = "")) +
    scale_colour_lancet() + 
    stat_summary(fun = median, geom = "point", size = 1, col = "black") +  # median
    stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", size = 0.3, width = 0.4, col = "black") +  
    stat_summary(fun.data = function(x) median_hilow(x, 0.5), geom = "errorbar", width = 0.25, col = "black") +  
    theme_classic() + theme(legend.position = 0) +
    theme(plot.margin = unit(c(0.2, 0.2, -0.3, 0.2), "cm"))
# ���pֵ
comps = t(combn(levels(group), 2))
my_comparisons = split(comps, 1:nrow(comps))  # ���бȽ���
my_comparisons = list(c("Shallow Deletion","Gain"), c("Shallow Deletion", "Amplification"),c("Diploid", "Gain"))
step = 0.2

p1 + ggsignif::geom_signif(comparisons = my_comparisons,  # Ҳ���������ض��ıȽ���
                          test = "t.test",
                          step_increase = step,
                          map_signif_level = T,
                          margin_top = 0.25,
                          color = "black") + 
    ylim(ymin-0.1, ymax+0.1*(ymax-ymin)+step*(ymax-ymin)*length(my_comparisons)) +
    theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 12))
    
# ����ļ�
ggsave(file = "CNV_vil.pdf", width = 5, height = 3.5)

#### ��������׻��� ----
data = getProfileData(mycgds, "PHGDH", c("ucec_tcga_methylation_hm450", "ucec_tcga_rna_seq_v2_mrna"), "ucec_tcga_all")
data = data[complete.cases(data),]
data_for_plot = data.frame(exp = log2(data$ucec_tcga_rna_seq_v2_mrna + 1), methy = data$ucec_tcga_methylation_hm450)

# ��ͼ
p3 = ggscatter(data_for_plot, x = "methy", y = "exp",
              add = "reg.line",  # Add regressin line
              add.params = list(color = "blue", fill = "gray60"),
              conf.int = TRUE, # Add confidence interval
              ggtheme = theme_bw())
# Add correlation coefficient
p3 + labs(x = "PHGDH Methylation(HM450)", y = "PHGDH expression(log2)") +
    stat_cor(method = "spearman", cor.coef.name = "cor", label.x.npc = 0.5, label.y.npc = 0.9)
ggsave(file = "RNA_Methy_cor.pdf", width = 5, height = 3.5)
