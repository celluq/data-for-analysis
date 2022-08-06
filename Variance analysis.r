##### RNA������ݲ����Է���
rm(list = ls())

# ������Ҫ��R��
library(pacman)
p_load(data.table, tidyverse, edgeR, ggbeeswarm, ggsignif)

# ���ù���Ŀ¼
setwd("D:\\R\\tcga_R\\DEA_NT")

# ע: ����������в���Ҫ�������������ٴ���Ϣ

#### ����׼��----
exp_data = read.table("clinic_data_clean_counts.txt", header = T, sep = "\t") %>% 
column_to_rownames("ensembl_gene_id")     #��ensembl_gene_id����Ϊ����
exp_data_T = exp_data %>% dplyr::select(str_which(colnames(.), ".01A$")) # ƥ������
nT = ncol(exp_data_T) # 369
exp_data_N = exp_data %>% dplyr::select(str_which(colnames(.), ".11A$"))
nN = ncol(exp_data_N) # 50
exp_data_N[1:3,1:3]

# �ϲ��������
exp_data_final = cbind(exp_data_N, exp_data_T)
dim(exp_data_final) # [1] ����60483������  419������
SampleLabel = factor(c(rep("Normal", nN), rep("EC", nT)), levels = c("Normal", "EC"))


# �趨Ŀ�����
target = "ENSG00000062282"
symbol = "PHGDH"

#### edgeR���������----
# ����DGEList����
dge = DGEList(counts = exp_data_final, group = SampleLabel)

# filter
keep = rowSums(cpm(dge) > 0.5) > nN # ���������趨���ڽ��������������еı����cpmֵ����0.5,��������ٿ��Լ�С��ֵ���߲����ɸ��
#keep = filterByExpr(dge)
dge_filter = dge[keep, keep.lib.sizes = FALSE]
dim(dge_filter$counts) # [1] 18672   419

# TMM��׼��
dge_norm = calcNormFactors(dge_filter)  # Ĭ��method = "TMM"���˱�׼�����ݲ��������ڲ������ķ���

# ע��: ԭ������ȡpseudo.counts ���ڿ��ӻ��ͺ�������, ��edgeR�ٷ����Ƽ�
# ���https://mp.weixin.qq.com/s?__biz=MzIyNzk1NjUxOA==&mid=2247485070&idx=1&sn=9de874f235aa6e83ef952d71590a9d42&chksm=e85804c7df2f8dd115069af7f509b3723bbe345ce9d7ea5a09f9189a92c186ddc32ed4a98757&token=1957354108&lang=zh_CN#rd
# ��ȡlogCPM,���ڿ��ӻ�������������ͼ������ͼ��PCAͼ��
logcpm = cpm(dge_norm, log = TRUE, prior.count = 1)
# exp_only_target = t(as.matrix(logcpm[target,]))
# write.table(data.frame(Symbol = symbol, exp_only_target, check.names = F), "logcpm_data.txt", row.names = F, sep = "\t", quote = F)

# Estimating the dispersion
design = model.matrix(~SampleLabel)
y = estimateDisp(dge_norm, design, robust = TRUE)

# DEA
fit = glmQLFit(y, design)
res = glmQLFTest(fit)
result = topTags(res, n = Inf)$table
result["GAL",] # Ŀ�����
result["GALR2",]
result["GALR1",]
result["GALR3",]
# ���
setwd(deaPath)
write_csv(data.frame(Symbol = rownames(result), result), "RNA_DEA_result_info.csv")
# ��������Ľ��
# sig_result = subset(result, abs(logFC) > 1 & FDR < 0.05)
# write_csv(data.frame(Symbol = rownames(sig_result), sig_result, check.names = F), "DEA_result_logFC1_fdr0.05.csv")

#### beeswarm plot ----
# ��ͼ����
DataForPlot = data.frame(Exp = logcpm["PHGDH",], Group = SampleLabel)
# y�᷶Χ
ymax = max(DataForPlot$Exp)
ymin = min(DataForPlot$Exp)
# p.adjֵ
fdrvalue = result["PHGDH",5] 
# ��ͼ
p = ggplot(DataForPlot, aes(x = Group, y = Exp, col = Group)) +
    geom_beeswarm(size = 1, cex = 1.5, alpha = 0.5) +  # cex���÷�ɢ�Ŀ�ȣ�alpha���õ��͸����
    ylim(1.1*ymin, 1.3*ymax) + 
    labs(x = "", y = bquote(italic(.(symbol))~expression~(log[2]))) +    #.(symbol) ������ʾ����
    scale_x_discrete(labels = paste(levels(SampleLabel), "\n", "(N=", table(SampleLabel), ")", sep = "")) +
    scale_colour_manual(values = c("#2700F6", "#C8181C")) + 
    stat_summary(fun = median, geom = "point", size = 1, col = "black") +  # median
    stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", size = 0.3, width = 0.4, col = "black") +  
    stat_summary(fun.data = function(x) median_hilow(x, 0.5), geom = "errorbar", width = 0.25, col = "black") +  
    theme_classic() + theme(legend.position = 0) +
    theme(plot.margin = unit(c(0.2, 0.2, -0.3, 0.2), "cm"))
# ���������pֵ
p + geom_signif(comparisons = list(c("NO", "YES")), y_position = 1.15*ymax, textsize = 3, annotation = paste0("list(~italic(p.adj)==", signif(fdrvalue, 3), ")"), parse = T, color = "black") 
ggsave(paste0(symbol, "_beeswarm.pdf"), width = 2, height = 3)

