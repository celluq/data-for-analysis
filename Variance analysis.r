##### RNA表达数据差异性分析
rm(list = ls())

# 加载需要的R包
library(pacman)
p_load(data.table, tidyverse, edgeR, ggbeeswarm, ggsignif)

# 设置工作目录
setwd("D:\\R\\tcga_R\\DEA_NT")

# 注: 差异表达分析中并不要求肿瘤样本的临床信息

#### 数据准备----
exp_data = read.table("clinic_data_clean_counts.txt", header = T, sep = "\t") %>% 
column_to_rownames("ensembl_gene_id")     #将ensembl_gene_id列作为行名
exp_data_T = exp_data %>% dplyr::select(str_which(colnames(.), ".01A$")) # 匹配列名
nT = ncol(exp_data_T) # 369
exp_data_N = exp_data %>% dplyr::select(str_which(colnames(.), ".11A$"))
nN = ncol(exp_data_N) # 50
exp_data_N[1:3,1:3]

# 合并表达数据
exp_data_final = cbind(exp_data_N, exp_data_T)
dim(exp_data_final) # [1] 共有60483个基因  419个样本
SampleLabel = factor(c(rep("Normal", nN), rep("EC", nT)), levels = c("Normal", "EC"))


# 设定目标基因
target = "ENSG00000062282"
symbol = "PHGDH"

#### edgeR差异表达分析----
# 构建DGEList对象
dge = DGEList(counts = exp_data_final, group = SampleLabel)

# filter
keep = rowSums(cpm(dge) > 0.5) > nN # 可以自由设定，在较少样本组样本中的表达量cpm值大于0.5,如果样本少可以减小此值或者不设此筛查
#keep = filterByExpr(dge)
dge_filter = dge[keep, keep.lib.sizes = FALSE]
dim(dge_filter$counts) # [1] 18672   419

# TMM标准化
dge_norm = calcNormFactors(dge_filter)  # 默认method = "TMM"，此标准化数据不能再用于差异基因的分析

# 注意: 原文是提取pseudo.counts 用于可视化和后续分析, 但edgeR官方不推荐
# 详见https://mp.weixin.qq.com/s?__biz=MzIyNzk1NjUxOA==&mid=2247485070&idx=1&sn=9de874f235aa6e83ef952d71590a9d42&chksm=e85804c7df2f8dd115069af7f509b3723bbe345ce9d7ea5a09f9189a92c186ddc32ed4a98757&token=1957354108&lang=zh_CN#rd
# 获取logCPM,用于可视化分析，包括热图、聚类图、PCA图等
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
result["GAL",] # 目标基因
result["GALR2",]
result["GALR1",]
result["GALR3",]
# 输出
setwd(deaPath)
write_csv(data.frame(Symbol = rownames(result), result), "RNA_DEA_result_info.csv")
# 显著差异的结果
# sig_result = subset(result, abs(logFC) > 1 & FDR < 0.05)
# write_csv(data.frame(Symbol = rownames(sig_result), sig_result, check.names = F), "DEA_result_logFC1_fdr0.05.csv")

#### beeswarm plot ----
# 绘图函数
DataForPlot = data.frame(Exp = logcpm["PHGDH",], Group = SampleLabel)
# y轴范围
ymax = max(DataForPlot$Exp)
ymin = min(DataForPlot$Exp)
# p.adj值
fdrvalue = result["PHGDH",5] 
# 绘图
p = ggplot(DataForPlot, aes(x = Group, y = Exp, col = Group)) +
    geom_beeswarm(size = 1, cex = 1.5, alpha = 0.5) +  # cex设置分散的宽度，alpha设置点的透明度
    ylim(1.1*ymin, 1.3*ymax) + 
    labs(x = "", y = bquote(italic(.(symbol))~expression~(log[2]))) +    #.(symbol) 才能提示变量
    scale_x_discrete(labels = paste(levels(SampleLabel), "\n", "(N=", table(SampleLabel), ")", sep = "")) +
    scale_colour_manual(values = c("#2700F6", "#C8181C")) + 
    stat_summary(fun = median, geom = "point", size = 1, col = "black") +  # median
    stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", size = 0.3, width = 0.4, col = "black") +  
    stat_summary(fun.data = function(x) median_hilow(x, 0.5), geom = "errorbar", width = 0.25, col = "black") +  
    theme_classic() + theme(legend.position = 0) +
    theme(plot.margin = unit(c(0.2, 0.2, -0.3, 0.2), "cm"))
# 添加显著性p值
p + geom_signif(comparisons = list(c("NO", "YES")), y_position = 1.15*ymax, textsize = 3, annotation = paste0("list(~italic(p.adj)==", signif(fdrvalue, 3), ")"), parse = T, color = "black") 
ggsave(paste0(symbol, "_beeswarm.pdf"), width = 2, height = 3)

