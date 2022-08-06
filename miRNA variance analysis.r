
rm(list = ls())
library(pacman)
p_load(data.table, tidyverse, edgeR, ggsci)

#### 读入数据 ----
exp_data = fread("UCEC_Portal_miRNA.mature_Counts_jioazheng.txt", h = T, sep = "\t", check.names = F) %>% column_to_rownames("Symbol")
#数据分组
exp_data_final = exp_data %>% dplyr::select(ends_with(c("-01A", "-11A")))
group = if_else(str_detect(colnames(exp_data_final), "-01A"), "EC", "NT") %>% 
    factor(levels = c("EC", "NT"))

# 目标因子
symbol = "hsa-miR-128-3p"
#### 差异表达分析----
setwd(deaPath)
dge = DGEList(counts = exp_data_final, group = group)

# filter
keep = filterByExpr(dge)
dge_filter = dge[keep, keep.lib.sizes = FALSE]
# 
dge_norm = calcNormFactors(dge_filter)  # 
logcpm = cpm(dge_norm, log = TRUE, prior.count = 1)
exp_only_target = t(as.matrix(logcpm[symbol,]))

# Estimating the dispersion
design = model.matrix(~0+group)
colnames(design) = levels(group)
y = estimateDisp(dge_norm, design, robust = TRUE)
fit = glmQLFit(y, design)

res = glmQLFTest(fit, contrast = makeContrasts(EC-NT, levels = design))
result = topTags(res, n = Inf)$table
# 全部因子的结果(用于后续靶基因分析)
write_csv(data.frame(miRNA = rownames(result), result), "miRNA_DEA_result_info.csv")
# # 目标因子的结果
# write_csv(data.frame(miRNA = rownames(result[symbol,]), result[symbol,]), paste0(symbol, "_DEA_result_info.csv"))

#### 画图 ----
data_for_plot = data.frame(exp = as.numeric(exp_only_target), group)

# y轴
ymax = max(data_for_plot$exp)
ymin = min(data_for_plot$exp)
p = ggplot(data_for_plot, aes(x = group, y = exp, col = group)) +
    geom_boxplot(width = 0.6, alpha = 0.6, outlier.shape = NA) + geom_jitter(width = 0.2) + 
    labs(x = "", y = bquote("hsa-miR-27b-3p"~expression~(log[2]))) +
    scale_x_discrete(labels = paste(levels(data_for_plot$group), "\n", "(N=", table(data_for_plot$group), ")", sep = "")) +
    scale_colour_manual(values = pal_lancet("lanonc", alpha = 0.5)(9)[c(3,2)]) +   
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
    theme(axis.title.y = element_text(size = 15), 
          axis.text.x = element_text(size = 14, vjust = 0.5), axis.text.y = element_text(size = 14), 
          panel.border = element_rect(colour = "black", size = 1), axis.ticks.length = unit(0.2, "cm")) + 
    theme(plot.margin = unit(c(0.3, 0.2, -0.3, 0.2), "cm"))

# 添加p
p + ggsignif::geom_signif(comparisons = list(c("EC", "NT")), # 
                          test = NULL, # 
                          annotations = paste0("list(~italic(p.adj)==", signif(result[symbol,5], 3), ")"), 
                          color = "black", # 
                          textsize = 4.5, #
                          parse = T) +
    ylim(ymin, ymax+0.2*(ymax-ymin)) +
    theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 12))

# 输出图
ggsave(file = paste0(symbol, "_boxplot.pdf"), width = 3, height = 4.5)
