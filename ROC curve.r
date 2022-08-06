##### 绘制ROC曲线 -----
rm(list = ls())

library(pacman)
p_load(tidyverse, ROCR, verification)

setwd(dataPath)

#### 数据准备----
# 读入数据
exp_data = read.table("Target_fpkm-uq_allsamples.txt", h = T, sep = "\t", check.names = F, row.names = 1)

# RNA差异表达分析的简易操作
exp_data = exp_data[,str_detect(colnames(exp_data), ".01A|.11A")]
# 反应变量和预测变量（用基因表达来预测肿瘤/正常）
labels = ifelse(grepl(".01A", colnames(exp_data)), 1, 0) # 设置样本标签

#### ROC计算 ----
# 构建prediction对象
predictions = as.numeric(exp_data)
pred = prediction(predictions, labels)
# 计算真阳性率和假阳性率
perf = performance(pred, "tpr", "fpr")
# 计算AUC相关统计量
res = roc.area(labels, predictions)
res$A   # AUC曲线下面积
res$p.value     # 统计量

#### ROC绘图 ----
setwd(rocPath)
pText = ifelse(res$p.value < 0.0001, "< 0.0001", paste0("= ",round(res$p.value, 3))) # p值显示: https://stackoverflow.com/questions/23969726/italics-and-normal-text-in-a-main-plot-title
pdf("ROC.pdf", width = 5, height = 4)
par(mar = c(3.5, 3.5, 1, 1), mgp = c(2, 0.6, 0), las = 1)
plot(perf, col = "#FF0000", lwd = 3, xlab = "1 - Specificity", ylab = "Sensitivity", cex.lab = 1.2, xaxs = "i", yaxs = "i") # xaxs = "i", yaxs = "i"起始点
plot(perf, col = "#FF0000", lwd = 3, xlab = "1 - Specificity", ylab = "Sensitivity", cex.lab = 1.2)
# 绘制对角线
abline(0, 1, lty = 2, lwd = 2, col = "#0026FA")
# 标注文字
text(0.65, 0.2, labels = paste0("AUC = ", round(res$A,3)), cex = 1.1, pos = 4)
text(0.65, 0.1, labels = bquote(italic("P")~.(pText)), cex = 1.1, pos = 4) # 
dev.off()
