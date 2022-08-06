##### ����ROC���� -----
rm(list = ls())

library(pacman)
p_load(tidyverse, ROCR, verification)

setwd(dataPath)

#### ����׼��----
# ��������
exp_data = read.table("Target_fpkm-uq_allsamples.txt", h = T, sep = "\t", check.names = F, row.names = 1)

# RNA����������ļ��ײ���
exp_data = exp_data[,str_detect(colnames(exp_data), ".01A|.11A")]
# ��Ӧ������Ԥ��������û�������Ԥ������/������
labels = ifelse(grepl(".01A", colnames(exp_data)), 1, 0) # ����������ǩ

#### ROC���� ----
# ����prediction����
predictions = as.numeric(exp_data)
pred = prediction(predictions, labels)
# �����������ʺͼ�������
perf = performance(pred, "tpr", "fpr")
# ����AUC���ͳ����
res = roc.area(labels, predictions)
res$A   # AUC���������
res$p.value     # ͳ����

#### ROC��ͼ ----
setwd(rocPath)
pText = ifelse(res$p.value < 0.0001, "< 0.0001", paste0("= ",round(res$p.value, 3))) # pֵ��ʾ: https://stackoverflow.com/questions/23969726/italics-and-normal-text-in-a-main-plot-title
pdf("ROC.pdf", width = 5, height = 4)
par(mar = c(3.5, 3.5, 1, 1), mgp = c(2, 0.6, 0), las = 1)
plot(perf, col = "#FF0000", lwd = 3, xlab = "1 - Specificity", ylab = "Sensitivity", cex.lab = 1.2, xaxs = "i", yaxs = "i") # xaxs = "i", yaxs = "i"��ʼ��
plot(perf, col = "#FF0000", lwd = 3, xlab = "1 - Specificity", ylab = "Sensitivity", cex.lab = 1.2)
# ���ƶԽ���
abline(0, 1, lty = 2, lwd = 2, col = "#0026FA")
# ��ע����
text(0.65, 0.2, labels = paste0("AUC = ", round(res$A,3)), cex = 1.1, pos = 4)
text(0.65, 0.1, labels = bquote(italic("P")~.(pText)), cex = 1.1, pos = 4) # 
dev.off()
