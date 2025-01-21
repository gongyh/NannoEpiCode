# 读取 featureCounts 输出
counts <- read.table("H3k9ac.consensus_peaks.featureCounts.txt", header = TRUE, row.names = 1, sep = "\t", skip = 1)

# 提取 counts 矩阵
count_matrix <- counts[, 6:ncol(counts)]  # 第 6 列开始是样本的 counts
colnames(count_matrix) <- c("H24_R2", "H0_R1", "H24_R1", "H0_R2")  # 修改列名

library(edgeR)

# 创建 DGEList 对象
dge <- DGEList(counts = count_matrix)

# 计算标准化因子
dge <- calcNormFactors(dge, method = "TMM")

# 获取标准化后的信号
normalized_counts <- cpm(dge)

# 查看归一化后的矩阵
head(normalized_counts)

# 保存为数据矩阵
write.table(normalized_counts, file = "H3K9ac_normalized_matrix.txt", sep = "\t", quote = FALSE)
