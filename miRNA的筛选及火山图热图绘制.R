# ============================
# TCGA-CESC 宫颈癌 miRNA 差异分析
# ============================
setwd("D:\\shengxin\\TCGA")   # 请修改为您的实际工作目录

# ---------------------------
# 0. 安装/加载所需包
# ---------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
pkgs <- c("TCGAbiolinks", "limma", "edgeR", "ggplot2", "ggrepel", "pheatmap")
for (p in pkgs) {
  if (!require(p, character.only = TRUE)) {
    if (p %in% c("TCGAbiolinks", "limma", "edgeR"))
      BiocManager::install(p)
    else
      install.packages(p)
  }
  library(p, character.only = TRUE)
}

# ---------------------------
# 1. 设置网络超时（避免下载失败）
# ---------------------------
options(timeout = 300)

# ---------------------------
# 2. 下载 miRNA 表达数据（带错误处理）
# ---------------------------
project <- "TCGA-CESC"
getProjectSummary(project)

download_success <- tryCatch({
  query <- GDCquery(
    project = project,
    data.category = "Transcriptome Profiling",
    data.type = "miRNA Expression Quantification",
    experimental.strategy = "miRNA-Seq",
    workflow.type = "BCGSC miRNA Profiling"
  )
  GDCdownload(query)
  expData <- GDCprepare(query)
  TRUE
}, error = function(e) {
  message("❌ GDC下载失败: ", e$message)
  message("请尝试手动下载数据，然后运行下方备用代码读取本地文件。")
  return(FALSE)
})

if (!download_success) {
  # 备用方案：从本地文件读取（请根据实际文件路径修改）
  local_file <- "miRNA_expression_matrix.txt"
  if (file.exists(local_file)) {
    message("从本地文件读取数据...")
    expData <- read.table(local_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  } else {
    stop("无法下载数据且未找到本地文件，请检查网络或手动下载数据后重试。")
  }
}

# ---------------------------
# 3. 提取 read_count 矩阵
# ---------------------------
if (is(expData, "RangedSummarizedExperiment")) {
  expr <- assay(expData)
  if (is.null(rownames(expr))) {
    rownames(expr) <- rowData(expData)$miRNA_ID
  }
} else {
  count_cols <- grep("^read_count", colnames(expData), value = TRUE)
  expr <- expData[, count_cols]
  rownames(expr) <- expData$miRNA_ID
  colnames(expr) <- gsub("^read_count_", "", colnames(expr))
}
colnames(expr) <- substr(colnames(expr), 1, 15)
cat("样本ID示例：", head(colnames(expr), 5), "\n")
cat("原始表达矩阵维度：", dim(expr), "\n")

# ---------------------------
# 4. 样本分组（肿瘤 vs 正常）
# ---------------------------
sample_codes <- substr(colnames(expr), 14, 15)
group <- factor(
  ifelse(as.integer(sample_codes) < 10, "tumor", "normal"),
  levels = c("normal", "tumor")
)
cat("分组统计：\n")
print(table(group))
if (sum(group == "normal") < 2) {
  stop("❌ 正常样本数量不足2个，无法进行差异分析。")
}

# ---------------------------
# 5. 过滤低表达 miRNA（至少半数样本表达量 > 0）
# ---------------------------
before <- nrow(expr)
keep <- rowSums(expr > 0) >= ncol(expr) / 2
expr <- expr[keep, ]
after <- nrow(expr)
cat(sprintf("过滤前 miRNA 数量：%d\n", before))
cat(sprintf("过滤后剩余 miRNA 数量：%d (移除 %d 个)\n", after, before - after))

expr <- as.matrix(expr)
mode(expr) <- "numeric"

# ---------------------------
# 6. 保存表达矩阵和分组信息
# ---------------------------
out_matrix <- cbind(miRNA_ID = rownames(expr), expr)
write.table(out_matrix, file = "CESC-miRNA-matrix.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)
group_df <- data.frame(sample = colnames(expr), group = group)
write.table(group_df, file = "CESC-miRNA-group.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---------------------------
# 7. edgeR 差异分析
# ---------------------------
logFC_cutoff <- 1
padj_cutoff <- 0.05

dge <- DGEList(counts = expr, group = group)
dge <- calcNormFactors(dge)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
dge <- estimateDisp(dge, design)
fit <- glmFit(dge, design)
lrt <- glmLRT(fit, contrast = c(-1, 1))

# 注意：topTags 默认按 PValue 排序，我们稍后会自行按 FDR 排序
DEG <- as.data.frame(topTags(lrt, n = Inf))
DEG$change <- factor(
  ifelse(DEG$FDR < padj_cutoff & abs(DEG$logFC) > logFC_cutoff,
         ifelse(DEG$logFC > 0, "UP", "DOWN"), "NOT")
)

# 输出上下调个数
cat("\n========== 差异统计 ==========\n")
cat("上调 miRNA 数量：", sum(DEG$change == "UP"), "\n")
cat("下调 miRNA 数量：", sum(DEG$change == "DOWN"), "\n")
cat("不显著 miRNA 数量：", sum(DEG$change == "NOT"), "\n")

# 保存全部结果（按原始 PValue 排序，但包含所有信息）
write.csv(DEG, file = "CESC-miRNA-edgeR.csv", quote = FALSE)
sig <- DEG[DEG$change != "NOT", ]
write.table(sig, file = "CESC-significant.xls", sep = "\t", quote = FALSE)

# ---------------------------
# 8. 上调 miRNA 按多种方式排序并输出（重点：按 FDR 升序）
# ---------------------------
up_regulated <- DEG[DEG$change == "UP", ]
if (nrow(up_regulated) > 0) {
  up_by_logFC <- up_regulated[order(up_regulated$logFC, decreasing = TRUE), ]
  up_by_PValue <- up_regulated[order(up_regulated$PValue, decreasing = FALSE), ]
  up_by_FDR <- up_regulated[order(up_regulated$FDR, decreasing = FALSE), ]   # 主要推荐
  up_by_LR <- up_regulated[order(up_regulated$LR, decreasing = TRUE), ]
  
  write.csv(up_by_logFC, "CESC-upregulated-sorted-by-logFC.csv")
  write.csv(up_by_PValue, "CESC-upregulated-sorted-by-PValue.csv")
  write.csv(up_by_FDR, "CESC-upregulated-sorted-by-FDR.csv")
  write.csv(up_by_LR, "CESC-upregulated-sorted-by-LR.csv")
  
  cat("\n上调 miRNA 排序文件已保存：\n",
      "  - CESC-upregulated-sorted-by-logFC.csv (按logFC降序)\n",
      "  - CESC-upregulated-sorted-by-PValue.csv (按PValue升序)\n",
      "  - CESC-upregulated-sorted-by-FDR.csv (按FDR升序，推荐)\n",
      "  - CESC-upregulated-sorted-by-LR.csv (按LR降序)\n")
}

# ---------------------------
# 9. 火山图（按 FDR 升序取前20个最显著miRNA，标注名称）
# ---------------------------
DEG_plot <- DEG
DEG_plot$FDR[DEG_plot$FDR == 0] <- .Machine$double.eps

sig_deg <- DEG[DEG$change != "NOT", ]
# 按 FDR 升序排列，取前20（FDR越小越显著）
top20 <- head(sig_deg[order(sig_deg$FDR, decreasing = FALSE), ], 20)

p <- ggplot(DEG_plot, aes(x = logFC, y = -log10(FDR), color = change)) +
  geom_point(size = 1.2, alpha = 0.6) +
  scale_color_manual(values = c("UP" = "red", "DOWN" = "blue", "NOT" = "gray")) +
  labs(x = expression(log[2]~Fold~Change), 
       y = expression(-log[10]~FDR), 
       title = "Top 20 Significant DE miRNAs by FDR (Volcano Plot)") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.title = element_blank()) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray30") +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "gray30")

if (nrow(top20) > 0) {
  p <- p + geom_text_repel(data = top20,
                           aes(label = rownames(top20)),
                           size = 4, fontface = "bold",
                           box.padding = 0.6, point.padding = 0.3,
                           segment.color = "grey30", segment.size = 0.3,
                           max.overlaps = 30, force = 2)
}

ggsave("CESC-miRNA-volcano.png", plot = p, width = 12, height = 8, dpi = 300)
ggsave("CESC-miRNA-volcano.pdf", plot = p, width = 12, height = 8)
print(p)

# ---------------------------
# 10. 热图（按 FDR 升序取前20个最显著miRNA，PNG格式）
# ---------------------------
if (nrow(sig_deg) == 0) {
  cat("没有显著差异 miRNA，跳过热图绘制。\n")
} else {
  top20_heat <- head(sig_deg[order(sig_deg$FDR, decreasing = FALSE), ], 20)
  top20_heat <- top20_heat[rownames(top20_heat) %in% rownames(expr), ]
  
  if (nrow(top20_heat) == 0) {
    cat("没有可用的miRNA用于热图。\n")
  } else {
    heat_data <- expr[rownames(top20_heat), ]
    logcpm <- cpm(heat_data, log = TRUE, prior.count = 1)
    
    anno_col <- data.frame(Group = group)
    rownames(anno_col) <- colnames(heat_data)
    
    png("CESC-miRNA-heatmap.png", width = 10, height = 8, units = "in", res = 300)
    pheatmap(logcpm,
             annotation_col = anno_col,
             show_colnames = FALSE,
             scale = "row",
             color = colorRampPalette(c("blue", "white", "red"))(100),
             border_color = NA,
             main = "Top 20 Significant DE miRNAs by FDR",
             fontsize_row = 10)
    dev.off()
    cat("热图已保存为 CESC-miRNA-heatmap.png\n")
  }
}

cat("\n========== 全部分析完成 ==========\n")