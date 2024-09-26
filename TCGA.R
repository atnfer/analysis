.libPaths("C:/Users/lyr/AppData/Local/R/win-library/4.3")
# 设置工作目录为包含数据文件的目录
setwd("C:/Users/lyr/Desktop/TCGA/gene_exp/case/")
# 搜索所有子文件夹中的TSV文件
# 找到所有子文件夹中的TSV文件
file_list <- list.files(path = "C:/Users/lyr/Desktop/TCGA/gene_exp/", pattern = "*.tsv", full.names = TRUE, recursive = TRUE)
# 创建一个空列表用于存储处理后的数据
processed_data_list <- list()

# 遍历所有文件
# 创建一个空列表用于存储处理后的数据
processed_data_list <- list()

# 遍历所有文件
for (file in file_list) {
  
  # 读取TSV文件
  data <- read.table(file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  
  # 删除第1, 3, 4, 5, 6行
  data <- data[-c(1, 3, 4, 5, 6), ]
  
  # 提取第1, 2, 3, 7列
  data <- data[, c(1, 2, 3, 7)]
  
  # 获取文件名（不带路径和扩展名）
  file_name <- basename(file)
  file_name <- sub("\\.tsv$", "", file_name)
  
  # 将第7列重命名为文件名
  colnames(data)[4] <- file_name
  
  # 将处理后的数据存入列表
  processed_data_list[[file_name]] <- data
}

# 将所有数据框按前三列合并
final_data <- Reduce(function(x, y) merge(x, y, by = c("gene_id", "gene_name", "gene_type"), all = TRUE), processed_data_list)

# 查看合并后的数据
head(final_data)
# 保存合并后的数据为CSV文件
write.csv(final_data, "C:/Users/lyr/Desktop/TCGA/TCGA_gene_exp.csv", row.names = FALSE,quote = FALSE)

#提取分组信息
gro<-read.csv("C:/Users/lyr/Desktop/gdc_manifest.2024-09-12.txt",sep = "\t")
new_df <- gro[, c("filename"), drop = FALSE]

# 添加 group 列，并将其值设为 'early_stage'
new_df$group <- 'early_stage'

# 设置主文件夹路径
main_dir <- "C:/Users/lyr/Desktop/TCGA/gene_exp/"

# 获取所有case文件夹中的TSV文件名
case_files <- list.files(path = file.path(main_dir, "case"), 
                         pattern = "\\.tsv$", 
                         full.names = TRUE, 
                         recursive = TRUE)

# 获取所有control文件夹中的TSV文件名
control_files <- list.files(path = file.path(main_dir, "control"), 
                            pattern = "\\.tsv$", 
                            full.names = TRUE, 
                            recursive = TRUE)

# 提取文件名
#case_files <- basename(case_files)
#control_files <- basename(control_files)

# 添加前缀 'X' 给文件名以数字开头的文件
#case_files <- ifelse(grepl("^\\d", case_files), paste0("X", case_files), case_files)
#control_files <- ifelse(grepl("^\\d", control_files), paste0("X", control_files), control_files)

# 创建数据框
#case_df <- data.frame(File = case_files, 
#                      Group = "case", 
#                      stringsAsFactors = FALSE)

#control_df <- data.frame(File = control_files, 
#                         Group = "control", 
#                         stringsAsFactors = FALSE)

# 合并两个数据框
all_data <- rbind(case_df, control_df)

# 显示前几行数据
head(all_data)

write.csv(all_data,"C:/Users/lyr/Desktop/TCGA/sample_info.csv",row.names = FALSE)

library(limma)
library(edgeR)

# 导入表达矩阵，假设数据已标准化为FPKM/TPM或为counts数据
expression_data <- read.csv("C:/Users/lyr/Desktop/TCGA/TCGA_gene_exp.csv", row.names = 1)
all_data<- read.csv("C:/Users/lyr/Desktop/TCGA/sample_info.csv", row.names = 1)
protein_genes <- c("ACE2", "ACP5", "ACY1", "ADA2", "ADGRE2", "ADGRG1", "BST2", "C19ORF12", 
                   "CCL15", "CD163", "CDCP1", "CDH6", "CDHR2", "CLSTN2", "CNDP1", "COL4A1D", 
                   "CTSD", "CTSL", "DDR1", "DPP10", "EFEMP1", "ENG", "ENPP2", "ERBB2D", "FGFR2B", 
                   "FSTL3", "FUT3_FUT5", "GGT1", "GRN", "HAVCR1", "IGFBP3", "IGFBP7", "IL10RB", 
                   "IL18BP", "IL18R1", "IL4RD", "IL6STD", "ITGA5", "ITGB2D", "ITGB7D", "KRT14", 
                   "LAG3 D", "LGALS9", "LTBP2", "MMED", "MSR1", "NFASC", "NOMO1", "NRCAM", 
                   "NRP2", "NT5E", "PCDH17", "PDGFRAB", "PGFD", "PLXNB2", "PVR", "SDC1", 
                   "SEMA7A", "SEZ6L2", "SLAMF1", "SPINT1", "SPP1", "TGFBI", "TIMP1", "TNFRSF11B", 
                   "TNFRSF21", "TNFSF13BD", "VCAM1", "VTCN1")
# 提取这些基因的表达数据
selected_data <- expression_data[expression_data$gene_name %in% protein_genes,]


row_names <- selected_data[, 1]
rownames(selected_data) <- row_names
selected_data<-selected_data[,-c(1,2)]
# 导入分组信息
group_info <- read.csv("C:/Users/lyr/Desktop/TCGA/sample_info.csv")  # 包含样本名称和分组信息
group <- factor(group_info$Group)  # 将分组信息转换为factor
group <- factor(group)
group <- relevel(group, ref = "control")


# 定义分组变量


# 创建表达矩阵
expr_matrix <- as.matrix(selected_data)

# 创建设计矩阵
design <- model.matrix(~group)

# 对表达数据进行voom转换
voom_data <- voom(expr_matrix, design, plot = TRUE)

# 拟合线性模型
fit <- lmFit(voom_data, design)

# 进行差异表达分析
fit <- eBayes(fit)

# 查看差异表达的结果
topTable(fit, coef = 2)

# 提取所有结果
all_results <- topTable(fit, number = Inf, adjust.method = "BH")

# 查看所有显著的结果（例如，FDR < 0.05）
significant_results <- all_results[all_results$adj.P.Val < 0.05, ]



library(ggplot2)



# 创建火山图
figure1<-ggplot(significant_results, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = adj.P.Val < 0.05 & abs(logFC) > 1)) +  # 高亮显著基因
  scale_color_manual(values = c("grey", "red")) +  # 设定颜色
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 P-Value")
# 创建 MA 图
figure2<-plotMA(fit, coef = 2, main = "MA Plot", ylim = c(-5, 5))


#pca分析
library(ggfortify)

# 进行PCA分析
pca_results <- prcomp(t(expr_matrix))

# 绘制PCA图
autoplot(pca_results, data = as.data.frame(group), colour = 'group') +
  theme_minimal() +
  labs(title = "PCA of Samples")

#箱图
# 加载必要的包
library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggsignif)
# 假设 expr_matrix 是你的表达矩阵，行名是蛋白质名称，列名是样本名称
expr_matrix <-subset(expr_matrix,expr_matrix %in% rownames(significant_results))
# group 是一个包含样本分组信息的向量

# 将矩阵转换为长格式数据框
expr_long <- melt(expr_matrix)
colnames(expr_long) <- c("Protein", "Sample", "Expression")
expr_long <-subset(expr_long,expr_long$Protein %in% rownames(significant_results))
# 添加样本分组信息
expr_long$Group <- group[match(expr_long$Sample, colnames(expr_matrix))]
# 使用 ggplot2 绘制每个蛋白质的箱线图


# 创建一个ggplot对象
expr_long$LogExpression <- log(expr_long$Expression + 1)  # +1 避免 log(0) 的情况

# 使用对数转换后的数据绘制箱线图
p <- ggplot(expr_long, aes(x = Group, y = LogExpression, fill = Group)) +
  geom_boxplot() +
  facet_wrap(~ Protein, scales = "free_y") +
  geom_signif(comparisons = list(c("case", "control")),
              map_signif_level = TRUE,
              step_increase = 0.1) +
  theme_minimal() +
  labs(title = "Boxplots of Protein Expression with Significance", x = "Group", y = "Expression") +
  theme(
    # 增大标题字体
    plot.title = element_text(size = 25),
    # 增大 x 轴和 y 轴标签字体
    axis.title.x = element_text(size = 21),
    axis.title.y = element_text(size = 21),
    # 增大 x 轴和 y 轴刻度标签字体
    axis.text.x = element_text(size = 19),
    axis.text.y = element_text(size = 19),
    # 增大 facet 标签字体
    strip.text = element_text(size = 19)
  )

# 添加显著性标记
figure4<-p + geom_signif(comparisons = list(c("case", "control")),
                map_signif_level = TRUE,
                step_increase = 0.1)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "pathview", "KEGGREST","GENIE3","WGCNA","ComplexHeatmap","rjson"))
library(clusterProfiler)
library(org.Hs.eg.db)  # 用于人类基因注释
library(pathview)
library(KEGGREST)

genes1 <- rownames(significant_results)
deg_genes <- bitr(genes1, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
genes2 <- subset(expression_data, gene_name %in% genes)
deg_genes<- rownames(genes2)
deg_genes <- sub("\\..*", "", deg_genes)
# 进行 GO 富集分析
go_enrichment <- enrichGO(gene = deg_genes,
                          OrgDb = org.Hs.eg.db,      # 根据物种选择相应的数据库
                          keyType = "ENSEMBL",       # 根据基因ID类型调整
                          ont = "BP",                # BP: Biological Process, CC: Cellular Component, MF: Molecular Function
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2)

# 查看和可视化结果
head(go_enrichment)
dotplot(go_enrichment, showCategory=10) + theme_minimal()

gene_entrez <- bitr(genes1, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
# 进行 KEGG 富集分析
kegg_enrichment <- enrichKEGG(gene = gene_entrez$ENTREZID,
                              organism = 'hsa',          # 根据物种选择，如人类为 'hsa'
                              pAdjustMethod = "BH",
                              pvalueCutoff = 0.5)

# 查看和可视化结果
head(kegg_enrichment)
figure5<-dotplot(kegg_enrichment, showCategory=10) + theme_minimal()

# 使用 pathview 可视化特定通路
library(pathview)
data<-pathview(gene.data = genes1,
         pathway.id = 'hsa04110',      # 选择感兴趣的 KEGG 通路 ID
         species = 'hsa')
# 安装并加载 STRINGdb
BiocManager::install("STRINGdb")
library(STRINGdb)
library(igraph)

# 初始化 STRINGdb
string_db <- STRINGdb$new(version="11.0", species=9606, score_threshold=400, 
                          input_directory="")

# 将基因列表映射到 STRING ID
mapped_genes <- string_db$map(data.frame(gene = genes), "gene", removeUnmappedRows = TRUE)

# 获取 PPI 网络
ppi_network <- string_db$get_interactions(mapped_genes$STRING_id)

# 构建 igraph 对象
g <- graph_from_data_frame(ppi_network, directed = FALSE)

# 简单可视化

library(igraph)
library(ggraph)
library(tidygraph)

# 使用 ggraph 可视化
figure6<-ggraph(g, layout = 'fr') +
  geom_edge_link(aes(edge_alpha = 0.5, edge_width = 1), show.legend = FALSE) +
  geom_node_point(aes(size = degree(g)), color = "lightblue", show.legend = FALSE) +
  geom_node_text(aes(label = name), repel = TRUE, size = 4) +
  theme_void()

library(GENIE3)

# 构建调控网络
weight_matrix <- GENIE3(expr_matrix)

# 获取高权重的调控关系
top_edges <- getLinkList(weight_matrix, threshold = 0.01)

# 构建 igraph 对象
library(igraph)
g <- graph_from_data_frame(top_edges, directed = TRUE)

# 简单可视化
plot(g, vertex.size=5, vertex.label=NA, edge.arrow.size=0.5)

library(WGCNA)

# 准备数据
datExpr <- as.data.frame(t(expr_matrix))  # 样本为行，基因为列

# 检查和处理缺失值
gsg = goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
}

# 构建共表达网络
powers = c(1:10)
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# 选择合适的软阈值
softPower = 6

# 构建邻接矩阵和拓扑重叠矩阵 (TOM)
adjacency = adjacency(datExpr, power = softPower)
TOM = TOMsimilarity(adjacency)
dissTOM = 1 - TOM

# 聚类
geneTree = hclust(as.dist(dissTOM), method = "average")
plot(geneTree, main = "Gene Clustering on TOM-based Dissimilarity", 
     xlab = "", sub = "")

library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
# 示例的列注释
coul <- colorRampPalette(brewer.pal(9, "OrRd"))(50)
annotation_col <- data.frame(Group = group)
rownames(annotation_col) <- colnames(expr_significant)
expr_matrix_z <- t(scale(t(expr_matrix)))
heatmap <- Heatmap(
  expr_matrix_z,
  name = "Expression",
  col = coul,
  cluster_rows = TRUE,
  cluster_columns = FALSE,  # 关闭列聚类
  show_column_names = FALSE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 10),
  top_annotation = HeatmapAnnotation(
    df = annotation_col,
    col = list(Group = c(case = "lightblue", control = "lightgreen"))
  )
)

# 绘制热图
draw(heatmap, heatmap_legend_side = "right")
# 假设 group 是分组因子
group <- factor(group)
group <- relevel(group, ref = "control")

# 创建设计矩阵
design <- model.matrix(~ group)
voom_data <- voom(expr_matrix, design, plot = TRUE)

# 拟合线性模型
fit <- lmFit(voom_data, design)

# 进行差异表达分析
fit <- eBayes(fit)

# 查看差异表达的结果
results <- topTable(fit, coef = 2)
significant_results <- results[results$adj.P.Val < 0.05, ]
