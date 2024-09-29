#path_to_your_script.R
.libPaths("C:/Users/lyr/AppData/Local/R/win-library/4.3")
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]  # 获取传递的文件路径

# 读取输入数据
data <- read.table(input_file, sep = "\t", header = TRUE)

# 在这里进行你需要的操作，举个例子，这里删除第二到第六行
processed_data <- data <- data[ , -1]

# 保存处理后的数据
output_file <- gsub(".tsv", "_processed.tsv", input_file)
write.table(processed_data, file = output_file, sep = "\t", row.names = FALSE, col.names = TRUE)

print(paste("Processed file saved to:", output_file))
