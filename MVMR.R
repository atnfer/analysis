.libPaths("C:/Users/lyr/AppData/Local/R/win-library/4.3")
# 加载必要的包
library(MVMR)
library(dplyr)
# 加载 MVMR 包
library(MVMR)

# 获取从 Python 传递的文件路径
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]

# 读取传递过来的 CSV 文件
data <- read.csv(input_file)

# 准备 IVW_MVMR 分析所需的数据
r_input <- format_mvmr(
    BXGs = merged_data[,c("phemo_beta","circ_beta")],
    BYG = merged_data$beta,
    seBXGs = merged_data[,c("phemo_se","circ_se")],
    seBYG = merged_data$sebeta,
    RSID = merged_data$SNP)
# 执行 MVMR 分析
mvmr_result<-ivw_mvmr(r_input, gencov = 0)
  
# 获取文件名（不包括扩展名），用于命名输出文件
  file_name <- tools::file_path_sans_ext(basename(input_file))
# 保存 MVMR 分析结果，使用文件名
write.csv(mvmr_result, paste0("mvmr_result_", file_name, ".csv"), row.names = FALSE)