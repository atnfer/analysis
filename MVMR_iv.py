import pandas as pd
import os

# 读取第一个和第二个文件
file1 = "path_to_first_file.csv"
file2 = "path_to_second_file.csv"

df1 = pd.read_csv(file1)
df2 = pd.read_csv(file2)

# 提取 SNP 列并合并，保留不重复的 SNP
snp_df = pd.DataFrame(pd.concat([df1['SNP'], df2['SNP']]).unique(), columns=['SNP'])

# 读取第三个文件，基于 SNP 列提取匹配行
file3 = "path_to_third_file.csv"
df3 = pd.read_csv(file3)
df3_filtered = df3[df3['SNP'].isin(snp_df['SNP'])]

# 修改 df3_filtered 列名，添加文件名后缀（不包括 SNP 列）
file3_name = os.path.splitext(os.path.basename(file3))[0]
df3_filtered.columns = ['SNP'] + [f"{col}_{file3_name}" for col in df3_filtered.columns if col != 'SNP']

# 对第四和第五个文件进行相同操作
file4 = "path_to_fourth_file.csv"
file5 = "path_to_fifth_file.csv"

df4 = pd.read_csv(file4)
df4_filtered = df4[df4['SNP'].isin(snp_df['SNP'])]

file4_name = os.path.splitext(os.path.basename(file4))[0]
df4_filtered.columns = ['SNP'] + [f"{col}_{file4_name}" for col in df4_filtered.columns if col != 'SNP']

df5 = pd.read_csv(file5)
df5_filtered = df5[df5['SNP'].isin(snp_df['SNP'])]

file5_name = os.path.splitext(os.path.basename(file5))[0]
df5_filtered.columns = ['SNP'] + [f"{col}_{file5_name}" for col in df5_filtered.columns if col != 'SNP']

# 合并所有数据
final_result = snp_df \
    .merge(df3_filtered, on='SNP', how='left') \
    .merge(df4_filtered, on='SNP', how='left') \
    .merge(df5_filtered, on='SNP', how='left')

# 保存结果
final_result.to_csv("final_result.csv", index=False)

