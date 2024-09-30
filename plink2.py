import subprocess
import pandas as pd
import os


# 定义要 clump 的文件列表，每个文件对应一个 GWAS 结果文件
input_dir = 'C:/Users/lyr/pyscript'
gwas_files  = [f for f in os.listdir(input_dir) if f.endswith('.tsv')]  # 替换为你的文件名
clumped_files = []

# 依次对每个文件运行 plink2 clump
for gwas_file in gwas_files:
    output_prefix = os.path.splitext(gwas_file)[0] + '_clumped'
    clumped_output = f'{output_prefix}.clumped'
    
    # plink2 clump 命令
    plink_cmd = [
        'plink2',
        '--bfile', input_file,  # 输入的文件前缀
        '--clump', gwas_file,  # GWAS 结果文件
        '--clump-kb', '250',  # 设置 clump 的窗口大小
        '--clump-p1', '5e-8',  # clump 的显著性阈值
        '--clump-r2', '0.1',  # clump 的 LD 阈值
        '--out', output_prefix  # 输出文件前缀
    ]
    
    try:
        subprocess.run(plink_cmd, check=True)
        print(f'PLINK clump 运行成功，输出结果保存到 {clumped_output}')
        clumped_files.append(clumped_output)
    except subprocess.CalledProcessError as e:
        print(f"PLINK clump 运行失败: {e}")
        continue

# 合并所有 clumped 结果
merged_clumped_file = 'merged_clumped_results.csv'
merged_df = pd.DataFrame()

for clumped_file in clumped_files:
    try:
        df = pd.read_csv(clumped_file, delim_whitespace=True)
        df['Source_File'] = clumped_file  # 添加一列标识来源文件
        merged_df = pd.concat([merged_df, df])
    except FileNotFoundError:
        print(f'无法找到文件: {clumped_file}')
        continue

# 保存合并后的结果
merged_df.to_csv(merged_clumped_file, index=False)
print(f'合并后的 clumped 结果保存到 {merged_clumped_file}')
