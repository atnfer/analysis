#annotation
import os
import pandas as pd
import subprocess

# 输入文件夹路径和输出文件夹路径
input_dir = 'C:/Users/lyr/pyscript'
output_dir = 'C:/Users/lyr/pyscript'
r_script_path = 'C:/Users/lyr/pyscript/script.R'  # R脚本路径

# 如果输出目录不存在，则创建
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# 获取文件夹中的所有TSV文件
files = [f for f in os.listdir(input_dir) if f.endswith('.tsv')]

# 遍历每个文件并处理
for file in files:
    input_file = os.path.join(input_dir, file)
    
    # Step 1: 读取文件
    df = pd.read_csv(input_file, sep='\t')
    
    # Step 2: 根据染色体编号拆分文件
    for chrom in df['hm_chrom'].unique():
        chrom_df = df[df['hm_chrom'] == chrom]
        
        # Step 3: 保存每个染色体的临时文件
        temp_file = os.path.join(output_dir, f'{file}_chrom_{chrom}.tsv')
        chrom_df.to_csv(temp_file, sep='\t', index=False)
        
        # Step 4: 调用Rscript处理文件
        result = subprocess.run(['Rscript', r_script_path, temp_file], capture_output=True, text=True)
        if result.returncode != 0:
            print(f'Error processing {temp_file}: {result.stderr}')
        else:
            print(f'Successfully processed {temp_file}')

# Step 5: 合并处理后的文件
processed_files = [os.path.join(output_dir, f) for f in os.listdir(output_dir) if f.endswith('_processed.tsv')]
dfs = [pd.read_csv(f, sep='\t') for f in processed_files]
merged_df = pd.concat(dfs)

# 保存合并后的文件
first_input_file = processed_files[0]
file_name, file_extension = os.path.splitext(os.path.basename(first_input_file))

# 创建新的文件名，添加 "_anno" 后缀，保留原始扩展名
merged_output_file = os.path.join(output_dir, f"{file_name}_anno{file_extension}")

merged_df.to_csv(merged_output_file, sep='\t', index=False)
print(f'Saved merged output file: {merged_output_file}')