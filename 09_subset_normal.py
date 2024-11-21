import os
import pandas as pd
import scanpy as sc
import scipy
from scipy.sparse import csr_matrix
import anndata as ad

# 1. Load the CSV file
file_list = pd.read_csv('./08_split_list.csv', sep = '\t')

# 2. Create 'subset' directory if it doesn't exist
os.makedirs('09_subset', exist_ok=True)

# 3. Loop through each row in the file list and process .h5ad files
file_info_list = []

for _, row in file_list.iterrows():
    absolute_path = row['Absolute_Path']
    download_id = row['Download_ID']
    print(f"Processing file Download ID : {download_id} -----------------------------------------------------")
    
    try:
        adata = sc.read_h5ad(absolute_path)
    
        # ------- 1.筛选 human normal 数据  -------------
        adata_filtered = adata[(adata.obs['organism'] == 'Homo sapiens') & 
                               (adata.obs['disease'].isin(['normal', 'Healthy'])), :]
        adata_info = str(adata_filtered).split('\n')
    
        # ------- 2. 创建 dict 存储信息  -------------
        data_dict = {}
        for line in adata_info:
            key_value = line.split(': ')
            if len(key_value) == 2:
                key, value = key_value[0].strip(), key_value[1].strip()
                data_dict[key] = [value]
    
        data_dict['Download_ID'] = download_id
        data_dict['sub_cell_number'] = adata_filtered.n_obs
        data_dict['sub_gene_number'] = adata_filtered.n_vars
        data_dict['sub_ENS_count'] = sum(['ENS' in gene for gene in adata_filtered.var_names])
        data_dict['sub_organism'] = ";".join(adata_filtered.obs['organism'].astype(str).unique()) if 'organism' in adata_filtered.obs else "Unknown"
    
        # ------- 3. 统计 X 信息  -------------
        if hasattr(adata_filtered, 'X'):
            data_dict.update({
                'sub_dim_X': ';'.join(map(str, adata_filtered.X.shape)),
                'sub_X_max': adata_filtered.X.max(),
                'sub_X_min': adata_filtered.X.min(),
                'sub_type_X': 'csr_matrix' if isinstance(adata_filtered.X, csr_matrix) else 'dense',
                'sub_feature_name_X': 'True' if 'feature_name' in adata_filtered.var else 'False'
            })
    
        ## Stats for adata_filtered.raw.X if it exists
        if adata_filtered.raw is not None and hasattr(adata_filtered.raw, 'X'):
            data_dict.update({
                'sub_dim_raw_X': ';'.join(map(str, adata_filtered.raw.X.shape)),
                'sub_raw_X_max': adata_filtered.raw.X.max() if adata_filtered.raw.X.size > 0 else '',
                'sub_raw_X_min': adata_filtered.raw.X.min() if adata_filtered.raw.X.size > 0 else '',
                'sub_type_raw': 'csr_matrix' if isinstance(adata_filtered.raw.X, csr_matrix) else 'dense',
                'sub_feature_name_raw': 'True' if 'feature_name' in adata_filtered.raw.var else 'False',
                'sub_X_size_same_with_raw_X': 'True' if adata_filtered.X.shape == adata_filtered.raw.X.shape else 'False'
            })
        else:
            data_dict.update({
                'sub_dim_raw_X': '',
                'sub_raw_X_max': '',
                'sub_raw_X_min': '',
                'sub_type_raw': '',
                'sub_feature_name_raw': '',
                'sub_X_size_same_with_raw_X': ''
            })
            
        # ------- 4. 保存 .h5ad 数据  -------------
        output_path = os.path.join('09_subset', f"{download_id}_human_normal.h5ad")
        data_dict['sub_path'] = os.path.abspath(output_path)
        
        adata_filtered.write(output_path)
        
        # ------- 5. 保存 meta 数据  -------------
        ## 将信息添加到 DataFrame 列表中
        # 删除 DataFrame 中的 'varm' 列（如果存在）
        if 'varm' in data_dict:
            del data_dict['varm']
        
        # 将 data_dict 转换为 DataFrame 后检查是否仍包含 'varm'
        df = pd.DataFrame(data_dict)
        if 'varm' in df.columns:
            df = df.drop(columns=['varm'])
        
        file_info_list.append(df)
        
    except Exception as e:
        print(f"Error processing {download_id}: {e}")

# 4. 合并所有文件的信息为一个 DataFrame
combined_df = pd.concat(file_info_list, ignore_index=True)
## 找到重复的列
common_columns = set(combined_df.columns) & set(file_list.columns)

## 移除 Download_ID，因为它是用于匹配的列，不需要重命名
common_columns.discard('Download_ID')

## 给 combined_df 中的重复列加上 'sub_' 前缀
combined_df = combined_df.rename(columns={col: f'sub_{col}' for col in common_columns})

## 使用 merge 函数根据 'Download_ID' 进行内连接合并
merged_df = pd.merge(file_list, combined_df, on='Download_ID', how='inner')

## 保存
merged_df.to_csv('./09_subset_normal_meta.csv', sep="\t", index=False)

# Load both CSV files
xx1 = pd.read_csv("./08_no_split.csv", sep='\t')
xx2 = pd.read_csv("./09_subset_normal_meta.csv", sep='\t')
xx2_filtered = xx2.iloc[:, [4, 20] + list(range(22, 43))]
# Remove 'sub_' prefix from column names in xx2_filtered
xx2_filtered.columns = [col.replace('sub_', '') if col.startswith('sub_') else col for col in xx2_filtered.columns]
# Rename 'path' to 'Absolute_Path' in xx2_filtered
xx2_filtered = xx2_filtered.rename(columns={'path': 'Absolute_Path'})
# Merge xx1 with xx2_filtered on the 'Download_ID' column
merged_df = pd.concat([xx1, xx2_filtered], ignore_index=True)

# Save the combined data to a new CSV file
output_file = "./09_hu_normal_meta.csv"
merged_df.to_csv(output_file, sep="\t", index=False)

print("Data processing completed. Files saved to '09_subset_normal_meta.csv', '09_subset/Download_ID_human_normal.h5ad','09_hu_normal_meta.csv'.")
