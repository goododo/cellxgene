import pandas as pd
import scanpy as sc
import os
import scipy
from scipy.sparse import csr_matrix
import anndata as ad
import numpy as np

# 1. 读取并处理数据
data_path = "./05_merged_data.csv"
meta = pd.read_csv(data_path, sep='\t')

## 修改 Assay 列
meta.loc[meta['Download_ID'] == 773, 'Assay'] = "10x scATAC-seq"

## 筛选 Organism 为 Homo sapiens 的数据
meta_human = meta[meta['Organism'].str.contains("Homo sapiens", na=False)]  # 1214, 24

## 筛选 Disease 为 normal 的数据
meta_human_nor = meta_human[meta_human['Disease'].str.contains("normal", na=False)]  # 996, 24

## 排除指定的 Assay 类型
exclude_assays = [
    "methylation profiling by high throughput sequencing",
    "snmC-Seq2",
    "mCT-seq",
    "BD Rhapsody Targeted mRNA",
    "10x scATAC-seq",
    "Abseq",
    "snm3C-seq"
]
hu_n_filtered = meta_human_nor[~meta_human_nor['Assay'].isin(exclude_assays)]  # 975, 24

## 统计细胞数
meta_cell_count = meta['Cell_count'].sum()
meta_human_cell_count = meta_human['Cell_count'].sum()
meta_human_nor_cell_count = meta_human_nor['Cell_count'].sum()
hu_n_filtered_cell_count = hu_n_filtered['Cell_count'].sum()

print(f"Total Cell Count: {meta_cell_count}")  # 176310638
print(f"Human Cell Count: {meta_human_cell_count}")  # 102043248
print(f"Human Normal Cell Count: {meta_human_nor_cell_count}")  # 94128533
print(f"Human Normal RNA Cell Count: {hu_n_filtered_cell_count}")  # 93097391

# 2. 根据 Download_ID 和 Dataset_id 构造文件路径，并存储到新列 location 中
base_dir = "/data/zgzheng/database_download/cellxgene/06_downloads"
meta['location'] = meta.apply(
    lambda row: os.path.join(base_dir, str(row['Download_ID']), f"{row['Download_ID']}.h5ad"), axis=1
)

# 3. 读取并提取信息，直接将新信息添加到 meta 数据框中

# 定义新列及其数据类型
new_columns_dtypes = {
    'cell_number': 'Int64',  # 可为空的整数类型
    'gene_number': 'Int64',
    'ENS_count': 'Int64',
    'organism_data': 'object',
    'dim_X': 'object',
    'X_max': 'float',
    'X_min': 'float',
    'type_X': 'object',
    'feature_name_X': 'object',
    'dim_raw_X': 'object',
    'raw_X_max': 'float',
    'raw_X_min': 'float',
    'type_raw': 'object',
    'feature_name_raw': 'object',
    'X_size_same_with_raw_X': 'object',
    'disease_data': 'object',
    'Absolute_Path': 'object',
    'Cell_number_match': 'boolean',
    'Gene_number_match': 'boolean',
    'Assay_match': 'boolean',
    'Tissue_match': 'boolean',
    'Cell_type_match': 'boolean',
    'Sex_match': 'boolean',
    'Disease_match': 'boolean',
    'Development_stage_match': 'boolean',
    'Organism_match': 'boolean',
    'Suspension_type_match': 'boolean',
    'processing_error': 'object'
}

# 在 meta 数据框中添加新列，指定数据类型
for col, dtype in new_columns_dtypes.items():
    if dtype == 'float':
        # 对于浮点类型的列，使用 np.nan 作为缺失值
        meta[col] = pd.Series([np.nan]*len(meta), dtype=dtype)
    else:
        # 其他类型的列，使用 pd.NA 作为缺失值
        meta[col] = pd.Series([pd.NA]*len(meta), dtype=dtype)

# 定义匹配函数
def match_column(obs_values, meta_value):
    # 将 meta_value 按逗号或中文逗号分割，并去除空白
    meta_list = [item.strip() for item in str(meta_value).split('，')]
    # 比较 obs_values 的唯一值与 meta_list
    return sorted(obs_values.unique()) == sorted(meta_list)

# 遍历 meta 数据框的每一行
for index, row in meta.iterrows():
    file = row['location']
    download_id = row['Download_ID']
    print(f"Processing file: {file} -----------------------------------------------------")
    
    try:
        adata = sc.read(file)
        adata_info = str(adata).split('\n')
        
        # ------- 1. 创建 dict 存储信息  -------------
        data_dict = {}
        for line in adata_info:
            key_value = line.split(': ')
            if len(key_value) == 2:
                key, value = key_value[0].strip(), key_value[1].strip()
                data_dict[key] = value

        # 提取信息
        data_dict['Download_ID'] = download_id
        data_dict['cell_number'] = adata.n_obs
        data_dict['gene_number'] = adata.n_vars
        data_dict['ENS_count'] = sum(['ENS' in gene for gene in adata.var_names])
        
        ## Organism 信息
        if 'organism' in adata.obs:
            organism = ";".join(adata.obs['organism'].astype(str).unique())
        else:
            organism = "Unknown"
        data_dict['organism_data'] = organism  # 修改键为 'organism_data'
        
        # ------- 2. 统计 X 信息  -------------
        ## adata.X 信息
        if hasattr(adata, 'X'):
            data_dict['dim_X'] = ';'.join(map(str, adata.X.shape))
            data_dict['X_max'] = adata.X.max()
            data_dict['X_min'] = adata.X.min()
            data_dict['type_X'] = 'csr_matrix' if isinstance(adata.X, csr_matrix) else 'dense'
            data_dict['feature_name_X'] = 'True' if 'feature_name' in adata.var else 'False'

        ## adata.raw.X 信息
        data_dict['dim_raw_X'] = pd.NA
        data_dict['raw_X_max'] = pd.NA
        data_dict['raw_X_min'] = pd.NA
        data_dict['type_raw'] = pd.NA
        data_dict['feature_name_raw'] = pd.NA
        data_dict['X_size_same_with_raw_X'] = pd.NA
        
        # 检查 adata.raw 和 adata.raw.X 是否存在且包含数据
        if adata.raw is not None and hasattr(adata.raw, 'X'):
            data_dict['dim_raw_X'] = ';'.join(map(str, adata.raw.X.shape)) if adata.raw.X.shape else pd.NA
            data_dict['raw_X_max'] = adata.raw.X.max() if adata.raw.X.size > 0 else pd.NA
            data_dict['raw_X_min'] = adata.raw.X.min() if adata.raw.X.size > 0 else pd.NA
            data_dict['type_raw'] = 'csr_matrix' if isinstance(adata.raw.X, csr_matrix) else 'dense'
            data_dict['feature_name_raw'] = 'True' if 'feature_name' in adata.raw.var else 'False'
            
            # 检查 X 和 raw.X 是否具有相同的形状
            data_dict['X_size_same_with_raw_X'] = 'True' if adata.X.shape == adata.raw.X.shape else 'False'
        
        ## Disease 信息
        if 'disease' in adata.obs:
            disease = ";".join(adata.obs['disease'].astype(str).unique())
        else:
            disease = "Unknown"
        data_dict['disease_data'] = disease  # 修改键为 'disease_data'
        
        ## 绝对路径
        data_dict['Absolute_Path'] = os.path.abspath(file)
        
        # ------- 3. 检查实际数据与 meta 表格信息是否一致 -------------
        # 准备匹配结果
        data_dict_all = data_dict.copy()
        
        ## 提取细胞数量和基因数量
        adata_cell_count = adata.n_obs
        adata_gene_count = adata.n_vars
        
        ## 匹配 cell_number 和 gene_number
        cell_match = adata_cell_count == row['Cell_count']
        gene_match = adata_gene_count == row.get('Gene_count', adata_gene_count)
        
        ## 匹配其他信息
        assay_match = match_column(adata.obs['assay'], row['Assay'])
        tissue_match = match_column(adata.obs['tissue'], row['Tissue'])
        cell_type_match = match_column(adata.obs['cell_type'], row['Cell_type'])
        sex_match = match_column(adata.obs['sex'], row['Sex'])
        disease_match = match_column(adata.obs['disease'], row['Disease'])
        development_stage_match = match_column(adata.obs['development_stage'], row['Development_stage'])
        organism_match = match_column(adata.obs['organism'], row['Organism'])
        suspension_type_match = match_column(adata.obs['suspension_type'], row['Suspension_type'])
        
        ## 储存匹配结果
        match_results = {
            "Cell_number_match": cell_match,
            "Gene_number_match": gene_match,
            "Assay_match": assay_match,
            "Tissue_match": tissue_match,
            "Cell_type_match": cell_type_match,
            "Sex_match": sex_match,
            "Disease_match": disease_match,
            "Development_stage_match": development_stage_match,
            "Organism_match": organism_match,
            "Suspension_type_match": suspension_type_match
        }
        
        ## 识别不匹配项
        mismatches = [name for name, matched in match_results.items() if not matched]
        
        ## 打印结果
        if not mismatches:
            print("All meta info matches with 05_merged_data.")
        else:
            print("The following information does not match:", ", ".join(mismatches))
        
        ## 将匹配结果添加到 data_dict_all
        data_dict_all.update(match_results)
        
        # ------- 4. 将数据赋值到 meta 的当前行 -------------
        for key, value in data_dict_all.items():
            if key in new_columns_dtypes:
                dtype = new_columns_dtypes[key]
                if pd.isna(value):
                    meta.at[index, key] = pd.NA if dtype != 'float' else np.nan
                else:
                    if dtype == 'Int64':
                        meta.at[index, key] = int(value)
                    elif dtype == 'float':
                        meta.at[index, key] = float(value)
                    elif dtype == 'boolean':
                        meta.at[index, key] = bool(value)
                    else:
                        meta.at[index, key] = value
            else:
                meta.at[index, key] = value

    except Exception as e:
        print(f"Error processing {file}: {e}")
        meta.at[index, 'processing_error'] = str(e)
        continue  # 继续处理下一个文件

# 4. 保存更新后的 meta 数据框
meta.to_csv("./08_check_info_all.csv", sep="\t", index=False)

# 选择需要的列保存为 08_check_info.csv
selected_columns = meta.columns  # 如果需要特定的列，可以在此修改
meta[selected_columns].to_csv("./08_check_info.csv", sep="\t", index=False)

# 5. 生成需要拆分和不需要拆分的列表
## 不需要拆分的数据
no_split_list = meta[
    (meta['Organism'].str.contains("Homo sapiens", na=False)) &
    (meta['Disease'].str.contains("normal", na=False))
]
no_split_list.to_csv("./08_no_split.csv", sep="\t", index=False)

## 需要拆分的数据
split_list = meta[~(
    (meta['Organism'].str.contains("Homo sapiens", na=False)) &
    (meta['Disease'].str.contains("normal", na=False))
)]
split_list.to_csv("./08_split_list.csv", sep="\t", index=False)

print("Data processing complete. Files saved to './08_check_info_all.csv', './08_check_info.csv', './08_no_split.csv', and './08_split_list.csv'")

