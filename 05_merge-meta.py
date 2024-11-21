import pandas as pd
import sys
import csv

anntation= pd.read_csv("04_download-anntation-new.csv",delimiter="\t")
anntation

collection=pd.read_csv("03_collection-author.csv",delimiter="\t")
collection

merge1=anntation.merge(collection, left_on='Collection_id',right_on='Collection_id')
merge1

if( len(merge1['Collection_id']) != len(anntation['Collection_id']) ):
    print("anntation and collection merged something wrong1")
    sys.exit(1)

dd= pd.read_csv("02_cellxgene-collections.csv",delimiter="\t")
dd= dd[['Dataset_id','Raw_data','Lab_website','Detail_url']]

merge2 = merge1.merge(dd, left_on=['Dataset_id'], right_on=['Dataset_id'] ,how='left')
merge2

# if( len(merge2['Dataset_id']) != len(anntation['Dataset_id']) ):
    # print("anntation and collection merged something wrong2")
    # sys.exit(1)

merge2.to_csv('05_merged_data.csv',encoding='utf-8', sep='\t', index=False)

# 创建子表，包含 Download_ID 和 下载链接
sub_table = merge2[['Download_ID', 'Dataset_id']].copy()
sub_table['h5ad'] = "https://datasets.cellxgene.cziscience.com/" + sub_table['Dataset_id'] + ".h5ad"

#保存子表为 txt 文件，制表符分隔
sub_table[['Download_ID', 'h5ad']].rename(columns={'Download_ID': 'download_id'}).to_csv(
    '05_download_urls.txt', sep='\t', index=False, header=True
)
