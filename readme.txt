

#第一次下载数据使用01_download_file.py得到01_download_anntation.csv总览表格有重复的Collection_id是数据集id，和01_data_url.txt所有数据下载地址
python 01_download_file.py     


#使用02_crawl_detail.py获取数据的原始数据来源GSE等得到02_cellxgene-collections.csv有重复的Collection_id是数据集id  现在不需要参数
python 02_crawl_detail.py      


#使用03_get_collection_author.py获取数据的文献相关信息得到03_collection-author.csv没有重复的Collection_id是数据下载id
python 03_get_collection_author.py   


#使用04_supp-download.py生成后续更新数据时需要下载的新数据表 需要参数给定已下载最新数据表格的绝对路径
#输出04_difference_list.csv是需要下载的数据信息
#输出04_download-anntation-new.csv 包含了全部数据的信息 可作为下次下载的参数,每次的Download_ID都一致
#输出04_supp_data_url.txt是需要补充下载的数据下载地址
python 04_supp-download.py ../download_0401/download-anntation-new.csv    


#合并04_download-anntation-new.csv 02_cellxgene-collections.csv 03_collection-author.csv 得到05_merged_data.csv 包含了全部数据的信息,每次的Download_ID都一致
#同时输出05_download_urls.txt是全部数据的下载地址
python 05_merge-meta.py

#使用06_download_data.pl下载补充数据04_supp_data_url.txt或者全部数据05_download_urls.txt，根据需要修改脚本               

#检查下载完整并删掉下载好的log
nohup ./07_check_download.sh > 07_check_download_output.log 2>&1 & 

#对下载的数据集进行质量检查和统计分析，通过匹配 05_merged_data.csv 数据并验证 .h5ad 格式数据文件的结构来确保数据完整性。
#输出文件：
#1. 08_check_h5ad_data_info_all.csv：包含每个 .h5ad 文件的统计和结构信息，如细胞数、基因数以及 所有meta数据一致性检查结果。
#2. 08_check_h5ad_data_info.csv：包含每个 .h5ad 文件的统计和结构信息，如细胞数、基因数以及 关键meta数据一致性检查结果。
#3. 08_no_split_list_combine.csv：列出所有可以直接用的human normal data。
#4. 08_split_list_combine.csv：列出需要基于Organism 和 Disease 元数据拆分出human normal的数据集。
python 08_Check_and_Statistical_Analysis.py

#基于08_split_list_combine.csv表格提取Organism 和 Disease中混有其他的物种或疾病的human normal data保存到09_subset文件夹
python 09_subset_normal.py

#目前这一版下载是基于20241101的meta信息，因为匹配之前版本（有数据删除或修改）的Download_ID导致Download_ID不连续有中断
#1622个数据集，其中human数据集1214个，human normal数据集996个，human normal表达谱数据集975个
#176M cell，102M human cell, 94M human normal cell, 92M human normal cell 表达谱数据
#最后所有数据的meta和path：08_check_h5ad_data_info.csv，所有 human normal cell 表达谱数据的meta和path：09_hu_normal_meta.csv
