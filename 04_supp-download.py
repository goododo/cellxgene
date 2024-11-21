import pandas as pd
import os
import time
from threading import Thread
import csv
import requests
from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.support.ui import WebDriverWait
from concurrent.futures import ThreadPoolExecutor
import datetime
import json
import sys
args = sys.argv
if len(args) > 1:
    dataset_id_list = args[1]
    print(f"The dataset_id_list is: {dataset_id_list}")
else:
    print("No arguments provided .")



def chrome_file(url, file_path, color, file_size, count):
    try:
        file_name = url.split('?')[0].split('/')[-1]
        # 注意需要和 Chrome 浏览器设置的下载路径一致
        download_file_path = os.path.join(file_path, file_name)
        if not os.path.exists(download_file_path):
            options = Options()
            options.add_argument('--no-sandbox')
            options.add_argument('--disable-dev-shm-usage')
            options.add_argument('--headless')
            options.add_argument('blink-settings=imagesEnabled=false')
            web = webdriver.Chrome(options=options)
            # 配置浏览器
            web.command_executor._commands["send_command"] = ("POST", '/session/$sessionId/chromium/send_command')
            params = {'cmd': 'Page.setDownloadBehavior',
                      'params': {'behavior': 'allow', 'downloadPath': file_path}}
            web.execute("send_command", params=params)

            web.get(url)

            # 等待下载完成
            # 通过文件大小变化(在上面)或者定时等待的方式均可判断
            # 这里用定时等待的方式
            print(color % ('%s/%s 开始下载 -> ' % (count, file_name)))
            while not os.path.exists(download_file_path):
                time.sleep(8)
                try:
                    print(color % (str(os.path.getsize(download_file_path + '.crdownload')) + ' kb-(%s kb)' % file_size))
                except Exception as e:
                    pass
            web.quit()
            print('\033[31m%s/%s下载完成\033[0m' % (count, file_name))
        else:
            print('\033[31m%s/%s已存在\033[0m' % (count, file_name))
        return True, '成功'
    except Exception as e:
        pass
    return False, "下载文件失败"


def download_file(h5ad_url, count, color):
    for item in range(5):
        try:
            resp = requests.post(h5ad_url)
            data = resp.json()
            url = data.get('presigned_url', "")
            file_size = data.get('file_size', "")

            if url:
                file_path = os.path.join(os.getcwd(), str(count))
                if not os.path.exists(file_path):
                    os.mkdir(file_path)  # 创建文件夹
                status, msg = chrome_file(url, file_path, color, file_size, count)
                if not status:
                    with open('log.txt', 'a', encoding='utf8') as f:
                        f.write("%d -> %s" % (count, msg))
                        f.write('\n')
            return
        except Exception as e:
            print(e)
    return


def parse_detail(url):
    resp = requests.get(url)
    data_ls = resp.json()
    count = anntation['Download_ID'].max()
    dataframes = []
    dataframes2 = []
    with ThreadPoolExecutor(2) as t:  # 添加线程池
        for item in data_ls:
            dataset = item['name']
            collection_id = item['collection_id']
            data_key = dataset + '_' + collection_id
            #dataset_id = item['dataset_assets'][0]['dataset_id']
            if (data_key in difference_list ):
                count += 1
                first_id = item['id']
                collection_id = item['collection_id']
                dataset = item['name']
                dataset_id = item['dataset_assets'][0]['dataset_id']
                Assay = "，".join([item['label'] for item in item['assay']])
                cell_count = item['cell_count']
                cell_type = "，".join([item['label'] for item in item['cell_type']])
                development_stage = "，".join([item['label'] for item in item['development_stage']])
                disease = "，".join([item['label'] for item in item['disease']])
                organism = "，".join([item['label'] for item in item['organism']])
                published_at = item['published_at']
                date = datetime.datetime.fromtimestamp(published_at)
                self_reported_ethnicity = "，".join([item['label'] for item in item['self_reported_ethnicity']])
                sex = "，".join([item['label'] for item in item['sex']])
                suspension_type = item['suspension_type'][0]
                tissue = "，".join([item['label'] for item in item['tissue']])
                data = {'Download_ID': [count], 'data_key':[data_key] ,'Dataset' : [dataset] ,'Dataset_id': [dataset_id], 'Assay':[Assay] , 'Tissue':[tissue], 'Cell_type': [cell_type], 'Cell_count':[cell_count], 'Development_stage':[development_stage], 'Disease':[disease],'Organism':[organism],'Published':[date],'Self_reported_ethnicity':[self_reported_ethnicity],'Sex':[sex],'Suspension_type':[suspension_type], 'Collection_id': [collection_id],'Collection_id': [collection_id]  }
                temp_df = pd.DataFrame(data)
                dataframes.append(temp_df)
                dataset_assets = item['dataset_assets']
                for row in dataset_assets:
                    if row['filetype'] == 'H5AD':
                        h5ad_id = row['id']
                    if row['filetype'] == 'RDS':
                        rds_id = row['id']
                h5ad_url = 'https://datasets.cellxgene.cziscience.com/' + first_id + '.h5ad'
                #print(h5ad_url)
                rds_url = 'https://datasets.cellxgene.cziscience.com/' + first_id + '.rds' 
                #t.submit(download_file, h5ad_url, count, "\033[32m%s\033[0m")
                #t.submit(download_file, rds_url, count, "\033[34m%s\033[0m")
                data2 = { 'download_id' : [count], 'h5ad' : [h5ad_url], 'rds' : [rds_url] }
                temp_df2 = pd.DataFrame(data2)
                dataframes2.append(temp_df2)
                
            else:
                download_id = anntation['Download_ID'][anntation['data_key'] == data_key ].to_string(index=False)
                #print(download_id)
                first_id = item['id']
                collection_id = item['collection_id']
                dataset = item['name']
                dataset_id = item['dataset_assets'][0]['dataset_id']
                Assay = "，".join([item['label'] for item in item['assay']])
                cell_count = item['cell_count']
                cell_type = "，".join([item['label'] for item in item['cell_type']])
                development_stage = "，".join([item['label'] for item in item['development_stage']])
                disease = "，".join([item['label'] for item in item['disease']])
                organism = "，".join([item['label'] for item in item['organism']])
                published_at = item['published_at']
                date = datetime.datetime.fromtimestamp(published_at)
                self_reported_ethnicity = "，".join([item['label'] for item in item['self_reported_ethnicity']])
                sex = "，".join([item['label'] for item in item['sex']])
                suspension_type = item['suspension_type'][0]
                tissue = "，".join([item['label'] for item in item['tissue']])
                data = {'Download_ID': [download_id], 'data_key':[data_key] , 'Dataset' : [dataset] ,'Dataset_id': [dataset_id], 'Assay':[Assay] , 'Tissue':[tissue], 'Cell_type': [cell_type], 'Cell_count':[cell_count], 'Development_stage':[development_stage], 'Disease':[disease],'Organism':[organism],'Published':[date],'Self_reported_ethnicity':[self_reported_ethnicity],'Sex':[sex],'Suspension_type':[suspension_type], 'Collection_id': [collection_id],'Collection_id': [collection_id]  }
                temp_df = pd.DataFrame(data)
                dataframes.append(temp_df)
        df = pd.concat(dataframes, ignore_index=True)
        df.to_csv('04_download-anntation-new.csv', index=False,sep="\t",quoting=csv.QUOTE_NONE)
        df2 = pd.concat(dataframes2, ignore_index=True)
        df2.to_csv('04_supp_data_url.txt', index=False,sep="\t",quoting=csv.QUOTE_NONE)
        #return df

def main():
    url = "https://api.cellxgene.cziscience.com/dp/v1/datasets/index"
    parse_detail(url)
    #df=parse_detail(url)
    #df2 = pd.concat([anntation,df], ignore_index=True)
    #df2.to_csv('download-anntation-new.csv', index=False,sep="\t",quoting=csv.QUOTE_NONE,encoding='utf-8')
############################################################################################################
############################################################################################################
dataset_id_list = args[1]
anntation= pd.read_csv( dataset_id_list ,delimiter="\t")
anntation['data_key'] = anntation['Dataset']+'_'+ anntation['Collection_id']
dataset_id1 = set(anntation['data_key'])
line_count = len(anntation['data_key'])

print(f"已下载数据个数为: {line_count}")
url = "https://api.cellxgene.cziscience.com/dp/v1/datasets/index"
resp = requests.get(url)
data_ls = resp.json()
new_dataset = []
for item in data_ls:
    #dataset_id = item['dataset_assets'][0]['dataset_id']
    Dataset = item['name']
    collection_id = item['collection_id']
    data_key = Dataset + '_' + collection_id
    new_dataset.append(data_key)

line_count2 = len(new_dataset)
print(f"现在cellxgene数据库数据个数为: {line_count2}")


new_dataset_set = set(new_dataset)
difference = new_dataset_set - dataset_id1
difference_list = list(difference)
df = pd.DataFrame({"Col1": difference_list})
df[['Name', 'Collection_id']] = df['Col1'].str.extract(r'^(.*)_([^_]*)$')
df.to_csv("04_difference_list.csv", index=False,sep="\t",quoting=csv.QUOTE_NONE)
difference_length = len(difference_list)
difference_length

if __name__ == '__main__':
    if (difference_length>0) :
        print(f"共有{difference_length}个新增数据需要下载")
        main()
    else :
        print(f"没有新增数据需要下载")
    

