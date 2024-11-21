# -*- encoding:utf-8 -*-
# @time: 2023/5/7 22:00
# @author: ifeng
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
    count = 0
    dataframes = []
    dataframes2 = []
    with ThreadPoolExecutor(10) as t:  # 添加线程池
        for item in data_ls:
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
            data = {'Download_ID': [count], 'Dataset' : [dataset] ,'Dataset_id': [dataset_id], 'Assay':[Assay] , 'Tissue':[tissue], 'Cell_type': [cell_type], 'Cell_count':[cell_count], 'Development_stage':[development_stage], 'Disease':[disease],'Organism':[organism],'Published':[date],'Self_reported_ethnicity':[self_reported_ethnicity],'Sex':[sex],'Suspension_type':[suspension_type], 'Collection_id': [collection_id],'Collection_id': [collection_id]  }
            temp_df = pd.DataFrame(data)
            dataframes.append(temp_df)
            
            dataset_assets = item['dataset_assets']
            for row in dataset_assets:
                if row['filetype'] == 'H5AD':
                    h5ad_id = row['id']
                if row['filetype'] == 'RDS':
                    rds_id = row['id']
            h5ad_url = 'https://api.cellxgene.cziscience.com/dp/v1/datasets/' + first_id + '/asset/' + h5ad_id
            rds_url = 'https://api.cellxgene.cziscience.com/dp/v1/datasets/' + first_id + '/asset/' + rds_id
            #t.submit(download_file, h5ad_url, count, "\033[32m%s\033[0m")
            #t.submit(download_file, rds_url, count, "\033[34m%s\033[0m")
            h5ad_url = 'https://datasets.cellxgene.cziscience.com/' + first_id + '.h5ad'
            rds_url = 'https://datasets.cellxgene.cziscience.com/' + first_id + '.rds' 
            data2 = { 'download_id' : [count], 'h5ad' : [h5ad_url], 'rds' : [rds_url] }
            temp_df2 = pd.DataFrame(data2)
            dataframes2.append(temp_df2)
            
            
        df = pd.concat(dataframes, ignore_index=True)
        df.to_csv('01_download_anntation.csv', index=False,sep="\t",quoting=csv.QUOTE_NONE,encoding='utf-8')
        df2 = pd.concat(dataframes2, ignore_index=True)
        df2.to_csv('01_data_url.txt', index=False,sep="\t",quoting=csv.QUOTE_NONE)

def main():
    url = "https://api.cellxgene.cziscience.com/dp/v1/datasets/index"
    parse_detail(url)


if __name__ == '__main__':
    main()



