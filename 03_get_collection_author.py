# -*- encoding:utf-8 -*-
# @time: 2023/5/12 9:23
# @author: ifeng
import requests
import csv
import pandas as pd

def main():
	url = 'https://api.cellxgene.cziscience.com/dp/v1/collections/index'
	resp = requests.get(url)
	data_list = resp.json()
	dataframes = []
	for item in data_list:
		collection = item['name']
		collection_id = item['id']
		consortia = item['consortia']
		#journal = item['publisher_metadata']['journal']
		try:
			journal = item['publisher_metadata']['journal']
		except Exception as e:
			journal = 'NA'
		try:
			published_year = item['publisher_metadata']['published_year']
		except Exception as e:
			published_year ='NA'
		try:
			author = "，".join([item['family'] for item in item['publisher_metadata']['authors']])
		except Exception as e:
			author='NA'
			#continue
		data = {'Collection': [collection], 'Collection_id' : [collection_id] , 'Consortia': [consortia], 'Journal' : [journal],'Published_year':[published_year], 'Author' : [author]}
		temp_df = pd.DataFrame(data)
		dataframes.append(temp_df)
	df = pd.concat(dataframes, ignore_index=True)
	df.to_csv('03_collection-author.csv', index=False,sep="\t",quoting=csv.QUOTE_NONE,encoding='utf-8')



if __name__ == '__main__':
	main()


# collection=pd.read_csv("collection-author.csv",delimiter="\t")
# anntation= pd.read_csv("download-anntation.csv",delimiter="\t")
# merge1=anntation.merge(collection, left_on='Collection_id',right_on='Collection_id')
# merge1.to_csv('merged_data2222.csv', sep='\t', index=False)
