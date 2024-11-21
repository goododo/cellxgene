import csv
import json
import requests
from concurrent.futures import ThreadPoolExecutor

CATE_DICT = {
    1: "collections",
    2: "datasets"
}


def write2csv(data_list, cate):
    with open('02_cellxgene-%s.csv' % cate, 'a', encoding='utf8', newline="") as f:
        writer = csv.writer(f,delimiter="\t")
        writer.writerow(data_list)


def parse_detail(p_url, p_name, p_id, cate, publisher):
    detail_url = "https://cellxgene.cziscience.com/collections/%s" % p_url.split('/')[-1]
    for item in range(3):
        try:
            headers = {
                "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/116.0.0.0 Safari/537.36",
                "Referer": "https://cellxgene.cziscience.com/"
            }
            resp = requests.get(p_url, headers=headers)
            pub_data = resp.json()
            consortia = ""
            if pub_data.get('consortia'):
                consortia = "，".join(pub_data.get('consortia', "")).replace(",", "，")
            description = pub_data.get('description', "").replace(",", "，").replace('\n', "")
            contact_name = pub_data.get('contact_name', "").replace(",", "，")
            links = pub_data.get('links', "")
            raw_data = ""
            law_website = ""
            other = []
            raw_data_list = []
            for link in links:
                link_type = link.get("link_type", "")
                if link_type == "LAB_WEBSITE":
                    law_website = link.get("link_url")
                elif link_type == "RAW_DATA":
                    raw_data = link.get("link_url")
                    raw_data_list.append(raw_data)
                else:
                    other.append(link)
            # 可能有很多个datasets
            raw_data_str = ';'.join(raw_data_list)
            for dataset in pub_data.get('datasets'):
                datasets = dataset.get('name', "").replace(",", "，")
                dataset_id = dataset['dataset_assets'][0]['dataset_id']
                tissue = ""
                if dataset.get('tissue'):
                    tissue = "，".join(
                        [item.get("label", "").replace(",", "，") for item in dataset.get('tissue')])
                disease = ""
                if dataset.get('disease'):
                    disease = "，".join(
                        [item.get("label", "").replace(",", "，") for item in dataset.get('disease')])
                assay = ""
                if dataset.get('assay'):
                    assay = "，".join([item.get("label", "").replace(",", "，") for item in dataset.get('assay')])
                organism = ""
                if dataset.get('organism'):
                    organism = "，".join([item.get("label", "").replace(",", "，") for item in dataset.get('organism')])
                cell_count = str(dataset.get('cell_count', "")).replace(",", "，")
                # print(consortia, description, contact_name, links, datasets, tissue)
                data_list = [p_name, p_id,consortia, description, publisher, contact_name, raw_data_str, law_website,
                             json.dumps(other), datasets,dataset_id, tissue, disease, assay, organism, cell_count, detail_url]
                write2csv(data_list, cate  )
            print(detail_url + " crawl over!")
            break
        except Exception as e:
            pass
    print(detail_url + " error -> %s" % str(e))


def main():
    # cate = CATE_DICT.get(
        # int(input(
            # """1. collections
# 2. Datasets
# 请输入序号: """)))
    cate= 'collections'
    with open('02_cellxgene-%s.csv' % cate, 'w', encoding='utf8', newline="") as f:
        writer = csv.writer(f,delimiter="\t")
        writer.writerow(
            ["Collection",'Collection_id', "Consortia", "Description", "Publication", "Contact_name", "Raw_data", "Lab_website", "Other",
             "Dataset","Dataset_id", "Tissue", "Disease", "Assay", "Organism", "Cell_count", "Detail_url"])
    print('start...')
    url = "https://api.cellxgene.cziscience.com/dp/v1/%s/index" % cate
    resp = requests.get(url)
    data_list = resp.json()
    with ThreadPoolExecutor(5) as t:
        for item in data_list:
            if cate == "datasets":
                p_id = item.get('collection_id')
            pub_obj = item.get("publisher_metadata", "")
            publisher = ""
            if pub_obj:
                journal = pub_obj.get('journal')
                published_year = pub_obj.get('published_year')
                family = pub_obj.get('authors')[0].get('family')
                publisher = f"{family} et al. ({published_year}) {journal}"
            p_name = item.get("name")
            p_id = item.get("id")
            p_url = 'https://api.cellxgene.cziscience.com/dp/v1/collections/%s' % p_id
            t.submit(parse_detail, p_url, p_name, p_id, cate, publisher)


if __name__ == '__main__':
    main()
