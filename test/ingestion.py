import requests
import json
import pprint
import re
import patoolib

from pathlib import WindowsPath as Path
import shutil
import os
import time


with open("config.json", "r") as f:
    cfg = json.load(f)
    num_cases = cfg['num_cases']


fields = [
    'files.data_category',
    'files.data_format',
    'files.file_id',
    'files.file_name',
    'case_id'
    ]

fields = ",".join(fields)

endpt = "https://api.gdc.cancer.gov/cases"

filters = {
    "op": "and",
    "content":[
        {
            'op': 'in',
            'content': {
                'field': 'project.program.name',  # program
                'value': ['TCGA']
            }
        },
        {
            'op': 'in',
            'content': {
                'field': 'project.name',  # project
                # either 'Glioblastoma Multiforme' or 'Brain Lower Grade Glioma'
                'value': ['Glioblastoma Multiforme']
            }
        },
        {
            "op": "in",
            "content":{
                "field": "cases.project.primary_site",  # primary site
                "value": ["Brain"]
            }
        },
        {
            "op": "in",
            "content":{
                "field": "demographic.race",  # race
                "value": ["White"]
            }
        }
    ]
}

os.chdir("data/download")

params = {
    "filters": filters,
    "fields": fields,
    'sort': 'submitter_id:asc',
    'size': 1000  # should be >= number of hits outputted
    }
js = requests.post(endpt, headers = {"Content-Type": "application/json"}, json = params).json()
print('warnings: ', js['warnings'])
print('number of hits: ', js['data']['pagination']['total'])
# print(json.dumps(js, indent=2))
case_files = {}

for case in js['data']['hits']:
    cnvs = []
    for d in case['files']:
        if d['data_category'] == 'Copy Number Variation' and \
            not 'focal_score' in d['file_name'] and not '.tsv' in d['file_name']:
            cnvs.append(d['file_id'])
    if not cnvs:
        continue
    case_files[case['id']] = {'cnvs': cnvs}


data_endpt = "https://api.gdc.cancer.gov/data"

num_files = 0

ids = []
for case in list(case_files.keys())[:num_cases]:
    ids.append((case, case_files[case]['cnvs']))
    num_files += len(case_files[case]['cnvs'])
    
print("num cases: ", num_cases)
print("num files: ", num_files)
    
for id in ids:
    params = {"ids": id[1]}
    response = requests.post(data_endpt, data = json.dumps(params),
                            headers={"Content-Type": "application/json"})
                    
    response_head_cd = response.headers["Content-Disposition"]

    file_name = re.findall("filename=(.+)", response_head_cd)[0]
    
    with open(file_name, "wb") as output_file:
        output_file.write(response.content)
        
    patoolib.extract_archive(file_name)
    
    for i, d in enumerate(Path(file_name[:-3]).iterdir()):
        if d.is_dir():
            for f in d.iterdir():
                new_name = str(i) + '.txt'
                try:
                    os.rename(str(f), new_name)
                except FileExistsError:
                    pass
                Path("../raw/" + id[0]).mkdir(parents=True, exist_ok=True)
                try:
                    shutil.move(new_name, "../raw/" + id[0])
                except shutil.Error:
                    pass

import pandas as pd
from pathlib import WindowsPath as Path
import os

os.chdir("../raw")

# FUTURE TODO: don't create entire df to check column names, just do a one-line read
# FUTURE TODO: find a better way to filter the 'Major_Copy_Number' files before downloading

for d in Path().iterdir():
    # iter through case files
    dfs = []
    for f in d.iterdir():
        df = pd.read_csv(str(f), delimiter="\t")
        if not 'Major_Copy_Number' in df.columns:
            dfs.append(df)
    concatenated_df = pd.concat(dfs, axis=0)
    concatenated_df.to_csv("../merged/" + str(d) + ".txt", index=False)

                