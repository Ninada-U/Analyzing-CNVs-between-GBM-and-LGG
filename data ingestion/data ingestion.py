import requests
import json
import pprint

fields = [
    'files.data_category',
    'files.data_format',
    'files.file_id'
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
#     pprint.pprint(case)
    vcfs = []
    if any((d['data_format'] == 'VCF' and \
            d['data_category'] == 'Simple Nucleotide Variation') for d in case['files']):
        for d in case['files']:
            if d['data_format'] == 'VCF':
                vcfs.append(d['file_id'])
        case_files[case['id']] = {'vcfs': vcfs}


import requests
import json
import re

data_endpt = "https://api.gdc.cancer.gov/data"

for case in list(case_files.keys())[:1]:
    ids = case_files[case]['vcfs']
    print(ids)


params = {"ids": ids}

response = requests.post(data_endpt,
                        data = json.dumps(params),
                        headers={
                            "Content-Type": "application/json"
                            })
                            
                            response_head_cd = response.headers["Content-Disposition"]

file_name = re.findall("filename=(.+)", response_head_cd)[0]

with open(file_name, "wb") as output_file:
    output_file.write(response.content)