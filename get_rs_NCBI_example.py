#!/usr/bin/env python
# _*_ coding: utf-8 _*_

'''
获得 rs 列表内所有 rs 的 hg19 位置信息。
'''

import requests
import json
import time
import sys
import re
from lxml import etree
from requests import HTTPError, ConnectionError


from getRsInfo import getRsInfo
rs_list = sys.argv[1]

key = '5bb720c52b556ce0164b88740f53b56a4d08'
work_dir = 'get_rs_test'
ref_scaffold = 'config/ref_scaffold.json'
ref_version = 'hg19'

exsits_rs = []
with open('%s/chr_pos.txt'%(work_dir), 'r') as f:
    for each in f:
        rs_id = each.strip().split()[1]
        exsits_rs.append(rs_id)

with open(rs_list, 'r') as f:
    for each in f:
        rs_id = each.strip()
        if re.match(r'rs(\d+){1,10}', rs_id):
            print (f'{rs_id} not standard rs name, next!!')
            continue
        print(f'{rs_id} starting')
        if rs_id not in exsits_rs:
            getRsInfo(rs_id, work_dir,'%s/chr_pos.txt'%(work_dir), ref_scaffold, ref_version, apikey=key)
        else:
            print(f"{rs_id} exists in results file, next!!")
    print(f"All finished")


