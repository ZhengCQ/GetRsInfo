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


rs_list = sys.argv[1]

class getRsInfo(object):
    def __init__(self, rs_id, work_dir, outfile, ref_config, ref_version, apikey='5bb720c52b556ce0164b88740f53b56a4d08'):
        self._rs_id = rs_id
        self._work_dir = work_dir
        self._ref_version = ref_version
        self._outfile = outfile
        self._key = apikey
        self._hg_scaf_info = self.readscaf(ref_config)
        self.main()

    def readscaf(self, file):
        '''
        读取配置文件，得到hg19和hg38的scaffold信息
        '''
        with open(file, 'r') as fh:
            hg_scaf_info = json.load(fh)
        return hg_scaf_info

    def requests_info(self, rs_ID, key):
        '''
        从NCBI E-utils 获取RS信息，如果因为链接不稳定失败，会重复最多5次。
        url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=snp&id=10745332&api_key=5bb720c52b556ce0164b88740f53b56a4d08'
        :return: NCBI返回的文本或者None
        '''
        url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=snp&id={rs_ID}&api_key={key}'
        t = 5
        while t > 0:
            try:
                response = requests.get(url, timeout=100)
                code = response.status_code
                if code == requests.codes.ok:
                    return response.text
                else:
                    self._logger.warning(f'URL: {url}; Response Code: {code}')
            except (ConnectionError, TimeoutError):
                time.sleep(0.5)
                continue
            except HTTPError as e:
                print(f"URL: {url}; Error: {e}")
                #break
                continue
            finally:
                t -= 1
        return None

    def parse_hgvs(self, hgvs_info):
        """
        ### 存在基因等信息 HGVS=NC_000008.11:g.18400285C>T,NC_000008.10:g.18257795C>T,NG_012246.1:g.14041C>T,NM_000015.2:c.282C>T,XM_017012938.1:c.282C>T|SEQ=[C/T]|GENE=NAT2:10
        ### 基因间区 HGVS=NC_000010.11:g.8983232T>C,NC_000010.10:g.9025195T>C|SEQ=[T/C]
        ### HGVS=NC_000003.12:g.169364845G>A,NC_000003.11:g.169082633G>A,NG_028279.1:g.303931C>T|SEQ=[G/A]|GENE=MECOM:2122
        ### 多个转录本？HGVS=NC_000020.11:g.54154632A>G,NC_000020.10:g.52771171A>G,NG_008334.1:g.24346T>C,NM_000782.4:c.*140T>C,NM_001128915.1:c.*140T>C,XM_005260304.5:c.*312T>C,XM_005260304.1:c.*312T>C,XM_017027691.2:c.*160T>C,XM_017027693.2:c.*312T>C|SEQ=[A/G]|GENE=CYP24A1:1591
        ### 具有多种突变NC_000009.12:g.97793827A>G,NC_000009.12:g.97793827A>T,NC_000009.11:g.100556109A>G,NC_000009.11:g.100556109A>T|SEQ=[A/G/T]|GENE=PTCSC2:101928337
        ### 具有多个基因symbol:NC_000011.10:g.111511840T>C,NC_000011.9:g.111382565T>C|SEQ=[T/C]|GENE=MIR34C:407042,MIR34B:407041,BTG4:54766,LOC728196:728196
        """
        def add_dict(dic_info, hg_ref, info_lst, category):
            if category == 'hgvs':
                scaf, mark, pos = info_lst
                if scaf in self._hg_scaf_info[hg_ref].keys():
                    dic_info.setdefault(hg_ref, {}).setdefault('chr', self._hg_scaf_info[hg_ref][scaf])
                    dic_info.setdefault(hg_ref, {}).setdefault('pos', pos)
            elif category == 'seq':
                ref = info_lst[0]
                alts = info_lst[1:]
                dic_info.setdefault(hg_ref, {}).setdefault('ref', ref)
                dic_info.setdefault(hg_ref, {}).setdefault('alt', alts)
            elif category == 'gene':
                dic_info.setdefault(hg_ref, {}).setdefault('gene', info_lst)

        dic_info = {}
        hgvs_parse = re.search(r'HGVS=(\S+)',hgvs_info).groups()[0].split('|')
        #parse hgvs突变，获取坐标
        #NC_000020.11:g.54154632A>G,NC_000020.10:g.52771171A>G,NG_008334.1:g.24346T>C,NM_000782.4:c.*140T>C,NM_001128915.1:c.*140T>C,XM_005260304.5:c.*312T>C,XM_005260304.1:c.*312T>C,XM_017027691.2:c.*160T>C,XM_017027693.2:c.*312T>C
        for each in hgvs_parse[0].split(','):
            if re.search(r':(g|m)\.', each):
                info_lst = re.match(r'(\S+):(\w).(\*{0,1}\d+)\w+', each).groups()[:]
                add_dict(dic_info, 'hg38', info_lst, 'hgvs')
                add_dict(dic_info, 'hg19', info_lst, 'hgvs')

        ## 取序列，ref和alt
        #SEQ=[A/G/T]
        #SEQ=[T/C]
        #SEQ=[TAATTATTAA/TAA/TAATTATTAATTATTAA]
        seq_lst = re.match(r'SEQ=\[(\S+)\]', hgvs_parse[1]).groups()[0].split('/')
        add_dict(dic_info, 'hg38', seq_lst, 'seq')
        add_dict(dic_info, 'hg19', seq_lst, 'seq')

        ## gene信息,不一定存在
        genelst = ['.']
        if len(hgvs_parse) > 2:
            genelst = re.match(r'GENE=(\S+)', hgvs_parse[2]).groups()[0].split(',')
            genelst = [i.split(':')[0] for i in genelst if not re.match(r'(^MIR|^LOC)', i)] #MIR和LOC开头的基因不要
        add_dict(dic_info, 'hg38', genelst, 'gene')
        add_dict(dic_info, 'hg19', genelst, 'gene')

        #将整个hgvs添加
        dic_info.setdefault('hg19',{}).setdefault('hgvs', hgvs_info)
        dic_info.setdefault('hg38',{}).setdefault('hgvs', hgvs_info)
        return dic_info

    def parse_xml(self, text):
        '''
        解析从网上获取的xml,或得hgvs信息，并转化为坐标
        '''
        try:
            text = text.encode('utf-8')
        except:
            pass
        xml = etree.XML(text)
        #position = xml.xpath('//DocSum/Item[@Name="CHRPOS"]/text()')[0]
        hgvs = xml.xpath('//DocSum/Item[@Name="DOCSUM"]/text()')[0]
        #hgvs = 'HGVS=NC_000010.11:g.8983232T>C,NC_000010.10:g.9025195T>C|SEQ=[T/C]'
        return hgvs

    def main(self):
        outf = open(self._outfile, 'a')
        url_info = self.requests_info(self._rs_id.replace('rs', ''), self._key)
        if url_info:
            hgvs = self.parse_xml(url_info)
            dic_info = self.parse_hgvs(hgvs)
            outf.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(self._ref_version, self._rs_id, dic_info[self._ref_version]['chr'],dic_info[self._ref_version]['pos'],\
                       dic_info[self._ref_version]['ref'],','.join(dic_info[self._ref_version]['alt']),\
                       ','.join(dic_info[self._ref_version]['gene']),dic_info[self._ref_version]['hgvs']))
            print(f'Got rs {self._rs_id} info into {self._outfile}')
        else:
            print(f"Error: {self._rs_id}")
        outf.close()