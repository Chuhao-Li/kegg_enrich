#!/usr/bin/env python3
# 1. 给定属名，获取该属下所有物种的KEGG organism代号。
# 2. 下载这些organism的通路。
# /mnt/sdb1/home/lch/database/KEGG/htext/Organisms.tab
# http://rest.kegg.jp/list/pathway/pae

import sys
import os

script_dir = os.path.dirname(os.path.abspath(__file__))
outdir = sys.argv[1]
genus = sys.argv[2]
outdir = os.path.abspath(outdir)

orgs = []
for line in open(f'{script_dir}/Organisms.tab'):
    line = line.strip().split('\t')
    # print(line)
    if genus in line[-2]:
        orgs.append(line[-1].split()[0])

if not os.path.exists(outdir):
    os.mkdir(outdir)

out = set()
for org in orgs:
    if not os.path.exists(f'{outdir}/{org}.pathway'):
        os.system(f'wget -O - "http://rest.kegg.jp/list/pathway/{org}" >{outdir}/{org}.pathway')
    for line in open(f'{outdir}/{org}.pathway'):
        term = line.split('\t')[0].split(org)[1]
        out.add(term)

with open(f'{outdir}/{genus}_pathway.list', 'w') as f:
    f.write('\n'.join(out))
