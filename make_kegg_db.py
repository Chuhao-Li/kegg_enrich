#!/usr/bin/env python3

import sys
import os

script_dir = os.path.dirname(os.path.abspath(__file__))

def build_kegg_db(outdir, genus=""):
    outdir = os.path.abspath(outdir)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    
    # 这个文件会下得比较慢，开日本节点的VPN可以加速。
    if not os.path.exists(f"{outdir}/ko00001.json"):
        os.system(f"wget -O {outdir}/ko00001.json 'https://www.kegg.jp/kegg-bin/download_htext?htext=ko00001&format=json&filedir='")
    
    if not os.path.exists(f"{outdir}/ko00001.tab"):
        os.system(f"python {script_dir}/json_to_db.py {outdir}/ko00001.json >{outdir}/ko00001.tab")

    if not (os.path.exists(f"{outdir}/TERM2GENE.csv") and os.path.exists(f"{outdir}/TERM2NAME.csv")):
        os.system(f"python {script_dir}/TERM2GENE.py {outdir}")
    
    if genus:
        os.system(f"python {script_dir}/kegg_slim.py {outdir}/slim {genus}")
