#!/usr/bin/env python3

import os
import argparse
import pandas as pd
from make_kegg_db import build_kegg_db

script_dir = os.path.dirname(os.path.abspath(__file__))

def annotate(args):
    args.outdir = os.path.abspath(args.outdir)
    if not os.path.exists(args.outdir):
        os.mkdir(args.outdir)
    os.system(f"emapper.py  -i {args.seq} --itype proteins -o {args.outdir}/out --cpu 32")
    print("done. ")

def gene2x(args):
    args.outdir = os.path.abspath(args.outdir)
    if not os.path.exists(args.outdir):
        os.mkdir(args.outdir)
    data = pd.read_csv(args.table, sep="\t", skiprows=4, skipfooter=3)
    gene2ko = data[["#query", "KEGG_ko"]]
    gene2go = data[["#query", "GOs"]]
    gene2ko.loc[:, "KEGG_ko"] = gene2ko["KEGG_ko"].apply(lambda x: x.replace("ko:", ""))
    gene2ko.loc[:, "KEGG_ko"] = gene2ko["KEGG_ko"].apply(lambda x: x.replace("-", ""))
    gene2ko.to_csv(f"{args.outdir}/gene2ko.xls", index=False, header=False, sep="\t")
    gene2go.to_csv(f"{args.outdir}/gene2go.xls", index=False, header=False, sep="\t")
    print("done. ")

def kegg_db(args):
    '''
    构建kegg数据库。
    首先下载ko00001.json，然后转化成ko00001.tab，再生成
    TERM2GENE.csv和TERM2NAME.csv。
    如果不希望细菌的结果中出现人类相关的通路，就要slim一下。
    具体就是使用kegg_slim.py指定一个属名，然后该脚本会下载属对应的通路。
    '''
    build_kegg_db(args.outdir, genus=args.genus)
    print("done. ")

def kegg_enrich(args):
    ''' 执行富集分析
    调用kegg.py
    '''
    gene_list = args.list
    gene2ko = args.gene2ko
    keggdb = args.db
    slimdb = args.slim
    outfile = args.outfile
    os.system(f"python {script_dir}/kegg.py -l {gene_list} -a {gene2ko} -d {keggdb} -s {slimdb} -o {outfile}")
    print("done. ")
    

parser = argparse.ArgumentParser()
parser.description = "kegg和go富集分析相关"
subparsers = parser.add_subparsers(help='')
t1 = subparsers.add_parser('annotate', help='user emapper to annotate. Please ensure the Id consistency')
t1.add_argument('seq', help='input protein sequence')
t1.add_argument('-o', '--outdir', help='输出目录')
t1.set_defaults(func=annotate)

t2 = subparsers.add_parser('gene2x', help='extract gene2go and gene2ko map table from emapper output')
t2.add_argument('table', help='emapper output file')
t2.add_argument('-o', '--outdir', help='输出目录')
t2.set_defaults(func=gene2x)

t3 = subparsers.add_parser('kegg_db', help='build kegg database. ')
t3.add_argument('-g', '--genus', help='genus of target species')
t3.add_argument('-o', '--outdir', help='output directory')
t3.set_defaults(func=kegg_db)

t3 = subparsers.add_parser('kegg_enrich', help='kegg over representative analysis. ')
t3.add_argument('-l', '--list', help='待富集的基因列表')
t3.add_argument('-a', '--gene2ko', help='gene2ko对应表，tab分隔，一对多的情况的话多个ko间用逗号分隔。')
t3.add_argument('-d', '--db', help='KEGG数据库')
t3.add_argument('-s', '--slim', help='KEGG slim数据库')
t3.add_argument('-o', '--outfile', help='输出文件')
t3.set_defaults(func=kegg_enrich)

args = parser.parse_args()
args.func(args)

