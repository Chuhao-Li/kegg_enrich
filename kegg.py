#!/usr/bin/env python3
# 参考文献：https://doi.org/10.1093/bioinformatics/bth456. # 这篇文献里面的超几何分布公式错了，分母中的i应该改为n。
# 2022-02-22: 增加命令行接口，考虑一个基因对多个K号的情况。

import sys
import os
import argparse
from math import factorial as fac
parser = argparse.ArgumentParser()
parser.description = "KEGG富集分析"
parser.add_argument('-l', '--list', help='待富集的基因列表')
parser.add_argument('-a', '--gene2ko', help='gene2ko对应表，tab分隔，一对多的情况的话多个ko间用逗号分隔。')
parser.add_argument('-d', '--db', help='KEGG数据库')
parser.add_argument('-s', '--slim', help='KEGG slim数据库')
parser.add_argument('-o', '--outfile', help='输出文件')
args = parser.parse_args()

keggdb = os.path.abspath(args.db)
slimdb = os.path.abspath(args.slim)
gene_list = open(args.list).read().strip().splitlines()
f = open(args.outfile, "w")

def com(m, n):
    # 组合公式
    out = (fac(n) // (fac(m)*fac(n-m)))
    return out

# term在这里的意思是通路Id。
print("# N represent total num of annotated genes in the whole genome", file=f)
print("# n represent total num of interested genes(Here is the differential expressed genes. )", file=f)
print("# k represent num interested genes that belong to a specific term", file=f)
print("# M represent num annotated genes that belong to a specific term", file=f)

gene2ko = {}
ko2term = {}
term2num = {} # term: M
term2name = {}

# 读取TERM2NAME文件
for line in open(f'{keggdb}/TERM2NAME.csv'):
    line = line.strip().split(',')
    term2name[line[0]] = line[1]

# 读取TERM2GENE文件
for line in open(f'{keggdb}/TERM2GENE.csv'):
    line = line.strip().split(',')
    if line[1] in ko2term:
        ko2term[line[1]].append(line[0])
    else:
        ko2term[line[1]] = [line[0]]

# 读取gene2ko文件，同时计算每个通路对应的背景基因的数量。
n_gene = 0
N = 0 # total num of genes that have a term annotation in the whole genome. 
gene2type = {"no_ko": set(), 
             "has_ko_but_no_map": set(), 
             "has_ko_and_has_map": set()
             }

for line in open(args.gene2ko):
    n_gene += 1
    line = line.strip('\n').split('\t')
    gene_id, kos = line
    if kos == '':
        gene2type["no_ko"].add(gene_id)
        continue
    kos = kos.split(',')
    kos1 = [] # 没有term的K号
    kos2 = [] # 有term的K号。
    for ko in kos:
        if ko in ko2term:
            kos2.append(ko)
        else:
            kos1.append(ko)
    if not kos2:
        gene2type["has_ko_but_no_map"].add(gene_id)
        continue
    else:
        gene2type["has_ko_and_has_map"].add(gene_id)
        N += 1
    # 一个基因可以对应多个K号，一个K号又可以对应多个TERM（通路）。
    gene2ko[gene_id] = kos2
    for ko in kos2:
        for term in ko2term[ko]:
            if term in term2num:
                term2num[term].add(gene_id)
            else:
                term2num[term] = set([gene_id])

term2num = {k: len(v) for k,v in term2num.items()}

n_with_no_KO_1 = len(gene2type["no_ko"])
n_with_no_map = len(gene2type["has_ko_but_no_map"])
print(f'# backgroupd: {n_gene} genes in total, while {n_with_no_KO_1} without KO annotation, {n_with_no_map} has no map, the left {N} genes were used for KEGG over representation analysis. ', file=f)

# 读取感兴趣的基因列表，同时计算每个通路对应的instreated gene的数量。
term2num_2 = {}
n_with_no_KO = 0
n_with_no_map2 = 0
n_gene2 = 0
n = 0 # total num of interested genes(Here is the differential expressed genes. )
for i in gene_list:
    n_gene2 += 1
    if i in gene2type["no_ko"]:
        n_with_no_KO += 1
        continue
    elif i in gene2type["has_ko_but_no_map"]:
        n_with_no_map2 += 1
        continue
    elif i in gene2type["has_ko_and_has_map"]:
        n += 1
        kos = gene2ko[i]
    for ko in kos:
        for term in ko2term[ko]:
            if term in term2num_2:
                term2num_2[term].add(i)
            else:
                term2num_2[term] = set([i])

print(f'# interested gene list: {n_gene2} genes in total, while {n_with_no_KO} without KO annotation, {n_with_no_map2} without map, the left {n} genes were used for KEGG over representation analysis.', file=f)

organism_path = set(open(args.slim).read().strip().splitlines())

print("term_name\tterm\tk/n\tM/N\tp-value\tincluded genes", file=f)
for term, genes in term2num_2.items():
    k = len(genes)
    included_genes = ';'.join(genes)
    term_name = term2name[term]
    M = term2num[term]
    p0 = 0
    for i in range(k):
        # print("k, n, M, N", k, n, M, N)
        p0 += (com(i, M)*com(n-i, N-M))/(com(n, N))
    p = 1-p0
    if p < 0.05 and term_name.split()[0] in organism_path:
        print(f"{term_name}\t{term}\t{k}/{n}\t{M}/{N}\t{p}\t{included_genes}", file=f)
