# 使用clusterProfiler进行KEGG富集分析时，使用enrichKEGG方法会遇到ID无法匹配的问题。
# 这时候使用通用的enricher可能会更加方便。
# 这个方法要求提供四个参数：
# 1. 感兴趣的基因列表
# 2. 背景基因列表
# 3. TERM to GENE dataframe
# 4. TERM to NAME dataframe

# 这个脚本用于得到3和4。
import re
import os
import sys

outdir = os.path.abspath(sys.argv[1])

htext = open(f"{outdir}/ko00001.tab")
term_to_name = open(f'{outdir}/TERM2NAME.csv', 'w')
term_to_gene = open(f'{outdir}/TERM2GENE.csv', 'w')
added_term = set()
no_map_id = 0
unexpected_cols = 0
for line in htext:
    line = line.strip().split('\t')
    level3 = line[2]
    try:
        name = re.search(r'.*?(?= \[)', level3).group()
        term = re.search(r'(?<=\[PATH:)[^\]]+', level3).group()
    except:
        # print(line)
        no_map_id += 1
        continue
    try:
        gene = line[3].split('  ')[0]
    except(IndexError):
        unexpected_cols += 1
        continue
    if term not in added_term:
        term_to_name.write(term + ',' + name + '\n')
        added_term.add(term)
    term_to_gene.write(term + ',' + gene + '\n')

print('line num without map id: ', no_map_id)
print('line num with unexpected columns: ', unexpected_cols)
