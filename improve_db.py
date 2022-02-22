import pandas as pd
import os
#bash: python3 json_to_db.py >ko00001.tab
#bash: grep -v "09190 Not Included in Pathway or Brite" ko00001.tab |grep -v "09180 Brite Hierarchies" >ko00001.1_6.tab 
# pathway + brite 
a = pd.read_csv('htext/ko00001.tab', sep='\t', header=None)
a = a.dropna() # Brite Hierarchies的第一行只有3列，没有对应的K号。
a.index = a.iloc[:, 3].apply(lambda x:x.split(' ')[0])
a.index.name = 'K number'
a.columns = ['level_1', 'level_2', 'level_3', 'level_4']

# brite.enzyme
b = pd.read_csv('htext/ko01000.tab', sep='\t', header=None)
b = b.dropna() # 去掉没有K号的行
b.index = b.iloc[:, 4].apply(lambda x:x.split(' ')[0]) # 总共5列
b.index.name = "K number"
b.columns = ['level_1', 'level_2', 'level_3', 'level_4', 'level_5']

all_K_nums = set(a.index) | set(b.index)

def arrange():
    indir = '/mnt/sdb1/home/lch/project/huming/jzl8/pathway_compare/kegg_annot'
    suffix = '.kegg.tab'
    infiles = [i for i in os.listdir(indir) if i.endswith(suffix)]
    series = []
    
    for fi in infiles:
        print('正在整理{}的注释信息...'.format(fi))
        k_nums = [i.split('\t')[1].rstrip() for i in open(indir + '/' + fi) if '\tK' in i]
        unknow_K_nums = set(k_nums) - all_K_nums
        # assert unknow_K_nums == set(), '数据库中的K号不全，{}个K号不在数据库中：{}'.format(str(len(unknow_K_nums)), ' '.join(unknow_K_nums))
        if unknow_K_nums != set():
            print('WARNING: 数据库中的K号不全，{}个K号不在数据库中：{}'.format(str(len(unknow_K_nums)), ' '.join(unknow_K_nums)))
        k_pathway = list(set(k_nums) & set(a.index))
        freq = a.loc[k_pathway,].level_2.value_counts()
        freq.name = fi.replace(suffix, '')
        series.append(freq)
    
    result = pd.concat(series, axis=1)
    result.to_csv('result.csv')

def explore(fi):
    indir = '/mnt/sdb1/home/lch/project/huming/jzl8/pathway_compare/kegg_annot'
    k_nums = [i.split('\t')[1].rstrip() for i in open(indir + '/' + fi) if '\tK' in i]
    k_pathway = list(set(k_nums) & set(a.index))
    return a.loc[k_pathway,]

x = explore('jzl8.kegg.tab')
