#!/usr/bin/env python3
# 从kegg数据库下载ko00001.json文件，该文件包含kegg通路图的层次结构。把该文件转化成容易处理的tab文件。
import json
import sys

infile = sys.argv[1]
a = json.load(open(infile))

def bo(c, index): # c是大洋葱，index如果为[3, 2, 1]表示子任务是要剥第1层的第3个，第二层的第二个，第三层的第1个小洋葱。
    if index:
        this = eval("c" + ''.join(["['children'][{}]".format(str(i)) for i in index]))
    else:
        this = c
    if 'children' in this:
        for i in range(len(this['children'])):
            # print(index + [i])
            bo(c, index + [i])
    else:
        out = []
        for i in range(len(index)):
            express = "c" + ''.join(["['children'][{}]".format(str(i)) for i in index[:i+1]]) + "['name']"
            out.append(eval(express))
        print('\t'.join(out))

bo(a, [])
