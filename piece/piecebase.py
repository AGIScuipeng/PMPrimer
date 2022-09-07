'''
创建人员: Nerium
创建日期: 2022/08/31
更改人员: Nerium
更改日期: 2022/09/05
'''

from piece.piecedefine import *

from collections import Counter
import math

#通过流程主类中的原始、比对序列信息，根据id、比对序列坐标找到原始序列坐标
def originpos(mainc, id, finalpos) :
    x, y, originstr, finalstr = 0, 0, mainc._origindata[id], mainc._comparedata[id]
    while y < finalpos :
        if originstr[x] == finalstr[y] : x += 1; y += 1
        else : y += 1
    return x

#计算对比序列区间的保守度：延续法
def calc_conserve_continue(seqdict, posl, posr, mem) :
    seqcnt, mem[1] = len(seqdict.values()), mem[0]
    mem[0] += Counter([seq[posr-1] for seq in seqdict.values()]).most_common(1)[0][1]

    return (mem[0]/seqcnt)/(posr-posl+1)

#计算对比序列区间的保守度：中断法（不到阈值立刻终止）
def calc_conserve_termina(seqdict, posl, posr, mem, rate) :
    seqcnt = len(seqdict.values())
    mem[0] += Counter([seq[posr-1] for seq in seqdict.values()]).most_common(1)[0][1]

    return (mem[0]/seqcnt) >= rate

#计算列表的标准差
def calc_std(num_list) :
    list_len = len(num_list)
    return math.sqrt(sum([(x - sum(num_list) / list_len) ** 2 for x in num_list]) / list_len)

def calc_cv_value(num_list) :
    return calc_std(num_list)/(sum(num_list)/len(num_list))

def calc_n_divide_average(num_list) :
    return len(num_list)**2/sum(num_list)

#根据列表1对两个列表进行排序
def rank_lists_byfirst(list1, list2, reverse=False) :
    return [list(x) for x in zip(*(sorted(zip(list1, list2), key=lambda x: (x[0], x[1]), reverse=reverse)))]

#基础模块，log等功能都在其中
class piecebase() :
    def __init__(self, level=BASE_DEBUG_LEVEL0) -> None:
        self._level = int(level)

    def baselog(self, setlevel, msg, ends='\n') :
        if setlevel & self._level : print(msg, end=ends)

    def debuglog(self, setlevel, msg, ends='\n') :
        if setlevel & self._level : print('\033[0;33;40m{}\033[0m'.format(msg), end=ends)