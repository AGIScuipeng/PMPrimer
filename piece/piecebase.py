'''
创建人员: Nerium
创建日期: 2022/08/31
更改人员: Nerium
更改日期: 2022/09/02
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
    #allseq = list(seqdict.values())
    seqcnt = len(seqdict.values())
    mem[0] += next(iter(Counter([seq[posr-1] for seq in seqdict.values()]).values()))
    '''
    for i in range(posl, posl+1 if posr-posl == 0 else posr+1) :
        frac += (next(iter(Counter([seq[i-1] for seq in allseq]).values())) / seqcnt)
    return frac/(posr-posl+1)
    '''
    return (mem[0]/seqcnt)/(posr-posl+1)

#计算对比序列区间的保守度：中断法（不到阈值立刻终止）
def calc_conserve_termina(seqdict, posl, posr, mem, rate) :
    #allseq = list(seqdict.values())
    seqcnt = len(seqdict.values())
    mem[0] += next(iter(Counter([seq[posr-1] for seq in seqdict.values()]).values()))
    '''
    for i in range(posl, posl+1 if posr-posl == 0 else posr+1) :
        frac += (next(iter(Counter([seq[i-1] for seq in allseq]).values())) / seqcnt)
    return frac/(posr-posl+1)
    '''
    return (mem[0]/seqcnt) >= rate

def calc_std(num_list) :
    return math.sqrt(sum([(x - sum(num_list) / len(num_list)) ** 2 for x in num_list]) / len(num_list))

#基础模块，log等功能都在其中
class piecebase() :
    def __init__(self, level=BASE_DEBUG_LEVEL0) -> None:
        self._level = int(level)

    def baselog(self, setlevel, msg, ends='\n') :
        if setlevel & self._level : print(msg, end=ends)

    def debuglog(self, setlevel, msg, ends='\n') :
        if setlevel & self._level : print('\033[0;33;40m{}\033[0m'.format(msg), end=ends)