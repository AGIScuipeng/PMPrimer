'''
创建人员: Nerium
创建日期: 2022/08/31
更改人员: Nerium
更改日期: 2022/09/01
'''

from piece.piecedefine import *

from collections import Counter

#通过流程主类中的原始、比对序列信息，根据id、比对序列坐标找到原始序列坐标
def originpos(mainc, id, finalpos) :
    x, y, originstr, finalstr = 0, 0, mainc._origindata[id], mainc._comparedata[id]
    while y < finalpos :
        if originstr[x] == finalstr[y] : x += 1; y += 1
        else : y += 1
    return x

def evaluate(seqdict, posl, posr) :
    allseq, frac = list(seqdict.values()), 0.0
    seqcnt = len(allseq)
    for i in range(posl, posl+1 if posr-posl == 0 else posr+1) :
        frac += (next(iter(Counter([seq[i-1] for seq in allseq]).values())) / seqcnt)
    return frac/(posr-posl+1)

#基础模块，log等功能都在其中
class piecebase() :
    def __init__(self, level=BASE_DEBUG_LEVEL0) -> None:
        self._level = level

    def baselog(self, setlevel, msg, ends='\n') :
        if setlevel & self._level : print(msg, end=ends)

    def debuglog(self, setlevel, msg, ends='\n') :
        if setlevel & self._level : print('\033[0;33;40m{}\033[0m'.format(msg), end=ends)