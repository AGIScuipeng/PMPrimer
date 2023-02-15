'''
创建人员: Nerium
创建日期: 2022/08/31
更改人员: Nerium
更改日期: 2023/02/15
'''

from .piecedefine import *

from collections import Counter
import primer3
import math
import json
import time

'''
创建人员: Nerium
创建日期: 2022/12/09
更改人员: Nerium
更改日期: 2022/12/09
'''
#写入JSON文件
def write_json(fname, data) :
    json.dump(data, open(fname, 'w'), indent=4, ensure_ascii=True)

'''
创建人员: Nerium
创建日期: 2023/02/15
更改人员: Nerium
更改日期: 2023/02/15
'''
#在字符串中统计list各项
#使用此种方式最快，是O(len(seq)*len(l))，否则是O(len(seq)*(len(l)^2))
def list_count(seq, l) :
    return {i : seq.count(i) for i in l}

'''
创建人员: Nerium
创建日期: 2022/12/09
更改人员: Nerium
更改日期: 2022/12/09
'''
#根据序列，计算TM和同源二聚体
def calc_tm_hairpin_homod(seq) :
    tm = primer3.calcTm(seq)
    hpin = primer3.calcHairpin(seq)
    homod = primer3.calcHomodimer(seq)

    return tm, False if hpin.tm < 25 else True, False if homod.tm < 25 else True

'''
创建人员: Nerium
创建日期: 2022/08/31
更改人员: Nerium
更改日期: 2023/02/10
'''
#通过流程主类中的原始、比对序列信息，根据id、比对序列坐标找到原始序列坐标
#primer3的索引是从0开始的，返回的x是从1开始的
def pos_translate(align_seq, ppos) :
    origin_seq = align_seq.replace('-', '')
    x, y = 0, 0
    while y <= ppos :
        if align_seq[x] == origin_seq[y] : x += 1; y += 1
        else : x += 1
    return x

'''
创建人员: Nerium
创建日期: 2022/08/31
更改人员: Nerium
更改日期: 2022/09/08
'''
#计算对比序列区间的保守度：延续法
def calc_conserve_continue(seqdict, posl, posr, mem, rate) :
    #if 1 <= posr <= 5 :print(posr, end=' | ')
    seqcnt, mem[1], tmp = len(seqdict.values()), mem[0], 0.0
    mem[0] += Counter([seq[posr-1] for seq in seqdict.values()]).most_common(1)[0][1]

    for ibp in range(posr, max(posl,posr-5) if posr != posl else posr+1, -1 if posr != posl else 1) :
        try : tmp += Counter([seq[ibp-1] for seq in seqdict.values()]).most_common(1)[0][1]
        except : break
    #if 1 <= posr <= 5 :print(posr, end=' | ')
    #if 1 <= posr <= 5 :print(ibp, (tmp/seqcnt)/((posr-ibp+1) if posr != posl else (ibp-posr+1)) >= rate)

    return (mem[0]/seqcnt)/(posr-posl+1) >= rate and (tmp/seqcnt)/((posr-ibp+1) if posr != posl else (ibp-posr+1)) >= rate

'''
创建人员: Nerium
创建日期: 2022/08/31
更改人员: Nerium
更改日期: 2022/09/18
'''
#计算对比序列区间的保守度：中断法（不到阈值立刻终止）
def calc_conserve_termina(seqdict, posl, posr, mem, rate) :
    seqcnt = len(seqdict.values())
    mem[0] = Counter([seq[posr-1] for seq in seqdict.values()]).most_common(1)[0][1]

    return (mem[0]/seqcnt) >= rate

'''
创建人员: Nerium
创建日期: 2022/08/31
更改人员: Nerium
更改日期: 2022/09/04
'''
#计算列表的标准差
def calc_std(num_list) :
    list_len = len(num_list)
    return math.sqrt(sum([(x - sum(num_list) / list_len) ** 2 for x in num_list]) / list_len)

'''
创建人员: Nerium
创建日期: 2022/08/31
更改人员: Nerium
更改日期: 2022/09/08
'''
def calc_cv_value(num_list) :
    return calc_std(num_list)/(sum(num_list)/len(num_list))

'''
创建人员: Nerium
创建日期: 2022/08/31
更改人员: Nerium
更改日期: 2022/09/08
'''
def calc_n_divide_average(num_list) :
    return len(num_list)**2/sum(num_list)

'''
创建人员: Nerium
创建日期: 2022/08/31
更改人员: Nerium
更改日期: 2022/09/20
'''
#列表作为一个系统，计算其香农熵
def calc_shannon_entropy(num_list) :
    return -sum([p*math.log(p, len(num_list)) for p in num_list if p])

'''
创建人员: Nerium
创建日期: 2022/08/31
更改人员: Nerium
更改日期: 2022/09/28
'''
#计算对比序列区间的保守度：香农熵终端和窗口延续法
def calc_conserve_termina_shannon(seqdict, posl, posr, mem, rate, window=1) :
    seqcnt, mem[1], mem[0] = len(seqdict), mem[0], sum(seqdict[max(posl-1, posr-window):posr])/(max(posr-max(posl, posr-window), 1))
    return mem[0] <= rate

'''
创建人员: Nerium
创建日期: 2022/08/31
更改人员: Nerium
更改日期: 2022/09/30
'''
#生成阈值最低香农熵
def generate_shannon_bynum(threshold) :
    #return round(calc_shannon_entropy([threshold, 1-threshold, 0, 0, 0]), 1)
    return calc_shannon_entropy([threshold, 1-threshold, 0, 0, 0])

'''
创建人员: Nerium
创建日期: 2022/08/31
更改人员: Nerium
更改日期: 2022/09/18
'''
#根据列表1对两个列表进行排序
def rank_lists_byfirst(list1, list2, reverse=False) :
    return [list(x) for x in zip(*(sorted(zip(list1, list2), key=lambda x: (x[0], x[1]), reverse=reverse)))]

'''
创建人员: Nerium
创建日期: 2022/11/18
更改人员: Nerium
更改日期: 2022/11/18
'''
#对两个对齐的序列计算差异分数
def mismatch_btw_2seq(seq1, seq2) :
    if len(seq1) != len(seq2) : return -1

    ret = [1 for i in range(len(seq1)) if seq1[i] != seq2[i]]
    return len(ret)

'''
创建人员: Nerium
创建日期: 2023/02/10
更改人员: Nerium
更改日期: 2023/02/10
'''
#计算区间的占位符个数
def calc_gaps_in_region(seq, ppos, llen, reverse=False) :
    tlen, gaps = llen, 0
    if reverse :
        for bpi in seq[:ppos][::-1] :
            if bpi != '-' : tlen -= 1
            else : gaps += 1

            if tlen == 0 : return gaps
        return 0

    for bpi in seq[ppos-1:] :
        if bpi != '-' : tlen -= 1
        else : gaps += 1

        if tlen == 0 : return gaps
    return 0

'''
创建人员: Nerium
创建日期: 2022/10/27
更改人员: Nerium
更改日期: 2022/10/27
'''
#从字符串中分离出id，属，种，亚种
def split_all_from_str(string) :
    strsplit = string.split(' ')

    if strsplit[1] == 'Unclassified' or len(strsplit) < 4 : return None

    if strsplit[3] == 'subsp.' or strsplit[3] == 'variant' :
        return strsplit[0], strsplit[1], strsplit[2], strsplit[4]

    return strsplit[0], strsplit[1], strsplit[2], strsplit[2]

'''
创建人员: Nerium
创建日期: 2022/08/31
更改人员: Nerium
更改日期: 2022/12/09
'''
#基础模块，log等功能都在其中
class piecebase() :
    def __init__(self, level=BASE_DEBUG_LEVEL0) -> None:
        self._level = max(BASE_DEBUG_LEVEL0, min(BASE_DEBUG_LEVEL3, level))
        self._time = str(int(time.time()))

    '''
    创建人员: Nerium
    创建日期: 2022/08/31
    更改人员: Nerium
    更改日期: 2022/11/16
    '''
    def baselog(self, msg, ends='\n') :
        print(msg, end=ends)

    def debuglog(self, setlevel, msg, ends='\n') :
        if setlevel & self._level : print('\033[0;33;40m{}\033[0m'.format(msg), end=ends)

    def successlog(self, msg, ends='\n') :
        print('\033[0;32;40m{}\033[0m'.format(msg), end=ends)

    def warnlog(self, msg, ends='\n') :
        print('\033[0;36;40m{}\033[0m'.format(msg), end=ends)

    def errorlog(self, msg, ends='\n') :
        raise SystemExit('\033[0;31;40m{}\033[0m'.format(msg))