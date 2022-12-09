'''
创建人员: Nerium
创建日期: 2022/09/29
更改人员: Nerium
更改日期: 2022/12/09
'''

from .piecedefine import *
from .piecebase import calc_tm_hairpin_homod, split_all_from_str
import primer3

'''
创建人员: Nerium
创建日期: 2022/09/29
更改人员: Nerium
更改日期: 2022/11/29
'''
class pieceevaluate() :
    def __init__(self, pbase, nonconser_sort, conser, primer_dict, seqdict, evaluate_opt) -> None:
        self._nonconser_sort = nonconser_sort
        self._conser = conser
        self._primer_dict = primer_dict
        self._seqdict = seqdict
        self.__evaluate_opt = evaluate_opt

        self._base = pbase

    '''
    创建人员: Nerium
    创建日期: 2022/10/08
    更改人员: Nerium
    更改日期: 2022/11/30
    '''
    #根据条件从区间中过滤出合适的保守区间和非保守区间(conser1)nonconser(conser2)
    #nonconser_sort是根据多样性由高到低排序的，conser是根据位置排序的
    def filter_area(self, minlen=None, hpcnt=None) :
        if minlen is None : minlen = self.__evaluate_opt['minlen']
        if hpcnt is None : hpcnt = self.__evaluate_opt['hpcnt']

        self._base.baselog('\n正在生成扩增子候选...')
        posmem = []
        for area in self._nonconser_sort :
            #if area[1]-area[0] < minlen : continue
            for idx, rang in enumerate(self._conser[1:], start=1) :
                #判断保守区间是否包含住了当前非保守区间
                if rang[0] > area[1] and self._conser[idx-1][1] < area[0] :
                    #判断前保守区间的F引物和后保守区间的R引物的haplotype个数是否小于最大值
                    if (len(self._primer_dict[self._conser[idx-1][0]][0].keys()) < hpcnt and len(self._primer_dict[rang[0]][1].keys()) < hpcnt) :
                        #判断扩增子长度是否满足
                        if rang[1]-self._conser[idx-1][0] < minlen : break

                        #判断前保守区间的F引物和后保守区间的R引物是否存在（不存在的可能性较低）（可以和上面if合并，但是为了方便调试信息的输出，所以分开）
                        if len(self._primer_dict.get(self._conser[idx-1][0], ({}, {}))[0]) and len(self._primer_dict.get(rang[0], ({}, {}))[1]) :
                            posmem.append((self._conser[idx-1], self._conser[idx]))
                        else : self._base.debuglog(BASE_DEBUG_LEVEL1, '{0} 或 {1} 没有引物/{0} Or {1} No Primer.'.format(rang, self._conser[idx-1]))
                        break
                    else :
                        self._base.debuglog(BASE_DEBUG_LEVEL1, '{0} 或 {1} 多样性超过阈值/ {0} Or {1} Haplotype Overtake Threshold'.format(rang, self._conser[idx-1]))

        self._base.successlog('扩增子候选生成完毕 共{}个'.format(len(posmem)))
        self._base.debuglog(BASE_DEBUG_LEVEL1, posmem)
        self._posmem = posmem
        return posmem

    '''
    创建人员: Nerium
    创建日期: 2022/11/30
    更改人员: Nerium
    更改日期: 2022/11/30
    '''
    #全排序获取合适的区间
    def full_permutation(self, minlen=None, maxlen=None, hpcnt=None) :
        if minlen is None : minlen = self.__evaluate_opt['minlen']
        if maxlen is None : maxlen = self.__evaluate_opt['maxlen']
        if hpcnt is None : hpcnt = self.__evaluate_opt['hpcnt']

        self._base.baselog('\n正在生成扩增子候选...')
        posmem = []
        for idi, area in enumerate(self._conser) :
            for idx, rang in enumerate(self._conser[idi+1:], start=idi+1) :
                if not(minlen < rang[1] - area[0] < maxlen) : continue

                #先判断前保守区间的F引物和后保守区间的R引物的haplotype个数是否小于最大值
                if len(self._primer_dict[area[0]][0].keys()) < hpcnt and len(self._primer_dict[rang[0]][1].keys()) < hpcnt :
                    #再判断前保守区间的F引物和后保守区间的R引物是否存在（可以和上面if合并，但是为了方便调试信息的输出，所以分开）
                    if len(self._primer_dict.get(area[0], ({}, {}))[0]) and len(self._primer_dict.get(rang[0], ({}, {}))[1]) :
                        posmem.append((self._conser[idi], self._conser[idx]))
                    else : self._base.debuglog(BASE_DEBUG_LEVEL1, '{0} 或 {1} 没有引物/{0} Or {1} No Primer.'.format(rang, area))
                else :
                    self._base.debuglog(BASE_DEBUG_LEVEL1, '{0} 或 {1} 多样性超过阈值/ {0} Or {1} Haplotype Overtake Threshold'.format(rang, area))

        self._base.successlog('扩增子候选生成完毕 共{}个'.format(len(posmem)))
        self._base.debuglog(BASE_DEBUG_LEVEL1, posmem)
        self._posmem = posmem
        return posmem

    '''
    创建人员: Nerium
    创建日期: 2022/10/11
    更改人员: Nerium
    更改日期: 2022/12/09
    '''
    #计算扩增子的覆盖度（目前是F、R计算亚种并取交集）
    def evaluate_cover_rate(self) :
        self._base.baselog('\n扩增子覆盖度为/ Cover Rate Of Amplicon：')
        rates = {}
        for amp in self._posmem :
            spdf, spdr, seqcnt = self._primer_dict[amp[0][0]], self._primer_dict[amp[1][0]], len(self._seqdict.values())

            if self.__evaluate_opt['rmlow'] :
                #根据TM值来排除不符合条件的引物，如果是一次设计的话，不要此参数否则会有错误
                for pri in list(spdf[0].keys()) :
                    if primer3.calcTm(pri) < self.__evaluate_opt['tm'] : spdf[0].pop(pri)
                for pri in list(spdr[1].keys()) :
                    if primer3.calcTm(pri) < self.__evaluate_opt['tm'] : spdr[1].pop(pri)

            #print([(len({'_'.join(split_all_from_str(s)[1:]) for v in spdf[0].values() for s in v if g in s}), len({'_'.join(split_all_from_str(s)[1:]) for v in spdr[1].values() for s in v if g in s})) for g in self._statistic_cnt[0]])
            #rates.setdefault('[{},{}]'.format(amp[0][0], amp[1][1]), [min(len({'_'.join(split_all_from_str(s)[1:]) for v in spdf[0].values() for s in v if g in s}), len({'_'.join(split_all_from_str(s)[1:]) for v in spdr[1].values() for s in v if g in s})) / len([True for s in self._statistic_cnt[2] if g in s]) for g in self._statistic_cnt[0]])

            self._base.debuglog(BASE_DEBUG_LEVEL1, ([amp[0][0], amp[1][1]], len({'_'.join(split_all_from_str(s)[1:]) for v in spdf[0].values() for s in v}), len({'_'.join(split_all_from_str(s)[1:]) for v in spdr[1].values() for s in v})))
            #rates.setdefault('[{},{}]'.format(amp[0][0], amp[1][1]), min(len({'_'.join(split_all_from_str(s)[1:]) for v in spdf[0].values() for s in v}), len({'_'.join(split_all_from_str(s)[1:]) for v in spdr[1].values() for s in v})))
            rates.setdefault('[{},{}]'.format(amp[0][0], amp[1][1]), len({'_'.join(split_all_from_str(s)[1:]) for v in spdf[0].values() for s in v} & {'_'.join(split_all_from_str(s)[1:]) for v in spdr[1].values() for s in v}))

        self._base.baselog('\n'.join(['{} : 亚种{:.2f}%'.format(k, v*100/len(self._statistic_cnt[2])) for k, v in rates.items()]))
        self._cover_rates = rates
        return rates

    '''
    创建人员: Nerium
    创建日期: 2022/10/12
    更改人员: Nerium
    更改日期: 2022/11/28
    '''
    #评估扩增子的分辨能力seq:set(species)
    def evaluate_resolution(self) :
        self._base.baselog('\n扩增子分辨力为/ Resolution Of Amplicon：')
        self._diversedict = {}
        reso, genuset, speset, subset = {}, set(), set(), set()
        merge_spe_set, merge_sub_set = set(), set()
        for amp in self._posmem :
            diverse1, diverse2, resdict = amp[0][1], amp[1][0], dict()
            for idallstr, seq in self._seqdict.items() : 
                #遍历过程中统计所有的种和亚种
                idsplit = split_all_from_str(idallstr)
                genuset.add(idsplit[1]); speset.add('_'.join(idsplit[1:3])); subset.add('_'.join(idsplit[1:4]))

                if seq[diverse1:diverse2] in resdict : resdict[seq[diverse1:diverse2]].add('_'.join(idsplit[1:]))
                else : resdict.setdefault(seq[diverse1:diverse2], {'_'.join(idsplit[1:]),})

            #请注意dict赋值存在直接指向、浅拷贝、深拷贝问题，如果后续resdict修改了，则此处需要改为深拷贝
            self._diversedict.setdefault('[{},{}]'.format(amp[0][0], amp[1][1]), resdict)
            self._base.debuglog(BASE_DEBUG_LEVEL3, resdict, ends='\n\n\n')

            #seq作为key，value只有一个种才计入分辨能力，亚种同理
            reso.setdefault('[{},{}]'.format(amp[0][0], amp[1][1]), ({list(v)[0].split('_')[0] for v in resdict.values() if len({s.split('_')[0] for s in list(v)}) == 1}, {'_'.join(list(v)[0].split('_')[:2]) for v in resdict.values() if len({'_'.join(s.split('_')[:2]) for s in list(v)}) == 1}, {list(v)[0] for v in resdict.values() if len(v) == 1}))

            if self.__evaluate_opt['merge'] == False : continue
            for v in resdict.values() :
                if len({'_'.join(s.split('_')[:2]) for s in list(v)}) == 1 :
                    merge_spe_set.add('_'.join(list(v)[0].split('_')[:2]))
                if len(v) == 1 :
                    merge_sub_set.add(list(v)[0])

        #在种和亚种的层次上都要统计分辨能力
        genuscnt, specnt, subcnt = len(genuset), len(speset), len(subset)
        self._base.debuglog(BASE_DEBUG_LEVEL1, '属数量：{0}；物种数量：{1}；亚种数量：{2}/ Genus Number: {0}; Species Number: {1}; Subspecies Number: {2}'.format(genuscnt, specnt, subcnt))

        self._base.baselog('\n'.join(['{} : 属 {}%; 种{:.2f}%; 亚种 {:.2f}%'.format(k, (len(v[0])/genuscnt)*100, (len(v[1])/specnt)*100, (len(v[2])/subcnt)*100) for k, v in reso.items()]))
        self._resolution = reso; self._statistic_cnt = (genuset, speset, subset)

        if self.__evaluate_opt['merge'] :
            self._base.baselog('合并后物种分辨力为/ Resolution Of Amplicon After Merged：{:.2f}%'.format(len(merge_spe_set)*100/len(speset)))
            self._base.baselog('合并后亚种分辨力为/ Resolution Of Amplicon After Merged：{:.2f}%'.format(len(merge_sub_set)*100/len(subset)))

        return reso

    '''
    创建人员: Nerium
    创建日期: 2022/12/09
    更改人员: Nerium
    更改日期: 2022/12/09
    '''
    #待选扩增子的引物
    def recommend_area_primer(self) :
        tmp_data = {}
        for amp in self._posmem :
            spdf, spdr = self._primer_dict[amp[0][0]], self._primer_dict[amp[1][0]]
            tmp_data.setdefault(str(amp), ({pri:calc_tm_hairpin_homod(pri) for pri in spdf[0].keys()}, {pri:calc_tm_hairpin_homod(pri) for pri in spdr[1].keys()}))

        self._base.warnlog('请注意引物设计和引物评估模块的的熔解温度可能存在差异 / Please Attention TM Of Primer Design Module & Primer Evaluate Maybe Exist Diff')

        return tmp_data