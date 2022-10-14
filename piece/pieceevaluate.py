'''
创建人员: Nerium
创建日期: 2022/09/29
更改人员: Nerium
更改日期: 2022/10/14
'''

from piece.piecedefine import *

'''
创建人员: Nerium
创建日期: 2022/09/29
更改人员: Nerium
更改日期: 2022/10/14
'''
class pieceevaluate() :
    def __init__(self, pbase, nonconser_sort, conser, primer_dict, seqdict) -> None:
        self._nonconser_sort = nonconser_sort
        self._conser = conser
        self._primer_dict = primer_dict
        self._seqdict = seqdict

        self._base = pbase

    '''
    创建人员: Nerium
    创建日期: 2022/10/08
    更改人员: Nerium
    更改日期: 2022/10/14
    '''
    #根据条件从区间中过滤出合适的保守区间和非保守区间(conser1)nonconser(conser2)
    #nonconser_sort是根据多样性由高到低排序的，conser是根据位置排序的
    def filter_area(self, minlen=80, hpcnt=10) :
        self._base.baselog('\n正在进行扩增子选定...')
        posmem = []
        for area in self._nonconser_sort :
            if area[1]-area[0] < minlen : continue

            for idx, rang in enumerate(self._conser) :
                #先判断保守区间是否包含住了当前非保守区间 以及 前保守区间的F引物和后保守区间的R引物的hypertpe个数是否小于最大值
                if rang[0] > area[1] and self._conser[idx-1][1] < area[0] and (len(self._primer_dict[self._conser[idx-1][0]][0].keys()) < hpcnt and len(self._primer_dict[rang[0]][1].keys()) < hpcnt) :
                    #再判断前保守区间的F引物和后保守区间的R引物是否存在（可以和上面if合并，但是为了方便调试信息的输出，所以分开）
                    if len(self._primer_dict.get(self._conser[idx-1][0], ({}, {}))[0]) and len(self._primer_dict.get(rang[0], ({}, {}))[1]) :
                        posmem.append((self._conser[idx-1], self._conser[idx]))
                    else : self._base.debuglog(BASE_DEBUG_LEVEL1, '{0} 或 {1} 没有引物/{0} Or {1} No Primer.'.format(rang, self._conser[idx-1]))
                    break

        self._base.successlog('已经完成扩增子选定 共{}个'.format(len(posmem)))
        self._base.debuglog(BASE_DEBUG_LEVEL1, posmem)
        self._posmem = posmem
        return posmem

    '''
    创建人员: Nerium
    创建日期: 2022/10/11
    更改人员: Nerium
    更改日期: 2022/10/14
    '''
    #计算扩增子的覆盖度（目前是计算所有输入序列的）
    def evaluate_cover_rate(self) :
        self._base.baselog('\n扩增子覆盖度为/ Cover Rate Of Amplicon：')
        rates = {}
        for amp in self._posmem :
            spdf, spdr, seqcnt = self._primer_dict[amp[0][0]], self._primer_dict[amp[1][0]], len(self._seqdict.values())
            rates.setdefault('[{},{}]'.format(amp[0][0], amp[1][1]), min(sum([len(v) for v in spdf[0].values()]), sum([len(v) for v in spdr[1].values()]))/seqcnt)

        self._base.baselog('\n'.join(['{} : {}%'.format(k, v*100) for k, v in rates.items()]))
        self._cover_rates = rates
        return rates

    '''
    创建人员: Nerium
    创建日期: 2022/10/12
    更改人员: Nerium
    更改日期: 2022/10/14
    '''
    #评估扩增子的分辨能力seq:set(species)
    def evaluate_resolution(self) :
        self._base.baselog('\n扩增子分辨力为/ Resolution Of Amplicon：')
        self._diversedict = {}
        reso, speset, subset = {}, set(), set()
        for amp in self._posmem :
            diverse1, diverse2, resdict = amp[0][1], amp[1][0], dict()
            for spe, seq in self._seqdict.items() : 
                #遍历过程中统计所有的种和亚种
                idsplit = spe.split('\t')[-1].split('_')
                speset.add('_'.join(idsplit[:2])); subset.add('_'.join(idsplit))

                if seq[diverse1:diverse2] in resdict : resdict[seq[diverse1:diverse2]].add('_'.join(idsplit))
                else : resdict.setdefault(seq[diverse1:diverse2], {'_'.join(idsplit),})

            #请注意dict赋值存在直接指向、浅拷贝、深拷贝问题，如果后续resdict修改了，则此处需要改为深拷贝
            self._diversedict.setdefault('[{},{}]'.format(amp[0][0], amp[1][1]), resdict)
            self._base.debuglog(BASE_DEBUG_LEVEL3, resdict, ends='\n\n\n')

            #seq作为key，value只有一个种才计入分辨能力，亚种同理
            reso.setdefault('[{},{}]'.format(amp[0][0], amp[1][1]), ({'_'.join(list(v)[0].split('_')[:2]) for v in resdict.values() if len({'_'.join(s.split('_')[:2]) for s in list(v)}) == 1}, {list(v)[0] for v in resdict.values() if len(v) == 1}))

        #在种和亚种的层次上都要统计分辨能力
        specnt, subcnt = len(speset), len(subset)
        self._base.debuglog(BASE_DEBUG_LEVEL1, '物种数量：{0}；亚种数量：{1}/ Species Number: {0}; Subspecies Number: {1}'.format(specnt, subcnt))

        self._base.baselog('\n'.join(['{} : 种{}% 亚种 {}%'.format(k, (len(v[0])/specnt)*100, (len(v[1])/subcnt)*100) for k, v in reso.items()]))
        self._resolution = reso; self._statistic_cnt = (specnt, subcnt)

        return reso