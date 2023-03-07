'''
创建人员: Nerium
创建日期: 2022/09/29
更改人员: Nerium
更改日期: 2023/03/03
'''

from .piecedefine import *
from .piecebase import calc_tm_hairpin_homod, generate_degene, generate_rep, split_all_from_str, write_json

'''
创建人员: Nerium
创建日期: 2022/12/14
更改人员: Nerium
更改日期: 2022/12/14
'''
def useret(fname) :
    from Bio.Blast import NCBIXML

    dup_dict = {}
    for record in NCBIXML.parse(open(fname)) :
        dup_dict.setdefault(record.query, {})
        for align in record.alignments :
            hsps_len = len(align.hsps)
            if hsps_len in dup_dict[record.query] : dup_dict[record.query][hsps_len].add(align.title)
            else : dup_dict[record.query].setdefault(hsps_len, {align.title, })
    return dup_dict

'''
创建人员: Nerium
创建日期: 2022/12/12
更改人员: Nerium
更改日期: 2023/03/02
'''
def json2fasta(data, fname) :
    with open('{}_tmp.fasta'.format(fname), 'w') as ff:
        for area, fr in data.items() :
            for idx, lv in enumerate([fr[1], fr[3]]) :
                for idy, pri in enumerate(list(lv.keys())) :
                    ff.write('>{}_{}_{}_{}\n'.format(area, DATA2JSON_REFELCT[idx], idy, pri))
                    if idx ==0 :ff.write(pri+'\n')
                    else : ff.write(''.join([DEFAULT_DNA_REFLECT_DICT.get(bp, '-') for bp in pri][::-1])+'\n')
    return '{}_tmp.fasta'.format(fname)

'''
创建人员: Nerium
创建日期: 2022/09/29
更改人员: Nerium
更改日期: 2023/03/03
'''
class pieceevaluate() :
    def __init__(self, pbase, nonconser_sort, conser, primer_dict, seqdict, evaluate_opt) -> None:
        self._nonconser_sort = nonconser_sort
        self._conser = conser
        self._primer_dict = primer_dict
        self._seqdict = seqdict
        self.__evaluate_opt = evaluate_opt
        self._effective_len = {}
        self._primer_degene = {'F': dict(), 'R': dict()}

        self._base = pbase

    '''
    创建人员: Nerium
    创建日期: 2023/02/24
    更改人员: Nerium
    更改日期: 2023/02/24
    '''
    #根据左右区间，算出平均有效长度，并保存待后使用
    def calc_aver_effec_len(self, posl, posr) :
        seq_len_list = []
        for seq in self._seqdict.values() :
            tseq = seq[posl-1:posr].replace('-', '')
            seq_len_list.append(len(tseq))
        self._effective_len['[{},{}]'.format(posl, posr)] = sum(seq_len_list)//len(seq_len_list)
        return self._effective_len['[{},{}]'.format(posl, posr)]

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
    更改日期: 2023/02/24
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
                if not(minlen < self.calc_aver_effec_len(area[0], rang[1]) < maxlen) : continue

                #先判断前保守区间的F引物和后保守区间的R引物的haplotype个数是否小于最大值
                if len(self._primer_dict[area[0]][0].keys()) < hpcnt and len(self._primer_dict[rang[0]][1].keys()) < hpcnt :
                    #再判断前保守区间的F引物和后保守区间的R引物是否存在（可以和上面if合并，但是为了方便调试信息的输出，所以分开）
                    if len(self._primer_dict.get(area[0], ({}, {}))[0]) and len(self._primer_dict.get(rang[0], ({}, {}))[1]) :
                        posmem.append((self._conser[idi], self._conser[idx]))
                    else : self._base.debuglog(BASE_DEBUG_LEVEL1, '{0} 或 {1} 没有引物/{0} Or {1} No Primer.'.format(area, rang))
                else :
                    self._base.debuglog(BASE_DEBUG_LEVEL1, '{0} 或 {1} 多样性{2}、 {3}超过阈值/ {0} Or {1} Haplotype Is {2}, {3} Overtake Threshold'.format(rang, area, len(self._primer_dict[area[0]][0].keys()), len(self._primer_dict[rang[0]][1].keys())))

        self._base.successlog('扩增子候选生成完毕 共{}个'.format(len(posmem)))
        self._base.debuglog(BASE_DEBUG_LEVEL1, posmem)
        self._posmem = posmem
        return posmem

    '''
    创建人员: Nerium
    创建日期: 2022/10/11
    更改人员: Nerium
    更改日期: 2023/03/07
    '''
    #计算扩增子的覆盖度（目前是F、R计算亚种并取交集）
    def evaluate_cover_rate(self) :
        #根据简并引物生成所有引物
        tmp_degene = {'F' : dict(), 'R': dict()}
        cover_record = {'F' : dict(), 'R' : dict()}
        for posi, (degene, haplo) in self._primer_degene['F'].items() :
            tmp_degene['F'].setdefault(posi, haplo)
            cover_record['F'].setdefault(posi, set())
        for posi, (degene, haplo) in self._primer_degene['R'].items() :
            tmp_degene['R'].setdefault(posi, [''.join([DEFAULT_DNA_REFLECT_DICT[bp] for bp in h[::-1]]) for h in haplo])
            cover_record['R'].setdefault(posi, set())

        #记录每条序列可以被覆盖的保守区位点
        for idstr, seq in self._seqdict.items() :
            tseq = seq.replace('-', '')
            for posi, pris in tmp_degene['F'].items() :
                for pri in pris :
                    if pri in tseq : 
                        temp = split_all_from_str(idstr)
                        if temp is None : break
                        cover_record['F'][posi].add(' '.join(temp[:4])); break
            for posi, pris in tmp_degene['R'].items() :
                for pri in pris :
                    if pri in tseq :
                        temp = split_all_from_str(idstr)
                        if temp is None : break
                        cover_record['R'][posi].add(' '.join(temp[:4])); break

        self._base.baselog('\n扩增子覆盖度为/ Cover Rate Of Amplicon：')
        rates = {}
        for amp in self._posmem :
            spdf, spdr, seqcnt = self._primer_dict[amp[0][0]], self._primer_dict[amp[1][0]], len(self._seqdict.values())

            if self.__evaluate_opt['rmlow'] :
                import primer3

                #根据TM值来排除不符合条件的引物，如果是一次设计的话，不要此参数否则会有错误
                for pri in list(spdf[0].keys()) :
                    if primer3.calcTm(pri) < self.__evaluate_opt['tm'] : spdf[0].pop(pri)
                for pri in list(spdr[1].keys()) :
                    if primer3.calcTm(pri) < self.__evaluate_opt['tm'] : spdr[1].pop(pri)

            #根据保守区位点获取交集（存在误差，相同三级的序列有未覆盖的情况下会统计成100%）
            cover_f, cover_r = cover_record['F'][amp[0][0]], cover_record['R'][amp[1][0]]
            cover_f = set([' '.join(split_all_from_str(f)[1:4]) for f in cover_f if split_all_from_str(f) is not None])
            cover_r = set([' '.join(split_all_from_str(r)[1:4]) for r in cover_r if split_all_from_str(r) is not None])
            coveres = cover_f & cover_r

            self._base.debuglog(BASE_DEBUG_LEVEL1, ([amp[0][0], amp[1][1]], len(cover_f), len(cover_r), len(coveres), len(self._statistic_cnt[2])))
            rates.setdefault('[{},{}]'.format(amp[0][0], amp[1][1]), len(coveres))

        self._base.baselog('\n'.join(['{}, 有效长度 {} bp : 亚种{:.2f}%'.format(k, self._effective_len[k], v*100/len(self._statistic_cnt[2])) for k, v in rates.items()]))
        self._cover_rates = rates
        return rates

    '''
    创建人员: Nerium
    创建日期: 2022/10/12
    更改人员: Nerium
    更改日期: 2023/02/23
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
                if idsplit is None : continue
                genuset.add(idsplit[1]); speset.add('_'.join(idsplit[1:3])); subset.add('_'.join(idsplit[1:4]))

                tseq = seq[diverse1-1:diverse2].replace('-', '')
                if tseq in resdict : resdict[tseq].add('_'.join(idsplit[1:]))
                else : resdict.setdefault(tseq, {'_'.join(idsplit[1:]),})

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

        self._base.baselog('\n'.join(['{}, 有效长度 {} bp : 属 {}%; 种{:.2f}%; 亚种 {:.2f}%'.format(k, self._effective_len[k], (len(v[0])/genuscnt)*100, (len(v[1])/specnt)*100, (len(v[2])/subcnt)*100) for k, v in reso.items()]))
        self._resolution = reso; self._statistic_cnt = (genuset, speset, subset)

        if self.__evaluate_opt['merge'] :
            self._base.baselog('合并后物种分辨力为/ Resolution Of Amplicon After Merged：{:.2f}%'.format(len(merge_spe_set)*100/len(speset)))
            self._base.baselog('合并后亚种分辨力为/ Resolution Of Amplicon After Merged：{:.2f}%'.format(len(merge_sub_set)*100/len(subset)))

        return reso

    '''
    创建人员: Nerium
    创建日期: 2022/12/09
    更改人员: Nerium
    更改日期: 2023/03/07
    '''
    #待选扩增子的引物
    def recommend_area_primer(self, dct=None) :
        tmp_data = {}
        for amp in self._posmem :
            spdf, spdr = self._primer_dict[amp[0][0]], self._primer_dict[amp[1][0]]

            #获取简并的结果
            fdegene, fhaplo = generate_degene(spdf[0], dct, self.__evaluate_opt['degene'])
            rdegene, rhaplo = generate_degene(spdr[1], dct, self.__evaluate_opt['degene'])
            self._primer_degene['F'].update({amp[0][0]: [fdegene, fhaplo]})
            self._primer_degene['R'].update({amp[1][0]: [rdegene, rhaplo]})
            #根据引物的原始数量进行排序
            fvalue = sorted({pri: (sum([dct.get(split_all_from_str(p)[0], 0)+1 for p in pset]) if dct is not None else len(pset), calc_tm_hairpin_homod(pri)) for pri, pset in spdf[0].items()}.items(), key=lambda z : z[1][0], reverse=True)
            rvalue = sorted({pri: (sum([dct.get(split_all_from_str(p)[0], 0)+1 for p in pset]) if dct is not None else len(pset), calc_tm_hairpin_homod(pri)) for pri, pset in spdr[1].items()}.items(), key=lambda z : z[1][0], reverse=True)
            tmp_data.setdefault(str(amp), (fdegene, {pri : info for pri, info in fvalue}, rdegene, {pri : info for pri, info in rvalue}))

        self._base.warnlog('请注意引物设计和引物评估模块的的熔解温度可能存在差异 / Please Attention TM Of Primer Design Module & Primer Evaluate Maybe Exist Diff')

        self._recommend = tmp_data
        return tmp_data

    '''
    创建人员: Nerium
    创建日期: 2022/12/12
    更改人员: Nerium
    更改日期: 2022/12/12
    '''
    def blast_db_search(self, jsonpath) :
        from base64 import b64encode as bencode
        import subprocess
        import json
        import os 

        #将结果转换为fasta作为查询序列
        data = json.load(open(jsonpath))
        query_name = json2fasta(data, self._base._time)

        #查看文件是否存在
        for pathx in self.__evaluate_opt['blast'] : 
            if os.path.exists(pathx) == False :
                self._base.errorlog('建库序列集合文件不存在/ Blast Use File Not Exsits')

        #文件名称算哈希作为新文件夹
        final_path_name = bencode(''.join(self.__evaluate_opt['blast']).encode()).decode()[:50]
        final_file_name = [bencode(f.encode()).decode() for f in self.__evaluate_opt['blast']]

        #blast建库
        self._base.baselog('\nBLAST建库中...')
        if os.path.exists('{}'.format(final_path_name)) == False :
            os.system('mkdir {0}'.format(final_path_name))

            for idx, fx in enumerate(self.__evaluate_opt['blast']) :
                cm = subprocess.Popen('makeblastdb -in {} -parse_seqids -dbtype nucl -out {}/{}'.format\
                    (fx, final_path_name, final_file_name[idx]),\
                    shell=True, stderr=subprocess.PIPE)
                cm.communicate()
        self._base.successlog('BLAST建库完毕')

        #搜索相似序列
        self._base.baselog('\nBLAST搜索中...')
        #blastn -db blastdb/ncbi_myco_all -query jsontmp.fasta -out pri_res.xml -outfmt 5 -num_alignments 100000 -task blastn-short -word_size 15 -evalue 1
        for idx, fx in enumerate(self.__evaluate_opt['blast']) :
            cm = subprocess.Popen('blastn -db {}/{} -query {} -out {}_blast_res_{}.xml -outfmt 5 -num_alignments 100000 -task blastn-short -word_size 15 -evalue 1'.format\
                (final_path_name, final_file_name[idx], query_name, self._base._time, idx),\
                shell=True, stderr=subprocess.PIPE)
            cm.communicate()
        self._base.successlog('BLAST搜索完毕')

        #根据blast的结果统计分布并保存
        analysis_rate = []
        for idx in range(len(self.__evaluate_opt['blast'])) :
            analysis_rate.append({})
            for id, distri in useret('{}_blast_res_{}.xml'.format(self._base._time, idx)).items() :
                all_cnt = sum([len(v) for v in distri.values()])
                analysis_rate[idx].setdefault(id, {i[0]:i[1] for i in sorted({hl : len(aset)/all_cnt for hl, aset in distri.items()}.items())})
        if self.__evaluate_opt['save'] : write_json('{}_blast_analysis.json'.format(self._base._time), analysis_rate)

        #根据唯一命中的概率是否小于50%来删除不达标的引物
        for idx in range(len(analysis_rate)) :
            for id, statis in analysis_rate[idx].items() :
                tmp = id.split('_')
                del_flag = True if statis.get(1, 0) < 0.5 and len(statis.keys()) else False

                if del_flag == False : continue
                pri = tmp[-1] if tmp[1] == 'F' else ''.join([DEFAULT_DNA_REFLECT_DICT.get(bp, '-') for bp in tmp[-1]][::-1])
                if pri in self._recommend[tmp[0]][JSON2DATA_REFELCT[tmp[1]]] : self._recommend[tmp[0]][JSON2DATA_REFELCT[tmp[1]]].pop(pri)

        return self._recommend