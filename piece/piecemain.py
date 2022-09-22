'''
创建人员: Nerium
创建日期: 2022/08/31
更改人员: Nerium
更改日期: 2022/09/21
'''

from piece.piecebase import *
from piece.primerdesign import piecedesign

import os

#流程主类
class piecemain() :
    def __init__(self, args, pbase) -> None:
        self.args = args

        #获取到包的绝对路径，方便让subprocess调用muscle
        self._todo_path = '{}/static'.format(os.path.dirname(os.path.abspath(__file__)))

        #原始、对比数据的保存
        self._origindata = {}
        self._origindata_shannon = []
        self._comparedata = {}
        self._comparedata_shannon = []

        #基础模块的获取，log等功能都在其中
        self._base = pbase

    #原始数据的保存
    def getorigin(self) :
        if self.args.file is not None :
            with open(self.args.file, 'r') as ff :
                slt_id, slt_seq = '', ''
                for line in ff :
                    if '>' in line :  
                        if slt_id != '' and slt_seq != '' : self._origindata.update({slt_id: slt_seq})
                        slt_id = line[1:].strip(); slt_seq = ''
                    else : slt_seq += line.strip()
                if slt_id != '' and slt_seq != '' : self._origindata.update({slt_id: slt_seq})

            #如果不需要MUSCLE对齐，那么认为原始内容就是对齐的，所以直接计算香农熵
            if 'muscle' not in self.args.alldesign :
                seqlen, seqcnt = len(next(iter(self._origindata.values()))), len(self._origindata.values())
                for bp in range(seqlen) : self._origindata_shannon.append(calc_shannon_entropy([Counter([seq[bp] for seq in self._origindata.values()]).get(slg, 0)/seqcnt for slg in DEFAULT_DNA_SINGLE_LIST]))

    #对比数据的保存
    def aftercmp(self, pcds) :
        if pcds.tmpfile_path is not None :
            with open(pcds.tmpfile_path, 'r') as ff :
                slt_id, slt_seq = '', ''
                for line in ff :
                    if '>' in line :
                        if slt_id != '' and slt_seq != '' : self._comparedata.update({slt_id: slt_seq})
                        slt_id = line[1:].strip(); slt_seq = ''
                    else : slt_seq += line.strip()
                if slt_id != '' and slt_seq != '' : self._comparedata.update({slt_id: slt_seq})

            seqlen, seqcnt = len(next(iter(self._comparedata.values()))), len(self._comparedata.values())
            for bp in range(seqlen) : self._comparedata_shannon.append(calc_shannon_entropy([Counter([seq[bp] for seq in self._comparedata.values()]).get(slg, 0)/seqcnt for slg in DEFAULT_DNA_SINGLE_LIST]))

    #主流程函数
    def maintrunk(self) :
        #先保存原始数据
        self.getorigin()

        #如果开启，则调用muscle进行多序列比对；探查保守区间；循环论证最佳引物
        if self.args.alldesign is not None :
            pcds = piecedesign(self._base, self._todo_path, self.args.file)
            if 'muscle' in self.args.alldesign :
                pcds.callmuscle()

                #保存对比后的数据
                self.aftercmp(pcds)

            #挖掘出所有符合条件的保守区间
            conser = pcds.detect_conser_area_shannon(self._comparedata_shannon if 'muscle' in self.args.alldesign else self._origindata_shannon)
            self._base.baselog(BASE_DEBUG_LEVEL1, '保守区间列表 / List Of Conservative Area is : \n{0}'.format(conser))

            #挖掘出所有符合条件的非保守区间
            nonconser = pcds.detect_non_conser_area(self._comparedata_shannon if 'muscle' in self.args.alldesign else self._origindata_shannon, conser)
            self._base.baselog(BASE_DEBUG_LEVEL1, '非保守区间列表 / List Of Non Conservative Area is : \n{0}'.format(nonconser))

            if 'rank1' in self.args.alldesign :
                #非保守区间标准差
                allshannon = []
                for rang in nonconser :
                    areashannon = pcds.calc_area_diverse(self._comparedata_shannon if 'muscle' in self.args.alldesign else self._origindata_shannon, rang[0], rang[1])
                    allshannon.append(areashannon)

                #非保守区间的多样性进行排名
                stdlist, arealist = rank_lists_byfirst(allshannon, nonconser, reverse=True)
                self._base.baselog(BASE_DEBUG_LEVEL1, '\n非保守区间多样性排名为：/ Non Conservative Area Rank Is :')
                for idx, std in enumerate(stdlist) : self._base.baselog(BASE_DEBUG_LEVEL1, '[{}]\tScore : {};\tArea : {}'.format(idx+1, std, arealist[idx]))

            if 'rank2' in self.args.alldesign :
                #非保守区间标准差
                allshannon = []
                for rang in conser :
                    areashannon = pcds.calc_area_diverse(self._comparedata_shannon if 'muscle' in self.args.alldesign else self._origindata_shannon, rang[0], rang[1])
                    allshannon.append(areashannon)

                #非保守区间的多样性进行排名
                stdlist, arealist = rank_lists_byfirst(allshannon, conser, reverse=True)
                self._base.baselog(BASE_DEBUG_LEVEL1, '\n保守区间多样性排名为：/ Conservative Area Rank Is :')
                for idx, std in enumerate(stdlist) : self._base.baselog(BASE_DEBUG_LEVEL1, '[{}]\tScore : {};\tArea : {}'.format(idx+1, std, arealist[idx]))

            if 'rankall' in self.args.alldesign :
                #所有区间的多样性排名
                allarea = conser+nonconser
                allshannon = []
                for rang in allarea : 
                    areashannon = pcds.calc_area_diverse(self._comparedata_shannon if 'muscle' in self.args.alldesign else self._origindata_shannon, rang[0], rang[1])
                    allshannon.append(areashannon)
                stdlist, arealist = rank_lists_byfirst(allshannon, allarea, reverse=True)
                self._base.baselog(BASE_DEBUG_LEVEL1, '\n所有区间多样性排名为：/ All Area Rank Is :')
                for idx, std in enumerate(stdlist) : self._base.baselog(BASE_DEBUG_LEVEL1, '[{}]\tScore : {};\tArea : {}'.format(idx+1, std, arealist[idx]))

            #根据hypertype进行分析和后续的引物设计
            for rang in conser :
                alltype = pcds.detect_hypertype(self._comparedata if 'muscle' in self.args.alldesign else self._origindata, rang[0], rang[1])
                self._base.baselog(BASE_DEBUG_LEVEL1, 'Area {} : {}; \tLen : {}'.format(rang, len(alltype), rang[1]-rang[0]+1))
                #for t in alltype : self._base.debuglog(BASE_DEBUG_LEVEL1, t)

            '''
            for rang in conser[1:2] :
                if 'musle' in self.args.alldesign : data = self._comparedata
                else : data = self._origindata

                primer_list = []
                for seq in data.values() :
                    print(seq[rang[0]-1:rang[1]])
                    seq_args = {
                        'SEQUENCE_ID': '{}-{}'.format(rang[0], rang[1]),
                        'SEQUENCE_TEMPLATE': seq[rang[0]-1:rang[1]],
                        'SEQUENCE_INCLUDED_REGION': [0, rang[1]-rang[0]],
                        }
                    opt_args = {
                        'PRIMER_OPT_SIZE': rang[1]-rang[0],
                        'PRIMER_PICK_INTERNAL_OLIGO': 1,
                        'PRIMER_INTERNAL_MAX_SELF_END': 8,
                        'PRIMER_MIN_SIZE': 15,
                        'PRIMER_MAX_SIZE': rang[1]-rang[0],
                        'PRIMER_OPT_TM': 60.0,
                        'PRIMER_MIN_TM': 57.0,
                        'PRIMER_MAX_TM': 63.0,
                        'PRIMER_MIN_GC': 20.0,
                        'PRIMER_MAX_GC': 80.0,
                        'PRIMER_MAX_POLY_X': 100,
                        'PRIMER_INTERNAL_MAX_POLY_X': 100,
                        'PRIMER_SALT_MONOVALENT': 50.0,
                        'PRIMER_DNA_CONC': 50.0,
                        'PRIMER_MAX_NS_ACCEPTED': 0,
                        'PRIMER_MAX_SELF_ANY': 12,
                        'PRIMER_MAX_SELF_END': 8,
                        'PRIMER_PAIR_MAX_COMPL_ANY': 12,
                        'PRIMER_PAIR_MAX_COMPL_END': 8,
                        'PRIMER_PRODUCT_SIZE_RANGE': [0, rang[1]-rang[0]],
                    }
                    primer_list.append(pcds.callprimer(target=seq_args))#, opt=opt_args))

            print(primer_list)
            '''