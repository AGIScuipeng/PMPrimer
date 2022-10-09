'''
创建人员: Nerium
创建日期: 2022/08/31
更改人员: Nerium
更改日期: 2022/10/09
'''

from piece.piecedefine import *
from piece.piecebase import calc_shannon_entropy, rank_lists_byfirst, generate_shannon_bynum
from piece.primerdesign import piecedesign
from piece.pieceevaluate import pieceevaluate

import os
from collections import Counter

#流程主类
class piecemain() :
    def __init__(self, args, pbase) -> None:
        self.args = args

        #获取到包的绝对路径，方便让subprocess调用muscle
        self._todo_path = '{}/static'.format(os.path.dirname(os.path.abspath(__file__)))

        #原始、对比数据的保存
        self._origindata = {}; self._origindata_shannon = []
        self._comparedata = {}; self._comparedata_shannon = []

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

    #根据多样性进行排名
    def rank_by_diverse(self, pcds, area, msg, logsw=False) :
        allshannon = []
        #计算区间的多样性
        for rang in area :
            areashannon = pcds.calc_area_diverse(self._comparedata_shannon if 'muscle' in self.args.alldesign else self._origindata_shannon, rang[0], rang[1])
            allshannon.append(areashannon)

        #根据多样性进行排名
        stdlist, arealist = rank_lists_byfirst(allshannon, area, reverse=True)
        if logsw :
            self._base.baselog('\n{}多样性排名为：/ Non Conservative Area Rank Is :'.format(msg))
            for idx, std in enumerate(stdlist) : self._base.baselog('[{}]\tScore : {};\tArea : {}'.format(idx+1, std, arealist[idx]))

        return arealist

    #调用primer3-py进行引物设计
    #可以先将序列去重再进行引物设计，但是保存原始信息会比较麻烦，快速开发先走流程
    def primer_design(self, pcds, area) :
        primer_dict, areacnt = {}, len(area)
        self._base.baselog('\n正在根据保守区间进行引物设计...', ends='')
        for numi, rang in enumerate(area, start=1) :
            self._base.baselog('\r正在根据保守区间进行引物设计... {}/{}'.format(numi, areacnt), ends='')

            data = self._comparedata if 'muscle' in self.args.alldesign else self._origindata
            #最后引物是dict中key:set()的形式，去重和保留原始样本名称信息
            primer_dict.setdefault(rang[0], (dict(), dict()))
            #如果是一个-区间取消，则使用temp + (temp=None)break处理后添加到primer_dict中

            for spe, seq in data.items() :
                if seq[rang[0]-1:rang[1]].count('-') : continue

                seq_args = {
                    'SEQUENCE_ID': '{}-{}'.format(rang[0], rang[1]),
                    'SEQUENCE_TEMPLATE': seq.replace('-', ''),
                    'SEQUENCE_INCLUDED_REGION': [rang[0]-seq[:rang[0]].count('-')-1, rang[1]-rang[0]],
                }
                opt_args = {
                    'PRIMER_MIN_SIZE':15,
                    'PRIMER_OPT_SIZE':((15+rang[1]-rang[0])//2) if rang[1]-rang[0]<35 else 25,
                    'PRIMER_MAX_SIZE':(rang[1]-rang[0]) if rang[1]-rang[0]<35 else 35,
                    'PRIMER_PRODUCT_SIZE_RANGE':[rang[1]-rang[0], 100],
                    'PRIMER_MIN_TM': 50.0,
                    'PRIMER_PICK_LEFT_PRIMER': 1,
                    'PRIMER_PICK_RIGHT_PRIMER': 1,
                    'PRIMER_NUM_RETURN': 1,
                }
                pair_primer = pcds.callprimer(target=seq_args, opt=opt_args)

                #将原始样本名称对应，保留原始信息
                for pri in pair_primer[0] :
                    if pri in primer_dict[rang[0]][0] : primer_dict[rang[0]][0][pri].add(spe)
                    else : primer_dict[rang[0]][0].setdefault(pri, {spe})
                for pri in pair_primer[1] :
                    if pri in primer_dict[rang[0]][1] : primer_dict[rang[0]][1][pri].add(spe)
                    else : primer_dict[rang[0]][1].setdefault(pri, {spe})

        self._base.successlog('\r\n已经根据保守区间完成引物设计')
        [self._base.debuglog(BASE_DEBUG_LEVEL2, '{} : {}, {}'.format(k,{kk:len(vv) for kk,vv in v[0].items()},{kk:len(vv) for kk,vv in v[1].items()})) if len(v[0])+len(v[1]) else self._base.debuglog(BASE_DEBUG_LEVEL2, '{} : None'.format(k)) for k, v in primer_dict.items()]
        [self._base.debuglog(BASE_DEBUG_LEVEL3, '{} : {}'.format(k,v)) if len(v[0])+len(v[1]) else self._base.debuglog(BASE_DEBUG_LEVEL3, '{} : None'.format(k)) for k, v in primer_dict.items()]

        return primer_dict

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
            conser = pcds.detect_conser_area_shannon(self._comparedata_shannon if 'muscle' in self.args.alldesign else self._origindata_shannon,
                                                    self._comparedata if 'muscle' in self.args.alldesign else self._origindata)
            self._base.baselog('保守区间列表 / List Of Conservative Area is : \n{0}'.format(conser))

            #挖掘出所有符合条件的非保守区间
            nonconser = pcds.detect_non_conser_area(self._comparedata_shannon if 'muscle' in self.args.alldesign else self._origindata_shannon, conser)
            self._base.baselog('非保守区间列表 / List Of Non Conservative Area is : \n{0}'.format(nonconser))

            #非保守区间多样性
            nonconser_sort = self.rank_by_diverse(pcds, nonconser, '非保守区间', 'rank1' in self.args.alldesign)

            #保守区间多样性
            conser_sort = self.rank_by_diverse(pcds, conser, '保守区间', 'rank2' in self.args.alldesign)

            #所有区间的多样性排名，保守和非保守分别排名是必须的，但是全排序不是必须的
            if 'rankall' in self.args.alldesign :
                self.rank_by_diverse(pcds, conser+nonconser, '所有区间', True)

            #根据hypertype进行分析和后续的引物设计
            self._alltype = {}
            self._base.baselog('\n保守区间的HyperType情况如下：')
            for rang in conser :
                alltype = pcds.detect_hypertype(self._comparedata if 'muscle' in self.args.alldesign else self._origindata, rang[0], rang[1])
                self._base.baselog('Area {}; \tLen : {}; \t {}'.format(rang, rang[1]-rang[0]+1, len(alltype)))
                self._alltype.setdefault(str(rang), alltype)

            if 'primer' in self.args.alldesign :
                self._primer_dict = self.primer_design(pcds, conser)

        #如果开启，则进行最终区域选择和评估等
        if self.args.evaluate is not None :
            try : pcel = pieceevaluate(self._base, nonconser_sort, conser, self._primer_dict)
            except : self._base.errorlog('\n未进行引物设计/Cannot Find Designed Primer')

            area_res = pcel.filter_area()
            self._base.baselog(area_res)