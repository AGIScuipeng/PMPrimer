'''
创建人员: Nerium
创建日期: 2022/08/31
更改人员: Nerium
更改日期: 2023/03/03
'''

from .piecedefine import *
from .piecebase import calc_shannon_entropy, list_count, rank_lists_byfirst, generate_shannon_bynum, calc_tm_hairpin_homod, write_json, pos_translate
from .piecedesign import piecedesign
from .pieceevaluate import pieceevaluate
from .piecedataprogress import piecedataprogress

from collections import Counter
import os

'''
创建人员: Nerium
创建日期: 2022/08/31
更改人员: Nerium
更改日期: 2023/03/03
'''
#流程主类
class piecemain() :
    def __init__(self, args, pbase) -> None:
        self.args = args

        #获取到包的绝对路径，方便让subprocess调用muscle
        self._todo_path = '{}/static'.format(os.path.dirname(os.path.abspath(__file__)))

        #原始、对比数据的保存
        self._origindata = {}; self._origindata_shannon = []
        self._comparedata = {}; self._comparedata_shannon = []

        '''
        扩增子设计相关参数配置
        '''
        self.__design_opt = {'threshold' : generate_shannon_bynum(0.95), 'minlen' : 15, 'merge': False, 'pdetail': False, 'pdetail2' : False,\
             'primer' : False, 'primer2' : False, 'gaps' : 0.1, 'tm' : 50.0, 'window' : 1}
        if self.args.alldesign is not None :
            #遍历alldesign找到threshold:0.xx等参数
            for x in self.args.alldesign[::-1] :
                if 'threshold' in x : 
                    try : self.__design_opt['threshold'] = generate_shannon_bynum(float(x.split(':')[-1]))
                    except : self._base.warnlong('阈值参数解析错误/ Threshold Patameter Parse Error')
            #遍历alldesign找到minlen:15等参数
            for x in self.args.alldesign[::-1] :
                if 'minlen' in x : 
                    try : self.__design_opt['minlen'] = int(x.split(':')[-1])
                    except : self._base.warnlong('设计最小值参数解析错误/ All Design Minlen Patameter Parse Error')
            #遍历alldesign找到gaps:0.1等参数
            for x in self.args.alldesign[::-1] :
                if 'gaps' in x : 
                    try : self.__design_opt['gaps'] = float(x.split(':')[-1])
                    except : self._base.warnlong('空白符占比参数解析错误/ Gaps Rate Patameter Parse Error')
            #遍历alldesign找到tm:50等参数
            for x in self.args.alldesign[::-1] :
                if 'tm' in x : 
                    try : self.__design_opt['tm'] = float(x.split(':')[-1])
                    except : self._base.warnlong('熔解温度参数解析错误/ DNA Primer TM Patameter Parse Error')
            #遍历alldesign找到window:1等参数
            for x in self.args.alldesign[::-1] :
                if 'window' in x : 
                    try : self.__design_opt['window'] = int(x.split(':')[-1])
                    except : self._base.warnlong('窗口参数解析错误/ Window Patameter Parse Error')
            self.__design_opt['merge'] = True if 'merge' in self.args.alldesign else False
            self.__design_opt['pdetail'] = True if 'pdetail' in self.args.alldesign else False
            self.__design_opt['pdetail2'] = True if 'pdetail2' in self.args.alldesign else False
            self.__design_opt['primer'] = True if 'primer' in self.args.alldesign else False
            self.__design_opt['primer2'] = True if 'primer2' in self.args.alldesign else False

        '''
        扩增子评估相关参数配置
        '''
        self.__evaluate_opt = {'minlen' : 150, 'maxlen' : 1500, 'hpcnt' : 10, 'merge': False, 'fullp' : True, 'save': False, 'tm' : 50.0, \
            'rmlow' : False, 'blast' : None, 'degene' : 12}
        #遍历evaluate找到hpcnt:10等参数
        for x in self.args.evaluate[::-1] :
            if 'hpcnt' in x : 
                try : self.__evaluate_opt['hpcnt'] = int(x.split(':')[-1])
                except : self._base.warnlong('引物特异参数解析错误/ Haplotype Number Patameter Parse Error')
        #遍历evaluate找到minlen:150等参数
        for x in self.args.evaluate[::-1] :
            if 'minlen' in x : 
                try : self.__evaluate_opt['minlen'] = int(x.split(':')[-1])
                except : self._base.warnlong('扩增子最小值参数解析错误/ Amplicon Minlen Patameter Parse Error')
        #遍历evaluate找到maxlen:1500等参数
        for x in self.args.evaluate[::-1] :
            if 'maxlen' in x : 
                try : self.__evaluate_opt['maxlen'] = int(x.split(':')[-1])
                except : self._base.warnlong('扩增子最大值参数解析错误/ Amplicon Maxlen Patameter Parse Error')
        #遍历evaluate找到blast:file1,filex等参数
        for x in self.args.evaluate[::-1] :
            if 'blast' in x : 
                try : self.__evaluate_opt['blast'] = x.split(':')[-1].split(',')
                except : self._base.warnlong('序列集合文件路径解析出错/ Blast File Parameter Parse Error')
        #遍历evaluate找到degene:12等参数
        for x in self.args.evaluate[::-1] :
            if 'degene' in x : 
                try : self.__evaluate_opt['degene'] = int(x.split(':')[-1])
                except : self._base.warnlong('简并参数解析出错/ Degenerate Parameter Parse Error')
        self.__evaluate_opt['merge'] = True if 'merge' in self.args.evaluate else False
        #self.__evaluate_opt['fullp'] = True if 'fullp' in self.args.evaluate else False
        self.__evaluate_opt['save'] = True if 'save' in self.args.evaluate else False
        self.__evaluate_opt['tm'] = self.__design_opt['tm']
        self.__evaluate_opt['rmlow'] = True if 'rmlow' in self.args.evaluate else False

        '''
        数据清洗相关参数配置
        '''
        #默认数据清洗相关参数，以及遍历alldesign找到相关参数
        self.__data_filt = {'len' : True, 'sameseq' : True, 'matrix' : False}
        if self.args.progress is not None and 'notlen' in self.args.progress : self.__data_filt.update({'len' : False})
        if self.args.progress is not None and 'notsameseq' in self.args.progress : self.__data_filt.update({'sameseq' : False})
        if self.args.progress is not None and 'matrix' in self.args.progress : self.__data_filt.update({'matrix' : True})

        #基础模块的获取，log等功能都在其中
        self._base = pbase

    '''
    创建人员: Nerium
    创建日期: 2022/08/31
    更改人员: Nerium
    更改日期: 2023/02/27
    '''
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
            if self.args.alldesign is not None and 'muscle' not in self.args.alldesign :
                try :
                    seqlen, seqcnt = len(next(iter(self._origindata.values()))), len(self._origindata.values())
                    for bp in range(seqlen) : self._origindata_shannon.append(calc_shannon_entropy([Counter([seq[bp] for seq in self._origindata.values()]).get(slg, 0)/seqcnt for slg in DEFAULT_DNA_SINGLE_LIST]))
                except : self._base.errorlog('序列并未对齐/ Sequences Not Align Yet')

    '''
    创建人员: Nerium
    创建日期: 2022/08/31
    更改人员: Nerium
    更改日期: 2022/09/12
    '''
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

    '''
    创建人员: Nerium
    创建日期: 2022/08/31
    更改人员: Nerium
    更改日期: 2022/09/22
    '''
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

    '''
    创建人员: Nerium
    创建日期: 2022/08/31
    更改人员: Nerium
    更改日期: 2023/03/01
    '''
    #调用primer3-py进行引物设计
    #可以先将序列去重再进行引物设计，但是保存原始信息会比较麻烦，快速开发先走流程
    def primer_design(self, pcds, area) :
        primer_dict, areacnt = {}, len(area)
        area_statistic = {}
        self._base.baselog('\n正在根据保守区间进行引物设计...', ends='')
        for numi, rang in enumerate(area, start=1) :
            self._base.baselog('\r正在根据保守区间进行引物设计... {}/{}'.format(numi, areacnt), ends='')

            data = self._comparedata if 'muscle' in self.args.alldesign else self._origindata
            #最后引物是dict中key:set()的形式，去重和保留原始样本名称信息
            primer_dict.setdefault(rang[0], (dict(), dict()))
            area_statistic.setdefault(rang[0], {'F': {}, 'R' : {}})

            #如果是一个-区间取消，则使用temp + (temp=None)break处理后添加到primer_dict中
            for spe, seq in data.items() :
                tmp_seq = seq[rang[0]-1:rang[1]]
                if len(tmp_seq.replace('-', '')) < 15 : continue
                #if tmp_seq.count('-') : continue

                #简并符号不算
                if len(tmp_seq) - sum(list_count(tmp_seq, ['A','T','C','G','-']).values()) : 
                    self._base.warnlog('\n{} 含有简并符号故跳过/ Skip Because Have Special BP'.format(spe)); continue

                seq_rm_gap_len = rang[1]-rang[0]-tmp_seq.count('-')+1
                seq_args = {
                    'SEQUENCE_ID': '{}-{}'.format(rang[0], rang[1]),
                    'SEQUENCE_TEMPLATE': seq.replace('-', ''),
                    'SEQUENCE_INCLUDED_REGION': [max(0, rang[0]-seq[:rang[0]-1].count('-')-1), seq_rm_gap_len],
                }
                opt_args = {
                    'PRIMER_MIN_SIZE':15,
                    'PRIMER_OPT_SIZE':((15+seq_rm_gap_len)//2) if seq_rm_gap_len<25 else 20,
                    'PRIMER_MAX_SIZE':min(25, seq_rm_gap_len),# if seq_rm_gap_len<35 else 35,
                    'PRIMER_PRODUCT_SIZE_RANGE':[150, 1500],
                    'PRIMER_MIN_TM': self.__design_opt['tm'],
                    'PRIMER_PICK_LEFT_PRIMER': 1,
                    'PRIMER_PICK_RIGHT_PRIMER': 1,
                    'PRIMER_NUM_RETURN': 1,
                }
                #pair_primer = pcds.callprimer(target=seq_args, opt=opt_args)
                try : pair_primer = pcds.callprimer(target=seq_args, opt=opt_args)
                except : self._base.errorlog([max(0, rang[0]-seq[:rang[0]-1].count('-')-1), seq_rm_gap_len])

                #将原始样本名称对应，保留原始信息
                for pri in pair_primer[0] :
                    if pri in primer_dict[rang[0]][0] : primer_dict[rang[0]][0][pri].add(spe)
                    else : primer_dict[rang[0]][0].setdefault(pri, {spe})
                for pri in pair_primer[1] :
                    if pri in primer_dict[rang[0]][1] : primer_dict[rang[0]][1][pri].add(spe)
                    else : primer_dict[rang[0]][1].setdefault(pri, {spe})

                if self.__design_opt['pdetail'] == False and self.__design_opt['primer2'] == False : continue
                if pair_primer[2] is not None : 
                    ppos = pos_translate(seq, pair_primer[2][0])
                    ass = '[{},{}]'.format(ppos, pos_translate(seq, pair_primer[2][0]+pair_primer[2][1]-1))
                    area_statistic[rang[0]]['F'].update({ass: area_statistic[rang[0]]['F'].get(ass, 0)+1})
                if pair_primer[3] is not None : 
                    ppos = pos_translate(seq, pair_primer[3][0])
                    ass = '[{},{}]'.format(pos_translate(seq, pair_primer[3][0]-pair_primer[3][1]+1), ppos)
                    area_statistic[rang[0]]['R'].update({ass: area_statistic[rang[0]]['R'].get(ass, 0)+1})

        self._base.successlog('\r\n已经根据保守区间完成引物设计')
        if self.__design_opt['pdetail'] : [self._base.baselog('{} : \nF{}\nR{}\n'.format(k,{kk:len(vv) for kk,vv in v[0].items()},{kk:len(vv) for kk,vv in v[1].items()})) if len(v[0])+len(v[1]) else self._base.baselog('{} : None'.format(k)) for k, v in primer_dict.items()]
        [self._base.debuglog(BASE_DEBUG_LEVEL3, '{} : {}'.format(k,v)) if len(v[0])+len(v[1]) else self._base.debuglog(BASE_DEBUG_LEVEL3, '{} : None'.format(k)) for k, v in primer_dict.items()]

        if self.__design_opt['pdetail'] : self._base.baselog(area_statistic)
        #if self.__evaluate_opt['save'] : write_json('{}_all_primer1.json'.format(self._base._time), {k:({kk:(len(vv), calc_tm_hairpin_homod(kk)) for kk, vv in v[0].items()}, {kk:(len(vv), calc_tm_hairpin_homod(kk)) for kk, vv in v[1].items()}) for k, v in primer_dict.items() if len(v[0])+len(v[1])})
        return primer_dict, area_statistic

    '''
    创建人员: Nerium
    创建日期: 2022/12/07
    更改人员: Nerium
    更改日期: 2023/02/15
    '''
    #调用primer3-py进行引物设计
    #通过第一次引物设计得到的区间进行二次引物设计
    def primer_design_2(self, pcds, area) :
        primer_dict, areacnt = {}, len(area)
        area_statistic = {}
        self._base.baselog('\n正在根据保守区间进行二次引物设计...', ends='')
        for numi, rang in enumerate(area, start=1) :
            self._base.baselog('\r正在根据保守区间进行二次引物设计... {}/{}'.format(numi, areacnt), ends='')

            data = self._comparedata if 'muscle' in self.args.alldesign else self._origindata
            #最后引物是dict中key:set()的形式，去重和保留原始样本名称信息
            primer_dict.setdefault(rang[0], (dict(), dict()))
            area_statistic.setdefault(rang[0], {'F': {}, 'R' : {}})

            #根据一次设计的结果获取到二次设计的区间
            tmp_rangef = sorted(self._area_statistic[rang[0]]['F'].items(), key=lambda z : z[1], reverse=True)[0][0].replace('[', '').replace(']', '').split(',')
            tmp_ranger = sorted(self._area_statistic[rang[0]]['R'].items(), key=lambda z : z[1], reverse=True)[0][0].replace('[', '').replace(']', '').split(',')
            tmp_rangef, tmp_ranger = [int(i) for i in tmp_rangef], [int(i) for i in tmp_ranger]
            self._base.debuglog(BASE_DEBUG_LEVEL1, (tmp_rangef, tmp_ranger))

            tmp_rang = rang
            #如果是一个-区间取消，则使用temp + (temp=None)break处理后添加到primer_dict中
            for spe, seq in data.items() :
                #二次设计F引物
                rang = tmp_rangef
                tmp_seq = seq[rang[0]-1:rang[1]]
                if tmp_seq.count('-') : continue

                #简并符号不算
                if len(tmp_seq) - sum(list_count(tmp_seq, ['A','T','C','G','-']).values()) : 
                    self._base.warnlog('\n{} 含有简并符号故跳过/ Skip Because Have Special BP'.format(spe)); continue

                seq_args = {
                    'SEQUENCE_ID': '{}-{}'.format(rang[0], rang[1]),
                    'SEQUENCE_TEMPLATE': seq.replace('-', ''),
                    'SEQUENCE_INCLUDED_REGION': [rang[0]-seq[:rang[0]].count('-')-1, rang[1]-rang[0]+1],
                }
                opt_args = {
                    'PRIMER_MIN_SIZE':15,
                    'PRIMER_OPT_SIZE':((15+rang[1]-rang[0]+1)//2) if rang[1]-rang[0]+1<35 else 25,
                    'PRIMER_MAX_SIZE':rang[1]-rang[0]+1 if rang[1]-rang[0]+1<35 else 35,
                    'PRIMER_PRODUCT_SIZE_RANGE':[rang[1]-rang[0]+1, 100],
                    'PRIMER_MIN_TM': self.__design_opt['tm'],
                    'PRIMER_PICK_LEFT_PRIMER': 1,
                    'PRIMER_PICK_RIGHT_PRIMER': 1,
                    'PRIMER_NUM_RETURN': 1,
                }
                p_left = pcds.callprimer_left(target=seq_args, opt=opt_args)

                #二次设计右引物
                rang = tmp_ranger
                tmp_seq = seq[rang[0]-1:rang[1]]
                if tmp_seq.count('-') : continue

                #简并符号不算
                if len(tmp_seq) - sum(list_count(tmp_seq, ['A','T','C','G','-']).values()) : 
                    self._base.warnlog('\n{} 含有简并符号故跳过/ Skip Because Have Special BP'.format(spe)); continue

                seq_args = {
                    'SEQUENCE_ID': '{}-{}'.format(rang[0], rang[1]),
                    'SEQUENCE_TEMPLATE': seq.replace('-', ''),
                    'SEQUENCE_INCLUDED_REGION': [rang[0]-seq[:rang[0]].count('-')-1, rang[1]-rang[0]+1],
                }
                opt_args = {
                    'PRIMER_MIN_SIZE':15,
                    'PRIMER_OPT_SIZE':((15+rang[1]-rang[0]+1)//2) if rang[1]-rang[0]+1<35 else 25,
                    'PRIMER_MAX_SIZE':rang[1]-rang[0]+1 if rang[1]-rang[0]+1<35 else 35,
                    'PRIMER_PRODUCT_SIZE_RANGE':[rang[1]-rang[0]+1, 100],
                    'PRIMER_MIN_TM': self.__design_opt['tm'],
                    'PRIMER_PICK_LEFT_PRIMER': 1,
                    'PRIMER_PICK_RIGHT_PRIMER': 1,
                    'PRIMER_NUM_RETURN': 1,
                }
                p_right = pcds.callprimer_right(target=seq_args, opt=opt_args)

                pair_primer = (p_left[0], p_right[0], p_left[1], p_right[1])
                rang = tmp_rang
                #将原始样本名称对应，保留原始信息
                for pri in pair_primer[0] :
                    if pri in primer_dict[rang[0]][0] : primer_dict[rang[0]][0][pri].add(spe)
                    else : primer_dict[rang[0]][0].setdefault(pri, {spe})
                for pri in pair_primer[1] :
                    if pri in primer_dict[rang[0]][1] : primer_dict[rang[0]][1][pri].add(spe)
                    else : primer_dict[rang[0]][1].setdefault(pri, {spe})

                if self.__design_opt['pdetail2'] == False : continue
                if pair_primer[2] is not None : 
                    ppos = pair_primer[2][0] + seq[:rang[0]].count('-') + 1
                    ass = '[{},{}]'.format(ppos, ppos+pair_primer[2][1]-1)
                    area_statistic[rang[0]]['F'].update({ass: area_statistic[rang[0]]['F'].get(ass, 0)+1})
                if pair_primer[3] is not None : 
                    ppos = pair_primer[3][0] + seq[:rang[0]].count('-') + 1
                    ass = '[{},{}]'.format(ppos-pair_primer[3][1]+1, ppos)
                    area_statistic[rang[0]]['R'].update({ass: area_statistic[rang[0]]['R'].get(ass, 0)+1})

        self._base.successlog('\r\n已经根据保守区间完成引物设计')
        if self.__design_opt['pdetail2'] : [self._base.baselog('{} : \nF{}\nR{}\n'.format(k,{kk:len(vv) for kk,vv in v[0].items()},{kk:len(vv) for kk,vv in v[1].items()})) if len(v[0])+len(v[1]) else self._base.baselog('{} : None'.format(k)) for k, v in primer_dict.items()]
        [self._base.debuglog(BASE_DEBUG_LEVEL3, '{} : {}'.format(k,v)) if len(v[0])+len(v[1]) else self._base.debuglog(BASE_DEBUG_LEVEL3, '{} : None'.format(k)) for k, v in primer_dict.items()]

        if self.__design_opt['pdetail2'] : self._base.baselog(area_statistic)
        return primer_dict, area_statistic

    '''
    创建人员: Nerium
    创建日期: 2022/12/08
    更改人员: Nerium
    更改日期: 2023/02/27
    '''
    #调用primer3-py进行引物设计
    #通过第一次引物设计得到的区间直接从序列中提取
    def primer_design_extract(self, pcds, area) :
        primer_dict, areacnt = {}, len(area)
        self._base.baselog('\n正在根据保守区间进行二次引物提取...', ends='')
        for numi, rang in enumerate(area, start=1) :
            self._base.baselog('\r正在根据保守区间进行二次引物提取... {}/{}'.format(numi, areacnt), ends='')

            data = self._comparedata if 'muscle' in self.args.alldesign else self._origindata
            #最后引物是dict中key:set()的形式，去重和保留原始样本名称信息
            primer_dict.setdefault(rang[0], (dict(), dict()))

            #根据一次设计的结果获取到二次设计的区间
            if len(self._area_statistic[rang[0]]['F']) == 0 or len(self._area_statistic[rang[0]]['R']) == 0 :
                self._base.warnlog('区间{0}无引物/ Area {0} No Primer'.format(rang)); continue
            tmp_rangef = sorted(self._area_statistic[rang[0]]['F'].items(), key=lambda z : z[1], reverse=True)[0][0].replace('[', '').replace(']', '').split(',')
            tmp_ranger = sorted(self._area_statistic[rang[0]]['R'].items(), key=lambda z : z[1], reverse=True)[0][0].replace('[', '').replace(']', '').split(',')
            tmp_rangef, tmp_ranger = [int(i) for i in tmp_rangef], [int(i) for i in tmp_ranger]
            self._base.debuglog(BASE_DEBUG_LEVEL1, (tmp_rangef, tmp_ranger))

            tmp_rang = rang
            len_stat = {'F': {}, 'R': {}}
            #如果是一个-区间取消，则使用temp + (temp=None)break处理后添加到primer_dict中
            for spe, seq in data.items() :
                #直接提取F引物
                rang = tmp_rangef
                tmp_seq = seq[rang[0]-1:rang[1]].replace('-', '')
                len_stat['F'].update({len(tmp_seq) : len_stat['F'].get(len(tmp_seq), 0)+1})

                #简并符号不算
                if len(tmp_seq) - sum(list_count(tmp_seq, ['A','T','C','G','-']).values()) : 
                    self._base.warnlog('\n{} 含有简并符号故跳过/ Skip Because Have Special BP'.format(spe)); continue
                p_left = ([tmp_seq, ], None)

                #直接提取R引物
                rang = tmp_ranger
                tmp_seq = seq[rang[0]-1:rang[1]].replace('-', '')
                len_stat['R'].update({len(tmp_seq) : len_stat['F'].get(len(tmp_seq), 0)+1})

                #简并符号不算
                if len(tmp_seq) - sum(list_count(tmp_seq, ['A','T','C','G','-']).values()) : 
                    self._base.warnlog('\n{} 含有简并符号故跳过/ Skip Because Have Special BP'.format(spe)); continue
                p_right = ([''.join([DEFAULT_DNA_REFLECT_DICT.get(bp, '-') for bp in tmp_seq])[::-1], ], None)

                pair_primer = (p_left[0], p_right[0], p_left[1], p_right[1])
                rang = tmp_rang
                #将原始样本名称对应，保留原始信息
                for pri in pair_primer[0] :
                    if pri in primer_dict[rang[0]][0] : primer_dict[rang[0]][0][pri].add(spe)
                    else : primer_dict[rang[0]][0].setdefault(pri, {spe})
                for pri in pair_primer[1] :
                    if pri in primer_dict[rang[0]][1] : primer_dict[rang[0]][1][pri].add(spe)
                    else : primer_dict[rang[0]][1].setdefault(pri, {spe})

            fplen = sorted(len_stat['F'].items(), key=lambda z:z[1], reverse=True)[0][0]
            rplen = sorted(len_stat['R'].items(), key=lambda z:z[1], reverse=True)[0][0]
            for pri in list(primer_dict[rang[0]][0].keys()) : 
                if len(pri) != fplen or len(pri) < 15 : primer_dict[rang[0]][0].pop(pri)
            for pri in list(primer_dict[rang[0]][1].keys()) : 
                if len(pri) != rplen or len(pri) < 15 : primer_dict[rang[0]][1].pop(pri)

        self._base.successlog('\r\n已经根据保守区间完成二次引物提取')
        if self.__design_opt['pdetail2'] : [self._base.baselog('{} : \nF{}\nR{}\n'.format(k,{kk:(len(vv), calc_tm_hairpin_homod(kk)) for kk,vv in v[0].items()},{kk:(len(vv), calc_tm_hairpin_homod(kk)) for kk,vv in v[1].items()})) if len(v[0])+len(v[1]) else self._base.baselog('{} : None'.format(k)) for k, v in primer_dict.items()]
        [self._base.debuglog(BASE_DEBUG_LEVEL3, '{} : {}'.format(k,v)) if len(v[0])+len(v[1]) else self._base.debuglog(BASE_DEBUG_LEVEL3, '{} : None'.format(k)) for k, v in primer_dict.items()]

        #if self.__evaluate_opt['save'] : write_json('{}_all_primer2.json'.format(self._base._time), {k:({kk:(len(vv), calc_tm_hairpin_homod(kk)) for kk, vv in v[0].items()}, {kk:(len(vv), calc_tm_hairpin_homod(kk)) for kk, vv in v[1].items()}) for k, v in primer_dict.items() if len(v[0])+len(v[1])})
        return primer_dict, self._area_statistic

    '''
    创建人员: Nerium
    创建日期: 2022/08/31
    更改人员: Nerium
    更改日期: 2023/03/03
    '''
    #主流程函数
    def maintrunk(self) :

        if self.args.file is not None :
            #先保存原始数据
            self.getorigin()

        #如果开启则进行数据清洗，最后将xx.fasta清洗存为xx.filt.fasta
        if self.args.progress is not None :
            pcdp = piecedataprogress(self._base, self._origindata, self.__data_filt)

            pcdp.check_info()

            self._origindata = pcdp.filt_data()

            self.args.file = self.args.file.replace('.fasta', '.filt.fasta')
            pcdp.write_in_file(self.args.file)

        #如果开启，则调用muscle进行多序列比对；探查保守区间；循环论证最佳引物
        if self.args.alldesign is not None :
            pcds = piecedesign(self._base, self._todo_path, self.args.file, self.__design_opt)
            if 'muscle' in self.args.alldesign :
                pcds.callmuscle()

                #保存对比后的数据
                self.aftercmp(pcds)

            self._base.debuglog(BASE_DEBUG_LEVEL2, self._comparedata_shannon if 'muscle' in self.args.alldesign else self._origindata_shannon)
            #挖掘出所有符合条件的保守区间
            conser = pcds.detect_conser_area_shannon(self._comparedata_shannon if 'muscle' in self.args.alldesign else self._origindata_shannon,
                                                    self._comparedata if 'muscle' in self.args.alldesign else self._origindata)
            self._base.baselog('保守区间列表 / List Of Conservative Area is : \n{0}'.format(conser))

            #挖掘出所有符合条件的非保守区间
            nonconser = pcds.detect_non_conser_area(self._comparedata_shannon if 'muscle' in self.args.alldesign else self._origindata_shannon, conser)
            self._base.baselog('非保守区间列表 / List Of Non Conservative Area is : \n{0}'.format(nonconser))

            #if self.__evaluate_opt['save'] : write_json('{}_conser_nonconser.json'.format(self._base._time), {'conser' : conser, 'nonconser' : nonconser})

            #非保守区间多样性
            nonconser_sort = self.rank_by_diverse(pcds, nonconser, '非保守区间', 'rank1' in self.args.alldesign)

            #保守区间多样性
            if 'rank2' in self.args.alldesign : conser_sort = self.rank_by_diverse(pcds, conser, '保守区间', True)

            #所有区间的多样性排名，保守和非保守分别排名是必须的，但是全排序不是必须的
            if 'rankall' in self.args.alldesign : self.rank_by_diverse(pcds, conser+nonconser, '所有区间', True)

            #根据haplotype进行分析和后续的引物设计
            self._alltype = {}
            if 'haplo' in self.args.alldesign : self._base.baselog('\n保守区间的haplotype情况如下：')
            for rang in conser :
                alltype = pcds.detect_haplotype(self._comparedata if 'muscle' in self.args.alldesign else self._origindata, rang[0], rang[1])
                if 'haplo' in self.args.alldesign : self._base.baselog('Area {}; \t {}'.format(rang, len(alltype)))
                self._alltype.setdefault(str(rang), alltype)

            if 'primer' in self.args.alldesign or 'primer2' in self.args.alldesign :
                self._primer_dict, self._area_statistic = self.primer_design(pcds, conser)

            if 'primer2' in self.args.alldesign :
                self._primer_dict, self._area_statistic = self.primer_design_extract(pcds, conser)

        #如果开启，则进行最终区域选择和评估等
        if self.args.evaluate is not None :
            try : pcel = pieceevaluate(self._base, nonconser_sort, conser, self._primer_dict, self._comparedata if 'muscle' in self.args.alldesign else self._origindata, self.__evaluate_opt)
            except : self._base.errorlog('\n未进行引物设计/Cannot Find Designed Primer')

            if self.__evaluate_opt['fullp'] : area_res = pcel.full_permutation()
            #根据条件从区间中过滤出合适的保守区间和非保守区间(conser1)nonconser(conser2)
            else : area_res = pcel.filter_area()

            if len(pcel._posmem) == 0 : self._base.errorlog('没有合适的扩增子区间/ No Right Amplicon')

            #对扩增子信息进行整合
            amplicon_info = pcel.recommend_area_primer(dct=pcdp._same_cnt if 'pcdp' in locals().keys() is not None else None)

            #评估扩增子分辨力
            reso = pcel.evaluate_resolution()

            #评估扩增子覆盖度
            cover_rate = pcel.evaluate_cover_rate()

            if self.__evaluate_opt['save'] : write_json('{}_recommand_area_primer.json'.format(self._base._time), amplicon_info)

            if self.__evaluate_opt['blast'] : write_json('{}_final_recommand_area_primer.json'.format(self._base._time), pcel.blast_db_search('{}_recommand_area_primer.json'.format(self._base._time)))