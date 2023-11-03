'''
创建人员: Nerium
创建日期: 2022/08/31
更改人员: Nerium
更改日期: 2023/10/14
'''

from .piecedefine import *
from .piecebase import calc_conserve_continue, calc_conserve_termina_shannon, generate_shannon_bynum

from primer3 import bindings
import subprocess, platform

'''
创建人员: Nerium
创建日期: 2022/08/31
更改人员: Nerium
更改日期: 2023/03/06
'''
#多序列比对、保守区间遍历、PCR设计等
class piecedesign() :
    def __init__(self, pbase, todo_path, filepath, design_opt) -> None:
        import re
        self._todo_path = todo_path
        self.__design_opt = design_opt
        self.__platform = platform.system()

        self.filepath = filepath
        self.tmpfile_path = re.sub(r"\.fasta|\.afa|\.fa", ".mc.fasta", self.filepath) if self.filepath is not None and '.fasta' in self.filepath else pbase.errorlog('文件路径为空，或命名错误')
        self._base = pbase

    '''
    创建人员: Nerium
    创建日期: 2022/08/31
    更改人员: Nerium
    更改日期: 2023/03/01
    '''
    #subprocess调用muscle进行多序列比对
    def callmuscle(self) :
        self._base.baselog('MUSCLE 多序列对比中...')
        #!!!因为MSUCLE输出很多，所以Popen()中的stdout=PIPE和wait()不能结合使用，因为操作系统的管道大小有上限
        #很久之前就解决过的问题，但是因为对MUSCLE不了解，再次出现
        cm = subprocess.Popen('{}/{} --super5 {} --output {}'.format\
            (self._todo_path, PLATFORM_TODO[self.__platform], self.filepath, self.tmpfile_path),\
            shell=True, stderr=subprocess.PIPE)
        cm.communicate()

        self._base.successlog('MUSCLE 多序列比对完成' if cm.returncode == 0 else self._base.errorlog('MUSCLE 多序列比对异常'))

    '''
    创建人员: Nerium
    创建日期: 2022/08/31
    更改人员: Nerium
    更改日期: 2023/03/06
    '''
    #调用primer3-py设计引物
    def callprimer(self, target, opt=None, tops=99) :
        #分别设计F引物，R引物
        opt['PRIMER_PICK_RIGHT_PRIMER']=0
        f_result = bindings.design_primers(target, opt)
        if f_result['PRIMER_LEFT_NUM_RETURNED'] == 0 : self._base.debuglog(BASE_DEBUG_LEVEL3, f_result)

        opt['PRIMER_PICK_RIGHT_PRIMER']=1; opt['PRIMER_PICK_LEFT_PRIMER']=0
        r_result = bindings.design_primers(target, opt)
        if r_result['PRIMER_RIGHT_NUM_RETURNED'] == 0 : self._base.debuglog(BASE_DEBUG_LEVEL3, r_result)

        #TOPS参数暂且无用，后续可能会有用，故暂时保留
        fcnt, rcnt = min(f_result['PRIMER_LEFT_NUM_RETURNED'], tops), min(r_result['PRIMER_RIGHT_NUM_RETURNED'], tops)
        return ([f_result['PRIMER_LEFT_{}_SEQUENCE'.format(i)] for i in range(fcnt)], [r_result['PRIMER_RIGHT_{}_SEQUENCE'.format(i)] for i in range(rcnt)], f_result['PRIMER_LEFT_0'] if fcnt else None, r_result['PRIMER_RIGHT_0'] if rcnt else None)

    '''
    创建人员: Nerium
    创建日期: 2022/12/07
    更改人员: Nerium
    更改日期: 2022/12/07
    '''
    #调用primer3-py设计引物
    def callprimer_left(self, target, opt=None, tops=99) :
        #分别设计F引物，R引物
        opt['PRIMER_PICK_RIGHT_PRIMER']=0
        f_result = bindings.designPrimers(target, opt)
        if f_result['PRIMER_LEFT_NUM_RETURNED'] == 0 : self._base.debuglog(BASE_DEBUG_LEVEL3, f_result)

        #TOPS参数暂且无用，后续可能会有用，故暂时保留
        fcnt = min(f_result['PRIMER_LEFT_NUM_RETURNED'], tops)
        return ([f_result['PRIMER_LEFT_{}_SEQUENCE'.format(i)] for i in range(fcnt)], f_result['PRIMER_LEFT_0'] if fcnt else None)

    '''
    创建人员: Nerium
    创建日期: 2022/12/07
    更改人员: Nerium
    更改日期: 2022/12/07
    '''
    #调用primer3-py设计引物
    def callprimer_right(self, target, opt=None, tops=99) :
        #分别设计F引物，R引物
        opt['PRIMER_PICK_LEFT_PRIMER']=0
        r_result = bindings.designPrimers(target, opt)
        if r_result['PRIMER_RIGHT_NUM_RETURNED'] == 0 : self._base.debuglog(BASE_DEBUG_LEVEL3, r_result)

        #TOPS参数暂且无用，后续可能会有用，故暂时保留
        rcnt = min(r_result['PRIMER_RIGHT_NUM_RETURNED'], tops)
        return ([r_result['PRIMER_RIGHT_{}_SEQUENCE'.format(i)] for i in range(rcnt)], r_result['PRIMER_RIGHT_0'] if rcnt else None)

    '''
    创建人员: Nerium
    创建日期: 2022/08/31
    更改人员: Nerium
    更改日期: 2022/09/30
    '''
    #挖掘所有保守区域
    def detect_conser_area(self, seqdict, threshold=0.95, minlen=15) :
        self._base.baselog('探测比对后序列的所有保守区域...', ends='')
        posl, posr, seqlen, posmem, mem = 1, 1, len(next(iter(seqdict.values()))), [], [0.0]*2

        while posr < seqlen :
            while calc_conserve_continue(seqdict, posl, posr, mem, threshold) and posr < seqlen : posr += 1
            if posr - posl + (1 if posr-posl == 0 else 0) >= minlen and ((mem[1]/len(seqdict.values()))/(1 if posr-posl == 0 else posr-posl)) >= threshold : posmem.append([posl, posr-1 if posr < seqlen else posr])
            posr += 1; posl = posr; mem = [0.0]*2
        self._base.baselog('\r比对后序列的所有保守区域探测完毕')
        return posmem

    '''
    创建人员: Nerium
    创建日期: 2022/11/22
    更改人员: Nerium
    更改日期: 2023/03/01
    '''
    #根据间隔的香农熵合并相近的保守区间
    def merge_conser_area(self, shannons, posmem, seqdict) :
        if len(posmem) == 0 : return []

        final_pos = [posmem[0]]
        for idx in range(1, len(posmem)) :
            #和前一个合并
            non_conserl, non_conserr = posmem[idx-1][1]+1, posmem[idx][0]-1
            part_shannons = shannons[non_conserl-1:non_conserr]

            self._base.debuglog(BASE_DEBUG_LEVEL1, posmem[idx], ends=' | ')
            self._base.debuglog(BASE_DEBUG_LEVEL1, sum(part_shannons), ends=' | ')

            #如果间隔香农熵和<2，则看final合并
            if sum(part_shannons) < 2 :
                self._base.debuglog(BASE_DEBUG_LEVEL1, final_pos, ends= ' | ')
                if final_pos[-1][1] == posmem[idx-1][1] : 
                    self._base.debuglog(BASE_DEBUG_LEVEL1, len(self.detect_haplotype(seqdict, final_pos[-1][0], posmem[idx][1])), ends=' | ')
                    #还要看合并后的haplotype，现行写死10
                    if len(self.detect_haplotype(seqdict, final_pos[-1][0], posmem[idx][1])) < 10 : final_pos[-1][1] = posmem[idx][1]
                    else : final_pos.append(posmem[idx])
                else : final_pos.append(posmem[idx])
            else : final_pos.append(posmem[idx])
            self._base.debuglog(BASE_DEBUG_LEVEL1, '')

        return final_pos

    '''
    创建人员: Nerium
    创建日期: 2022/08/31
    更改人员: Nerium
    更改日期: 2023/03/01
    '''
    #通过香农熵挖掘所有保守区域
    def detect_conser_area_shannon(self, shannons, seqdict, threshold=None, minlen=None, pwindow=None, merge=None) :
        if threshold is None : threshold = self.__design_opt['threshold']
        if minlen is None : minlen = self.__design_opt['minlen']
        if merge is None : merge = self.__design_opt['merge']
        if pwindow is None : pwindow = self.__design_opt['window']

        self._base.baselog('探测比对后序列的所有保守区域...')
        posl, posr, seqlen, posmem, mem, window = 1, 1, len(shannons), [], [1.0]*2, pwindow

        #如果保守区域个数<2 且 窗口<4 则扩大窗口继续尝试
        while len(posmem) < 1 and window < 4 :
            while posr < seqlen :
                while calc_conserve_termina_shannon(shannons, posl, posr, mem, threshold, window=window) and posr < seqlen : posr += 1
                if posr - posl + (1 if posr-posl == 0 else 0) >= minlen and mem[1] <= threshold : posmem.append([posl, posr-1 if posr < seqlen and posl != posr else posr])
                posr += 1; posl = posr; mem = [1.0]*2
            window += 1; posl, posr, seqlen, mem = 1, 1, len(shannons), [1.0]*2

        if merge : posmem = self.merge_conser_area(shannons, posmem, seqdict)

        seqcnt = len(seqdict)
        #正序删除list会导致元素迁移从而删除错误，故倒序遍历
        for rang in posmem[::-1] :
            tcnt, rmflag = 0, False
            for ibp in range(max(0, rang[0]-1), rang[1]) :
                tmp_rate = [seq[ibp] for seq in seqdict.values()].count('-')/seqcnt
                if tmp_rate >= 0.9 : tcnt += 1
                if tmp_rate >= self.__design_opt['gaps'] : posmem.remove(rang); rmflag = True; self._base.debuglog(BASE_DEBUG_LEVEL1, '{0} {2}空白符占比 {1:.2f}% 区间删除/{0} Deleted For {2} Gaps Rate {1:.2f}%'.format(rang, tmp_rate*100, ibp+1)); break
            if rang[1] - rang[0] + 1 - tcnt < 15 and rmflag is not True : posmem.remove(rang); self._base.debuglog(BASE_DEBUG_LEVEL1, '{0} 有效长度不足15/ Region {0} Truly Length Less Than 15 bp'.format(rang))


        #没有足够的保守区间，程序直接退出
        if len(posmem) < 2 : self._base.errorlog('\n香农熵中断和延续法无法探测到足够的保守区域/ Shannon Terminate Or Continue Cannot Detect Enough Conserved Region')
        self._base.successlog('\r比对后序列的所有保守区域探测完毕')
        #[[57,77], [165,184], [339,360], [453,474], [698,717], [835,852]]
        #[[148,168], [569,589]]
        return posmem

    '''
    创建人员: Nerium
    创建日期: 2022/08/31
    更改人员: Nerium
    更改日期: 2023/02/17
    '''
    #挖掘所有非保守区域
    def detect_non_conser_area(self, shannons, posmem, minlen=1) :
        posmem_len, seqlen = len(posmem), len(shannons)
        ret = [[rang[1]+1, posmem[idx+1][0]-1] for idx, rang in enumerate(posmem) if idx != posmem_len-1 and posmem[idx+1][0] - rang[1] > minlen]
        #if posmem[-1][-1] != seqlen : ret.append([posmem[-1][-1]+1, seqlen])
        return ret

    '''
    创建人员: Nerium
    创建日期: 2022/08/31
    更改人员: Nerium
    更改日期: 2023/06/06
    '''
    #计算区域的多样性，结果越接近1，多样性越高
    def calc_area_diverse(self, shannons, posl, posr) :
        self._base.debuglog(BASE_DEBUG_LEVEL3, shannons[posl-1:posr])
        allshannon = shannons[posl-1:posr]
        return round(sum(allshannon)/len(allshannon), 8)

    '''
    创建人员: Nerium
    创建日期: 2022/08/31
    更改人员: Nerium
    更改日期: 2022/10/17
    '''
    #通过将列表集合化，可以得到有多少个不同的序列，从而把相同的排除掉
    def detect_haplotype(self, seqdict, posl, posr) :
        return set([seq[posl-1:posr] for seq in seqdict.values()])