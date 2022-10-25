'''
创建人员: Nerium
创建日期: 2022/08/31
更改人员: Nerium
更改日期: 2022/10/25
'''

from piece.piecedefine import *
from piece.piecebase import calc_conserve_continue, calc_conserve_termina_shannon, generate_shannon_bynum

from primer3 import bindings
import subprocess, platform

'''
创建人员: Nerium
创建日期: 2022/08/31
更改人员: Nerium
更改日期: 2022/10/25
'''
#多序列比对、保守区间遍历、PCR设计等
class piecedesign() :
    def __init__(self, pbase, todo_path, filepath, model='simple') -> None:
        self._todo_path = todo_path
        self.model = model
        self.__platform = platform.system()

        self.filepath = filepath
        self.tmpfile_path = self.filepath.replace('.fasta', '.mc.fasta') if self.filepath is not None and '.fasta' in self.filepath else pbase.errorlog('文件路径为空，或命名错误')

        self._base = pbase

    '''
    创建人员: Nerium
    创建日期: 2022/08/31
    更改人员: Nerium
    更改日期: 2022/10/25
    '''
    #subprocess调用muscle进行多序列比对
    def callmuscle(self) :
        self._base.baselog('MUSCLE 多序列对比中...', ends='')
        #!!!因为MSUCLE输出很多，所以Popen()中的stdout=PIPE和wait()不能结合使用，因为操作系统的管道大小有上限
        #很久之前就解决过的问题，但是因为对MUSCLE不了解，再次出现
        cm = subprocess.Popen('{}/{} --align {} --output {}'.format\
            (self._todo_path, PLATFORM_TODO[self.__platform], self.filepath, self.tmpfile_path),\
            shell=True, stderr=subprocess.PIPE)
        cm.communicate()

        self._base.baselog('\rMUSCLE 多序列比对完成' if cm.returncode == 0 else self._base.errorlog('MUSCLE 多序列比对异常'))

    '''
    创建人员: Nerium
    创建日期: 2022/08/31
    更改人员: Nerium
    更改日期: 2022/10/14
    '''
    #调用primer3-py设计引物
    def callprimer(self, target, opt=None, tops=99) :
        #分别设计F引物，R引物
        opt['PRIMER_PICK_RIGHT_PRIMER']=0
        f_result = bindings.designPrimers(target, opt)
        if f_result['PRIMER_LEFT_NUM_RETURNED'] == 0 : self._base.debuglog(BASE_DEBUG_LEVEL3, f_result)

        opt['PRIMER_PICK_RIGHT_PRIMER']=1; opt['PRIMER_PICK_LEFT_PRIMER']=0
        r_result = bindings.designPrimers(target, opt)
        if r_result['PRIMER_RIGHT_NUM_RETURNED'] == 0 : self._base.debuglog(BASE_DEBUG_LEVEL3, r_result)

        #TOPS参数暂且无用，后续可能会有用，故暂时保留
        fcnt, rcnt = min(f_result['PRIMER_LEFT_NUM_RETURNED'], tops), min(r_result['PRIMER_RIGHT_NUM_RETURNED'], tops)
        return ([f_result['PRIMER_LEFT_{}_SEQUENCE'.format(i)] for i in range(fcnt)], [r_result['PRIMER_RIGHT_{}_SEQUENCE'.format(i)] for i in range(rcnt)])

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
    创建日期: 2022/08/31
    更改人员: Nerium
    更改日期: 2022/10/10
    '''
    #通过香农熵挖掘所有保守区域
    def detect_conser_area_shannon(self, shannons, seqdict, threshold=generate_shannon_bynum(0.95), minlen=15) :
        self._base.baselog('探测比对后序列的所有保守区域...')
        posl, posr, seqlen, posmem, mem, window = 1, 1, len(shannons), [], [1.0]*2, 1

        #如果保守区域个数<2 且 窗口<4 则扩大窗口继续尝试
        while len(posmem) < 2 and window < 4 :
            while posr < seqlen :
                while calc_conserve_termina_shannon(shannons, posl, posr, mem, threshold) and posr < seqlen : posr += 1
                if posr - posl + (1 if posr-posl == 0 else 0) >= minlen and mem[1] <= threshold : posmem.append([posl, posr-1 if posr < seqlen and posl != posr else posr])
                posr += 1; posl = posr; mem = [1.0]*2
            window += 1; posl, posr, seqlen, mem = 1, 1, len(shannons), [1.0]*2

        seqcnt = len(seqdict)
        #正序删除list会导致元素迁移从而删除错误，故倒序遍历
        for rang in posmem[::-1] :
            for ibp in range(max(0, rang[0]-1), rang[1]) :
                if [seq[ibp] for seq in seqdict.values()].count('-')/seqcnt >= 0.1 : posmem.remove(rang); self._base.debuglog(BASE_DEBUG_LEVEL1, '{0} 空白符过多区间删除/{0} Deleted For Gaps'.format(rang)); break

        #没有足够的保守区间，程序直接退出
        if len(posmem) < 2 : self._base.errorlog('\n香农熵中断和延续法无法探测到足够的保守区域/ Shannon Terminate Or Continue Cannot Detect Enough Conservate Area')
        self._base.successlog('\r比对后序列的所有保守区域探测完毕')
        return posmem

    '''
    创建人员: Nerium
    创建日期: 2022/08/31
    更改人员: Nerium
    更改日期: 2022/09/10
    '''
    #挖掘所有非保守区域
    def detect_non_conser_area(self, shannons, posmem, minlen=1) :
        posmem_len, seqlen = len(posmem), len(shannons)
        ret = [[rang[1]+1, posmem[idx+1][0]-1] for idx, rang in enumerate(posmem) if idx != posmem_len-1 and posmem[idx+1][0] - rang[1] > minlen]
        if posmem[-1][-1] != seqlen : ret.append([posmem[-1][-1]+1, seqlen])
        return ret

    '''
    创建人员: Nerium
    创建日期: 2022/08/31
    更改人员: Nerium
    更改日期: 2022/10/14
    '''
    #计算区域的多样性，结果越接近1，多样性越高
    def calc_area_diverse(self, shannons, posl, posr) :
        self._base.debuglog(BASE_DEBUG_LEVEL3, shannons[posl-1:posr])
        allshannon = shannons[posl-1:posr]
        return round(sum(allshannon)/len([x for x in allshannon if x]), 8)

    '''
    创建人员: Nerium
    创建日期: 2022/08/31
    更改人员: Nerium
    更改日期: 2022/10/17
    '''
    #通过将列表集合化，可以得到有多少个不同的序列，从而把相同的排除掉
    def detect_haplotype(self, seqdict, posl, posr) :
        return set([seq[posl-1:posr] for seq in seqdict.values()])