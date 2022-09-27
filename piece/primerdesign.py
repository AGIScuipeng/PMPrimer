'''
创建人员: Nerium
创建日期: 2022/08/31
更改人员: Nerium
更改日期: 2022/09/27
'''

from piece.piecedefine import *
from piece.piecebase import calc_conserve_continue, calc_conserve_termina_shannon, calc_shannon_entropy, generate_shannon_bynum

from primer3 import bindings
from collections import Counter
import subprocess, platform

#多序列比对、保守区间遍历、PCR设计等
class piecedesign() :
    def __init__(self, pbase, todo_path, filepath, model='simple') -> None:
        self._todo_path = todo_path
        self.model = model
        self.__platform = platform.system()

        self.filepath = filepath
        self.tmpfile_path = self.filepath.replace('.fasta', '.mc.fasta')

        self._base = pbase

    #subprocess调用muscle进行多序列比对
    def callmuscle(self) :
        self._base.baselog(BASE_DEBUG_LEVEL1, 'MUSCLE 多序列对比中...', ends='')
        cm = subprocess.Popen('{}/{} --align {} --output {}'.format\
            (self._todo_path, PLATFORM_TODO[self.__platform], self.filepath, self.tmpfile_path),\
            shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        cm.wait()
        #self._base.debuglog(BASE_DEBUG_LEVEL1, cm.communicate()[1].decode())
        self._base.baselog(BASE_DEBUG_LEVEL1, '\rMUSCLE 多序列比对完成' if cm.returncode == 0 else 'MUSCLE 多序列比对异常')

    #调用primer3-py设计引物
    def callprimer(self, target, opt=None) :
        primer3_result = bindings.designPrimers(target, opt)
        #print(primer3_result)

        rescnt = primer3_result['PRIMER_LEFT_NUM_RETURNED']
        return ([primer3_result['PRIMER_LEFT_{}_SEQUENCE'.format(i)] for i in range(rescnt)], [primer3_result['PRIMER_RIGHT_{}_SEQUENCE'.format(i)] for i in range(rescnt)])

    #挖掘所有保守区域
    def detect_conser_area(self, seqdict, threshold=0.95, minlen=15) :
        self._base.baselog(BASE_DEBUG_LEVEL1, '探测比对后序列的所有保守区域...', ends='')
        posl, posr, seqlen, posmem, mem = 1, 1, len(next(iter(seqdict.values()))), [], [0.0]*2

        while posr < seqlen :
            while calc_conserve_continue(seqdict, posl, posr, mem, threshold) and posr < seqlen : posr += 1
            if posr - posl + (1 if posr-posl == 0 else 0) >= minlen and ((mem[1]/len(seqdict.values()))/(1 if posr-posl == 0 else posr-posl)) >= threshold : posmem.append([posl, posr-1 if posr < seqlen else posr])
            posr += 1; posl = posr; mem = [0.0]*2
        self._base.baselog(BASE_DEBUG_LEVEL1, '\r比对后序列的所有保守区域探测完毕')
        return posmem

    #通过香农熵挖掘所有保守区域
    def detect_conser_area_shannon(self, seqdict, threshold=generate_shannon_bynum(0.95), minlen=15) :
        self._base.baselog(BASE_DEBUG_LEVEL1, '探测比对后序列的所有保守区域...', ends='')
        posl, posr, seqlen, posmem, mem = 1, 1, len(seqdict), [], [1.0]*2

        while posr < seqlen :
            while calc_conserve_termina_shannon(seqdict, posl, posr, mem, threshold) and posr < seqlen : posr += 1
            if posr - posl + (1 if posr-posl == 0 else 0) >= minlen and mem[1] <= threshold : posmem.append([posl, posr-1 if posr < seqlen and posl != posr else posr])
            posr += 1; posl = posr; mem = [1.0]*2
        self._base.baselog(BASE_DEBUG_LEVEL1, '\r比对后序列的所有保守区域探测完毕')
        return posmem

    #挖掘所有非保守区域
    def detect_non_conser_area(self, seqdict, posmem, minlen=1) :
        posmem_len, seqlen = len(posmem), len(seqdict)
        ret = [[rang[1]+1, posmem[idx+1][0]-1] for idx, rang in enumerate(posmem) if idx != posmem_len-1 and posmem[idx+1][0] - rang[1] > minlen]
        if posmem[-1][-1] != seqlen : ret.append([posmem[-1][-1]+1, seqlen])
        return ret

    #计算非保守区域的多样性，结果越接近1，多样性越高
    def calc_area_diverse(self, seqdict, posl, posr) :
        seqcnt, allshannon = len(seqdict), []
        self._base.debuglog(BASE_DEBUG_LEVEL3, 'A\tT\tC\tG\t-\t{0},{1}'.format(posl, posr))

        self._base.debuglog(BASE_DEBUG_LEVEL3, seqdict[posl:posr])
        allshannon.extend(seqdict[posl-1:posr])
        return round(sum(allshannon)/len(allshannon), 8)

    #通过将列表集合化，可以得到有多少个不同的序列，从而把相同的排除掉
    def detect_hypertype(self, seqdict, posl, posr) :
        return set([seq[posl:posr+1] for seq in seqdict.values()])