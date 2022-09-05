'''
创建人员: Nerium
创建日期: 2022/08/31
更改人员: Nerium
更改日期: 2022/09/05
'''

from piece.piecedefine import *
from piece.piecebase import calc_conserve_continue, calc_std

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
    def callprimer(self, target, opt) :
        return bindings.designPrimers(target, opt)

    #挖掘所有保守区域
    def detectarea(self, seqdict, threshold=0.95, minlen=15) :
        self._base.baselog(BASE_DEBUG_LEVEL1, '探测比对后序列的所有保守区域...', ends='')
        posl, posr, seqlen, posmem, mem = 1, 1, len(next(iter(seqdict.values()))), [], [0.0]

        while posr < seqlen :
            while calc_conserve_continue(seqdict, posl, posr, mem) >= threshold and posr < seqlen : posr += 1
            #后续可以指定保守区域的最小长度再判断是否添加
            if posr - posl + 1 >= minlen : posmem.append([posl, posr])
            posr += 1; posl = posr; mem[0] = 0.0
        self._base.baselog(BASE_DEBUG_LEVEL1, '\r比对后序列的所有保守区域探测完毕')
        return posmem

    #挖掘所有非保守区域
    def detect_non_conser_area(self, posmem) :
        posmem_len = len(posmem)
        return [[rang[1]+1, posmem[idx+1][0]-1] for idx, rang in enumerate(posmem) if idx != posmem_len-1 and posmem[idx+1][0] - rang[1] > 1]

    #计算非保守区域的标准差，结果越接近0，多样性越高
    def calc_non_conserve_area_std(self, seqdict, posl, posr) :
        seqcnt, allstd = len(seqdict.values()), []
        self._base.debuglog(BASE_DEBUG_LEVEL3, 'A\tT\tC\tG\t-')

        for i in range(posl, posl+1 if posr-posl == 0 else posr+1) :
            self._base.debuglog(BASE_DEBUG_LEVEL3, [Counter([seq[i-1] for seq in seqdict.values()]).get(slg, 0)/seqcnt for slg in DEFAULT_DNA_SINGLE_LIST])
            allstd.append(calc_std([Counter([seq[i-1] for seq in seqdict.values()]).get(slg, 0)/seqcnt for slg in DEFAULT_DNA_SINGLE_LIST]))
        return round(calc_std(allstd), 8)