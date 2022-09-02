'''
创建人员: Nerium
创建日期: 2022/8/31
更改人员: Nerium
更改日期: 2022/09/02
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
    def detectarea(self, seqdict) :
        self._base.baselog(BASE_DEBUG_LEVEL1, '探测比对后序列的所有保守区域...', ends='')
        posl, posr, seqlen, posmem, mem = 1, 1, len(next(iter(seqdict.values()))), [], [0.0]
        while posr < seqlen :
            while calc_conserve_continue(seqdict, posl, posr, mem) >= 0.8 and posr < seqlen : posr += 1
            #后续可以指定保守区域的最小长度再判断是否添加
            #if posr - posl + 1 >= 15 : posmem.append((posl, posr))
            posmem.append((posl, posr))
            posr += 1; posl = posr; mem[0] = 0.0
        self._base.baselog(BASE_DEBUG_LEVEL1, '\r比对后序列的所有保守区域探测完毕')
        return posmem

    #计算非保守区域的标准差，结果越接近0，多样性越高
    def calc_non_conserve_area_std(self, seqdict, posl, posr) :
        seqcnt, allstd = len(seqdict.values()), []
        for i in range(posl, posl+1 if posr-posl == 0 else posr+1) :
            self._base.debuglog(BASE_DEBUG_LEVEL3, [Counter([seq[i-1] for seq in seqdict.values()]).get(slg, 0)/seqcnt for slg in DEFAULT_DNA_SINGLE_LIST])
            allstd.append(calc_std([Counter([seq[i-1] for seq in seqdict.values()]).get(slg, 0)/seqcnt for slg in DEFAULT_DNA_SINGLE_LIST]))
        return calc_std(allstd)