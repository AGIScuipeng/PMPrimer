'''
创建人员: Nerium
创建日期: 2022/8/31
更改人员: Nerium
更改日期: 2022/09/02
'''

from piece.piecedefine import *
from piece.piecebase import evaluate

from primer3 import bindings
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
        cm = subprocess.Popen('{}/{} --align {} --output {}'.format\
            (self._todo_path, PLATFORM_TODO[self.__platform], self.filepath, self.tmpfile_path),\
            shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        cm.wait()
        #self._base.debuglog(BASE_DEBUG_LEVEL1, cm.communicate()[1].decode())
        self._base.baselog(BASE_DEBUG_LEVEL1, 'MUSCLE 多序列比对完成' if cm.returncode == 0 else 'MUSCLE 多序列比对异常')

    #调用primer3-py设计引物
    def callprimer(self, target, opt) :
        return bindings.designPrimers(target, opt)

    #挖掘所有保守区域
    def detectarea(self, seqdict) :
        self._base.baselog(BASE_DEBUG_LEVEL1, '探测比对后序列的所有保守区域...', ends='')
        posl, posr, seqlen, posmem = 1, 1, len(next(iter(seqdict.values()))), []
        while posr < seqlen :
            while evaluate(seqdict, posl, posr) > 0.6 and posr < seqlen : posr += 1
            posmem.append((posl, posr)); posr += 1; posl = posr
        self._base.baselog(BASE_DEBUG_LEVEL1, '\r比对后序列的所有保守区域探测完毕')
        return posmem