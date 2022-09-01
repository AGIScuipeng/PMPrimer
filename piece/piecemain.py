'''
创建人员: Nerium
创建日期: 2022/8/31
更改人员: Nerium
更改日期: 2022/8/31
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
        self._comparedata = {}

        #基础模块的获取，log等功能都在其中
        self._base = pbase

    #原始数据的保存
    def getorigin(self) :
        if self.args.file is not None :
            with open(self.args.file, 'r') as ff :
                slt_id = ''
                for line in ff :
                    if '>' in line : slt_id = line[1:].strip()
                    else : self._origindata.update({slt_id: line.strip()})
        #print(self._origindata)

    #对比数据的保存
    def aftercmp(self, pcds) :
        if pcds.tmpfile_path is not None :
            with open(pcds.tmpfile_path, 'r') as ff :
                slt_id = ''
                for line in ff :
                    if '>' in line : slt_id = line[1:].strip()
                    else : self._comparedata.update({slt_id: line.strip()})
        #print(self._comparedata)

    #主流程函数
    def maintrunk(self) :
        #先保存原始数据
        self.getorigin()

        #调用muscle进行多序列比对
        #if self.args.alldesign is not None :
        pcds = piecedesign(self._base, self._todo_path, self.args.file)
        pcds.callmuscle()

        #保存对比后的数据
        self.aftercmp(pcds)

        #print(originpos(self, 'sth2', 13))