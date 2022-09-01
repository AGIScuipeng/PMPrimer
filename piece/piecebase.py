'''
创建人员: Nerium
创建日期: 2022/8/31
更改人员: Nerium
更改日期: 2022/8/31
'''

from piece.piecedefine import *

def originpos(mainc, id, finalpos) :
    x, y, originstr, finalstr = 0, 0, mainc._origindata[id], mainc._comparedata[id]
    while y < finalpos :
        if originstr[x] == finalstr[y] : x += 1; y += 1
        else : y += 1
    return x

#基础模块，log等功能都在其中
class piecebase() :
    def __init__(self, level=BASE_DEBUG_LEVEL0) -> None:
        self._level = level

    def baselog(self, setlevel, msg) :
        if setlevel & self._level : print(msg)

    def debuglog(self, setlevel, msg) :
        if setlevel & self._level : print('\033[0;33;40m{}\033[0m'.format(msg))