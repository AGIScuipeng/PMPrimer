'''
创建人员: Nerium
创建日期: 2022/08/31
更改人员: Nerium
更改日期: 2022/08/31
'''

from piece.__auxiliary__ import __version__
from piece.piecedefine import *
from piece.piecebase import piecebase
from piece.piecemain import piecemain

from argparse import ArgumentParser

#封装包程序的入口
def entry() :
    #命令行程序
    paramparse = ArgumentParser(description='primer piece', epilog = __version__)
    paramparse.add_argument('--file', '-f', help='初始序列文件')
    paramparse.add_argument('--alldesign', '-a', help='多序列引物设计')
    paramparse.add_argument('--debuglevel', '-d', help='调试等级')
    paramparse.add_argument('--evaluate', '-e', help='引物标准评估')
    paramparse.add_argument('--nextuse', '-n', help='后续SNP挖掘和物种鉴定')

    #解析命令行
    args = paramparse.parse_args()

    #生成流程主类
    pc = piecemain(args, piecebase(BASE_DEBUG_LEVEL1 if args.debuglevel is None else args.debuglevel))
    #开始进入主流程
    pc.maintrunk()