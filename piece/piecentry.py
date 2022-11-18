'''
创建人员: Nerium
创建日期: 2022/08/31
更改人员: Nerium
更改日期: 2022/11/18
'''

from .__auxiliary__ import __version__
from .piecedefine import *
from .piecebase import piecebase
from .piecemain import piecemain

from argparse import ArgumentParser
import sys

'''
创建人员: Nerium
创建日期: 2022/08/31
更改人员: Nerium
更改日期: 2022/11/18
'''
#封装包程序的入口
def entry() :
    #命令行程序
    paramparse = ArgumentParser(description='primer piece', epilog = __version__)
    paramparse.add_argument('--file', '-f', help='初始序列文件')
    paramparse.add_argument('--progress', '-p', nargs='+', help='基础数据清洗')
    paramparse.add_argument('--alldesign', '-a', nargs='+', help='多序列引物设计', default=['default'])
    paramparse.add_argument('--debuglevel', '-d', help='调试等级', type=int, default=0)
    paramparse.add_argument('--evaluate', '-e', nargs='+', help='引物标准评估', default=['default'])
    paramparse.add_argument('--nextuse', '-n', help='后续SNP挖掘和物种鉴定')

    #解析命令行
    args = paramparse.parse_args()
    #命令行为空时，直接打印使用信息并退出
    if len(sys.argv) == 1 : paramparse.print_help(); sys.exit()

    #使用流程主类生成对象
    pc = piecemain(args, piecebase(BASE_DEBUG_LEVEL0 if args.debuglevel is None else args.debuglevel))
    #开始进入主流程
    pc.maintrunk()