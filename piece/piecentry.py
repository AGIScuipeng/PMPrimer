'''
创建人员: Nerium
创建日期: 2022/08/31
更改人员: Nerium
更改日期: 2023/02/27
'''

from .__auxiliary__ import __version__
from .piecedefine import *
from .piecebase import piecebase
from .piecemain import piecemain

from argparse import ArgumentParser
from argparse import RawTextHelpFormatter
import sys

'''
创建人员: Nerium
创建日期: 2022/08/31
更改人员: Nerium
更改日期: 2023/02/27
'''
#封装包程序的入口
def entry() :
    #命令行程序
    paramparse = ArgumentParser(description='primer piece', epilog = __version__,formatter_class=RawTextHelpFormatter)
    paramparse.add_argument('--file', '-f', help='初始序列文件')
    paramparse.add_argument('--progress', '-p', nargs='+', 
                            help='''基础数据清洗，参数有：\
                            \nnotlen 不根据长度分布清洗\
                            \nnotsameseq 不根据相同序列合并\
                            \nmatrix 根据序列相似度画出聚类热图''')
    paramparse.add_argument('--alldesign', '-a', nargs='+',
                            help='''多序列引物设计，参数有：\
                            \nthreshold 保守区间判定的阈值，默认为0.95，使用： threshold:0.95\
                            \nminlen 保守区间连续长度最小值，默认为15bp，使用： minlen:15\
                            \ngaps 空白符占比最高比例，默认为0.1，使用： gaps:0.1\
                            \ntm 熔解温度默认值，默认为50，使用： tm:50.0\
                            \nmuscle 使用muscle进行多序列比对并另存文件\
                            \nmerge 保守区间合并\
                            \nprimer 第一次引物设计\
                            \npdetail 第一次引物设计详细信息输出\
                            \nprimer2 第二次引物提取\
                            \npdetail2 第二次引物提取详细信息输出''')
    paramparse.add_argument('--debuglevel', '-d', help='调试等级', type=int, default=0)
    paramparse.add_argument('--evaluate', '-e', nargs='+', default=['default'], 
                            help='''引物标准评估，参数有：\
                            \nhpcnt HaploType最大值，默认为10，使用： hpcnt:10\
                            \nminlen 扩增子区间最小值，默认为150，使用： minlen:150\
                            \nmaxlen 扩增子区间最大值，默认为1500，使用： maxlen:1500\
                            \nblast 使用文件建库和查询，多个fasta文件使用,来分隔，使用： blast:../taxo1.fasta,homo.fasta\
                            \nmerge 多个扩增子合并\
                            \nfullp 使用全排列生成待选扩增子\
                            \nrmlow 去除引物提取熔解温度低于设计模块配置值的引物\
                            \nsave 保存所有中间文件''')
    #paramparse.add_argument('--nextuse', '-n', help='后续SNP挖掘和物种鉴定')

    #解析命令行
    args = paramparse.parse_args()
    #命令行为空时，直接打印使用信息并退出
    if len(sys.argv) == 1 : paramparse.print_help(); sys.exit()

    #使用流程主类生成对象
    pc = piecemain(args, piecebase(BASE_DEBUG_LEVEL0 if args.debuglevel is None else args.debuglevel))
    #开始进入主流程
    pc.maintrunk()