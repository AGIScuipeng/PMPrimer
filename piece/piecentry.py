'''
创建人员: Nerium
创建日期: 2022/08/31
更改人员: Nerium
更改日期: 2023/09/06
'''

from .__auxiliary__ import __version__
from .piecedefine import *
from .piecebase import piecebase
from .piecemain import piecemain

from argparse import ArgumentParser
from argparse import RawTextHelpFormatter
import sys, locale

'''
创建人员: Nerium
创建日期: 2022/08/31
更改人员: Nerium
更改日期: 2023/09/06
'''
#封装包程序的入口
def entry() :
    #命令行程序
    paramparse = ArgumentParser(description='PMPrimer is a Python-based tool for automated design and evaluation of multiplex PCR primer pairs for diverse templates.',\
         epilog = __version__,formatter_class=RawTextHelpFormatter)
    if 'zh_CN' == locale.getdefaultlocale()[0] : 
        paramparse.add_argument('--file', '-f', help='FASTA文件')
        paramparse.add_argument('--progress', '-p', nargs='+', 
                                help='''基础数据清洗，参数有：\
                                \nnotlen \t\t不根据长度分布清洗\
                                \nnotsameseq \t不根据最后一级分类合并相同序列\
                                \nmatrix \t\t根据序列距离画出聚类热图''')
        paramparse.add_argument('--alldesign', '-a', nargs='+',
                                help='''多序列引物设计，参数有：\
                                \nmuscle \t\t使用MUSCLE进行多序列比对并另存文件\
                                \nthreshold \t保守区间判定的阈值，默认为0.95，使用 threshold:0.95\
                                \nminlen \t\t保守区间连续长度最小值，默认为15bp，使用 minlen:15\
                                \ngaps \t\t空白符占比最高比例，默认为0.1，使用 gaps:0.1\
                                \nmerge \t\t保守区间合并\
                                \nrank1 \t\t展示非保守区间的多样性分数值排名\
                                \nrank2 \t\t展示保守区间的多样性分数值排名\
                                \nhaplo \t\t展示保守区间的haploType个数\
                                \ntm \t\t熔解温度最小值，默认为50，使用 tm:50.0\
                                \nprimer2 \t引物设计''')
        paramparse.add_argument('--evaluate', '-e', nargs='+', default=['default'], 
                                help='''扩增子设计和评估，参数有：\
                                \nhpcnt \t\tHaploType个数最大值，默认为10，使用： hpcnt:10\
                                \nminlen \t\t扩增子区间最小值，默认为150，使用： minlen:150\
                                \nmaxlen \t\t扩增子区间最大值，默认为1500，使用： maxlen:1500\
                                \nblast \t\t使用文件建库和查询，多个fasta文件使用,来分隔，使用： blast:../taxo1.fasta,homo.fasta\
                                \nrmlow \t\t去除引物提取熔解温度低于设计模块配置值的引物\
                                \nsave \t\t保存结果文件''')
        paramparse.add_argument('--debuglevel', '-d', type=int, default=0,
                                help='调试等级 \t1 基础调试信息，2 模块调试信息，4 复杂调试信息')
        #paramparse.add_argument('--nextuse', '-n', help='后续SNP挖掘和物种鉴定')
    else :
        paramparse.add_argument('--file', '-f', help='Input FASTA file. For example, seqs.fasta')
        paramparse.add_argument('--progress', '-p', nargs='+', 
                                help='''Basic data processing module, has the following parameters: \
                                \nnotlen \t\tDo not processing sequences according to distribution of length\
                                \nnotsameseq \tDo not remove redundancy sequences in terminal taxa\
                                \nmatrix \t\tDraw heatmap according to distances of multiple templates''')
        paramparse.add_argument('--alldesign', '-a', nargs='+',
                                help='''Primer design module, has the following parameters: \
                                \nmuscle \t\tUse MUSCLE to align multiple sequences, the alignment sequences will be saved as *.mc.fasta\
                                \nthreshold \tThreshold for identify conservative region, default is 0.95 frequency for major allele, use like threshold:0.95\
                                \nminlen \t\tMinimum length for identifying conservative region, default is 15, use like minlen:15\
                                \ngaps \t\tMaximum rate of gaps, default is 0.1, use like gaps:0.1\
                                \nmerge \t\tUse merge module to merge conservative regions\
                                \nrank1 \t\tRank according to diversity score of non conservative regions and display them\
                                \nrank2 \t\tRank according to diversity score of conservative regions and display them\
                                \nhaplo \t\tDisplay haploType sequence number of conservative regions\
                                \ntm \t\tMinimum melting temperature, default is 50, use like tm:50.0\
                                \nprimer2 \tPrimer design''')
        paramparse.add_argument('--evaluate', '-e', nargs='+', default=['default'], 
                                help='''Amplicon selection and evaluation module, has the following parameters: \
                                \nhpcnt \t\tMaximum count of haploType sequences, default is 10, use like hpcnt:10\
                                \nminlen \t\tMinimum length of amplicon, default is 150bp, use like minlen:150\
                                \nmaxlen \t\tMaximum length of amplicon, default is 1500bp, use like maxlen:1500\
                                \nblast \t\tUse blast to evaluate amplicon specificity by aligning amplicon primer pairs to target fasta files, use ',' to split multiple target files, use like blast:seqs.fasta,hosts.fasta\
                                \nrmlow \t\tRemove primers with melting temperature lower than parameter "tm" in primer design module\
                                \nsave \t\tSave final results''')
        paramparse.add_argument('--debuglevel', '-d', type=int, default=0,
                                help='''Debug level \t1 Basic debug info, 2 Module debug info, 4 Complex debug info''')

    #解析命令行
    args = paramparse.parse_args()
    #命令行为空时，直接打印使用信息并退出
    if len(sys.argv) == 1 : paramparse.print_help(); sys.exit()

    #使用流程主类生成对象
    pc = piecemain(args, piecebase(BASE_DEBUG_LEVEL0 if args.debuglevel is None else args.debuglevel))
    #开始进入主流程
    pc.maintrunk()