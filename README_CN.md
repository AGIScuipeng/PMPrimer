# 扩增子流程化
English ReadMe, [Click Here](README_EN.md)

## 0x0 模块解析
> 本软件需要依赖第三方：primer3-py、biopython、seaborn等；

> 使用BLAST需要本地安装有ncbi-blast-2.13.0+；

> 本软件官方只支持Linux，对于Windows默认不支持，但不代表Windows不能运行。
+ 输入文件
+ 基础数据清洗
+ 多序列引物设计
+ 扩增子设计和评估
+ 调试模式

## 0x1 命令行参数
1. 输入文件
   > --file -f

    使用如下：-f seqs.fasta 或 --file seqs.fasta

2. 基础数据清洗
   > --progress -p

    参数共有:
    + notlen 不根据长度分布清洗
    + notsameseq 不根据相同序列合并
    + matrix 根据序列差异度画出聚类热图

    >如果序列长度都相差不大，但是分布不太平均，则使用notlen防止去除占比少的序列。

    将需要的参数直接在 -p 后添加即可，示例如下：-p notlen notsameseq matrix，无所谓顺序。

3. 多序列引物设计
   >--alldesign -a

    参数共有：
    + threshold 保守区间判定的阈值，默认为0.95，使用： threshold:0.95
    + minlen 保守区间连续长度最小值，默认为15bp，使用： minlen:15
    + gaps 空白符占比最高比例，默认为0.1，使用： gaps:0.1
    + tm 熔解温度默认值，默认为50，使用： tm:50.0
    + muscle 使用muscle进行多序列比对并另存文件
    + merge 保守区间合并
    + primer2 第二次引物提取
    + pdetail2 第二次引物提取详细信息输出

    将需要的参数直接在 -a 后添加即可，示例如下：-a threshold:0.95 minlen:16 merge primer2，无所谓顺序。

4. 扩增子设计和评估
   > --evaluate -e

    参数共有：
    + hpcnt HaploType最大值，默认为10，使用： hpcnt:10
    + minlen 扩增子区间最小值，默认为150，使用： minlen:150
    + maxlen 扩增子区间最大值，默认为1500，使用： maxlen:1500
    + blast 使用文件建库和查询，多个fasta文件使用,来分隔，使用： blast:../taxo1.fasta,homo.fasta
    + merge 多个扩增子合并
    + rmlow 去除引物提取熔解温度低于设计模块配置值的引物
    + save 保存所有中间文件

    将需要的参数直接在 -e 后添加即可，示例如下：-e hpcnt:50 maxlen:1200 save，无所谓顺序。

5. 调试模式
   > --debuglevel -d

    参数为数字：
    + 1 代表着基础调试信息
    + 2 代表着模块调试信息
    + 4 代表着输出复杂的调试信息

    如需调试信息输出，一般情况下使用 -d 1即可，本模块采取的是二进制计算，如3就代表着1+2，7就是1+2+4全部调试信息输出等。

## 0x2 小建议
1. 如果输入文件是多序列比对后的文件，可以使用参数如下：

    pmprimer -f seqs.fasta -a merge primer2 -e hpcnt:70 save

2. 如果输入文件需要数据清洗，那么就需要先进行数据清洗：

    pmprimer -f seqs_need_filt.fasta -p default
    > 此命令会根据长度分布、相同序列保留亚种的方式对序列进行清洗

3. 如果输入文件需要对齐后再进行流程：

    pmprimer -f seqs_need_align.fasta -a muscle merge primer2 -e hpcnt:70 save
    > 此命令会将对齐后的文件保存为seqs_need_align.mc.fasta