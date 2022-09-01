# 扩增子流程化

## 0x0 预想流程
+ 第一步，使用primer3对输入的所有序列进行引物设计，并通过muscle多序列比对、清洗合并等方式完成最终的引物设计；
+ 第二步，通过评估标准对上一步的结果进行评估，从而为第三步的使用提供更好的数据、参数支撑；
+ 第三步，参考第二步的数据、参数，将物种鉴定、SNP挖掘等步骤进行关联。

## 0x1 第一步
1. primer3有python的调用接口，为primer3-py，但是都是2年以前的代码，有参考和使用价值，但是不知道是否改用primer5，primer5没有python接口。
    > primer3-py，地址为：[github of primer3-py](https://github.com/libnano/primer3-py)

2. muscle每次release都会发布源代码和Linux、Windows、MacOS下的可执行文件。可以直接将可执行文件打包进工具，根据执行平台直接运行。
    > muscle地址为：[muscle5](https://github.com/rcedgar/muscle)，release地址为：[release of muscle5](https://github.com/rcedgar/muscle/releases)。

    > 多序列比对后，muscle会把序列对齐，会在中间增加gap，但是物种鉴定等应用需要的是原始序列、坐标信息，故需要将对比前后的坐标有转换系统（原始序列可能也有gap）；

3. 清洗、合并等操作需要自行完成
    1. 保守区间的查找
        > 通过窗口和步长在muscle多序列对比之后，找到保守率最高的区间，未达到阈值则继续，达到阈值则保留
    2. 引物的设计15~25bp
        > 通过primer3设计合适的引物，并持续通过TM值等进行评估
    
    需要将3.1和3.2进行循环论证，找到最合适的保守区间和保守区间的引物设计

## 0x2 第二步
> 评估系统根据扩增子异常识别、覆盖度、分辨能力、差异分析等维度来进行

## 0x3 第三步
+ 物种鉴定已经使用Perl完成demo
+ SNP挖掘已经使用Perl完成demo
> 后续不知道是否改由Python重写，或者直接由Python调用Perl