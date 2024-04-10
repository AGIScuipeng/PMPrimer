'''
创建人员: Nerium
创建日期: 2022/10/25
更改人员: Nerium
更改日期: 2023/02/27
'''

from .piecedefine import *
from .piecebase import split_all_from_str, mismatch_btw_2seq
from multiprocessing import Pool, cpu_count

'''
创建人员: Nerium
创建日期: 2022/10/25
更改人员: Nerium
更改日期: 2022/10/25
'''
#获取序列集的平均长度 data = {id : seq}
def seq_everage_len(data) :
    return sum([len(x) for x in data.values()]) // len(data)

'''
创建人员: Nerium
创建日期: 2022/10/25
更改人员: Nerium
更改日期: 2022/11/16
'''
#获取 data = {id : seq}; ret = [id1, id2]
def seq_after_filt_len(data) :
    datacnt = {}
    for id, seq in data.items() : datacnt.update({round(len(seq), -1) : datacnt.get(round(len(seq), -1), set()).union({id, })})
    datacnt = {x:y for x, y in sorted(datacnt.items(), key=lambda z : len(z[1]), reverse=True)}

    key_list = list(datacnt.keys())

    rate_content = 0
    all_content = len(data)
    for idx, k in enumerate(key_list) :
        rate_content += len(datacnt[k])
        if rate_content / all_content > 0.9 : break

    return [id for k in key_list[:idx+1] for id in datacnt[k]]

'''
创建人员: Nerium
创建日期: 2022/10/25
更改人员: Nerium
更改日期: 2023/02/27
'''
#序列集去重data = {id: seq}; seqset = {seq : {id1, id2, ...}}（id1和id2保证了seq相同的情况下亚种层次的不同）
#same_cnt用以后续简并引物作为原始序列排序的记录
def seq_set(data, flag=False) : 
    seqset, same_cnt = dict(), dict()
    for id, seq in data.items() : 
        idsplit = split_all_from_str(id)
        if idsplit is None : continue

        spe_id_dict = {' '.join(split_all_from_str(x)[1:]) : split_all_from_str(x)[0] for x in seqset.get(seq, set())}
        t_spe = ' '.join(idsplit[1:])
        #如果当前三级结构已经在seqset里面了，则使用seqset中记录的那个三级结构的唯一id作为key来保存同样的三级和序列一共重复了多少次
        if t_spe in spe_id_dict :
            same_cnt.update({spe_id_dict[t_spe] : same_cnt.get(spe_id_dict[t_spe], 0)+1})
            continue
        seqset.update({seq : seqset.get(seq, set()).union({id, })})

    if flag : return seqset, same_cnt
    return seqset

'''
创建人员: Nerium
创建日期: 2022/10/25
更改人员: Nerium
更改日期: 2022/10/27
'''
#序列集去重后id调整 seqset = {seq1 : {id1, id2, ...}}; return {id1 :seq1, id2:seq1}
def seq_keep1id_after_set(seqset) :
    temp = {}
    for k, v in seqset.items() :
        for id in v :
            temp.setdefault(id, k)

    return temp

'''
创建人员: Nerium
创建日期: 2022/10/28
更改人员: Nerium
更改日期: 2022/11/10
'''
#去除sp.未分类的序列 seqset = {id1: seq1, id2: seqx}
def seq_remove_unclassfied(seqset) :
    temp = {}
    for id, seq in seqset.items() :
        idsplit = split_all_from_str(id)
        if idsplit is None : continue
        if idsplit[2] == 'sp.' or idsplit[1] == 'Uncultured' or 'UNVERIFIED' in idsplit[1] or idsplit[1].isalpha() is False : continue

        temp.update({id: seq})

    return temp

'''
创建人员: Nerium
创建日期: 2022/10/25
更改人员: Nerium
更改日期: 2023/02/27
'''
#数据清洗功能类
class piecedataprogress() :
    def __init__(self, pbase, data, filt_opt) -> None:
        self._data = data
        self.__filt_opt = filt_opt

        self._base = pbase

    '''
    创建人员: Nerium
    创建日期: 2022/10/25
    更改人员: Nerium
    更改日期: 2022/10/25
    '''
    #数据集基本信息
    def check_info(self) :
        self._base.baselog('序列集共有序列 {0} 条 / Number Of Sequece Data is {0}'.format(len(self._data)))
        self._base.baselog('序列集平均长度 {0} 位 / Everage Length Of Sequece Data is {0} bp'.format(seq_everage_len(self._data)))
        self._base.baselog('去重后序列集共 {0} 条 / Number After Duplicate Remove is {0}'.format(len(seq_set(self._data))))

    '''
    创建人员: Nerium
    创建日期: 2022/11/10
    更改人员: Nerium
    更改日期: 2022/11/25
    '''
    #数据集基本信息
    def check_gnu_spe_sub_info(self) :
        ids_list = list(self._data.keys())
        subset = set()
        speset = set()
        gnuset = set()
        for id in ids_list : 
            id_temp = split_all_from_str(id)
            if id_temp[2] == 'sp.' or id_temp[1] == 'Uncultured' or 'UNVERIFIED' in id_temp[1] or id_temp[1].isalpha() is False : continue

            subset.add(' '.join(id_temp[1:]))
            speset.add(' '.join(id_temp[1:3]))
            gnuset.add(id_temp[1])
        self._base.baselog('清洗后亚种共 {} 个'.format(len(subset)))
        self._base.baselog('清洗后物种共 {} 个'.format(len(speset)))
        self._base.baselog('清洗后属共 {} 个'.format(len(gnuset)))
        self._base.debuglog(BASE_DEBUG_LEVEL1, gnuset)
        self._base.debuglog(BASE_DEBUG_LEVEL3, speset)
        self._base.debuglog(BASE_DEBUG_LEVEL3, subset)

    '''
    创建人员: Nerium
    创建日期: 2022/11/18
    更改人员: Nerium
    更改日期: 2022/11/21
    '''
    #物种内计算差异矩阵
    def calc_mismath_matrix(self) :
        tmp_spes = {}
        if len(self._data) == 1 : self._base.errorlog('序列只有1条无法计算差异图/Cannot Calculate Matrix Because Only 1 Sequence.')
        #分成物种存储相应id, seq
        for id, seq in self._data.items() :
            #后续可能要根据配置的物种等名称进行统计，故保留注释代码
            spe = 'all'#' '.join(split_all_from_str(id)[1:3])
            if spe in tmp_spes : tmp_spes[spe].update({id: seq})
            else : tmp_spes.setdefault(spe, {id: seq})

        #进程池进行差异矩阵计算
        jobs = []
        pool = Pool(cpu_count()//2)
        for spe, v in tmp_spes.items() :
            for id1, seq1 in v.items() :
                for id2, seq2 in v.items() :
                    if id1 == id2 : jobs.append((spe, id1, id2, 0)); continue
                    jobs.append((spe, id1, id2, pool.apply_async(func=mismatch_btw_2seq, args=(seq1, seq2))))
        pool.close()
        pool.join()

        #存储为差异矩阵
        mismatch_matrix = {}
        for spe in tmp_spes.keys() : mismatch_matrix.update({spe : dict()})
        for job in jobs :
            spe, id1, id2, score = job[0], split_all_from_str(job[1])[0], split_all_from_str(job[2])[0], 0 if job[3] == 0 else job[3].get()
            if score == -1 : self._base.errorlog('序列并未对齐/ Sequences Not Alignment')

            if id1 in mismatch_matrix[spe] : mismatch_matrix[spe][id1].update({id2 : score})
            else : mismatch_matrix[spe].setdefault(id1, {id2: score})

        #绘制聚类热力图

        import seaborn as sns
        import pandas as pd
        import matplotlib.pyplot as plt

        tmp_spe = list(mismatch_matrix.keys())[0]
        tmp_dict = mismatch_matrix[tmp_spe]
        data = pd.DataFrame.from_dict(tmp_dict)

        cmap = sns.color_palette("coolwarm", 200)
        sns_res = sns.clustermap(data, cmap=cmap)
        plt.savefig('{}.png'.format('tmp'))

        #输出聚类后的横坐标
        tmp_key_list = list(mismatch_matrix[tmp_spe].keys())
        self._base.debuglog(BASE_DEBUG_LEVEL2, [tmp_key_list[idx] for idx in sns_res.dendrogram_row.reordered_ind])

    '''
    创建人员: Nerium
    创建日期: 2022/10/25
    更改人员: Nerium
    更改日期: 2023/02/27
    '''
    #先长度清洗，后去重，得到的是最后的结果
    def filt_data(self) : 
        ids_list = seq_after_filt_len(self._data) if self.__filt_opt.get('len', True) else list(self._data.keys())
        self._base.baselog('长度清洗后共 {0} 条 / Number After Keep Majority Length is {0}'.format(len(ids_list)))

        if self.__filt_opt.get('sameseq', True) :
            seqset, self._same_cnt = seq_set({id : self._data[id] for id in ids_list}, flag=True)
            resdata = seq_keep1id_after_set(seqset)
        else : resdata = {id : self._data[id] for id in ids_list}
        self._base.baselog('同序列保留亚种后共 {0} 条 / Number After Keep Different Subspecies When Same Sequece is {0}'.format(len(resdata)))

        self._data = seq_remove_unclassfied(resdata)
        self._base.baselog('未分类序列清洗后共 {0} 条 / Number After Remove Unclassfied is {0}'.format(len(self._data)))

        if self.__filt_opt.get('matrix', False) : self.calc_mismath_matrix()

        self.check_gnu_spe_sub_info()
        return self._data

    '''
    创建人员: Nerium
    创建日期: 2022/10/25
    更改人员: Nerium
    更改日期: 2023/02/27
    '''
    #清洗完成的数据，写入文件中
    def write_in_file(self, file_path) : 
        self._base.successlog('清洗完成后序列共 {0} 条 / Number After Progress is {0}'.format(len(self._data)))
        with open(file_path, 'w') as tf :
            for id, seq in self._data.items() :
                tf.write('>' + id)
                tf.write('\n' + seq + '\n')