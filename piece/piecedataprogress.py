'''
创建人员: Nerium
创建日期: 2022/10/25
更改人员: Nerium
更改日期: 2022/10/31
'''

from piece.piecedefine import *
from piece.piecebase import split_all_from_str

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
更改日期: 2022/10/25
'''
#获取 data = {id : seq}; ret = [id1, id2]
def seq_after_filt_len(data) :
    datacnt = {}
    for id, seq in data.items() : datacnt.update({round(len(seq), -1) : datacnt.get(round(len(seq), -1), set()).union({id, })})
    datacnt = {x:y for x, y in sorted(datacnt.items(), key=lambda x : x[1])}

    key_list = list(datacnt.keys())
    rate_content = 0
    all_content = len(data)
    for idx, k in enumerate(key_list) :
        rate_content += len(datacnt[k])
        if rate_content / all_content > 0.7 : break

    return [id for k in key_list[:idx+1] for id in datacnt[k]]

'''
创建人员: Nerium
创建日期: 2022/10/25
更改人员: Nerium
更改日期: 2022/10/27
'''
#序列集去重data = {id: seq}; seqset = {seq : {id1, id2, ...}}（id1和id2保证了seq相同的情况下亚种层次的不同）
def seq_set(data) : 
    seqset = dict()
    for id, seq in data.items() : 
        idsplit = split_all_from_str(id)
        if ' '.join(idsplit[1:]) in [' '.join(split_all_from_str(x)[1:]) for x in seqset.get(seq, set())] : continue
        seqset.update({seq : seqset.get(seq, set()).union({id, })})

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
更改日期: 2022/10/31
'''
#去除sp.未分类的序列 seqset = {id1: seq1, id2: seqx}
def seq_remove_unclassfied(seqset) :
    temp = {}
    for id, seq in seqset.items() :
        idsplit = split_all_from_str(id)
        if idsplit[2] == 'sp.' or idsplit[1] == 'Uncultured' : continue

        temp.update({id: seq})

    return temp

'''
创建人员: Nerium
创建日期: 2022/10/25
更改人员: Nerium
更改日期: 2022/10/27
'''
#数据清洗功能类
class piecedataprogress() :
    def __init__(self, pbase, data) -> None:
        self._data = data

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
        self._base.baselog('序列集平均长度 {0} / Everage Length Of Sequece Data is {0} bp'.format(seq_everage_len(self._data)))
        self._base.baselog('去重后序列集共 {0} 条 / Number After Duplicate Remove is {0}'.format(len(seq_set(self._data))))

    '''
    创建人员: Nerium
    创建日期: 2022/10/25
    更改人员: Nerium
    更改日期: 2022/10/28
    '''
    #先长度清洗，后去重，得到的是最后的结果
    def filt_data(self) : 
        ids_list = seq_after_filt_len(self._data)
        resdata = seq_keep1id_after_set(seq_set({id : self._data[id] for id in ids_list}))

        self._data = seq_remove_unclassfied(resdata)
        return self._data

    '''
    创建人员: Nerium
    创建日期: 2022/10/25
    更改人员: Nerium
    更改日期: 2022/10/25
    '''
    #清洗完成的数据，写入文件中
    def write_in_file(self, file_path) : 
        self._base.baselog('清洗完成后序列共 {0} 条 / Number After Progress is {0}'.format(len(self._data)))
        with open(file_path, 'w') as tf :
            for id, seq in self._data.items() :
                tf.write('>' + id)
                tf.write('\n' + seq + '\n')