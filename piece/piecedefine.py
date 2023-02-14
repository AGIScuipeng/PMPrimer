'''
创建人员: Nerium
创建日期: 2022/08/31
更改人员: Nerium
更改日期: 2023/02/14
'''

#全局变量都在此声明
BASE_DEBUG_LEVEL0, BASE_DEBUG_LEVEL1, BASE_DEBUG_LEVEL2, BASE_DEBUG_LEVEL3 = 0, 1, 2, 4
CMD_PARAMETER_ALLDESIGN = ['rank1', 'rank2', 'rankall', 'haplo', 'threshold:0.xx', 'primer']

DEFAULT_DNA_SINGLE_LIST = ['A', 'T', 'C', 'G', '-']
DEFAULT_DNA_REFLECT_DICT = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

GENE_RELEASE ={'R' : ['G', 'A'],
            'Y' : ['T', 'C'],
            'K' : ['G', 'T'],
            'M' : ['A', 'C'],
            'S' : ['G', 'C'],
            'W' : ['A', 'T'],
            'B' : ['G', 'T', 'C'],
            'D' : ['G', 'A', 'T'],
            'H' : ['A', 'C', 'T'],
            'V' : ['G', 'C', 'A'],
            'N' : ['A', 'T', 'C', 'G'],
        }

GENE_DEGENE = {'AG': 'R', 
            'CT': 'Y', 
            'GT': 'K', 
            'AC': 'M', 
            'CG': 'S', 
            'AT': 'W', 
            'CGT': 'B', 
            'AGT': 'D', 
            'ACT': 'H', 
            'ACG': 'V', 
            'ACGT': 'N'}

PLATFORM_WINDOWS, PLATFORM_LINUX = 'Windows', 'Linux'
PLATFORM_TODO = {PLATFORM_WINDOWS: 'muscle5.1.win64.exe', PLATFORM_LINUX: 'muscle5.1.linux_intel64'}

DATA2JSON_REFELCT = {0: 'F', 1 : 'R'}
JSON2DATA_REFELCT = {'F': 0, 'R': 1}