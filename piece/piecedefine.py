'''
创建人员: Nerium
创建日期: 2022/08/31
更改人员: Nerium
更改日期: 2022/12/08
'''

#全局变量都在此声明
BASE_DEBUG_LEVEL0, BASE_DEBUG_LEVEL1, BASE_DEBUG_LEVEL2, BASE_DEBUG_LEVEL3 = 0, 1, 2, 4
CMD_PARAMETER_ALLDESIGN = ['rank1', 'rank2', 'rankall', 'haplo', 'threshold:0.xx', 'primer']

DEFAULT_DNA_SINGLE_LIST = ['A', 'T', 'C', 'G', '-']
DEFAULT_DNA_REFLECT_DICT = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

PLATFORM_WINDOWS, PLATFORM_LINUX = 'Windows', 'Linux'
PLATFORM_TODO = {PLATFORM_WINDOWS: 'muscle5.1.win64.exe', PLATFORM_LINUX: 'muscle5.1.linux_intel64'}