'''
创建人员: Nerium
创建日期: 2022/09/29
更改人员: Nerium
更改日期: 2022/09/30
'''

from piece.piecedefine import *

class pieceevaluate() :
    def __init__(self, pbase, nonconser_sort, conser, primer_dict) -> None:
        self._nonconser_sort = nonconser_sort
        self._conser = conser
        self._primer_dict = primer_dict

        self._base = pbase

    #根据条件从区间中过滤出合适的保守区间和非保守区间(conser1)nonconser(conser2)
    #nonconser_sort是根据多样性由高到低排序的，conser是根据位置排序的
    def filter_area(self, minlen=100, hpcnt=5) :
        posmem = []
        for area in self._nonconser_sort :
            if area[1]-area[0] < minlen : continue

            for idx, rang in enumerate(self._conser) :
                if rang[0] > area[1] and self._conser[idx-1][1] < area[0] and (len(self._primer_dict[rang[0]][0].keys()) < hpcnt and len(self._primer_dict[rang[0]][1].keys()) < hpcnt) :
                    posmem.append((self._conser[idx-1], self._conser[idx]))
                    break
        return posmem