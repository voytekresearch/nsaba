# -*- coding: utf-8 -*-
"""
Nsaba Testing Framework
Created on Mon Jul 13 16:33:15 2015

@author: Simon Haxby
"""

from nsaba import Nsaba
import timeit

if __name__ == '__main__':

    # Simon Path IO
    data_dir = 'C:\Users\\John\\Documents\\Code\\Python\\nsaba\\data_dir'
    Nsaba.aba_load(data_dir)
    Nsaba.ns_load(data_dir)

    A = Nsaba()
    A.get_aba_ge([733, 33, 88])
    A.ge.keys()

    A.get_aba_ge_all()

    #A.coord_to_ge([10, 20, 30], [733, 33, 88], search_radii=20)

    # term = 'attention'
    # A.get_ns_act(term, thresh=.2, k=20, search_radii=3)
    # A.get_ns_act(term, thresh=.2, method='sphere')
    # print A.make_ge_ns_mat(term, [733, 33, 88])