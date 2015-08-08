# -*- coding: utf-8 -*-
"""
Nsaba Testing Framework
Created on Mon Jul 13 16:33:15 2015

@author: Simon Haxby
"""

from nsaba import Nsaba

if __name__ == '__main__':

    # Simon Path IO
    data_dir = 'C:\Users\\John\\Documents\\Code\\Python\\nsaba\\data_dir'
    Nsaba.aba_load(data_dir)
    Nsaba.ns_load(data_dir)

    A = Nsaba()
    A.get_aba_ge([733, 33, 88])
    A.ge.keys()

    term = 'attention'
    A.get_ns_act(term, thresh=.2)
    A.get_ns_act(term, thresh=.2, method='sphere')

    print A.make_ge_ns_mat(term, 733)