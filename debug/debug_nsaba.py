# -*- coding: utf-8 -*-
"""
Nsaba Main Module Debugging/Profiling Script
Created on Mon Jul 13 16:33:15 2015

@author: Simon Haxby
"""

from nsaba import Nsaba

if __name__ == '__main__':

    # Simon Path IO
    data_dir = 'C:\Users\\John\\Documents\\Code\\Python\\nsaba\\data_dir'
    Nsaba.aba_load(data_dir)
    Nsaba.ns_load(data_dir)

    # Getting gene expression for Entrez IDs: 733, 33, 88
    A = Nsaba()
    A.get_aba_ge([733, 33, 88])
    # A.ge.keys()

    # Getting NS activation vector for term 'attention'
    # term = 'attention'
    # A.get_ns_act(term, thresh=.3, k=20, search_radii=3)
    #
    # Make matrix combining activation and gene expression vectors
    # print A.matrix_builder(term, [733, 33, 88])

    # Coordinate to gene expression
    A.coords_to_ge([[10, 20, 30]], [733, 33, 88], search_radii=20)
