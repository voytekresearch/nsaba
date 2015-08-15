# -*- coding: utf-8 -*-
"""
builder.py:
Contains NsabaBuilder class designed single use,
heavy duty Nsaba operations such as creating a GE
matrix for all registered Entrez IDs.

Author: Simon Haxby
"""

from nsaba import Nsaba
import numpy as np


class NsabaBuilder(Nsaba):
    """ Nsaba heavy duty building tasks"""
    def __init__(self):
        Nsaba.__init__(self)

    def __proceed_check(self):
        warning_flag = True
        while warning_flag:
            y_n = raw_input("WARNING: this operation can take upwards of an hour, proceed? (Y/n): ")
            if y_n == 'Y':
                warning_flag = False
            elif y_n == 'n':
                return 1
            else:
                print "Invalid response: %s" % y_n
        return 0

    def get_aba_ge_all(self):
        """ Returns a dictionary with ABA gene expression coefficient across all genes
        at sampled locations"""

        if self.__check_static_members() == 1:
            return 1
        if self.__proceed_check() == 1:
            return 1

        entrez_ids = self.aba['probe_df']['entrez_id'][
            self.aba['probe_df']['entrez_id'].notnull()].unique().astype(int)

        self.get_aba_ge(entrez_ids)

    def build_sparse_ge_mat(self, mni_grid_size=(200,200,200)):
        """ Builds sparse 3D MNI numpy grid, and assigns a gene expression pointer to that coordinate"""

        mni_space = np.zeros(mni_grid_size)




