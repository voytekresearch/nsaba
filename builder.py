# -*- coding: utf-8 -*-
"""
builder.py: Contains NsabaBuilder, a class
designed for one-time-use, heavy duty, time consuming
operations with the Nsaba.

Author: Simon Haxby & Torben Noto
"""
from nsaba import Nsaba
import numpy as np

from testtools import not_operational

class NsabaBuilder(Nsaba):
    """Nsaba heavy duty building tasks"""
    def __init__(self):
        Nsaba.__init__(self)

    def _proceed_check(self):
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
        """Returns a dictionary with ABA gene expression coefficient across all genes
        at sampled locations"""

        if self._proceed_check() == 1:
            return 1

        entrez_ids = self.aba['probe_df']['entrez_id'][
            self.aba['probe_df']['entrez_id'].notnull()].unique().astype(int)

        self.get_aba_ge(entrez_ids)

    @not_operational
    def build_sparse_ge_mat(self, mni_grid_size=(200, 200, 200)):
        """Builds sparse 3D MNI numpy grid, and assigns a gene expression pointer to that coordinate"""

        if self._proceed_check() == 1:
            return 1

        # Talk to Torben about implementation
        # Development on hold 8/27/15
        grid_shift = [dim - 0.5*dim for dim in mni_grid_size]
        mni_space = np.zeros(mni_grid_size)
        for coord in self.aba['mni_coords'].data:
            for entrez_id in self.ge['aba']:
                mni_space[coord[0]+grid_shift[0], coord[1]+grid_shift[1], coord[2]+grid_shift[2]] = 0

    @not_operational
    def build_sparse_ns_mat(self, save_location='.'):
        """Builds a 4D matrix of the term heats where we have NS studies """

        if self._proceed_check() == 1:
            return 1

        matrix_size = 100
        ns_big_matrix = np.zeros((matrix_size*2, matrix_size*2, matrix_size*2, 3406))

        for x in xrange(matrix_size*2):
            for y in xrange(matrix_size*2):
                for z in xrange(matrix_size*2):
                    if self.is_coord((x-matrix_size, y-matrix_size, z-matrix_size)):
                        ns_big_matrix[x][y][z][:] = self.coord_to_terms((x-matrix_size, y-matrix_size, z-matrix_size))
            if np.mod(x, 2):
                print str(x) + ' percent complete'

        self.ns['ns_study_matrix'] = ns_big_matrix
        if save_location:
            np.save(save_location, ns_big_matrix)
        else:
            print 'You have the option to save this matrix by inputing a save location and filename'
        return ns_big_matrix
