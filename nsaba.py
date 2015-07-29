# -*- coding: utf-8 -*-
"""
Neurosynth & Allen Brain Atlas: nsaba
Authors: Simon Haxby, Scott Susi, Torben Noto

Last Updated: 7/26/2015
"""
import numpy as np
import pandas as pd
import os

class Nsaba(object):
    aba = {
        'exp_mat': None,
        'exp_df': None,
        'probe_mat': None,
        'probe_df': None,
        'si_mat': None,
        'si_df': None,
        'mni_coords': None
    }

    ns = {
        'mni_term_table': None,
        'terms': None,
        'study_ids': None,
        'unique_ids': None,
        'ids_x_features': None,
        'database_labels': None,
        'x_mni': None,
        'y_mni': None,
        'z_mni': None,
    }

    @classmethod
    def aba_load(cls, path='.', csv_names=None):

        if not csv_names:
            csv_names = [
                'MicroarrayExpression.csv',
                'SampleAnnot.csv',
                'Probes.csv']

        if len(csv_names) != 3:
            raise IndexError, "'csv_names' must a list of 3 'str' variables"

        print 'This may take a minute or two ...'

        csv_path = os.path.join(path, csv_names[0])
        cls.aba['exp_df'] = pd.read_csv(csv_path)
        cls.aba['exp_mat'] = cls.aba['exp_df'].as_matrix()
        print '%s loaded.' % csv_names[0]

        csv_path = os.path.join(path, csv_names[1])
        cls.aba['si_df'] = pd.read_csv(csv_path)
        cls.aba['si_mat'] = cls.aba['si_df'].as_matrix()
        print '%s loaded.' % csv_names[1]

        csv_path = os.path.join(path, csv_names[2])
        cls.aba['probe_df'] = pd.read_csv(csv_path)
        cls.aba['probe_mat'] = cls.aba['probe_df'].as_matrix()
        print '%s loaded.' % csv_names[2]

        cls.aba['mni_coords'] = cls.aba['si_mat'][1:, 10:].astype(float)
        print "Nsaba.aba['mni_coords'] initialized."
        return 0

    @classmethod
    def ns_load(cls, ns_path, ns_files=None):
        """ Initialization 'ns' dictionary """
        if not ns_files:
            ns_files = ['database.txt', 'features.txt']

        df = pd.read_table(os.path.join(ns_path, ns_files[0]))
        mni_space = df[df['space'].str.contains("MNI")]

        cls.ns['x_mni'] = list(mni_space['x'])
        cls.ns['y_mni'] = list(mni_space['y'])
        cls.ns['z_mni'] = list(mni_space['z'])
        cls.ns['study_ids'] = list(mni_space['id'])
        cls.ns['unique_ids'] = list(set(cls.ns['study_ids']))  # removing duplicates
        print '%s keys loaded.' % ns_files[0]

        term_table = pd.read_table(os.path.join(ns_path, ns_files[1])
        cls.ns['mni_term_table'] = term_table.loc[term_table['pmid'].isin(cls.ns['unique_ids'])]
        cls.ns['terms'] = term_table.columns.values
        cls.ns['id_x_features'] = np.array(cls.ns['mni_term_table'])
        cls.ns['database_labels'] = df.columns.values
        print '%s keys loaded.' % ns_files[1]
        return 0

    def __init__(self):

        self.ge = {}
        # ...
        # ...

    def get_ge(self, entrez_ids):

        if self.__check_static_members() == 1:
            return 1

        try:
            iter(entrez_ids)
        except TypeError:
            print "Invalid parameter form; please contain entrez ids as 'str' types in iterable container"
            return 1
        else:
            if isinstance(entrez_ids, str):
                print "Invalid parameter form; please contain entrez ids as 'str' types in iterable container"
                return 1

        for entrez_id in entrez_ids:
            probe_ids = pd.DataFrame(self.aba['probe_df'].loc[self.aba['probe_df'][5] == \
                                                              entrez_id]).index.tolist()

            if len(probe_ids) == 0:
                print 'Entrez ID: %s not registered with ABA database' % entrez_id
                continue

            ge_df = pd.DataFrame(self.aba['exp_df'].loc[probe_ids.pop(0) - 1, :])
            for probe_id in probe_ids:
                ge_df = ge_df.join(self.aba['exp_df'].loc[probe_id - 1, :])

            ge_mat = ge_df.as_matrix().astype(float)[1:, :]
            self.ge[entrez_id] = []

            self.ge[entrez_id].append([np.mean(row) for row in ge_mat])
            diff_arr = [np.amax(row) - np.amin(row) for row in ge_mat]
            self.ge[entrez_id].append(diff_arr)
            self.ge[entrez_id].append(np.amax(diff_arr))
            self.ge[entrez_id].append(len(ge_mat[0]))

        return 0

    def make_ge_ns_mat(self):
        if self.__check_static_members() == 1:
            return 1

    def get_ns_act(self, query, term_freq=.005):
        if self.__check_static_members() == 1:
            return 1

    def __check_static_members(self):
        for val in self.aba.itervalues():
            if val is None:
                print "Unassigned Nsaba 'aba' static variable: see Nsaba.aba_load(path)"
                return 1
        if self.ns_db is None:
            print "Unassigned Nsaba 'ns_db' static variable: see Nsaba.ns_load(path)"
            return 1
        return 0


# -----------------------------------------------------------

# Example Usage

'''
        
from nsaba import Nsaba

ns_path = 'blah/blah'
aba_path = 'blah/bluh'

Nsaba.aba_load(aba_path)
Nsaba.ns_load(ns_path)

A = Nsaba()

A.get_ge('7889')
A.get_ns_act('Alzhiemers')
mat = A.make_ge_ns_mat()

B = Nsaba()

'''

# Run science here

# -----------------------------------------------------------
