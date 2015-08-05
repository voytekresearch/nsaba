# -*- coding: utf-8 -*-
"""
Neurosynth & Allen Brain Atlas: nsaba
Authors: Simon Haxby, Scott Susi, Torben Noto

Last Updated: 7/26/2015
"""
import numpy as np
import pandas as pd
import os
import itertools

from scipy import spatial


class Nsaba(object):

    aba = {
        'exp_df': None,
        'probe_df': None,
        'si_df': None,
        'mni_coords': None
    }

    ns = {
        'database_df': None,
        'features_df': None,
    }

    @classmethod
    def aba_load(cls, aba_path=".", csv_names=None):
        """Initialization of 'aba' dictionary"""
        if not csv_names:
            csv_names = [
                'MicroarrayExpression.csv',
                'SampleAnnot.csv',
                'Probes.csv']

        if len(csv_names) != 3:
            raise IndexError, "'csv_names' must a list of 3 'str' variables"

        print 'This may take a minute or two ...'

        csv_path = os.path.join(aba_path, csv_names[1])
        cls.aba['si_df'] = pd.read_csv(csv_path)
        print '%s loaded.' % csv_names[1]

        csv_path = os.path.join(aba_path, csv_names[0])
        cls.aba['exp_df'] = pd.read_csv(csv_path, header=None)
        cls.aba['exp_df'].columns = list(
            itertools.chain.from_iterable(
                [['probe_id'], range(cls.aba['si_df'].shape[0])]))

        print '%s loaded.' % csv_names[0]

        csv_path = os.path.join(aba_path, csv_names[2])
        cls.aba['probe_df'] = pd.read_csv(csv_path)
        print '%s loaded.' % csv_names[2]

        cls.aba['mni_coords'] = cls.aba['si_df'].loc[:, 'mni_x':'mni_z'].as_matrix().astype(float)
        print "Nsaba.aba['mni_coords'] initialized."


    @classmethod
    def ns_load(cls, ns_path=".", ns_files=None):
        """Initialization of 'ns' dictionary"""
        if not ns_files:
            ns_files = ['database.txt', 'features.txt']

        df = pd.read_table(os.path.join(ns_path, ns_files[0]))
        cls.ns['database_df'] = df.loc[df.space == 'MNI', ['id', 'x', 'y', 'z']]
        print "%s loaded." % ns_files[0]

        cls.ns['features_df'] = pd.read_table(os.path.join(ns_path, ns_files[1]))
        print "%s loaded." % ns_files[1]

    def __init__(self):

        self.ge = {}
        self.term = {}

        # Will integrate these into self.term
        self.ns_coord_tree = None
        self.ns_coord_act_df = None
        self.ns_act_vector = None
        self.aba_void_indices = None

    def get_aba_ge(self, entrez_ids):

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
            probe_ids = self.aba['probe_df'].loc[self.aba['probe_df']['entrez_id']
                                                 == entrez_id]['probe_id'].tolist()

            if len(probe_ids) == 0:
                print 'Entrez ID: %s not registered with ABA database' % entrez_id
                continue

            ge_df = self.aba['exp_df'].loc[self.aba['exp_df']['probe_id'].isin(probe_ids)]
            ge_mat = ge_df.as_matrix().astype(float)[:, 1:].T

            #self.ge[entrez_id] = []

            self.ge[entrez_id] = np.mean(ge_mat, axis=1)

            # Additional Metrics
            # self.ge[entrez_id].append(np.mean(ge_mat, axis=1))
            # diff_arr = [np.amax(row) - np.amin(row) for row in ge_mat]
            # self.ge[entrez_id].append(diff_arr)
            # self.ge[entrez_id].append(np.amax(diff_arr))
            # self.ge[entrez_id].append(len(ge_mat[0]))
            # self.ge[entrez_id] = np.squeeze(self.ge[entrez_id])


    # NeuroSynth book-keeping methods

    def is_term(self, term):
        """ Checks if this term is in the neurosynth database """
        if term in self.ns['features_df'].columns:
            return True
        else:
            return False

    def is_id(self, ID):
        """ Checks if ID is registered """
        if (self.ns['features_df']['pmid'] == ID).any():
            return True
        else:
            return False

    def coord_to_ids(self, coord):
        """ Uses the study dictionary above to find study ids from x,y,z coordinates """
        # Use Later?

        return

    def id_to_terms(self, ID):
        """ Finds all of the term heat values of a given ID """
        # Use Later?

        return

    def __term_to_coords(self, term, thresh=0):
        """ Finds coordinates associated with a given term.
        Returns NS coordinate tree and ID/coordinate/activation DataFrame"""

        term_ids_act = self.ns['features_df'].loc[self.ns['features_df'][term] > thresh, ['pmid', term]]
        term_ids = term_ids_act['pmid'].tolist()
        term_coords = self.ns['database_df'].loc[self.ns['database_df']['id'].isin(term_ids)]
        ns_coord_tree = spatial.KDTree(term_coords.loc[:, 'x':'z'].as_matrix().astype(float))

        term_ids_act.rename(columns={'pmid':'id'}, inplace=True)
        return ns_coord_tree, term_coords.merge(term_ids_act)

    def __sphere(self, xyz):
        """ Returns 3D Array containing coordinates in each layer of the sphere """
        sphere_bucket = []
        set_bucket = []

        # Needs work; generalize
        for i, r in enumerate(range(4, 0, -1)):
            pts = self.ns_coord_tree.query_ball_point(xyz, r)
            set_bucket.append(set(map(tuple, self.ns_coord_tree.data[pts])))

        for i in range(0,3):
            sphere_bucket.append(list(set_bucket[i].difference(set_bucket[i+1])))
        sphere_bucket.append(list(set_bucket[3]))
        rev_iter = reversed(sphere_bucket)

        return [layer for layer in rev_iter]

    def __get_act_values(self, bucket, weight, term):
        """ Returns weighted NS activation """
        bucket_act_vec = []
        for coords in bucket:
            bucket_act_vec.append(self.ns_coord_act_df[(self.ns_coord_act_df['x'] == coords[0])
                                 & (self.ns_coord_act_df['y'] == coords[1]) &
                                 (self.ns_coord_act_df['z'] == coords[2])][term].mean())

        return np.mean(bucket_act_vec)*weight

    def get_ns_act(self, term, thresh=0):
        """ Generates NS activation vector about ABA MNI coordinates  """
        if self.__check_static_members() == 1:
            return 1
        if not self.is_term(term):
            print "'%s' is not a registered term."
            return 1

        self.aba_void_indices = []
        self.ns_act_vector = []
        self.ns_coord_tree, self.ns_coord_act_df = self.__term_to_coords(term, thresh)

        for irow, xyz in enumerate(self.aba['mni_coords']):
            sphere_bucket = self.__sphere(xyz)
            sphere_vals = [0, 0]
            for w, bucket in enumerate(sphere_bucket):
                weight = 1./(w+1)
                bucket_mean = self.__get_act_values(bucket, weight, term)
                if np.isnan(bucket_mean):
                    sphere_vals[0] += 0
                    sphere_vals[1] += 0
                else:
                    sphere_vals[0] += bucket_mean
                    sphere_vals[1] += weight
            if sphere_vals[1] == 0:
                self.aba_void_indices.append(irow)
            else:
                act_coeff = sphere_vals[0]/sphere_vals[1]
                self.ns_act_vector.append(act_coeff)


    def make_ge_ns_mat(self, ns_term, entrez):
        if self.__check_static_members() == 1:
            return 1

        # Once dictionary is made, term[ns_term].ns_act_vector
        aba_indices = np.array([i for i in range(len(self.aba['mni_coords'])) if i not in self.aba_void_indices])
        ge = self.ge[entrez][aba_indices]
        return np.vstack((ge, self.ns_act_vector)).T

    def __check_static_members(self):
        for val in self.aba.itervalues():
            if val is None:
                print "Unassigned Nsaba 'aba' static variable: see Nsaba.aba_load(path)"
                return 1
        for val in self.ns.itervalues():
            if val is None:
                print "Unassigned Nsaba 'ns' static variable: see Nsaba.ns_load(path)"
                return 1
        return 0