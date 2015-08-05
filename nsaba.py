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

        raw_coords = cls.aba['si_df'].loc[:, 'mni_x':'mni_z'].as_matrix().astype(float)
        cls.aba['mni_coords'] = spatial.KDTree(raw_coords)
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
        self.term_act = []

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
        """Checks if ID is registered"""
        if (self.ns['features_df']['pmid'] == ID).any():
            return True
        else:
            return False

    def coord_to_ids(self, coord):
        """Uses the study dictionary above to find study ids from x,y,z coordinates"""
        # Use
        else:
            return "These coordinates don't match any studies"

    def id_to_terms(self, ID):
        """Finds all of the term heat values of a given ID"""
        if self.is_id(ID):
            ind = int(np.squeeze(np.where(self.ns['id_x_features'] == ID)[0]))
            self.ns['temp_term'] = list(self.ns['mni_term_table'].iloc[ind][1:])
            return self.ns['temp_term']

        else:
            return 'Not an ID of a study in NMI space'

    def id_to_coords(self, ID):
        """ Finds coordinates associated with a given study ID """
        if self.is_id(ID):
            self.ns['temp_coords'] = self.ns['id_dict'][ID]
            return self.ns['temp_coords']
        else:
            return 'Not an ID of a study in NMI space'

    def coord_to_terms(self, coord):
        """Returns the vector of term heats for a given (x,y,z) coordinate set.
        If there are multiple studies that mention the same coordinates, the average is taken."""
        if self.is_location(coord):
            ids = self.coord_to_ids(coord)
            if len(ids) == 1:
                return self.id_to_terms(ids[0])
            else:
                temp = np.zeros((len(ids), 3406))
                for i in xrange(len(ids)):
                    temp[i, :] = self.id_to_terms(ids[i])
                self.ns['temp_term'] = list(np.mean(temp, 0))
                return self.ns['temp_term']
        else:
            return 'not valid location'

    def term_to_ids(self, term, thresh):
        """Matches a term to the IDs of studies that use that term above a given threshold"""
        if self.is_term(term):
            term_all_ids = np.array(self.ns['mni_term_table'][term])
            id_inds = np.squeeze(np.where(term_all_ids > thresh))
            self.ns['temp_IDs'] = self.ns['id_x_features'][id_inds]
            return self.ns['temp_IDs']
        else:
            return 'This is not a valid term'

    def term_to_coords(self, term, thresh):
        '''Finds the coordinates that are associated with a given term up to a given threshold'''
        ids = self.term_to_ids(term, thresh)
        self.ns['temp_coords'] = [self.id_to_coords(i) for i in ids]
        return self.ns['temp_coords']

    # estimating term weights of unknown location
    def term_vector_of_unknown_point(self, coord, maxDist):
        '''Estimates the terms of an unknown point by drawing from known points around it using a sphere of radius maxDist'''
        self.ns['termVect'] = []
        if self.is_location(coord) == True:
            self.ns['termVect'] = self.coord_to_terms(coord)
            print 'This point exists!'
            return self.ns['termVect']
        else:
            ind_sphere = self.sphere(coord, maxDist)
            weight_vect = self.assign_weights(ind_sphere)
            self.ns['termVect'] = []
            for layer in xrange(len(ind_sphere)):
                for c in ind_sphere[layer]:
                    if self.is_location(c) == True:
                        print 'found a nearby point'
                        print c
                        if ~np.isnan(np.sum(self.coord_to_terms(c))):
                            temp = self.coord_to_terms(c)
                            self.ns['termVect'].append([t * weight_vect[0, layer] for t in temp])
                            print weight_vect[0, layer]

            if len(self.ns['termVect']) > 1:
                # Need a better way to normalize
                self.ns['termVect'] = np.sum(self.ns['termVect'], 0)

            return self.ns['termVect']

    # NS/ABA methods
    def generate_ns_vector(self, term):

        self.ns['activation_vector'] = np.zeros((len(self.aba['xyz'])))
        c = 0
        for xyz_set in self.aba['mni_coords']:
            xyz_set = np.floor(xyz_set)
            temp_term_vector = self.term_vector_of_unknown_point(xyz_set, 5)
            if not isinstance(temp_term_vector, list):
                temp_term_vector = list(temp_term_vector)
            # print type(temp_term_vector)
            # print temp_term_vector
            if len(temp_term_vector) == 1:
                self.ns['activation_vector'][c] = temp_term_vector[0][
                    int(np.squeeze(np.where(self.ns['terms'] == term)))]
            elif len(temp_term_vector) > 1:
                self.ns['activation_vector'][c] = temp_term_vector[int(np.squeeze(np.where(self.ns['terms'] == term)))]
            else:
                self.ns['activation_vector'][c] = 0
                c += 1
        return self.ns['activation_vector']

    def correlate_ge_ns(self, term, entrez):
        self.generate_ns_vector(term)
        self.get_ge(entrez)
        correlation = np.corrcoeff(self.ns['activation_vector'], self.ge[entrez])
        print 'FUCK YEAAAAHHHHHH'
        return correlation

    def make_ge_ns_mat(self, term, entrez):
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
        for val in self.ns.itervalues():
            if val is None:
                print "Unassigned Nsaba 'ns' static variable: see Nsaba.ns_load(path)"
                return 1
        return 0
