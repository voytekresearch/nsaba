# -*- coding: utf-8 -*-
"""
nsaba.py: (N)euro(s)ynth, (A)llen (B)rain (A)tlas
Methods to analyze genome-scale gene expression data
from the Allen Human Brain Atlas in conjunction with
fMRI activation maps from Neurosynth

Authors: Simon Haxby & Torben Noto
"""
import pickle
import os
import itertools
import numpy as np
import pandas as pd
from scipy import spatial


class NsabaBase(object):
    """Contains essential base data structures and methods which derived
    Nsaba classes all depend upon"""

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
            raise IndexError("'csv_names' must a list of 3 'str' variables")

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

        mni_coords = cls.aba['si_df'].loc[:, 'mni_x':'mni_z'].as_matrix().astype(float)
        cls.aba['mni_coords'] = spatial.KDTree(mni_coords)
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

        mni_coords = cls.ns['database_df'].loc[:, 'x':'z'].as_matrix().astype(float)
        cls.ns['mni_coords'] = spatial.KDTree(mni_coords)
        print "Nsaba.ns['mni_coords'] initialized."

        cls.ns['id_dict'] = {}
        c = 0

        for i in cls.ns['database_df'].loc[:, 'id']:
            if i not in cls.ns['id_dict']:
                cls.ns['id_dict'][i] = [(np.floor(cls.ns['database_df']['x'].iloc[c]), np.floor(cls.ns['database_df']['y'].iloc[c]), np.floor(cls.ns['database_df']['z'].iloc[c]))]
                c += 1
            else:
                cls.ns['id_dict'][i].append((np.floor(cls.ns['database_df']['x'].iloc[c]), np.floor(cls.ns['database_df']['y'].iloc[c]), np.floor(cls.ns['database_df']['z'].iloc[c])))
                c += 1
        return 0

class Nsaba(NsabaBase):

    def __init__(self):
        self.ge = {}
        self.term = {}
        self.__ns_weight_f = lambda r: 1. / r ** 2

    def check_entrez_struct(self, entrez_ids):
        """Checks if 'entrez_ids' parameter is an non-str iterable"""
        try:
            iter(entrez_ids)
        except TypeError:
            print "Invalid parameter form; please contain entrez ids as 'str' types in iterable container"
            return 1
        else:
            if isinstance(entrez_ids, str):
                print "Invalid parameter form; please contain entrez ids as 'str' types in iterable container"
                return 1
            else:
                return 0

    def get_aba_ge(self, entrez_ids):
        """Retrieves an stores gene expression coefficients in ABA dictionary based on a
        a passed list of Entrez IDs"""

        if self._check_static_members() == 1:
            return 1
        if self.check_entrez_struct(entrez_ids) == 1:
            return 1

        for entrez_id in entrez_ids:
            probe_ids = self.aba['probe_df'].loc[self.aba['probe_df']['entrez_id']
                                                 == entrez_id]['probe_id'].tolist()

            if len(probe_ids) == 0:
                print 'Entrez ID: %s not registered with ABA database' % entrez_id
                continue

            ge_df = self.aba['exp_df'].loc[self.aba['exp_df']['probe_id'].isin(probe_ids)]
            ge_mat = ge_df.as_matrix().astype(float)[:, 1:].T
            self.ge[entrez_id] = np.mean(ge_mat, axis=1)

    def pickle_ge(self, pkl_file="Nsaba_ABA_ge.pkl", output_dir='.'):
        if self._check_static_members() == 1:
            return 1
        pickle.dump(self.ge, open(os.path.join(output_dir, pkl_file), 'wb'))
        print "%s successfully created" % pkl_file

    def load_ge_pickle(self, file_path="ABA_ge.pkl"):
        self.ge = pickle.load(open(file_path, 'rb'))
        print "'ge' dictionary successfully loaded"

    def is_term(self, term):
        """Checks if this term is in the neurosynth database """
        if self._check_static_members() == 1:
            return 1
        if term in self.ns['features_df'].columns:
            return True
        else:
            return False

    def is_id(self, ID):
        """Checks if ID is registered """
        if self._check_static_members() == 1:
            return 1
        if (self.ns['features_df']['pmid'] == ID).any():
            return True
        else:
            return False

    # def is_id(self, study_id):
    #     """Checks if ID is registered """
    #     if len(self.ns['features_df']['pmid'] == study_id) > 0:
    #         if len(self.ns['database_df']['id'] == study_id) > 0:  # added layer because motherfuckers are missing data
    #             return True
    #     else:
    #         return False

    def is_coord(self, coordinate):
        """Checks if an x,y,z coordinate in list form matches a NS data point"""
        zipped_coordinates = np.squeeze(zip(self.ns['mni_coords'].data))
        for this_coordinate in zipped_coordinates:
            if this_coordinate[0] == coordinate[0]:
                if this_coordinate[1] == coordinate[1]:
                    if this_coordinate[2] == coordinate[2]:
                        return True
        return False

    def coord_to_ids(self, coordinate):
        """Uses the study dictionary above to find study ids from x,y,z coordinates """
        ids = []
        for i, coords in self.ns['id_dict'].items():
            for this_coordinate in coords:
                if this_coordinate[0] == coordinate[0]:
                    if this_coordinate[1] == coordinate[1]:
                        if this_coordinate[2] == coordinate[2]:
                            if i not in ids:
                                if self.is_id(i):
                                    ids.append(i)
        return ids

    def _id_to_terms(self, study_id):
        """Finds all of the term heat values of a given ID """
        if self.is_id(study_id):
            term_vector_off_by_1 = np.squeeze(self.ns['features_df'].loc[self.ns['features_df']['pmid'] == study_id].as_matrix())
            # shifting to remove id index from vector
            return term_vector_off_by_1[1:]
        else:
            return 'Invalid study id'

    def coord_to_terms(self, coord):
        ids = self.coord_to_ids(coord)
        if len(ids) == 1:
            # one study
            terms = self._id_to_terms(ids)
        elif len(ids) > 1:
            # multiple studies
            temp = []
            for multiple_id in ids:
                temp.append(self._id_to_terms(multiple_id))
                terms = np.mean(temp, 0)
        else:
            # print 'No terms found for id' + str(ids) + 'using coordinates:' + str(coord)
            terms = []
        return terms

    def _term_to_coords(self, term, thresh=0):
        """Finds coordinates associated with a given term.
        Returns NS coordinate tree and ID/coordinate/activation DataFrame"""
        term_ids_act = self.ns['features_df'].loc[self.ns['features_df'][term] > thresh, ['pmid', term]]
        term_ids = term_ids_act['pmid'].tolist()
        term_coords = self.ns['database_df'].loc[self.ns['database_df']['id'].isin(term_ids)]
        try:
            ns_coord_tree = spatial.KDTree(term_coords.loc[:, 'x':'z'].as_matrix().astype(float))
        except ValueError:
            print "No studies with term: '%s' and threshold: %.2f found" % (term, thresh)
            return 1, 1
        else:
            term_ids_act.rename(columns={'pmid': 'id'}, inplace=True)
            return ns_coord_tree, term_coords.merge(term_ids_act)

    def _sphere(self, xyz, coord_tree, max_rad=5):
        """Returns 3D Array containing coordinates in each layer of the sphere """
        sphere_bucket = []
        set_bucket = []

        # Needs work; generalize
        for i, r in enumerate(range(max_rad, 0, -1)):
            pts = coord_tree.query_ball_point(xyz, r)
            set_bucket.append(set(map(tuple, coord_tree.data[pts])))

        for i in range(0, 3):
            sphere_bucket.append(list(set_bucket[i].difference(set_bucket[i + 1])))
        sphere_bucket.append(list(set_bucket[3]))
        rev_iter = reversed(sphere_bucket)

        return [layer for layer in rev_iter]

    def _knn_search(self, xyz, coord_tree, max_rad=5, k=20):
        """KNN search of NS coordinates about ABA coordinates """
        r, inds = coord_tree.query(xyz, k)
        return inds[r < max_rad], r[r < max_rad]

    def _get_act_values(self, bucket, weight, term, ns_coord_act_df):
        """Returns weighted NS activation """
        bucket_act_vec = []
        for coords in bucket:
            coord = ns_coord_act_df[(ns_coord_act_df['x'] == coords[0])
                                    & (ns_coord_act_df['y'] == coords[1])
                                    & (ns_coord_act_df['z'] == coords[2])][term]
            bucket_act_vec.append(coord.mean())

        return np.array(bucket_act_vec)*weight

    def _knn_method(self, term, ns_coord_act_df, ns_coord_tree, search_radii, k):
        """KNN method """
        for irow, xyz in enumerate(self.aba['mni_coords'].data):
            coord_inds, radii = self._knn_search(xyz, ns_coord_tree, search_radii, k)
            coords = ns_coord_tree.data[coord_inds]
            weight = self.__ns_weight_f(radii)
            weighted_means = self._get_act_values(coords, weight, term, ns_coord_act_df)
            if len(weighted_means) == 0:
                self.term[term]['aba_void_indices'].append(irow)
            else:
                act_coeff = np.sum(weighted_means) / np.sum(weight)
                self.term[term]['ns_act_vector'].append(act_coeff)

    def _sphere_method(self, term, ns_coord_act_df, ns_coord_tree, search_radii):
        """Sphere buckets method"""
        for irow, xyz in enumerate(self.aba['mni_coords'].data):
            sphere_bucket = self._sphere(xyz, ns_coord_tree, search_radii)
            sphere_vals = [0, 0]
            for w, bucket in enumerate(sphere_bucket):
                weight = self.__ns_weight_f(w + 1)
                bucket_mean = np.mean(self._get_act_values(bucket, weight, term, ns_coord_act_df))
                if np.isnan(bucket_mean):
                    sphere_vals[0] += 0
                    sphere_vals[1] += 0
                else:
                    sphere_vals[0] += bucket_mean
                    sphere_vals[1] += weight
            if sphere_vals[1] == 0:
                self.term[term]['aba_void_indices'].append(irow)
            else:
                act_coeff = sphere_vals[0] / sphere_vals[1]
                self.term[term]['ns_act_vector'].append(act_coeff)

    def get_ns_act(self, term, thresh=0, method='knn', search_radii=5, k=None):
        """Generates NS activation vector about ABA MNI coordinates  """
        if self._check_static_members() == 1:
            return 1
        if not self.is_term(term):
            print "'%s' is not a registered term." % term
            return 1

        ns_coord_tree, ns_coord_act_df = self._term_to_coords(term, thresh)
        if type(ns_coord_tree) == int:
            return 1

        self.term[term] = {}
        self.term[term]['ns_act_vector'] = []
        self.term[term]['aba_void_indices'] = []

        if method == 'knn':
            if k is None:
                k = 20
            self._knn_method(term, ns_coord_act_df, ns_coord_tree, search_radii, k)
        elif method == 'sphere':
            if k is not None:
                raise ValueError("'k' parameter cannot be used with 'sphere' method.")
            self._sphere_method(term, ns_coord_act_df, ns_coord_tree, search_radii)
        else:
            raise TypeError("'%s' is not a valid parameter value for 'method' parameter, use either 'knn' or 'sphere"
                            % method)

    def make_ge_ns_mat(self, ns_term, entrez_ids):
        if self._check_static_members() == 1:
            return 1
        if self.check_entrez_struct(entrez_ids) == 1:
            return 1

        if ns_term in self.term and all([key in self.ge for key in entrez_ids]):
            ge_ns_mat = []
            for entrez_id in entrez_ids:
                aba_indices = np.array([i for i in range(len(self.aba['mni_coords'].data))
                                        if i not in self.term[ns_term]['aba_void_indices']])
                ge_ns_mat.append(self.ge[entrez_id][aba_indices])
            ge_ns_mat.append(self.term[ns_term]['ns_act_vector'])
            return np.vstack(ge_ns_mat).T
        else:
            print "Either term['%s'] or one or more Entrez ID keys does not exist; please check arguments" \
                  % ns_term

    def _coord_to_ge(self, coord, entrez_ids, search_radii=10, k=20):
        """Returns weighted ABA gene expression mean for some MNI coordinate based
        on a list of passed Entrez IDs"""

        ge_for_coord = []
        for entrez_id in entrez_ids:
            coord_inds, radii = self._knn_search(coord, self.aba['mni_coords'], search_radii, k)
            if len(coord_inds) == 0:
                print "No ABA coordinates are within search radius of specified coordinate"
                break
            weight = self.__ns_weight_f(radii)
            local_ge = self.ge[entrez_id][coord_inds]
            weighted_ge_mean = np.sum(local_ge*weight)/np.sum(weight)
            ge_for_coord.append(weighted_ge_mean)

        return ge_for_coord

    def coords_to_ge(self, coords, entrez_ids, search_radii=10, k=20):
        """Returns Returns weighted ABA gene expression mean for a list MNI coordinate based
        on a list of passed Entrez IDs"""
        if self._check_static_members() == 1:
            return 1
        if self.check_entrez_struct(entrez_ids) == 1:
            return 1

        ge_for_coords = []
        for coord in coords:
            ge_for_coord = self._coord_to_ge(coord, entrez_ids, search_radii, k)
            ge_for_coords.append(ge_for_coord)

        return np.array(ge_for_coords)

    def set_ns_weight_f(self, f):
        try:
            print "Test: f(e) = %.2f" % f(np.e)
            self.__ns_weight_f = f
        except TypeError:
            print "'f' is improper, ensure 'f' receives only one parameter and returns a numeric type"

    def _check_static_members(self):
        for val in self.aba.itervalues():
            if val is None:
                print "Unassigned Nsaba 'aba' static variable: see Nsaba.aba_load(path)"
                return 1
        for val in self.ns.itervalues():
            if val is None:
                print "Unassigned Nsaba 'ns' static variable: see Nsaba.ns_load(path)"
                return 1
        return 0
