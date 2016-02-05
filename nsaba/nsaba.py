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
import collections
import numpy as np
import pandas as pd
from scipy import spatial

from sklearn.neighbors import KNeighborsRegressor
from nsabatools import not_operational, preprint


class NsabaBase(object):
    """
    Contains essential base data structures and methods which derived
    Nsaba classes all depend upon.

    Fields
    ------
    __aba : Contains panda.Dataframe objects representating the
        structure of MicroarrayExpression.CSV, SampleAnnot.csv and Probes.CSV
        (default names), as well as numpy.array representing the MNI coordinates
        of each location sampled by ABA.

    __ns : Contains pandas.DataFrame objects representating Neurosynth's
        database.txt and features.txt CSV-style fields.

    Methods
    -------
    aba_load()
    ns_load()
    ns_load_id_dict()

    WARNING: NsabaBase is not meant to instantiated explicitly, only Nsaba should
    be publicly interfaced.

    """
    _aba = {
        'exp_df': None,
        'probe_df': None,
        'si_df': None,
        'mni_coords': None
    }

    _ns = {
        'database_df': None,
        'features_df': None,
    }

    @classmethod
    @preprint('This may take a minute or two ...')
    def aba_load(cls, aba_path=".", csv_names=None):
        """
        Initialization of aba dictionary

        Parameters
        ----------
        aba_path : string, optional
                Root directory of ABA files (defaultly named) MicroarrayExpression.csv,
                SampleAnnot.csv and Probes.txt.
        csv_names: tuple-like, optional
                Tuple specifying alternative names for MicroarrayExpression.csv, SampleAnnot.txt
                and Probes.txt. (NOTE: order affects aba instantiation).
                Default = ('MicroarrayExpression.csv', 'SampleAnnot.csv', 'Probes.txt').
        """

        if not csv_names:
            csv_names = [
                'MicroarrayExpression.csv',
                'SampleAnnot.csv',
                'Probes.csv']

        if len(csv_names) != 3:
            raise IndexError("'csv_names' must a list of 3 'str' variables")
        if not isinstance(aba_path, str):
            raise ValueError("'aba_path' must be a str.")

        csv_path = os.path.join(aba_path, csv_names[1])
        cls._aba['si_df'] = pd.read_csv(csv_path)
        print '%s loaded.' % csv_names[1]

        csv_path = os.path.join(aba_path, csv_names[0])
        cls._aba['exp_df'] = pd.read_csv(csv_path, header=None)
        print '%s loaded.' % csv_names[0]

        cls._aba['exp_df'].columns = list(
            itertools.chain.from_iterable(
                [['probe_id'], range(cls._aba['si_df'].shape[0])]))

        csv_path = os.path.join(aba_path, csv_names[2])
        cls._aba['probe_df'] = pd.read_csv(csv_path)
        print '%s loaded.' % csv_names[2]

        mni_coords = cls._aba['si_df'].loc[:, 'mni_x':'mni_z'].as_matrix().astype(float)
        cls._aba['mni_coords'] = spatial.KDTree(mni_coords)
        print "Nsaba.aba['mni_coords'] initialized.\n"

    @classmethod
    @preprint('This may take a minute or two ...')
    def ns_load(cls, ns_path=".", ns_files=None):
        """
        Initialization of ns dictionary

        Parameters
        ----------
        ns_path : string, optional
                Root directory of neurosynth files (defaultly named) database.txt and features.txt.
        ns_files : tuple-like, optional
                Tuple specifying alternative names for database.txt and features.txt
                (NOTE: order affects ns instantiation). Default = ('database.txt', 'features.txt').
        """

        if not ns_files:
            ns_files = ('database.txt', 'features.txt')

        if len(ns_files) != 2:
            raise IndexError("'ns_files' must a list of 2 'str' variables")
        if not isinstance(ns_path, str):
            raise ValueError("'ns_path' must be a str.")

        df = pd.read_table(os.path.join(ns_path, ns_files[0]))
        cls._ns['database_df'] = df.loc[df.space == 'MNI', ['id', 'x', 'y', 'z']]
        print "%s loaded." % ns_files[0]

        cls._ns['features_df'] = pd.read_table(os.path.join(ns_path, ns_files[1]))
        print "%s loaded." % ns_files[1]

        mni_coords = cls._ns['database_df'].loc[:, 'x':'z'].as_matrix().astype(float)
        cls._ns['mni_coords'] = spatial.KDTree(mni_coords)
        print "Nsaba.ns['mni_coords'] initialized.\n"

    @classmethod
    @preprint('This may take a minute or two ...')
    def ns_load_id_dict(cls):
        """ID dictionary thing needed for doing some NS analyses"""
        cls._check_static_members()
        cls._ns['id_dict'] = {}
        c = 0
        for i in cls._ns['database_df'].loc[:, 'id']:
            if i not in cls._ns['id_dict']:
                cls._ns['id_dict'][i] = [(np.floor(cls._ns['database_df']['x'].iloc[c]),
                                         np.floor(cls._ns['database_df']['y'].iloc[c]),
                                         np.floor(cls._ns['database_df']['z'].iloc[c]))]
                c += 1
            else:
                cls._ns['id_dict'][i].append((np.floor(cls._ns['database_df']['x'].iloc[c]),
                                             np.floor(cls._ns['database_df']['y'].iloc[c]),
                                             np.floor(cls._ns['database_df']['z'].iloc[c])))
                c += 1

    @classmethod
    def _check_static_members(self):
        """ Ensures Nsaba class is not instantiated without initalizing NsabaBase.aba and NsabaBase.ns."""
        for val in self._aba.itervalues():
            if val is None:
                raise AttributeError("Unassigned Nsaba 'aba' static variable: see Nsaba.aba_load(path)")
        for val in self._ns.itervalues():
            if val is None:
                raise AttributeError("Unassigned Nsaba 'ns' static variable: see Nsaba.ns_load(path)")


class Nsaba(NsabaBase):
    """
    Principal Nsaba class.
    Contains methods for data fetching and estimation.
    """

    def __init__(self):
        """Nsaba init method; terminates instantiation if Nsaba.ns or Nsaba.aba are not loaded."""

        self._check_static_members()
        self.ge = {}
        self.term = {}
        self.ns_weight_f = lambda r: 1. / np.power(r, 2)

    def get_ns_struct(self, key=None):
        """
        Returns _ns internal Dictionary or specified sub-structure.
        See class documentation for more information.

        Parameters
        ----------
        key: string, optional
            Name of specified sub-structure of _ns if provided;
            else _ns dictionary is returned.
        """

        if not key:
            return self._ns
        else:
            try:
                return self._ns[key]
            except KeyError:
                opts = " / ".join(self._ns.keys())
                raise KeyError("'key' argument invalid; options are: %s" % opts)

    def get_aba_struct(self, key=None):
        """
        Returns _aba internal dictionary or specified sub-structure.
        See class documentation for more information.

        Parameters
        ----------
        key: string, optional
            Name of specified sub-structure of _aba if provided;
            else _ns dictionary is returned.
       """

        if not key:
            return self._aba
        else:
            try:
                return self._aba[key]
            except KeyError:
                opts = " / ".join(self._aba.keys())
                raise KeyError("'key' argument invalid; options are: %s" % opts)

    def _check_entrez_struct(self, entrez_ids):
        """
        Checks if 'entrez_ids' parameter is an non-str iterable; type-checking method.
        Raises errors if entrez_ids is not as specified above; ensures that methods and
        data structures use 'entrez_ids' are well-behaved.

        Parameters
        ----------
        entrez_ids: List-like
            list-like structure containing NIH Entrez IDs.
        """

        try:
            iter(entrez_ids)
        except TypeError:
            raise TypeError("Invalid parameter form; please contain entrez ids in iterable container")
        else:
            if isinstance(entrez_ids, str):
                raise TypeError("Invalid parameter form; please contain entrez ids in iterable container")

    def est_aba_ge(self, entrez_ids, coords=None, **kwargs):
        """
        Retrieves, estimates and stores gene expression coefficients in ABA dictionary based on a
        a passed list of NIH Entrez IDs.

        Parameters
        ----------
        entrez_ids: List-like
            list-like structure containing NIH Entrez IDs.
        """

        self._check_entrez_struct(entrez_ids)

        for entrez_id in entrez_ids:
            # Fetch probe IDs for Entrez ID
            probe_ids = self._aba['probe_df'].loc[self._aba['probe_df']['entrez_id']
                                                 == entrez_id]['probe_id'].tolist()

            if len(probe_ids) == 0:
                print 'Entrez ID: %s not registered with ABA database' % entrez_id
                continue

            # Return gene expression on given probes across sampled locations.
            ge_df = self._aba['exp_df'].loc[self._aba['exp_df']['probe_id'].isin(probe_ids)]
            ge_mat = ge_df.as_matrix().astype(float)[:, 1:].T

            # Take average gene expression across probes at a given sampled location.

            ge_vec = np.mean(ge_mat, axis=1)
            self.ge[entrez_id] = {}

            if coords is None:
                self.ge[entrez_id]['GE'] = ge_vec
                self.ge[entrez_id]['coord_type'] = 'ABA'

            # Estimate gene expression at custom coordinates
            else:
                if 'knn_args' in kwargs:
                    self.ge[entrez_id]['classifer'] = KNeighborsRegressor(**kwargs['knn_args'])
                else:
                    self.ge[entrez_id]['classifer'] = KNeighborsRegressor()

                X = self._aba['mni_coords']
                y = ge_vec

                self.ge[entrez_id]['classifer'].fit(X.data, y)
                if 'store_coords' in kwargs:
                    if kwargs['store_coords']:
                        self.ge[entrez_id]['coords'] = coords

                if 'coord_type' in kwargs:
                    self.ge[entrez_id]['coord_type'] = kwargs['coord_type']
                else:
                    self.ge[entrez_id]['coord_type'] = 'Custom'

                self.ge[entrez_id]["GE"] = self.ge[entrez_id]['classifer'].predict(coords)

    def ge_ratio(self, entrez_ids, coords=None, **kwargs):
        """
        Calculates the ratio of gene expression at each ABA sampled MNI coordinate
        or custom coordinates.

        NOTE: This methods overwrites previously stored gene expression coefficients.

        Parameters
        ----------
        entrez_ids: (tuple-like) 2
            Entrez IDs of genes whose expression ratio is to be calculated.

        Returns
        -------
        ratio: np.array() [1 x N]
            Array of ratios of gene expression at each ABA sampled MNI coordinate,
            or at custom MNI coordinates where N is the number of sampled locations
            or custom coordinates provided.
        """

        if len(entrez_ids) != 2:
            raise ValueError("Invalid parameter form: entrez_ids should be in a 2-tuple")
        self._check_entrez_struct(entrez_ids)

        self.est_aba_ge(entrez_ids, coords=coords, **kwargs)

        ei1, ei2 = entrez_ids

        ratio = self.ge[ei1]['GE']/self.ge[ei2]['GE']

        return ratio

    def pickle_ge(self, pkl_file="Nsaba_ABA_ge.pkl", output_dir='.'):
        """
        Stores Nsaba.ge as pickle named by 'pkl_file' in directory 'output_dir'.

        Parameters
        ----------
        pkl_file: string, optional
            Name of pickle file.
        output_dir: string, optional
            Name of directory the pickle is to be written to;
            '/' automatically added via os.path.join.
        """

        pickle.dump(self.ge, open(os.path.join(output_dir, pkl_file), 'wb'))
        print "%s successfully created" % pkl_file

    @preprint('This may take a minute or two ...')
    def load_ge_pickle(self, pkl_file="Nsaba_ABA_ge.pkl", path='.'):
        """
        Loads pickle named by 'pkl_file' in directory 'output_dir' into Nsaba.ge.

        Parameters
        ----------
        pkl_file: string, optional
            Name of pickle file.
        path: string, optional
            Path to directory the pickle is written to;
            '/' automatically added via os.path.join.
        """

        self.ge = pickle.load(open(os.path.join(path, pkl_file), 'rb'))
        print "'ge' dictionary successfully loaded"

    def is_gene(self, gene):
        """
        Parameters
        ----------
        gene: int
            Checks whether gene is a registered NIH Entrez ID within ABA.
        """

        if isinstance(gene, str):
            raise ValueError("%s is a string; please pass as a numeric." % gene)
        if gene in self._aba['probe_df']['entrez_id']:
            return True
        else:
            return False

    def is_term(self, term):
        """
        Parameters
        ----------
        term: string
            Checks if this term is in the NS term database.
        """

        if term in self._ns['features_df'].columns:
            return True
        else:
            return False

    def is_id(self, study_id):
        """
        Parameters
        ----------
        study_id: int
            Checks if study_id is a registered NS study ID.
        """
        if any(self._ns['features_df'].pmid == study_id) or \
                any(self._ns['database_df'].id == study_id):
                    return True
        else:
            return False

    def is_coord(self, coord):
        """
        Parameters
        ----------
        term: tuple-like (3)
            Checks if an (x,y,z) MNI coordinate matches an NS data point.
        """

        if len(coord) == 3 and not isinstance(coord, str):
            if self._ns['mni_coords'].query(coord, distance_upper_bound=1)[0] == 0:
                return True
            else:
                return False
        else:
            raise ValueError("MNI coordinate in improper form; must be 3-tuple-like")

    def coord_to_ids(self, coord):
        """
        Uses the study dictionary above to find NS study IDs from x,y,z coordinates.

        Parameters
        ----------
        coordinate: tuple-like (3)
            Checks if an MNI (x,y,z) coordinate matches an NS data point.

        Returns
        -------
        ids: list (int)
            NS study IDs that have a data point corresponding to coord.
        """

        try:
            self._ns['id_dict']
        except KeyError:
            raise NameError("Study ID dictionary not initialized; see/call NsabaBase.ns_load_id_dict()")

        ids = []
        if len(coord) == 3 and not isinstance(coord, str):
            for i, coords in self._ns['id_dict'].items():
                for this_coordinate in coords:
                    if this_coordinate == tuple(coord):
                        if i not in ids:
                            if self.is_id(i):
                                ids.append(i)
            return ids
        else:
            raise ValueError("Argument form improper; check function documentation.")

    # KNN
    @not_operational
    def _coord_to_ge(self, coord, entrez_ids, search_radii=3, k=20):
        """
        Returns weighted ABA gene expression statistic for some MNI coordinate based
        on a list of passed Entrez IDs.

        Parameters
        ----------
        coord: tuple-like (3)
            Reference MNI coordinate.
        entrez_ids: list-like
            Entrez IDs for gene expressions to be estimated around 'coord'.
        search_radii: numeric, optional
            Max search radii for KNN gene expression estimation technique.
        k: int, optional
            k parameter for KNN estimation.

        Returns
        -------
        ge_for_coord: list [N]
            Estimated gene expression at 'coord' for genes specified by 'entrez_ids'.
            Where N is the size of 'entrez_ids' (number of genes for expression
            to be estimated).

        """

        ge_for_coord = []
        if len(coord) == 3 and not isinstance(coord, str):
            for entrez_id in entrez_ids:
                coord_inds, radii = self._knn_search(coord, self._aba['mni_coords'], search_radii, k)
                if len(coord_inds) == 0:
                    # print "No ABA coordinates are within search radius of specified coordinate"
                    break
                weight = self.ns_weight_f(radii)
                local_ge = self.ge[entrez_id][coord_inds]
                weighted_ge_mean = np.sum(local_ge*weight)/np.sum(weight)
                ge_for_coord.append(weighted_ge_mean)
        else:
            raise ValueError("MNI coordinate in improper form; must be 3-tuple-like")
        return ge_for_coord

    # KNN
    @not_operational
    def coords_to_ge(self, coords, entrez_ids, search_radii=3, k=20):
        """
        Returns weighted ABA gene expression statistic for a list MNI coordinate based
        on a list of passed Entrez IDs.

        Parameters
        ----------
        coords: list-like
            List of reference MNI coordinates for gene expression to be estimated at.
        entrez_ids: list-like
            Entrez IDs for gene expressions to be estimated around 'coord'.
        search_radii: numeric, optional
            Max search radii for KNN gene expression estimation technique.
        k: int, optional
            k parameter for KNN estimation.

        Returns
        -------
        np.array(ge_coords): numpy.array [M x N]
            Returns as a matrix of M coordinates by N estimated gene expression
            coefficients. See _coord_to_ge() documentation for more information.

        """
        self._check_entrez_struct(entrez_ids)
        if not all([key in self.ge for key in entrez_ids]):
            raise ValueError("One or more Entrez IDs not loaded into ge.")

        ge_for_coords = []
        for coord in coords:
            ge_for_coord = self._coord_to_ge(coord, entrez_ids, search_radii, k)
            if ge_for_coords > 0:
                ge_for_coords.append(ge_for_coord)

        return np.array(ge_for_coords)

    def id_to_ns_act(self, study_id):
        """
        Returns activations for all terms for a given study.

        Parameters
        ----------
        study_id: int
            int representing a paper/study in the NS framework.

        Return
        ------
        term_vector_off_by_1[1:]: numpy.array [1 x 3406]
            Vector of term activations for all terms for a specified NS study.
        """

        if self.is_id(study_id):
            term_vector_off_by_1 = np.squeeze(self._ns['features_df'].loc[
                                   self._ns['features_df'].pmid == study_id].as_matrix())
            # Shifting to remove ID index from vector
            return term_vector_off_by_1[1:]
        else:
            raise ValueError("Invalid NS study ID; check 'study_id' parameter")

    def coord_to_ns_act(self, coord, return_type='list'):
        """
        -- LEGACY --
        Used to support visualize_ns_old(); itself unsupported.

        Returns list of terms activations for a MNI coordinate
        for all NS terms.

        Parameters
        ----------
        coord: tuple-like (3)
            Reference MNI coordinate.

        return_type: str
            OPTIONS:
                'dict': See Returns.
                'list': Returns list of activations for each term.
                OTHER: Raises ValueError.

        Returns
        -------
        terms: dict OR list
            A dictionary with (term: activation) key pairs for
            the specified MNI coordinate.
        """
        ids = self.coord_to_ids(coord)
        if len(ids) == 1:
            terms = self.id_to_ns_act(ids)
        elif len(ids) > 1:
            temp = []
            for multiple_id in ids:
                temp.append(self.id_to_ns_act(multiple_id))
                terms = np.mean(temp, 0)
        else:
            return []

        # [1:] to remove 'PMID' column header
        if return_type == 'dict':
            return {term: act for term, act in zip(self._ns['features_df'].columns[1:], terms)}
        elif return_type == 'list':
            return terms
        else:
            raise ValueError("Invalid return_type argument; use 'list' or 'dict'.")

    @not_operational
    def term_to_id_coords(self, term, no_ids=3):
        """
        Returns coordinates associated with studies that have the
        greatest term activation.

        Parameters
        ----------
        term : string
            NS term of interest

        no_ids : numeric, optional
            Number of studies to return coordinates for.

        Returns
        -------
        coords : list [ PMID_coord_pair (int, list [tuple(1x3)] ) ]
            Returns a list of len(no_ids) lists containing a namedtuple:
            "PMID_coord_pair". PMID_coord_pair contains two arguments: PMID
            and a list of coordinates (in tuple form) for that study.
        """
        id_coord_pair = collections.namedtuple("PMID_coord_pair", "pmid coords")
        if self.is_term(term):
            try:
                self._ns['id_dict'][24379394]
            except KeyError:
                self.ns_load_id_dict()
            heat = self._ns['features_df'][term]
            sorted_heat_vals = sorted(enumerate(heat), key=lambda x: x[1], reverse=True)[0:no_ids]
            inds = zip(*sorted_heat_vals)[0]
            pmids = [self._ns['features_df']['pmid'].ix[ind] for ind in inds]
            coords = []
            for pmid in pmids:
                if self.is_id(pmid):
                    coords.append(id_coord_pair(pmid, self._ns['id_dict'][pmid]))
            return coords
        else:
            raise ValueError("No previous estimation found for '%s'." % term)

    def _term_to_coords(self, term, thresh=0):
        """
        Finds coordinates associated with a given term above
        a specified threshold.

        Parameters
        ----------
        term : string
            NS term of interest

        thresh : numeric
            NS activation threshold.
            (activations < threshold: are discarded)

        Returns
        -------
        (spatial.KDTree, pandas.DataFrame)

        """
        term_ids_act = self._ns['features_df'].loc[self._ns['features_df'][term] >= thresh, ['pmid', term]]
        term_ids = term_ids_act['pmid'].tolist()
        term_coords = self._ns['database_df'].loc[self._ns['database_df']['id'].isin(term_ids)]
        try:
            ns_coord_tree = spatial.KDTree(term_coords.loc[:, 'x':'z'].as_matrix().astype(float))
        except ValueError:
            raise ValueError("No studies with term: '%s' and threshold: %.2f found" % (term, thresh))
        else:
            term_ids_act.rename(columns={'pmid': 'id'}, inplace=True)
            return ns_coord_tree, term_coords.merge(term_ids_act)

    def est_ns_act(self, term, coords=None, **kwargs):
        """Generates NS activation vector about ABA MNI coordinates; timed at 26.1 s"""
        if not self.is_term(term):
            raise ValueError("'%s' is not a registered term." % term)

        self.term[term] = {}

        if coords is None:
            coords = self._aba['mni_coords']
            self.term[term]['coord_type'] = 'ABA MNI'
        else:
            self.term[term]['coords'] = coords
            if 'coord_type' in kwargs:
                self.term[term]['coord_type'] = kwargs['coord_type']

        ns_coord_tree, ns_coord_act_df = self._term_to_coords(term, 0)

        if 'knn_args' in kwargs:
            self.term[term]['classifer'] = KNeighborsRegressor(**kwargs['knn_args'])
        else:
            self.term[term]['classifer'] = KNeighborsRegressor()

        X = ns_coord_tree.data
        y = ns_coord_act_df[term].as_matrix()

        self.term[term]['classifer'].fit(X,y)

        self.term[term]['act'] = self.term[term]['classifer'].predict(coords.data)

    def make_ge_ns_mat(self, ns_term, entrez_ids):
        self._check_entrez_struct(entrez_ids)

        if ns_term in self.term and all([key in self.ge for key in entrez_ids]):
            ge_ns_mat = []
            act_vec_len = len(self.term[ns_term]['act'])
            for entrez_id in entrez_ids:
                if len(self.ge[entrez_id]['GE']) == act_vec_len:
                    ge_ns_mat.append(self.ge[entrez_id]['GE'])
                else:
                    raise ValueError("Size mismatch between ge[%d] and term[%s]['act']"
                                     % (entrez_id, ns_term))
            ge_ns_mat.append(self.term[ns_term]['act'])
            return np.vstack(ge_ns_mat).T
        else:
            raise ValueError("Either term['%s'] or one or more Entrez ID keys does not exist in ge; "
                             "please check arguments" % ns_term)

