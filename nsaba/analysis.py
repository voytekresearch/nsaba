"""
analysis.py:
Statistical testing and analysis tools for Nsaba.
Author: Simon Haxby & Torben Noto
"""

from nsaba import Nsaba
from nsabatools import preprint, not_operational
from geneinfo import load_gene_file, get_local_gene_info, scrape_nih_gene_info

import random
import collections
import csv
import os
import warnings

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.cluster import KMeans
from sklearn import mixture
from sklearn.gaussian_process import GaussianProcess

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import seaborn as sns
    import matplotlib.pyplot as plt



def cohen_d(x1, x2, n1, n2):
    """
    Computes Cohen's d; an effect size statistic:
    Overview: https://en.wikipedia.org/wiki/Effect_size#Cohen.27s_d
    """
    return (np.mean(x2) - np.mean(x1)) / np.sqrt(((n1-1)*np.var(x1) + (n2-1)*np.var(x2)) / (n1 + n2-2))  #cohen's d is now positive


class NsabaAnalysis(object):
    """
    Provides a suite of statistical analysis and informatic tools for exploring
    gene to term relationships.
    """
    def __init__(self, nsaba_obj):
        if type(nsaba_obj) == Nsaba:
            self.no = nsaba_obj
        else:
            raise ValueError("NsabaAnalysis() parameter not a Nsaba instance")

        self.gene_rec = collections.namedtuple("gene_rec", "entrez cohen_d p_value")
        self.term_rec = collections.namedtuple("term_rec", "term cohen_d p_value")
        self.gene_rec_spearman = collections.namedtuple("gene_rec", "term r_value")
        self.term_rec_spearman = collections.namedtuple("term_rec", "entrez r_value")
        self.gene_rec_pearson = collections.namedtuple("gene_rec", "term r_value p_value")
        self.term_rec_pearson = collections.namedtuple("term_rec", "entrez r_value p_value")
        self.default_quant = 85
        print "To use inline plotting functionality in Jupyter, '%matplotlib inline' must be enabled"

    def _mask_check(self, mask, method):
        """
        Checks mask contains only zeros or ones.

        Parameters
        ----------
        mask : np.array([ ABA samples x 1 ] )
            Bool mask for splitting ABA samples into control and
            functional groups.

        method : str
            See _split_mask().
        """
        if sum(mask) == 0 or sum(mask) == len(mask):
            raise ValueError("No coordinates assigned to functional group; try another "
                             "splitting method other than '%s.'" % method)

    def _split_mask(self, term_vec, method='quant', **kwargs):
        """
        Vector of term activations is split in to two groups; for
        t-testing between them.
        This method constructs splitting mask for t-test via specified method.

        Parameters
        ----------
        term_vec : np.array()

        method : str, optional
            OPTIONS:
                'quant' : Quantile based split; i.e: default is 85/15 to control/functional.
                'kmeans' : k-means, k=2.
                'mog' : Mixture model (Gaussian), with 2 mixing components.

            Specifies splitting method.

        **kwargs : dict
            OPTIONS:
                'quant' : Quantile split to be made with 'quant' method.
                            i.e: 85 would split 85% of coordinates in control
                                and 15% of coordinates in activation.
        Returns
        -------
        mask : np.array([ ABA samples x 1 ] )
            Bool mask for splitting ABA samples into control and
            functional groups.

        """
        if method == 'quant':
            if 'quant' not in kwargs:
                kwargs['quant'] = self.default_quant
            thres = np.percentile(term_vec, kwargs['quant'])
            mask = np.array([True if coeff > thres else False for coeff in term_vec])
            self._mask_check(mask, method)
            return mask
        elif method == 'kmeans':
            X = [[x] for x in term_vec]
            init_clusters = np.array([[min(term_vec)], [max(term_vec)]])
            kmn = KMeans(n_clusters=2, init=init_clusters, n_init=1)
            mask = kmn.fit_predict(X).astype(bool)
            self._mask_check(mask, method)
            return mask
        elif method == 'mog':
            X = [[x] for x in term_vec]
            gm = mixture.GMM(n_components=2)
            gm.fit(X)
            if stats.mode(gm.predict(X))[0][0] == 0:
                mask = np.array([True if gm.predict([[coeff]]) == 1 else False for coeff in term_vec])
            else:
                mask = np.array([False if gm.predict([[coeff]]) == 1 else True for coeff in term_vec])
            self._mask_check(mask, method)
            return mask
        else:
            raise ValueError("method parameter '%s' not recognized" % method)

    def _split_groups(self, ge_vec, mask):
        """
        Splits gene expression samples into control and functional network
        based on mask generated by self._split_mask().
        """
        functional_grp = ge_vec[mask]
        diff = set(ge_vec) - set(functional_grp)
        control_grp = np.array(list(diff))
        return control_grp, functional_grp

    def term_ge_ttest(self, term, gene, split_method='quant', graphops='density', **kwargs):
        """
        Performs gene expression t-test between coordinates in control and functional
        network based on term activation, and generates associated plot.
        Intuitively, coordinates with high term activation are likely to be local to a
        functional network associated with that term. Split into control and functional
        network groups is dependent on term activation and splitting method (see _split_mask()).

        Parameters
        ----------
            term : str
                NS term of functional network of interest.

            gene : int
                Entrez ID of gene to be t-tested.

            split_method : str, optional

            graphops : str, optional
                OPTIONS:
                    'density' : Distribution density plot of control and functional network expression.
                    'box' : Box plot of control and functional network expression.
                    'violin' : Violin plot of control and functional network expression.

            kwargs : dict
                OPTIONS:
                    'log' : term activation values are log-spaced.

                PASSED:
                    _split_mask()

        """
        analymat = self.no.matrix_builder([term], [gene])

        non_nans = []
        for ind, row in enumerate(analymat):
            if not any(np.isnan(row)):
                non_nans.append(ind)

        analymat = analymat[non_nans]

        # Splitting groups

        mask = self._split_mask(analymat[:, 1], method=split_method, **kwargs)
        cont_grp, funct_grp = self._split_groups(analymat[:, 0], mask)
        if 'log' in kwargs:
            if kwargs['log']:
                funct_grp = np.log(funct_grp)
                cont_grp = np.log(cont_grp)

        # T-Test
        print "t-value: %.4f \np-value: %.3E" % stats.ttest_ind(cont_grp, funct_grp)
        print "Effect size: %.4f" % cohen_d(cont_grp, funct_grp, len(cont_grp), len(funct_grp))
        print "Control/Functional Split: %d/%d\n" % (len(mask)-sum(mask), sum(mask))
        # Histogram/KDE Plots
        if graphops == 'density':
            ax = plt.axes()
            ax.set_title('Gene Expression Distributions')
            ax.set_xlabel(str(gene))
            ax.set_ylabel('density')
            sns.distplot(funct_grp, ax=ax, label=term)
            sns.distplot(cont_grp, label='null')
            plt.legend()
        elif graphops == 'box':
            ax = plt.axes()
            ax.boxplot([cont_grp, funct_grp])
            ax.set_xticks([1, 2])
            ax.set_xticklabels(["Low '"+term+"'", "High '"+term+"'"])
            ax.set_ylabel('Gene Expression')
        elif graphops == 'violin':
            ax = plt.axes()
            ax.violinplot([cont_grp, funct_grp])
            ax.set_xticks([1, 2])
            ax.set_xticklabels(["Low '"+term+"'", "High '"+term+"'"])
            ax.set_ylabel('Gene Expression')
            ax.plot(np.ones(len(cont_grp)), cont_grp, 'b.')
            ax.plot(1, np.mean(cont_grp), 'bs')
            ax.plot(2*np.ones(len(funct_grp)), funct_grp, 'g.')
            ax.plot(2, np.mean(funct_grp), 'gs')
        else:
            raise ValueError("graphops parameter '%s' not recognized" % graphops)

    def ge_term_ttest(self, term, gene, split_method='quant', graphops='density', **kwargs):
            """
            Performs a t-test on the association with a given term between regions with
            low expression and high expression of a given gene.

            Parameters
            ----------
                term : str
                    NS term to be t-tested.

                gene : int
                    Entrez ID to index regions for t-test.

                split_method : str, optional

                graphops : str, optional
                    OPTIONS:
                        'density' : Distribution density plot of control and functional network expression.
                        'box' : Box plot of control and functional network expression.
                        'violin' : Violin plot of control and functional network expression.

                kwargs : dict
                    OPTIONS:
                        'log' : term activation values are log-spaced.

                    PASSED:
                        _split_mask()

            """
            analymat = self.no.matrix_builder([term], [gene])

            non_nans = []
            for ind, row in enumerate(analymat):
                if not any(np.isnan(row)):
                    non_nans.append(ind)

            analymat = analymat[non_nans]

            # Splitting groups

            mask = self._split_mask(analymat[:, 0], method=split_method, **kwargs) #1 to 0 switch
            cont_grp, funct_grp = self._split_groups(analymat[:, 1], mask) #0 to 1 switch
            if 'log' in kwargs:
                if kwargs['log']:
                    funct_grp = np.log(funct_grp)
                    cont_grp = np.log(cont_grp)

            # T-Test
            print "t-value: %.4f \np-value: %.3E" % stats.ttest_ind(cont_grp, funct_grp)
            print "Effect size: %.4f" % cohen_d(cont_grp, funct_grp, len(cont_grp), len(funct_grp))
            print "Control/Functional Split: %d/%d\n" % (len(mask)-sum(mask), sum(mask))
            # Histogram/KDE Plots
            if graphops == 'density':
                ax = plt.axes()
                ax.set_title('Term Distributions')
                ax.set_xlabel(str(term))
                ax.set_ylabel('density')
                sns.distplot(funct_grp, ax=ax, label=term)
                sns.distplot(cont_grp, label='null')
                plt.legend()
            elif graphops == 'box':
                ax = plt.axes()
                ax.boxplot([cont_grp, funct_grp])
                ax.set_xticks([1, 2])
                ax.set_xticklabels(["Low "+str(gene), "High "+str(gene)])
                ax.set_ylabel('Term Association')
            elif graphops == 'violin':
                ax = plt.axes()
                ax.violinplot([cont_grp, funct_grp])
                ax.set_xticks([1, 2])
                ax.set_xticklabels(["Low "+str(gene), "High "+str(gene)])
                ax.set_ylabel('Term Association')
                ax.plot(np.ones(len(cont_grp)), cont_grp, 'b.')
                ax.plot(1, np.mean(cont_grp), 'bs')
                ax.plot(2*np.ones(len(funct_grp)), funct_grp, 'g.')
                ax.plot(2, np.mean(funct_grp), 'gs')
            else:
                raise ValueError("graphops parameter '%s' not recognized" % graphops)



    @preprint('This may take a couple of minutes ...')
    def term_ge_spearman_rho(self, term, sample_num=None, *kwargs):
        """"
        Calculates Spearman's Rho for each across all genes for a given term.

        Parameters
        ----------
        term : str
            NS term of interest

        Returns
        -------
        list : [ len(ge.keys()) ]
            Returns a list of Spearman coefficients corresponding to the Spearman
            correlation between term activation for 'term' and gene expression for
            all genes loaded into the base Nsaba object.
        """
        if term not in self.no.term:
            raise ValueError("Term activation not generated for '%s" % term)

        if sample_num == None:
            sample_num = len(self.no.ge.keys())

        sam_ids = random.sample(self.no.ge.keys(), sample_num)

        spearman_metrics = {'term': term, "gene_sample_size": sample_num}

        # Use parameters nih_only=True and use gi_csv_path='..' to specify path to 'gene_info.csv'
        if 'nih_only' in kwargs:
            if kwargs['nih_only']:
                if 'gi_csv_path' in kwargs:
                    gi_path = kwargs['gi_csv_path']
                    df = load_gene_file(gi_path)
                else:
                    df = load_gene_file()
                nih_ids = df['Entrez'].as_matrix()
                sam_ids = [entrez_id for entrez_id in nih_ids if entrez_id in nih_ids]
                print "Using NIH described genes only; Entrez ID sample size now %d" % (len(sam_ids))

        ge_mat = self.no.matrix_builder([term], sam_ids)

        non_nans = []
        for ind, row in enumerate(ge_mat):
            if not any(np.isnan(row)):
                non_nans.append(ind)

        ge_mat = ge_mat[non_nans]
        # Calculates Spearman's Rho for each across all genes for a given term
        spearman_vals = [stats.spearmanr(ge_mat[:, ge_mat.shape[1]-1], ge_mat[:, r])[0] for r in xrange(len(sam_ids))]


        term_stats = []
        for entrez in xrange(len(sam_ids)):
            term_stats.append(self.term_rec_spearman(sam_ids[entrez], spearman_vals[entrez]))

        # Sort effect sizes from greatest to smallest in magnitude
        term_stats.sort(key=lambda rec: rec.r_value)
        spearman_metrics['results'] = term_stats

        return spearman_metrics

    @preprint('This may take a couple of minutes ...')
    def ge_term_spearman_rho(self, entrez):
        """"
        Calculates Spearman's Rho for each across all genes for a given term.

        Parameters
        ----------
        entrez : str
            entrez id of gene of interest

        Returns
        -------
        list : [ len(ge.keys()) ]
            Returns a list of Spearman coefficients corresponding to the Spearman
            correlation between term activation for 'term' and gene expression for
            all genes loaded into the base Nsaba object.
        """
        if entrez not in self.no.ge:
            raise ValueError("Gene estimation not generated for '%s" % entrez)
        ge_mat = self.no.matrix_builder(self.no.term.keys(), [entrez])

        non_nans = []
        for ind, row in enumerate(ge_mat):
            if not any(np.isnan(row)):
                non_nans.append(ind)

        ge_mat = ge_mat[non_nans]

        # Calculates Spearman's Rho on all terms for a given gene
        return [stats.spearmanr(ge_mat[:, 0], ge_mat[:, r+1])[0]
                for r in xrange(len(self.no.term.keys()))]

    @preprint('This may take a couple of minutes ...')
    def term_ge_pearson(self, term, sample_num=None, *kwargs):
        """"
        Calculates Pearson's r and p for each across all genes for a given term.

        Parameters
        ----------
        term : str
            NS term of interest

        Returns
        -------
        list : [ len(ge.keys()) ]
            Returns a list of Pearson coefficients corresponding to the Pearson
            correlation between term activation for 'term' and gene expression for
            all genes loaded into the base Nsaba object.
        """
        if term not in self.no.term:
            raise ValueError("Term activation not generated for '%s" % term)

        if sample_num == None:
            sample_num = len(self.no.ge.keys())

        sam_ids = random.sample(self.no.ge.keys(), sample_num)
        pearson_metrics = {'term': term, "gene_sample_size": sample_num}

        # Use parameters nih_only=True and use gi_csv_path='..' to specify path to 'gene_info.csv'
        if 'nih_only' in kwargs:
            if kwargs['nih_only']:
                if 'gi_csv_path' in kwargs:
                    gi_path = kwargs['gi_csv_path']
                    df = load_gene_file(gi_path)
                else:
                    df = load_gene_file()
                nih_ids = df['Entrez'].as_matrix()
                sam_ids = [entrez_id for entrez_id in nih_ids if entrez_id in nih_ids]
                print "Using NIH described genes only; Entrez ID sample size now %d" % (len(sam_ids))

        ge_mat = self.no.matrix_builder([term], sam_ids)

        non_nans = []
        for ind, row in enumerate(ge_mat):
            if not any(np.isnan(row)):
                non_nans.append(ind)

        ge_mat = ge_mat[non_nans]

        # Calculates Pearson's Rho for each across all genes for a given term
        term_stats = []
        for entrez in xrange(len(sam_ids)):
            r, p = stats.pearsonr(ge_mat[:, ge_mat.shape[1]-1], ge_mat[:, entrez])

            term_stats.append(self.term_rec_pearson(sam_ids[entrez], r, p))

        # Sort effect sizes from greatest to smallest in magnitude
        term_stats.sort(key=lambda rec: rec.r_value)
        pearson_metrics['results'] = term_stats

        return pearson_metrics

    @not_operational
    def t_test_custom_ge(self, coords1, coords2, gene, quant=None, log=False, graphops='density'):
        """ T-Test of gene expression between term and non-term coordinates"""
        if not quant:
            quant = 85
        ge1 = self.no.coords_to_ge(coords1, [gene])
        ge2 = self.no.coords_to_ge(coords2, [gene])

        if log:
            ge1 = np.log(ge1)
            ge2 = np.log(ge2)

        # T-Test
        print "t-value: %.4f \np-value: %.3E" % stats.ttest_ind(ge1, ge2)
        print "Effect size: %.4f \n" % cohen_d(ge1, ge2, len(ge1), len(ge2))
        # Histogram/KDE Plots
        if graphops == 'density':
            ax = plt.axes()
            ax.set_title('Gene Expression Distributions')
            ax.set_xlabel('gene expression')
            ax.set_ylabel('density')
            sns.distplot(ge1, ax=ax, label='coordinate list 1')
            sns.distplot(ge2, label='coordinate list 2')
            plt.legend()
        elif graphops == 'box':
            ax = plt.axes()
            ax.boxplot([ge1, ge2])
            ax.set_xticks([1, 2])
            ax.set_xticklabels(['coordinate list 1', 'coordinate list 2'])
        elif graphops == 'violin':
            ax = plt.axes()
            ax.violinplot([ge1, ge2])
            ax.set_xticks([1, 2])
            ax.set_xticklabels(['coordinate list 1', 'coordinate list 2'])
            ax.set_ylabel('Gene Expression: %d' % gene)
            ax.plot(np.ones(len(ge1)), ge1, 'b.')
            ax.plot(1, np.mean(ge1), 'bs')
            ax.plot(2*np.ones(len(ge2)), ge2, 'g.')
            ax.plot(2, np.mean(ge2), 'gs')
        else:
            raise ValueError("graphops parameter '%s' not recognized" % graphops)

    @not_operational
    def t_test_custom_term(self, coords1, coords2, term, quant=None, log=False, graphops='density'):
        """ T-Test of gene expression between term and non-term coordinates"""
        if not quant:
            quant = 85
        term1 = self.no.coords_to_ns_act(coords1, term)
        term2 = self.no.coords_to_ns_act(coords2, term)

        if log:
            term1 = np.log(term1)
            term2 = np.log(term2)

        # T-Test
        print "t-value: %.4f \np-value: %.3E" % stats.ttest_ind(term1, term2)
        print "Effect size: %.4f \n" % cohen_d(term1, term2, len(term1), len(term2))
        # Histogram/KDE Plots
        if graphops == 'density':
            ax = plt.axes()
            ax.set_title(term + ' Distributions')
            ax.set_xlabel(term)
            ax.set_ylabel('density')
            sns.distplot(term1, ax=ax, label='coordinate list 1')
            sns.distplot(term2, label='coordinate list 2')
            plt.legend()
        elif graphops == 'box':
            ax = plt.axes()
            ax.boxplot([term1, term2])
            ax.set_xticks([1, 2])
            ax.set_xticklabels(['coordinate list 1', 'coordinate list 2'])
        elif graphops == 'violin':
            ax = plt.axes()
            ax.violinplot([term1, term2])
            ax.set_xticks([1, 2])
            ax.set_xticklabels(['coordinate list 1', 'coordinate list 2'])
            ax.set_ylabel(term)
            ax.plot(np.ones(len(term1)), term1, 'b.')
            ax.plot(1, np.mean(term1), 'bs')
            ax.plot(2*np.ones(len(term2)), term2, 'g.')
            ax.plot(2, np.mean(term2), 'gs')
        else:
            raise ValueError("graphops parameter '%s' not recognized" % graphops)

    @preprint('This may take a couple of minutes ...')
    def term_ge_ttest_multi(self, term, split_method='quant', sample_num=None, **kwargs):
        """
        Performs t-test equivalent to term_ge_t_test() across a subsample of Entrez IDs
        loaded into self.no. Default is to use all loaded Entrez ID genes.

        Parameters
        ----------
        term : str

        split_method : str, optional
            OPTIONS:
                'quant' : Quantile based split; i.e: default is 85/15 to control/functional.
                'kmeans' : k-means, k=2.
                'mog' : Mixture model (Gaussian), with 2 mixing components.

        sample_num : int, optional

        kwargs : dict
            OPTIONS:
                'genes_of_interest' : list
                'nih_only' : bool
                'gi_csv_path' : str
            PASSED:
                _split_mask()


        Returns
        -------
        ttest_metrics : dict
            Contains meta-information and results of multi-gene t-tests.
            See below for construction.

        """
        # Setting contingency variables
        if term not in self.no.term:
            raise ValueError("Term activation not generated for '%s" % term)
        if sample_num == None:
            sample_num = len(self.no.ge.keys())
        elif sample_num <= 0:
            raise ValueError("'sample_num' parameter must be greater than 0")
        if 'gene_of_interest' not in kwargs:
            kwargs['genes_of_interest'] = []

        if len(self.no.ge) < sample_num:
            raise ValueError("Sample number exceeds stored number of Entrez IDs")

        # Use parameters nih_only=True and use gi_csv_path='..' to specify path to 'gene_info.csv'
        sam_ids = random.sample(self.no.ge.keys(), sample_num)
        if 'nih_only' in kwargs:
            if kwargs['nih_only']:
                if 'gi_csv_path' in kwargs:
                    gi_path = kwargs['gi_csv_path']
                    df = load_gene_file(gi_path)
                else:
                    df = load_gene_file()
                nih_ids = df['Entrez'].as_matrix()
                sam_ids = [entrez_id for entrez_id in nih_ids if entrez_id in nih_ids]
                print "Using NIH described genes only; Entrez ID sample size now %d" % (len(sam_ids))

        # Fetching GE/NS activation matrix
        matrix = self.no.matrix_builder([term], sam_ids)

        non_nans = []
        for ind, row in enumerate(matrix):
            if not any(np.isnan(row)):
                non_nans.append(ind)

        matrix = matrix[non_nans]

        ge_mat = matrix.T[:-1]
        term_act_vector = matrix.T[-1:][0]

        mask = self._split_mask(term_act_vector, method=split_method, **kwargs)

        # Prepping ttest_metrics results dictionary
        ttest_metrics = {'term': term, 'split_method': split_method, "gene_sample_size": sample_num}
        if split_method == 'quant':
            if 'quant' in kwargs:
                ttest_metrics['quant'] = kwargs['quant']
            else:
                ttest_metrics['quant'] = self.default_quant

        gene_stats = []
        for eid, ge in zip(sam_ids, ge_mat):
            # Split coordinates in to term and non-term groups
            cont_grp, funct_grp = self._split_groups(ge, mask)
            test_stats = stats.ttest_ind(cont_grp, funct_grp)
            d = cohen_d(cont_grp, funct_grp, len(cont_grp), len(funct_grp))
            # One-sided T-Test
            gene_stats.append(self.gene_rec(eid, d, test_stats[1]))

            if eid in kwargs['genes_of_interest']:
                print 'Gene: ' + str(eid) + '  Effect Size: '+str(d)
        # Sort effect sizes from greatest to smallest in magnitude
        gene_stats.sort(key=lambda rec: rec.cohen_d)
        ttest_metrics['results'] = gene_stats
        return ttest_metrics

    @preprint('This may take a couple of minutes ...')
    def ge_term_ttest_multi(self, entrez, split_method='quant', sample_num=None, **kwargs):
        """
        Performs t-test equivalent to ge_term_t_test() across a subsample of term IDs
        loaded into self.no. Default is to use all loaded Entrez ID genes.

        Parameters
        ----------
        term : str

        split_method : str, optional
            OPTIONS:
                'quant' : Quantile based split; i.e: default is 85/15 to control/functional.
                'kmeans' : k-means, k=2.
                'mog' : Mixture model (Gaussian), with 2 mixing components.

        sample_num : int, optional

        kwargs : dict
            OPTIONS:
                'genes_of_interest' : list
                'nih_only' : bool
                'gi_csv_path' : str
            PASSED:
                _split_mask()


        Returns
        -------
        ttest_metrics : dict
            Contains meta-information and results of multi-gene t-tests.
            See below for construction.

        """
        # Setting contingency variables
        if entrez not in self.no.ge:
            raise ValueError("Gene estimation not generated for '%s" % gene)
        if sample_num == None:
            sample_num = len(self.no.term.keys())
        elif sample_num <= 0:
            raise ValueError("'sample_num' parameter must be greater than 0")
        if 'terms_of_interest' not in kwargs:
            kwargs['terms_of_interest'] = []

        if len(self.no.term) < sample_num:
            raise ValueError("Sample number exceeds stored number of terms")

        # Use parameters nih_only=True and use gi_csv_path='..' to specify path to 'gene_info.csv'
        sam_ids = random.sample(self.no.term.keys(), sample_num)

        # Fetching GE/NS activation matrix
        matrix = self.no.matrix_builder(sam_ids, [entrez])

        non_nans = []
        for ind, row in enumerate(matrix):
            if not any(np.isnan(row)):
                non_nans.append(ind)

        matrix = matrix[non_nans]

        term_mat = matrix.T[1:]
        gene_act_vector = matrix.T[0]

        mask = self._split_mask(gene_act_vector, method=split_method, **kwargs)
        # Prepping ttest_metrics results dictionary
        ttest_metrics = {'Entrez': entrez, 'split_method': split_method, "term_sample_size": sample_num}
        if split_method == 'quant':
            if 'quant' in kwargs:
                ttest_metrics['quant'] = kwargs['quant']
            else:
                ttest_metrics['quant'] = self.default_quant

        term_stats = []
        for tid, term in zip(sam_ids, term_mat):
            # Split coordinates in to gene and non-gene groups
            cont_grp, funct_grp = self._split_groups(term, mask)
            test_stats = stats.ttest_ind(cont_grp, funct_grp)
            d = cohen_d(cont_grp, funct_grp, len(cont_grp), len(funct_grp))
            # One-sided T-Test
            term_stats.append(self.term_rec(tid, d, test_stats[1]))

            if tid in kwargs['terms_of_interest']:
                print 'Term: ' + tid + '  Effect Size: '+str(d)
        # Sort effect sizes from greatest to smallest in magnitude
        term_stats.sort(key=lambda rec: rec.cohen_d)
        ttest_metrics['results'] = term_stats
        return ttest_metrics

    def fetch_gene_descriptions(self, metrics, coeff='cohen', nih_fetch_num=20, alpha=.05, **kwargs):
        """
        Uses ttest_metrics dictionary returned by term_ge_ttest_multi() or list of
        Spearman coefficents generated by term_ge_spearman_rho(); fetches gene information
        for genes with largest effect sizes.

        Parameters
        ----------
        coeff : str
            OPTIONS:
                'cohen' : Cohen's d
                'spearman' : Spearman's Rho

        metrics : dict OR list
            See term_ge_ttest_multi() 'Returns' or term_ge_spearman_rho() 'Returns'.
            NOTE: method functionality is contingent on 'coeff'.

        nih_fetch_num : int, optional
            Number of top results to fetch from NIH website or pre-generated NIH csv file.

        alpha : float
            t-test significance level; used in Bonferroni correction.

        kwargs : dict
            OPTIONS:
                'nih_dl' : bool, optional
                    Download gene name and description from NIH website.
                'csv_path' : str
                    Specifies path to gene_info.csv
                'verbose' : bool, optional
                    Print most information about most significant genes.

        Returns
        -------
        top_genes : list
            Returns a list of records with largest effect sizes from ttest_metrics with
            additional gene description information.

        """
        if 'verbose' not in kwargs:
            kwargs['verbose'] = True
        if 'nih_dl' not in kwargs:
            kwargs['nih_dl'] = False

        if kwargs['nih_dl'] == False:
            if 'csv_path' not in kwargs:
                raise ValueError("'csv_path' argument in **kwargs missing; provide argument or use 'nih_dl':True .")

        top_genes = []
        if coeff == 'cohen':
            # Checking if users mixed up coeff/metric parameters
            if isinstance(metrics, list):
                raise ValueError("list passed with coeff='cohen'; if you want to use Spearman's Rho use coeff='spearman'.")

            for rec in metrics['results'][:nih_fetch_num]:
                try:
                    if kwargs['nih_dl']:
                        gene_name, gene_description = gene_info(str(rec.entrez))
                    else:
                        gene_dat = get_local_gene_info(kwargs['csv_path'], [rec.entrez])
                        gene_name = gene_dat[0].name
                        gene_description = gene_dat[0].description
                    top_genes.append((rec.entrez, rec.cohen_d, rec.p_value, gene_name, gene_description))
                except IndexError:
                    continue

            if kwargs['verbose']:
                print "\nCorrected Bonferroni Alpha: %.3E\n\n" % (alpha / float(metrics['gene_sample_size']))
                for eid, coh_d, p_val, gene_i, descr in top_genes:
                    if len(descr) == 1:
                        print "%d (p = %.3E; d = %.3f): < No description found >\n\n" % (eid, p_val, coh_d)
                    else:
                        print "%d (p = %.3E; d = %.3f): %s\n\n" % (eid, p_val, coh_d, descr)
        elif coeff == 'spearman':
            # Checking if users mixed up coeff/metric parameters
            if isinstance(metrics, dict):
                raise ValueError("dict passed with coeff='spearman'; if you want to use Cohen's d use coeff='cohen'.")

            top_ids = np.argsort(metrics)[:nih_fetch_num]

            top_rs = [r for r in reversed(np.sort(metrics))][:nih_fetch_num]
            top_genes = []
            for x in xrange(len(top_ids)):
                try:
                    if kwargs['nih_dl']:
                        gene_name, gene_description = gene_info(str(top_ids[x]))
                    else:
                        gene_dat = get_gene_info(kwargs['csv_path'], [top_ids[x]])
                        gene_name = gene_dat[0].name
                        gene_description = gene_dat[0].description
                    top_genes.append((int(top_ids[x]), top_rs[x], gene_name, gene_description))
                except IndexError:
                    continue

            if kwargs['verbose']:
                print "\nCorrected Bonferroni Alpha: %.3E\n\n" % (alpha/float(len(self.no.ge.keys())))
                for eid, rho, gene_i, descr in top_genes:
                    if len(descr) == 1:
                        print "%d (r = %.3f): < No description found >\n\n" % (eid, rho)
                    else:
                        print "%d (r = %.3f): %s\n %s\n\n" % (eid, rho, gene_i, descr)

        else:
            raise ValueError("Invalid parameter value for 'coeff'; use 'spearman' or 'cohen'.")

        return top_genes

    def p_val_distr(self, ttest_metrics):
        """Visualizing p-value distribution"""
        p_vals = [rec.p_value for rec in ttest_metrics['results']]
        sig = sum((p < .05/float(ttest_metrics['gene_sample_size']) for p in p_vals))
        print "Percent Significant (Bonferroni Correction; alpha = .05): %.3f %%" % \
              (100*sig/float(ttest_metrics['gene_sample_size']))
        plt.hist(p_vals, bins=75)
        plt.title('P-value Distribution')
        plt.xlabel('p-values')
        plt.ylabel('frequency')

    def term_rho_distr_ge(self, r_values, genes_of_interest=None):
        """Visualizing effect-size distribution (Spearman's Rho)"""

        if genes_of_interest is None:
            genes_of_interest = []
        ax = plt.axes()
        ax.hist(r_values, bins=75)
        ax.set_title("Effect Size Distribution (Spearman's Rho)")
        ax.set_xlabel("effect sizes")
        ax.set_ylabel('frequency')

        if genes_of_interest != []:
            offsetter = 500/len(genes_of_interest)
            for r in xrange(len(r_values)):
                if self.no.ge.keys()[r] in genes_of_interest:
                    plt.plot([r_values[r], r_values[r]], [0, offsetter])
                    plt.annotate('Gene:'+str(self.no.ge.keys()[r])+' rho='+str(r_values[r]),
                                 [r_values[r], offsetter])
                    if genes_of_interest != []:
                        offsetter += 500/len(genes_of_interest)

    def gene_rho_distr(self, r_values, terms_of_interest=None):
        """Visualizing effect-size distribution (Spearman's Rho)"""

        if terms_of_interest is None:
            terms_of_interest = []
        ax = plt.axes()
        ax.hist(r_values, bins=75)
        ax.set_title("Effect Size Distribution (Spearman's Rho)")
        ax.set_xlabel("effect sizes")
        ax.set_ylabel('frequency')

        if terms_of_interest != []:
            offsetter = 500/len(terms_of_interest)
            for r in xrange(len(r_values)):
                if self.no.term.keys()[r] in terms_of_interest:
                    plt.plot([r_values[r], r_values[r]], [0, offsetter])
                    plt.annotate(self.no.term.keys()[r]+' rho='+str(r_values[r]), [r_values[r], offsetter])
                    if terms_of_interest != []:
                        offsetter += 500/len(terms_of_interest)

    def gene_cohen_d_distr(self, ttest_metrics, genes_of_interest=None, return_fig=False):
        """Visualizing effect-size distribution (Cohen's d)"""

        if genes_of_interest is None:
            genes_of_interest = []

        d_vals = [rec.cohen_d for rec in ttest_metrics['results']]
        ax = plt.axes()
        ax.hist(d_vals, bins=75)
        ax.set_title("Effect Size Distribution (Cohen's d) of " + ttest_metrics['term'])
        ax.set_xlabel('effect sizes')
        ax.set_ylabel('frequency')

        if genes_of_interest != []:
            offsetter = 450/len(genes_of_interest)
        for rec in ttest_metrics['results']:
            if int(rec.entrez) in genes_of_interest:
                plt.plot([rec.cohen_d, rec.cohen_d], [0, offsetter])
                plt.annotate('Gene:'+str(int(rec.entrez))+' d='+str(rec.cohen_d),
                             [rec.cohen_d, offsetter])
                if genes_of_interest != []:
                    offsetter += 425/len(genes_of_interest)
        if return_fig:
            return ax

    def term_cohen_d_distr(self, ttest_metrics, terms_of_interest=None, return_fig=False):
        """Visualizing effect-size distribution (Cohen's d)"""

        if terms_of_interest is None:
            terms_of_interest = []

        d_vals = [rec.cohen_d for rec in ttest_metrics['results']]
        ax = plt.axes()
        ax.hist(d_vals, bins=75)
        ax.set_title("Effect Size Distribution (Cohen's d) of gene #" + str(ttest_metrics['Entrez']))
        ax.set_xlabel('effect sizes')
        ax.set_ylabel('frequency')

        if terms_of_interest != []:
            offsetter = 450/len(terms_of_interest)
        for rec in ttest_metrics['results']:
            if rec.term in terms_of_interest:
                plt.plot([rec.cohen_d, rec.cohen_d], [0, offsetter])
                plt.annotate(rec.term +' d='+str(rec.cohen_d),
                             [rec.cohen_d, offsetter])
                if terms_of_interest != []:
                    offsetter += 425/len(terms_of_interest)
        if return_fig:
            return ax

    def save_ttest_metrics(self, ttest_metrics, fname, no_genes=20):
        """
        Saves results of term_ge_ttest_multi() with gene descriptions in a csv file.

        Parameters
        ----------
        ttest_metrics : dict
            See term_ge_ttest_multi() 'Returns'.

        fname : str
            Name of csv. i.e: fname.csv

        no_genes : int
            Number of records to write to csv.
        """

        top_genes = self.fetch_gene_descriptions(ttest_metrics, nih_fetch_num=no_genes, printme=False)
        eids = [int(i[0]) for i in top_genes]
        myfig = self.effect_size_distr(ttest_metrics, genes_of_interest=eids[0:no_genes], return_fig=True)
        plt.savefig(fname+'.png')

        with open(fname+'.csv', 'wb') as csvfile:
            writer = csv.writer(csvfile)
            for i in top_genes:
                writer.writerow([i[0], i[3], i[1], i[2], i[4]])

    @not_operational
    def validate(self, term, genes, alt_genes=None, method='t_test', quant=85):

        real_gene_output = []
        random_gene_output = []

        if method == 't_test':
            ttest_metrics = self.t_test_multi(term, quant=quant)
            if alt_genes == None:
                alt_genes = random.randint(1, len(ttest_metrics['results']), len(genes))
            for i in ttest_metrics['results']:
                if int(i[0]) in genes:
                    real_gene_output.append(i[1])
                if int(i[0]) in alt_genes:
                    random_gene_output.append(i[1])
            #return stats.ttest_ind(real_gene_output, random_gene_output)[1]  # p value
            return real_gene_output, random_gene_output

        if method == 'pearson':
            alt_indices = random.randint(1, len(self.no.aba['probe_df']['entrez_id']), len(genes))
            alt_genes = self.no.aba['probe_df']['entrez_id'][alt_indices]
            for gene in genes:
                ge_ns_mat = self.no.matrix_builder([term], [gene])
                real_gene_output.append(np.corrcoef(ge_ns_mat[:, 0], np.log(ge_ns_mat[:, 1])))
            for gene in alt_genes:
                ge_ns_mat = self.no.matrix_builder([term], [gene])
                real_gene_output.append(np.corrcoef(ge_ns_mat[:, 0], np.log(ge_ns_mat[:, 1])))
            return real_gene_output, random_gene_output

        if method == 'spearman':
            alt_indices = random.randint(1, len(self.no.aba['probe_df']['entrez_id']), len(genes))
            alt_genes = self.no.aba['probe_df']['entrez_id'][alt_indices]
            for gene in genes:
                ge_ns_mat = self.no.matrix_builder([term], [gene])
                real_gene_output.append(stats.spearmanr(ge_ns_mat[:, 0], ge_ns_mat[:, 1]))
            for gene in alt_genes:
                ge_ns_mat = self.no.matrix_builder([term], [gene])
                real_gene_output.append(stats.spearmanr(ge_ns_mat[:, 0], ge_ns_mat[:, 1]))
            return real_gene_output, random_gene_output

        if method == 'regression':
            alt_indices = random.randint(1, len(self.no.aba['probe_df']['entrez_id']), len(genes))
            alt_genes = self.no.aba['probe_df']['entrez_id'][alt_indices]
            for gene in genes:
                ge_ns_mat = self.no.matrix_builder([term], [gene])
                X = np.vstack([ge_ns_mat[:, 0], np.ones(len(ge_ns_mat[:, 0]))]).T
                m, c = np.linalg.lstsq(X, np.log(ge_ns_mat[:, 1]))[0]
                real_gene_output.append(m, c)
            for gene in alt_genes:
                ge_ns_mat = self.no.matrix_builder([term], [gene])
                X = np.vstack([ge_ns_mat[:, 0], np.ones(len(ge_ns_mat[:, 0]))]).T
                m, c = np.linalg.lstsq(X, np.log(ge_ns_mat[:, 1]))[0]
                random_gene_output.append(m, c)
            return real_gene_output, random_gene_output

    def gp_ns_ge(self, term, gene, logy=False, logx=False, only_term=False, **kwargs):
        """
        Fits Gaussian process to gene expression predictor of term activation.
        Choice of covariance kernel and hyperparameters can be passed via **kwargs.

        Parameters
        ----------
        term : str
            NS term to be linearly correlated/predicted with a given gene's
            gene expression.

        gene : str
            Entrez ID of gene whose gene expression is to be used as
            a linear predictor of term activation.

        logy : bool
            log-space of y-axis.

        logx : bool
            log-space of x-axis.

        kwargs : dict
            See: http://scikit-learn.org/stable/modules/generated/sklearn.gaussian_process.GaussianProcess.html
            For valid hyperparameters and covariance function values.

        """
        if gene in self.no.ge:
            if term in self.no.term:
                ge_ns_mat = self.no.matrix_builder([term], [gene])
                if only_term:
                    if ge_ns_mat.shape[0] > 900:  # check this num later
                        print 'reinitializing ' + term + ' for hacky plotting method'
                        self.no.est_ns_act(term, radius=10)
                        ge_ns_mat = self.no.matrix_builder([term], gene)
                if logy:
                    plt.yscale('log')
                if logx:
                    plt.xscale('log')

                # remove nans
                valid1 = np.isfinite(ge_ns_mat[:, 0])
                valid2 = np.isfinite(ge_ns_mat[:, 1])
                valid_inds = np.logical_and(valid1, valid2)

                X = ge_ns_mat[:,0][valid_inds]
                y = ge_ns_mat[:,1][valid_inds]

                # GP parameters

                if 'corr' not in kwargs:
                    kwargs['corr']  = 'squared_exponential'
                if 'theta0' not in kwargs:
                    kwargs['theta0'] = 1e-3
                if 'thetaL' not in kwargs:
                    kwargs['thetaL'] = 1e-4
                if 'thetaU' not in kwargs:
                    kwargs['thetaU'] = 1e-3
                if 'random_start' not in kwargs:
                    kwargs['random_start'] = 100
                if 'nugget' in kwargs:
                    kwargs['nugget'] = kwargs['nugget'](y)
                else:
                    kwargs['nugget'] = np.var(y)

                # GP
                gp = GaussianProcess(**kwargs)

                gp.fit(X.reshape(-1, 1), y.reshape(-1, 1))

                x = np.linspace(np.min(X), np.max(X), 5000)
                y_pred, MSE = gp.predict(x.reshape(-1,1), eval_MSE=True)
                sigma = np.sqrt(MSE)

                plt.plot(x, y_pred, 'r-', label=u'$\hat{f}$', linewidth=2)
                plt.plot(X,y, 'bo', ms=4, alpha=.3)
                ci = np.squeeze(y_pred.reshape(1,-1))
                plt.fill_between(x, ci - 1.9600 * sigma, ci + 1.9600 * sigma, alpha=.5, facecolor='red');

                plt.xlabel(str(gene))
                plt.ylabel(term)


            else:
                raise ValueError("Term '%s' has not been initialized. Use get_ns_act('%s')" % term)
        else:
            raise ValueError("Gene %s has not been initialized. "
                             "Use self.no.get_aba_ge([%s])" % str(g))

    def gp_ns_ns(self, term1, term2, logy=False, logx=False, **kwargs):
        """
        Fits Gaussian process to term activation predictor of term activation.
        Choice of covariance kernel and hyperparameters can be passed via **kwargs.

        Parameters
        ----------
        term1 : str
            NS predictor term.

        term2 : str
            NS term to be linearly correlated/predicted with a predictor term's
            term activation.

        logy : bool
            log-space of y-axis.

        logx : bool
            log-space of x-axis.

        kwargs : dict
            See: http://scikit-learn.org/stable/modules/generated/sklearn.gaussian_process.GaussianProcess.html
            For valid hyperparameters and covariance function values.

        """
        if term1 in self.no.term:
            if term2 in self.no.term:
                if logy:
                    plt.yscale('log')
                if logx:
                    plt.xscale('log')
                # remove nans
                valid1 = np.isfinite(self.no.term[term1]['act'])
                valid2 = np.isfinite(self.no.term[term1]['act'])
                valid_inds = np.logical_and(valid1, valid2)

                X = self.no.term[term1]['act'][valid_inds]
                y = self.no.term[term2]['act'][valid_inds]

                # GP parameters

                if 'corr' not in kwargs:
                    kwargs['corr']  = 'squared_exponential'
                if 'theta0' not in kwargs:
                    kwargs['theta0'] = 1e-3
                if 'thetaL' not in kwargs:
                    kwargs['thetaL'] = 1e-4
                if 'thetaU' not in kwargs:
                    kwargs['thetaU'] = 1e-3
                if 'random_start' not in kwargs:
                    kwargs['random_start'] = 100
                if 'nugget' in kwargs:
                    kwargs['nugget'] = kwargs['nugget'](y)
                else:
                    kwargs['nugget'] = np.var(y)

                # GP
                gp = GaussianProcess(**kwargs)

                # Implementation hiccup; measure against ill-conditioned covariance matrices
                # See: https://github.com/scikit-learn/scikit-learn/issues/4916
                X_ = (X + np.random.normal(0,1e-8, len(X))).reshape(-1,1)
                gp.fit(X_, y)

                x = np.linspace(np.min(X), np.max(X), 5000)
                y_pred, MSE = gp.predict(x.reshape(-1,1), eval_MSE=True)
                sigma = np.sqrt(MSE)

                plt.plot(x, y_pred, 'r-', label=u'$\hat{f}$', linewidth=2)
                plt.plot(X,y, 'bo', ms=4, alpha=.3)
                ci = np.squeeze(y_pred.reshape(1,-1))
                plt.fill_between(x, ci - 1.9600 * sigma, ci + 1.9600 * sigma, alpha=.5, facecolor='red');
                cushion = (np.max(X)-np.min(X))/10.
                plt.xlim(np.min(X)-cushion, np.max(X)+cushion)
                plt.xlabel(term1)
                plt.ylabel(term2)

            else:
                raise ValueError("Term '%s' has not been initialized. Use get_ns_act('%s')" % term2)
        else:
            raise ValueError("Term '%s' has not been initialized. Use get_ns_act('%s')" % term1)

    def gp_ge_ge(self, gene1, gene2, **kwargs):
        """
        Fits Gaussian process to gene expression predictor of gene expression.
        Choice of covariance kernel and hyperparameters can be passed via **kwargs.

        Parameters
        ----------
        gene1 : str
            Entrez ID of ABA predictor gene.

        gene2 : str
            Entrez ID of ABA gene to be linearly correlated/predicted with a predictor gene's
            gene expression.

        logy : bool
            log-space of y-axis.

        logx : bool
            log-space of x-axis.

        kwargs : dict
            See: http://scikit-learn.org/stable/modules/generated/sklearn.gaussian_process.GaussianProcess.html
            For valid hyperparameters and covariance function values.

        """

        genes = [gene1, gene2]

        for gene in genes:
            if gene not in self.no.ge:
                raise ValueError("Gene %s has not been initialized. "
                                 "Use self.no.get_aba_ge([%s])" % str(gene))

        if len(self.no.ge[genes[0]]["mean"]['GE']) != len(self.no.ge[genes[1]]["mean"]['GE']):
            raise ValueError("'GE' size mismatched rerun Nsaba.estimate_aba_ge() again.")

        X = self.no.ge[genes[0]]["mean"]['GE']
        y = self.no.ge[genes[1]]["mean"]['GE']

        # GP parameters

        if 'corr' not in kwargs:
            kwargs['corr']  = 'squared_exponential'
        if 'theta0' not in kwargs:
            kwargs['theta0'] = 1e-3
        if 'thetaL' not in kwargs:
            kwargs['thetaL'] = 1e-4
        if 'thetaU' not in kwargs:
            kwargs['thetaU'] = 1e-3
        if 'random_start' not in kwargs:
            kwargs['random_start'] = 100
        if 'nugget' in kwargs:
            kwargs['nugget'] = kwargs['nugget'](y)
        else:
            kwargs['nugget'] = np.var(y)

        # GP
        gp = GaussianProcess(**kwargs)

        # Implementation hiccup; measure against ill-conditioned covariance matrices
        # See: https://github.com/scikit-learn/scikit-learn/issues/4916
        X_ = (X + np.random.normal(0,1e-8, len(X))).reshape(-1,1)
        gp.fit(X_, y)

        x = np.linspace(np.min(X), np.max(X), 5000)
        y_pred, MSE = gp.predict(x.reshape(-1,1), eval_MSE=True)
        sigma = np.sqrt(MSE)

        plt.plot(x, y_pred, 'r-', label=u'$\hat{f}$', linewidth=2)
        plt.plot(X,y, 'bo', ms=4, alpha=.3)
        ci = np.squeeze(y_pred.reshape(1,-1))
        plt.fill_between(x, ci - 1.9600 * sigma, ci + 1.9600 * sigma, alpha=.5, facecolor='red');
        cushion = (np.max(X)-np.min(X))/10.
        plt.xlim(np.min(X)-cushion, np.max(X)+cushion)
        plt.xlabel(gene1)
        plt.ylabel(gene2)


def load_gene_list(csv_path):
    """
    Prints a list NIH registered Entrez ID for which there exists a gene description.

    Parameters
    ----------
    csv_path : str
        path to gene_info.csv

    """
    if isinstance(csv_path, str):
        my_dataframe = pd.read_csv(os.path.join(csv_path, 'gene_info.csv'))
        my_genes = my_dataframe['Entrez'].as_matrix()
        return my_genes
    else:
        raise ValueError('csv_path parameter must be a str')


