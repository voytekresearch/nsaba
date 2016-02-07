"""
analysis.py:
Statistical testing and analysis tools for Nsaba.
Author: Simon Haxby & Torben Noto
"""

from nsaba import Nsaba
from nsabatools import preprint, not_operational
from geneinfo import load_gene_file, get_gene_info, gene_info

import random
import collections
import csv
import os

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.cluster import KMeans
from sklearn import mixture


def cohen_d(x1, x2, n1, n2):
    """
    Computes Cohen's d; an effect size statistic:
    Overview: https://en.wikipedia.org/wiki/Effect_size#Cohen.27s_d
    """
    return (np.mean(x1) - np.mean(x2)) / np.sqrt(((n1-1)*np.var(x1) + (n2-1)*np.var(x2)) / (n1 + n2-2))


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
                kwargs['quant'] = 85
                # NOTE!!: This value ^ has a dependence in term_ge_ttest_multi
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

    @preprint('This may take a couple of minutes ...')
    def term_ge_spearman_rho(self, term):
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
        ge_mat = self.no.matrix_builder([term], self.no.ge.keys())

        # Calculates Spearman's Rho for each across all genes for a given term
        return [stats.spearmanr(ge_mat[:, ge_mat.shape[1]-1], ge_mat[:, r])[0]
                for r in xrange(len(self.no.ge.keys()))]

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
            ax.set_ylabel('Gene Expression')
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
    def term_ge_ttest_multi(self, term, split_method='quant',sample_num=None, **kwargs):
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
                nih_ids = df['entrez'].as_matrix()
                sam_ids = [entrez_id for entrez_id in nih_ids if entrez_id in nih_ids]
                print "Using NIH described genes only; Entrez ID sample size now %d" % (len(sam_ids))

        # Fetching GE/NS activation matrix
        ge_mat = self.no.matrix_builder([term], sam_ids).T[:-1]
        term_act_vector = self.no.matrix_builder([term], sam_ids).T[-1:][0]

        mask = self._split_mask(term_act_vector, method=split_method, **kwargs)

        # Prepping ttest_metrics results dictionary
        ttest_metrics = {'term': term, 'split_method': split_method, "gene_sample_size": sample_num}
        if split_method == 'quant':
            if 'quant' in kwargs:
                ttest_metrics['quant'] = kwargs['quant']
            else:
                ttest_metrics['quant'] = 85
                # NOTE!!! This dependent on _mask_split()'s default quant splitting value

        gene_stats = []
        for eid, ge in zip(sam_ids, ge_mat):
            # Split coordinates in to term and non-term groups
            cont_grp, funct_grp = self._split_groups(ge, mask)
            test_stats = stats.ttest_ind(cont_grp, funct_grp)
            d = cohen_d(cont_grp, funct_grp, len(cont_grp), len(funct_grp))
            # One-sided T-Test
            if test_stats[0] <= 0:
                gene_stats.append(self.gene_rec(eid, d, test_stats[1]))
            else:
                continue
            if eid in kwargs['genes_of_interest']:
                print 'Gene: ' + str(eid) + '  Effect Size: '+str(d)
        # Sort effect sizes from greatest to smallest in magnitude
        gene_stats.sort(key=lambda rec: rec.cohen_d)
        ttest_metrics['results'] = gene_stats
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
                'printer' : bool, optional
                    Print most information about most significant genes.

        Returns
        -------
        top_genes : list
            Returns a list of records with largest effect sizes from ttest_metrics with
            additional gene description information.

        """
        if 'printer' not in kwargs:
            kwargs['printer'] = True
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
                        gene_dat = get_gene_info(kwargs['csv_path'], [rec.entrez])
                        gene_name = gene_dat[0].name
                        gene_description = gene_dat[0].description
                    top_genes.append((rec.entrez, rec.cohen_d, rec.p_value, gene_name, gene_description))
                except IndexError:
                    continue

            if kwargs['printer']:
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

            if kwargs['printer']:
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

    def rho_distr(self, r_values, genes_of_interest=None):
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

    def cohen_d_distr(self, ttest_metrics, genes_of_interest=None, return_fig=False):
        """Visualizing effect-size distribution (Cohen's d)"""

        if genes_of_interest is None:
            genes_of_interest = []

        d_vals = [rec.cohen_d for rec in ttest_metrics['results']]
        ax = plt.axes()
        ax.hist(d_vals, bins=75)
        ax.set_title("Effect Size Distribution (Cohen's d)")
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
                    offsetter += 450/len(genes_of_interest)
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


