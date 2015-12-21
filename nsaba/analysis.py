"""
analysis.py:
Statistical testing and analysis tools for Nsaba.
Author: Simon Haxby
"""

from nsaba import Nsaba
from nsabatools import preprint, not_operational
<<<<<<< HEAD
from geneinfo import load_gene_file, get_gene_info

import random
import collections
import csv

=======
import random
import collections
from scipy import stats
>>>>>>> f5e917ca8716185fddaa6d6197779d53278fb897
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn import mixture
import csv
import pandas as pd


def cohen_d(x1, x2, n1, n2):
    return (np.mean(x1) - np.mean(x2)) / np.sqrt(((n1-1)*np.var(x1) + (n2-1)*np.var(x2)) / (n1 + n2-2))


class NsabaAnalysis(object):

    def __init__(self, nsaba_obj):
        if type(nsaba_obj) == Nsaba:
            self.no = nsaba_obj
        else:
            raise ValueError("NsabaAnalysis() parameter not a Nsaba instance")

        self.gene_rec = collections.namedtuple("gene_rec", "entrez cohen_d p_value")
        print "To use inline plotting functionality in Jupyter, '%matplotlib inline' must be enabled"

    def _split_mask(self, term_vec, method='var', **kwargs):
        """ Constructs splitting mask via specified method"""
        if method == 'var':
            if 'quant' not in kwargs:
                kwargs['quant'] = 85
            thres = np.percentile(term_vec, kwargs['quant'])
            mask = np.array([True if coeff > thres else False for coeff in term_vec])
            return mask
        elif method == 'kmeans':
            X = [[x] for x in term_vec]
            init_clusters = np.array([[min(term_vec)], [max(term_vec)]])
            kmn = KMeans(n_clusters=2, init=init_clusters, n_init=1)
            mask = kmn.fit_predict(X)
            return mask.astype(bool)
        elif method == 'mog':
            X = [[x] for x in term_vec]
            gm = mixture.GMM(n_components=2)
            gm.fit(X)
            if stats.mode(gm.predict(X))[0][0] == 0:
                mask = np.array([True if gm.predict([coeff]) == 1 else False for coeff in term_vec])
            else:
                mask = np.array([False if gm.predict([coeff]) == 1 else True for coeff in term_vec])
            return mask
        else:
            raise ValueError("method parameter '%s' not recognized" % method)

    def _split_groups(self, ge_vec, mask):
        """Splits gene samples into control and functional network"""
        functional_grp = ge_vec[mask]
        diff = set(ge_vec) - set(functional_grp)
        control_grp = np.array(list(diff))
        return control_grp, functional_grp

    def t_test(self, term, gene, split_method='var', log=False, graphops='density', **kwargs):
        """ T-Test of gene expression between term and non-term coordinates"""
        analmat = self.no.make_ge_ns_mat(term, [gene])
        # Splitting groups
        mask = self._split_mask(analmat[:, 1], method=split_method, **kwargs)
        cont_grp, funct_grp = self._split_groups(analmat[:, 0], mask)
        if log:
            funct_grp = np.log(funct_grp)
            cont_grp = np.log(cont_grp)

        # T-Test
        print "t-value: %.4f \np-value: %.3E" % stats.ttest_ind(cont_grp, funct_grp)
        print "Effect size: %.4f \n" % cohen_d(cont_grp, funct_grp, len(cont_grp), len(funct_grp))
        # Histogram/KDE Plots
        if graphops == 'density':
            ax = plt.axes()
            ax.set_title('Gene Expression Distributions')
            ax.set_xlabel('gene expression')
            ax.set_ylabel('density')
            sns.distplot(funct_grp, ax=ax, label=term)
            sns.distplot(cont_grp, label='null')
            plt.legend()
        elif graphops == 'box':
            ax = plt.axes()
            ax.boxplot([cont_grp, funct_grp])
            ax.set_xticks([1, 2])
            ax.set_xticklabels(['Low '+term, 'High '+term])
        elif graphops == 'violin':
            ax = plt.axes()
            ax.violinplot([cont_grp, funct_grp])
            ax.set_xticks([1, 2])
            ax.set_xticklabels(['Low '+term, 'High '+term])
            ax.set_ylabel('Gene Expression')
            ax.plot(np.ones(len(cont_grp)), cont_grp, 'b.')
            ax.plot(1, np.mean(cont_grp), 'bs')
            ax.plot(2*np.ones(len(funct_grp)), funct_grp, 'g.')
            ax.plot(2, np.mean(funct_grp), 'gs')
        else:
            raise ValueError("graphops parameter '%s' not recognized" % graphops)

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
    def t_test_multi(self, term, quant=None, sample_num=None, split_method='var',
                     genes_of_interest=None, **kwargs):
        if term not in self.no.term:
            raise ValueError("Term activation not generated for '%s" % term)
        if sample_num == None:
            sample_num = len(self.no.ge.keys())
        elif sample_num <= 0:
            raise ValueError("'sample_num' parameter must be greater than 0")
        if genes_of_interest == None:
            genes_of_interest = []

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
                sam_ids = df['Entrez'].as_matrix()
                # does nothing nih_ids = df['Entrez'].as_matrix()
                # sam_ids = [entrez_id for entrez_id in nih_ids if entrez_id in nih_ids]
                print "Using NIH described genes only; Entrez ID sample size now %d" % (len(sam_ids))

        ge_mat = self.no.make_ge_ns_mat(term, sam_ids).T[:-1]
        term_act_vector = self.no.make_ge_ns_mat(term, sam_ids).T[-1:][0]

        mask = self._split_mask(term_act_vector, method=split_method, quant=quant)

        ttest_metrics = {'term': term, 'split_method': split_method, "quantile": quant, "gene_sample_size": sample_num}
        gene_stats = []
        for eid, ge in zip(sam_ids, ge_mat):
            # Split coordinates in to term and non-term groups
            cont_grp, funct_grp = self._split_groups(ge, mask)
            test_stats = stats.ttest_ind(cont_grp, funct_grp)
            d = cohen_d(cont_grp, funct_grp, len(cont_grp), len(funct_grp))
            # One-sided T-Test
            if 'full' in kwargs:
                gene_stats.append(self.gene_rec(int(eid), d, test_stats[1]))
            else:
                if test_stats[0] <= 0:
                    gene_stats.append(self.gene_rec(int(eid), d, test_stats[1]))
                else:
                    continue
            if eid in genes_of_interest:
                print 'Gene: ' + str(int(eid)) + '  Effect Size: '+str(d)
        # Sort effect sizes from greatest to smallest in magnitude
        gene_stats.sort(key=lambda rec: rec.cohen_d)
        ttest_metrics['results'] = gene_stats
        return ttest_metrics

    @preprint('Fetching NIH gene descriptions ...')
    def fetch_gene_descriptions(self, ttest_metrics, gene_path='.', nih_fetch_num=20, alpha=.05, printme=True):
        """Prints: ID, p-value, Cohen's d, gene description for genes with the largest effect sizes"""
        top_genes = []
        for rec in ttest_metrics['results'][:nih_fetch_num]:
            print rec
            try:
                gene_dat = get_gene_info(gene_path, [str(int(rec.entrez))])
                gene_name = gene_dat[0][1]
                gene_description = gene_dat[0][2]
                top_genes.append((rec.entrez, rec.cohen_d, rec.p_value, gene_name, gene_description))
            except IndexError:
                continue

        if printme:
            print "\nCorrected Bonferroni Alpha: %.3E\n\n" % (alpha/float(ttest_metrics['gene_sample_size']))
            for eid, coh_d, p_val, gene_i, descr in top_genes:
                if len(descr) == 1:
                    print "%d (p = %.3E; d = %.3f): < No description found >\n\n" % (eid, p_val, coh_d)
                else:
                    print "%d (p = %.3E; d = %.3f): %s\n\n" % (eid, p_val, coh_d, descr)
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

    def effect_size_distr(self, ttest_metrics, genes_of_interest=[], return_fig=False):
        """Visualizing effect-size distribution"""
        d_vals = [rec.cohen_d for rec in ttest_metrics['results']]
        ax = plt.axes()
        ax.hist(d_vals, bins=75)
        ax.set_title("Effect Size Distribution (Cohen's d)")
        ax.set_xlabel('effect sizes')
        ax.set_ylabel('frequency')

        offsetter = 450/len(genes_of_interest)
        for rec in ttest_metrics['results']:
            if int(rec.entrez) in genes_of_interest:
                plt.plot([rec.cohen_d, rec.cohen_d], [0, offsetter])
                plt.annotate('Gene:'+str(int(rec.entrez))+' d='+str(rec.cohen_d),
                             [rec.cohen_d, offsetter])
                offsetter += 450/len(genes_of_interest)
        if return_fig:
            return ax

    def save_results(self, ttest_metrics, direct, fname=None, no_genes=20):
        if fname is None:
            fname = 'nsaba_output_figure'

        top_genes = self.fetch_gene_descriptions(ttest_metrics, nih_fetch_num=no_genes, printme=False)
        eids = [int(i[0]) for i in top_genes]
        myfig = self.effect_size_distr(ttest_metrics, genes_of_interest=eids[0:20], return_fig=True)
        plt.savefig(direct+fname+'.png')

        with open(direct+fname+'.csv', 'wb') as csvfile:
            writer = csv.writer(csvfile)
            for i in top_genes:
                writer.writerow([i[0], i[3], i[1], i[2], i[4]])

    def term_to_genes(self, term, method='spearman', **kwargs):
        if term not in self.no.term:
            raise ValueError("Term activation not generated for '%s" % term)
        if 'nih_only' in kwargs:
            if kwargs['nih_only']:
                if 'gi_csv_path' in kwargs:
                    gi_path = kwargs['gi_csv_path']
                    df = load_gene_file(gi_path)
                else:
                    df = load_gene_file()
                nih_ids = df['Entrez'].as_matrix()
                nih_genes = [entrez_id for entrez_id in nih_ids if entrez_id in nih_ids]
                print "Using NIH described genes only; Entrez ID sample size now %d" % (len(nih_genes))

            ge_mat = self.no.make_ge_ns_mat(term, nih_genes)
        else:
            ge_mat = self.no.make_ge_ns_mat(term, self.no.aba.keys())

        if method == 'spearman':
            if nih_genes:
                r_vals = [stats.spearmanr(ge_mat[:, ge_mat.shape[1]-1], ge_mat[:, r])[0] for r in xrange(len(nih_genes))]
            else:
                r_vals = [stats.spearmanr(ge_mat[:, ge_mat.shape[1]-1], ge_mat[:, r])[0] for r in xrange(len(self.no.ge.keys()))]
            return r_vals
        if method == 'pearson':
            if nih_genes:
                r_vals = [np.corrcoef(ge_mat[:, ge_mat.shape[1]-1], ge_mat[:, r])[1, 0] for r in xrange(len(nih_genes))]
            else:
                r_vals = [np.corrcoef(ge_mat[:, ge_mat.shape[1]-1], ge_mat[:, r])[1, 0] for r in xrange(len(self.no.ge.keys()))]
            return r_vals
        if method == 'regression':
            m_vals = []
            if nih_genes:
                for gene in xrange(len(nih_genes)):
                    X = np.vstack([ge_mat[:, ge_mat.shape[1]-1], np.ones(len(ge_mat[:, 0]))]).T
                    m, c = np.linalg.lstsq(X, ge_mat[:, gene])[0]
                    m_vals.append(m)
            else:
                for gene in xrange(len(self.no.ge.keys())):
                    X = np.vstack([ge_mat[:, ge_mat.shape[1]-1], np.ones(len(ge_mat[:, 0]))]).T
                    m, c = np.linalg.lstsq(X, ge_mat[:, gene])[0]
                    m_vals.append(m)
            return m_vals

    def validate_with_t_test(self, term, genes, method='t_test', quant=None, split_method='var',
                             non_zero=False, gi_csv_path= '/Users/Torben/Code/nsaba/'):
        stat_output = []

        if method == 't_test':
            ttest_metrics = self.t_test_multi(term, quant=quant,
                                              split_method=split_method,
                                              full=True, nih_only=True,
                                              gi_csv_path=gi_csv_path)
            for i in ttest_metrics['results']:
                if int(i[0]) in genes:
                    stat_output.append(i[1])
            return stat_output

        if method == 'pearson':
            for gene in genes:
                ge_mat = self.no.make_ge_ns_mat(term, [gene])
                if non_zero:
                    min_val = np.min(ge_mat[:, ge_mat.shape[1]][np.nonzero(ge_mat[:, ge_mat.shape[1]])]) # min nonzero value
                    term_axis = ge_mat[:, ge_mat.shape[1]][ge_mat[:, ge_mat.shape[1]] > min_val]
                    ge_axis = ge_mat[:, 0][ge_mat[:, ge_mat[:, ge_mat.shape[1]]] > min_val]
                else:
                    term_axis = ge_mat[:, ge_mat.shape[1]]
                    ge_axis = ge_mat[:, 0]
                stat_output.append(np.corrcoef(ge_axis, term_axis)[1, 0])

            return stat_output

        if method == 'spearman':
            for gene in genes:
                ge_ns_mat = self.no.make_ge_ns_mat(term, [gene])
                if non_zero:
                    min_val = np.min(ge_mat[:, ge_mat.shape[1]][np.nonzero(ge_mat[:, ge_mat.shape[1]])]) # min nonzero value
                    term_axis = ge_mat[:, ge_mat.shape[1]][ge_mat[:, ge_mat.shape[1]] > min_val]
                    ge_axis = ge_mat[:, 0][ge_mat[:, ge_mat[:, ge_mat.shape[1]]] > min_val]
                else:
                    term_axis = ge_mat[:, ge_mat.shape[1]]
                    ge_axis = ge_mat[:, 0]
                stat_output.append(stats.spearmanr(ge_axis, term_axis)[0])
            return stat_output

        if method == 'regression':
            for gene in genes:
                ge_ns_mat = self.no.make_ge_ns_mat(term, gene)
                if non_zero:
                    min_val = np.min(ge_mat[:, ge_mat.shape[1]][np.nonzero(ge_mat[:, ge_mat.shape[1]])]) # min nonzero value
                    term_axis = ge_mat[:, ge_mat.shape[1]][ge_mat[:, ge_mat.shape[1]] > min_val]
                    ge_axis = ge_mat[:, 0][ge_mat[:, ge_mat[:, ge_mat.shape[1]]] > min_val]
                else:
                    term_axis = ge_mat[:, ge_mat.shape[1]]
                    ge_axis = ge_mat[:, 0]
                X = np.vstack([ge_axis, np.ones(len(ge_axis))]).T
                m, c = np.linalg.lstsq(X, term_axis)[0]
                stat_output.append(m)
            return stat_output


    def validate_by_alpha(self, term, genes, method='t_test', quant=85, alpha=0.05, nih_only=False, gi_csv_path='.'):
        genes_found = []
        if method == 't_test':
            ttest_metrics = self.t_test_multi(term, quant=quant, nih_only=nih_only, gi_csv_path=gi_csv_path)
            alpha_studies = int(alpha * len(ttest_metrics['results']))

            sorted_studies = sorted(ttest_metrics['results'], key=lambda x: x.cohen_d)[:alpha_studies]

            for i in sorted_studies:
                for gene in genes:
                    if int(gene) == int(i.entrez):
                        genes_found.append(gene)
            return genes_found

        if method == 'pearson':
            r_vals = self.term_to_genes(term=term, method='pearson', nih_only=nih_only, gi_csv_path=gi_csv_path)
            alpha_study_threshold = int(alpha * len(r_vals))
            alpha_index = [i for i in reversed(sorted(r_vals))][alpha_study_threshold]
            print "pearson's r must be > " + str(alpha_index)
            for gene in genes:
                ge_ns_mat = self.no.make_ge_ns_mat(term, [gene])
                r_val = np.corrcoef(ge_ns_mat[:, 0], ge_ns_mat[:, 1])[1, 0]
                if r_val > alpha_index:
                    genes_found.append(gene)
            return genes_found

        if method == 'spearman':
            r_vals = self.term_to_genes(term=term, method='spearman', nih_only=nih_only, gi_csv_path=gi_csv_path)
            alpha_study_threshold = int(alpha * len(r_vals))
            alpha_index = [i for i in reversed(sorted(r_vals))][alpha_study_threshold]
            print "spearman's r must be > " + str(alpha_index)
            for gene in genes:
                ge_ns_mat = self.no.make_ge_ns_mat(term, [gene])
                r_val = stats.spearmanr(ge_ns_mat[:, 0], ge_ns_mat[:, 1])[0]
                if r_val > alpha_index:
                    genes_found.append(gene)
            return genes_found

        if method == 'regression':
            m_vals = self.term_to_genes(term=term, method='regression',nih_only=nih_only, gi_csv_path=gi_csv_path)
            alpha_study_threshold = int(alpha * len(m_vals))
            alpha_index = [i for i in reversed(sorted(m_vals))][alpha_study_threshold]
            print "slope of linear regression must be > " + str(alpha_index)
            for gene in genes:
                ge_ns_mat = self.no.make_ge_ns_mat(term, [gene])
                X = np.vstack([ge_ns_mat[:, 0], np.ones(len(ge_ns_mat[:, 0]))]).T
                m, c = np.linalg.lstsq(X, ge_ns_mat[:, 1])[0]
                if m > alpha_index:
                    genes_found.append(gene)
            return genes_found
        else:
            print 'Invalid analysis method. Use pearson, spearman, regression, or t_test'
            return []


def load_gene_list(path, filename):
        if isinstance(path, str):
            if isinstance(filename, str):
                my_dataframe = pd.read_csv(path+filename)
                my_genes = my_dataframe['Entrez'].as_matrix()
                return my_genes
            else:
                print 'filename must be a string'
        else:
            print 'filename must be a string'

