"""
analysis.py:
Statistical testing and analysis tools for Nsaba.
Author: Simon Haxby
"""
from nsaba import Nsaba
from geneinfo import gene_info
from nsabatools import preprint
import collections
from scipy import stats
import numpy as np
import seaborn as sns
import random


def cohen_d(x1, x2, n1, n2):
    return (np.mean(x1) - np.mean(x2)) / np.sqrt(((n1-1)*np.var(x1) + (n2-1)*np.var(x2)) / (n1 + n2-2))


class NsabaAnalysis(object):

    def __init__(self, nsaba_obj):
        if type(nsaba_obj) == Nsaba:
            self.no = nsaba_obj
        else:
            raise ValueError("NsabaAnalysis() parameter not a Nsaba instance")

        self.gene_rec = collections.namedtuple("gene_rec", "entrez cohen_d p_value")
        print "To use seaborn plotting functionality in Jupyter, '%matplotlib inline' must be enabled"

    def t_test(self, term, gene, quant=None):
        """ T-Test of gene expression between term and non-term coordinates"""
        if not quant:
            quant = 85
        analmat = self.no.make_ge_ns_mat(term, [gene])
        thres = np.percentile(analmat[:, 1], quant)
        # Splitting groups
        gt_thres = [row[0] for row in analmat if row[1] > thres]
        lt_thres = [row[0] for row in analmat if row[1] <= thres]
        # T-Test
        print "t-value: %.4f \np-value: %.3E" % stats.ttest_ind(lt_thres, gt_thres)
        print "Effect size: %.4f \n" % cohen_d(lt_thres, gt_thres, len(lt_thres), len(gt_thres))
        # Histogram/KDE Plots
        sns.distplot(gt_thres);
        sns.distplot(lt_thres);

    @preprint('This may take a couple of minutes ...')
    def t_test_multi(self, term, quant=None, sample_num=None):
        if term not in self.no.term:
            raise ValueError("Term activation not generated for '%s" % term)
        if not sample_num:
            sample_num = len(self.no.ge.keys())
        elif sample_num <= 0:
            raise ValueError("'sample_num' parameter must be greater than 0")
        if not quant:
            quant = 85

        thres = np.percentile(self.no.term[term]['ns_act_vector'], quant)

        if len(self.no.ge) < sample_num:
            raise ValueError("Sample number exceeds stored number of Entrez IDs")


        sam_ids = random.sample(self.no.ge.keys(), sample_num)
        ge_mat = self.no.make_ge_ns_mat(term, sam_ids).T[:-1]
        aba_sam_num = len(ge_mat[0])

        ttest_metrics = {'term': term, "quantile": quant, "gene_sample_size": sample_num}
        gene_stats = []
        for eid, ge in zip(sam_ids, ge_mat):
            # Split coordinates in to term and non-term groups
            gt_thres = [ge[i] for i in xrange(aba_sam_num) if self.no.term[term]['ns_act_vector'][i] > thres]
            lt_thres = [ge[i] for i in xrange(aba_sam_num) if self.no.term[term]['ns_act_vector'][i] <= thres]
            test_stats = stats.ttest_ind(lt_thres, gt_thres)
            d = cohen_d(lt_thres, gt_thres, len(lt_thres), len(gt_thres))
            # One-sided T-Test
            if test_stats[0] <= 0:
                gene_stats.append(self.gene_rec(eid, d, test_stats[1]))
            else:
                continue

        # Sort effect sizes from greatest to smallest in magnitude
        gene_stats.sort(key=lambda rec: rec.cohen_d)
        ttest_metrics['results'] = gene_stats

        return ttest_metrics

    @preprint('Fetching NIH gene descriptions ...')
    def fetch_gene_descriptions(self, ttest_metrics, nih_fetch_num=20, alpha=.05):
        """Prints: ID, p-value, Cohen's d, gene description for genes with the largest effect sizes"""
        top_genes = []
        for rec in ttest_metrics['results'][:nih_fetch_num]:
            try:
                top_genes.append((rec.entrez, rec.cohen_d, rec.p_value, gene_info(str(rec.entrez))[0]))
            except TypeError:
                continue

        print "\nCorrected Bonferroni Alpha: %.3E\n\n" % (alpha/float(ttest_metrics['gene_sample_size']))
        for eid, coh_d, p_val, descr in top_genes:
            if len(descr) == 1:
                print "%d (p = %.3E; d = %.3f): < No description found >\n\n" % (eid, p_val, coh_d)
            else:
                print "%d (p = %.3E; d = %.3f): %s\n\n" % (eid, p_val, coh_d, descr)

    def p_val_distr(self, ttest_metrics):
        """Visualizing p-value distribution"""
        p_vals = [rec.p_value for rec in ttest_metrics['results']]
        sig = sum((p < .05/float(ttest_metrics['gene_sample_size']) for p in p_vals))
        print "Percent Significant (Bonferroni Correction; alpha = .05): %.3f %%" % \
              (100*sig/float(ttest_metrics['gene_sample_size']))
        sns.distplot(p_vals, norm_hist=False, bins=75, kde=False, axlabel='p-values');

    def effect_size_distr(self, ttest_metrics):
        """Visualizing effect-size distribution"""
        d_vals = [rec.cohen_d for rec in ttest_metrics['results']]
        sns.distplot(d_vals, norm_hist=False, bins=75, kde=False, axlabel='effect sizes');
