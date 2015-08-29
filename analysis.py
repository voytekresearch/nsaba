"""
analysis.py:
Statistical testing and analysis tools for Nsaba.
Author: Simon Haxby
"""
from nsaba import Nsaba
from geneinfo import gene_info
from scipy import stats
import seaborn as sns
import random
import numpy as np
import pandas as pd


def cohen_d(x1, x2, n1, n2):
    return (np.mean(x1) - np.mean(x2)) / np.sqrt( ((n1-1)*np.var(x1) + (n2-1)*np.var(x2)) / (n1 + n2 -2) )


class NsabaAnalysis(object):

    def __init__(self, nsaba_obj):
        if type(nsaba_obj) == Nsaba:
            self.no = nsaba_obj
        else:
            raise ValueError("NsabaAnalysis() parameter not a Nsaba instance")
        
        print "To use seaborn plotting functionality in Jupyter, '%matplotlib inline' must be enabled"

    def t_test(self, term, gene, quant):
        """ T-Test of gene expression between term and non-term coordinates"""
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

    def t_test_multi(self, term, quant=None, sample_num=None, nih_fetch_num=20):

        if term not in self.no.term:
            raise ValueError("Term activation not generated for '%s" % term)
        if not sample_num:
            sample_num = len(self.no.ge.keys())
        if not quant:
            thres = np.percentile(self.no.term[term]['ns_act_vector'], 85)
        else:
            thres = np.percentile(self.no.term[term]['ns_act_vector'], quant)

        if len(self.no.ge) < sample_num:
            raise ValueError("Sample number exceeds stored number of Entrez IDs")

        aba_sam_num = len(self.no.ge[random.choice(self.no.ge.keys())])
        sam_ids = random.sample(self.no.ge.keys(), sample_num)
        ge_mat = self.no.make_ge_ns_mat(term, sam_ids).T[:-1]

        gene_p = []
        for eid, ge in zip(sam_ids, ge_mat):
            gt_thres = [ge[i] for i in xrange(aba_sam_num) if self.no.term[term]['ns_act_vector'][i] > thres]
            lt_thres = [ge[i] for i in xrange(aba_sam_num) if self.no.term[term]['ns_act_vector'][i] <= thres]
            test_stats = stats.ttest_ind(lt_thres, gt_thres)
            d = cohen_d(lt_thres, gt_thres, len(lt_thres), len(gt_thres))
            if test_stats[0] <= 0:
                gene_p.append((eid, d, test_stats[1]))
            else:
                continue

        top_genes = []
        for i in xrange(nih_fetch_num):
            try:
                top_genes.append( (gene_p[i][0], gene_p[i][1], gene_p[i][2],
                                   gene_info(str(gene_p[i][0]))[0]) )
            except TypeError:
                continue
        alpha = .05

        print "Corrected Bonferoni Alpha: %.3E\n\n" % (alpha/float(sample_num))
        for eid, coh_d, p_val, descr in top_genes:
            if len(descr) == 1:
                print "%d (p = %.3E; d = %.3f): < No description found >\n\n" % (eid, p_val, coh_d)
            else:
                print "%d (p = %.3E; d = %.3f): %s\n\n" % (eid, p_val, coh_d, descr)

        # Visualizing p-value distribution
        p_vals = [p[2] for p in gene_p ]
        sig = sum([ p < .05/sample_num for p in p_vals])
        print "Percent Significant (Bonferroni Correction; alpha = .05): %.3f %%" % (100*sig/float(sample_num))
        sns.distplot(p_vals, norm_hist=False, bins=75, kde=False);

        # Visualizing effect-size distribution
        p_vals = [p[1] for p in gene_p ]
        sns.distplot(p_vals, norm_hist=False, bins=75, kde=False);
