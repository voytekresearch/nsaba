__author__ = 'Torben'

from nsaba import Nsaba
from geneinfo import gene_info
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt


class NsabaAnalysis2(object):

    def __init__(self, nsaba_obj):
        if type(nsaba_obj) == Nsaba:
            self.no = nsaba_obj
        else:
            raise ValueError("NsabaAnalysis() parameter not a Nsaba instance")

        print "To use inline plotting functionality in Jupyter, '%matplotlib inline' must be enabled"

    def term_to_genes(self, term):
        if term not in self.no.term:
            raise ValueError("Term activation not generated for '%s" % term)
        ge_mat = self.no.make_ge_ns_mat(term, self.no.ge.keys())

        return [stats.spearmanr(ge_mat[:, ge_mat.shape[1]-1], ge_mat[:, r])[0] for r in xrange(len(self.no.ge.keys()))]

    def plot_rho_distr(self, r_values, genes_of_interest=[]):
        """Visualizing effect-size distribution"""
        ax = plt.axes()
        ax.hist(r_values, bins=75)
        ax.set_title("Distribution of spearman r values")
        ax.set_xlabel('spearman r')
        ax.set_ylabel('count')

        if genes_of_interest:
            offsetter = 450/len(genes_of_interest)
            for r in xrange(len(r_values)):
                if self.no.ge.keys()[r] in genes_of_interest:
                    plt.plot([r_values[r], r_values[r]], [0, offsetter])
                    plt.annotate('Gene:'+str(self.no.ge.keys()[r])+' d='+str(r_values[r]),
                                 [r_values[r], offsetter])
                    offsetter += 450/len(genes_of_interest)

    def fetch_gene_descriptions(self, r_values, nih_fetch_num=20, alpha=.05, printme=True):
        """Prints: ID, spearman's r, gene name, and gene description for genes with the largest effect sizes"""
        top_ids = np.argsort(r_values)[:nih_fetch_num]

        top_rs = [r for r in reversed(np.sort(r_values))][:nih_fetch_num]
        top_genes = []
        for x in xrange(len(top_ids)):
            try:
                top_genes.append((int(top_ids[x]), top_rs[x], gene_info(str(int(top_ids[x])))[0],
                                  gene_info(str(int(top_ids[x])))[1]))
            except TypeError:
                continue

        if printme:
            print "\nCorrected Bonferroni Alpha: %.3E\n\n" % (alpha/float(len(self.no.ge.keys())))
            for eid, rho, gene_i, descr in top_genes:
                if len(descr) == 1:
                    print "%d (r = %.3f): < No description found >\n\n" % (eid, rho)
                else:
                    print "%d (r = %.3f): %s\n %s\n\n" % (eid, rho, gene_i, descr)
        return top_genes
