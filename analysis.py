"""
analysis.py:
Statistical testing and analysis tools for Nsaba.
Author: Simon Haxby
"""
from nsaba import Nsaba
from geneinfo import gene_info
from scipy import stats
import numpy as np


class NsabaAnalysis(object):

    def __init__(self, nsaba_obj):
        if type(nsaba_obj) == Nsaba:
            self.no = nsaba_obj
        else:
            raise ValueError("NsabaAnalysis() parameter not a Nsaba instance")

    def t_test(self, term, gene, quant):

    def t_test_multi(self, term, quant=None, sample_num=None, nih_fetch_num=20):
        self.no
        ## Parameters:
        term = term_
        sample_num = len(tnsaba.ge.keys())
        top_id_return = 25
        thres = pd.DataFrame(tnsaba.term[term]['ns_act_vector']).quantile(.85)[0]
        ##

        if len(tnsaba.ge) < sample_num:
            raise ValueError("Sample number exceeds stored number of Entrez IDs")

        aba_sam_num = len(tnsaba.ge[random.choice(tnsaba.ge.keys())])
        sam_ids = random.sample(tnsaba.ge.keys(), sample_num)
        ge_mat = tnsaba.make_ge_ns_mat(term, sam_ids).T[:-1]

        gene_p = []
        for eid, ge in zip(sam_ids, ge_mat):
            gt_thres = [ge[i] for i in xrange(aba_sam_num) if tnsaba.term[term]['ns_act_vector'][i] > thres]
            lt_thres = [ge[i] for i in xrange(aba_sam_num) if tnsaba.term[term]['ns_act_vector'][i] <= thres]
            test_stats = stats.ttest_ind(lt_thres, gt_thres)
            d = cohen_d(lt_thres, gt_thres, len(lt_thres), len(gt_thres))
            if test_stats[0] <= 0:
                gene_p.append( (eid, d, test_stats[1]) )
            else:
                continue
