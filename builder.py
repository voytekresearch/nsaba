from nsaba import Nsaba

class NsabaBuilder(Nsaba):
    def __init__(self):
        Nsaba.__init__(self)

    def get_aba_ge_all(self):
        """ Returns a dictionary with ABA gene expression coefficient across all genes
        at sampled locations"""

        if self.__check_static_members() == 1:
            return 1

        warning_flag = True
        while warning_flag:
            y_n = raw_input("WARNING: this operation can take upwards of an hour, proceed? (Y/n): ")
            if y_n == 'Y':
                warning_flag = False
            elif y_n == 'n':
                return 0
            else:
                print "Invalid response: %s" % y_n

        entrez_ids = self.aba['probe_df']['entrez_id'][
            self.aba['probe_df']['entrez_id'].notnull()].unique().astype(int)

        self.get_aba_ge(entrez_ids)


