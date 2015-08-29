import pandas as pd
import numpy as np
from nsaba import Nsaba


class NsabaTestMethods(Nsaba):
    def __init__(self):
        Nsaba.__init__(self)

    def cross_check_coords(self):
        if self._check_static_members() == 1:
             return 1

        print "NS MNI Total: %d" % self.ns['database_df'].shape[0]
        unc = self.ns['database_df'].loc[:, 'x':'z'].drop_duplicates()
        print "NS MNI Unique: %d " % unc.shape[0]

        count = 0
        for coord in self.aba['mni_coords'].data:
            if ((unc['x'] == coord[0]) & (unc['y'] == coord[1]) & (unc['z'] == coord[2])).any():
                count += 1

        print "NS / ABA Shared: %d" % count

