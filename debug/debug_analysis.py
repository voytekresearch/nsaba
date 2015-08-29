from analysis import NsabaAnalysis
from nsaba import Nsaba

data_dir = '../data_dir'
Nsaba.aba_load(data_dir)
Nsaba.ns_load(data_dir)

tsaba = Nsaba()
tsaba.load_ge_pickle(path=data_dir)
tsaba.get_ns_act('reward')
anal = NsabaAnalysis(tsaba)
anal.t_test_multi('reward')
