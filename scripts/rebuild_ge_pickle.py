from nsaba.nsaba import Nsaba
from nsaba.nsaba.builder import NsabaBuilder
import os

data_dir = '/Users/simonhaxby/Code/Python/nsaba//data_dir'
os.chdir(data_dir)
Nsaba.aba_load()
Nsaba.ns_load()

builder = NsabaBuilder()
builder.get_aba_ge_all()
builder.pickle_ge(output_dir=data_dir)
