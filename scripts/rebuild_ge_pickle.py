#!/bin/python

# Rebuilds/builds gene expression pickle
# for all genes; this script should be run
# after pulling from master to ensure that
# the pickled representation of the data 
# matches the methods assumptions about it's
# form.

# Author: Simon Haxby 

from nsaba.nsaba import Nsaba
from nsaba.nsaba.builder import NsabaBuilder
import os

aba_dir = <INSERT DIRECTORY HERE>
ns_dir = <INSERT DIRECTORY HERE>

Nsaba.aba_load(aba_dir)
Nsaba.ns_load(ns_dir)

builder = NsabaBuilder()
builder.get_aba_ge_all()

pkl_dir = <INSERT DIRECTORY HERE>

builder.pickle_ge(output_dir=pkl))
