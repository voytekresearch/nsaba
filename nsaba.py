# -*- coding: utf-8 -*-
"""
Neurosynth & Allen Brain Atlas: nsaba
Authors: Simon Haxby, Scott Susi, Torben Noto

Last Updated: 7/9/2015
"""
import nibabel
import numpy as np
import pandas as pd
import csv
import os
from neurosynth.base.dataset import Dataset
from neurosynth.analysis import meta, decode, network


def create_aba_npys(folder_path='', output_path='',
        npy_names=['ABA_MicroExpression.npy',
                   'ABA_SampleAnnot.npy',
                   'ABA_Probes.npy']):
                       
    exp_file = 'MicroarrayExpression.csv'
    si_file = 'SampleAnnot.csv'    
    probe_file = 'Probes.csv'
    
    print 'Output Directory: %s' % os.getcwd()

    with open(os.path.join(folder_path, exp_file), 'rb') as f:
        reader = csv.reader(f)
        expression_mat = np.array(list(reader))
     
    np.save(os.path.join(output_path, npy_names[0]), expression_mat)
    expression_mat = None
    print '%s created.' % npy_names[0]

    
    with open(os.path.join(folder_path, si_file), 'rb') as f:
        reader = csv.reader(f)
        sample_info = np.array(list(reader))
        
    np.save(os.path.join(output_path, npy_names[1]), sample_info)
    sample_info = None
    print '%s created.' % npy_names[1]

    
    with open(os.path.join(folder_path, probe_file), 'rb') as f:
        reader = csv.reader(f)
        probes_mat = np.array(list(reader))
        
    np.save(os.path.join(output_path, npy_names[2]), probes_mat)
    probes_mat = None
    print '%s created.' % npy_names[2]
    return 0
    
    
def create_ns_pkl(folder_path='.', output_path='.'):
    db_path = os.path.join(folder_path, 'database.txt')
    dataset = Dataset(db_path)
    dataset.add_features('features.txt')
    pkl_path = os.path.join(output_path, 'ns_dataset.pkl')
    print 'Output Directory: %s' % output_path
    dataset.save(pkl_path)
    print 'ns_dataset.pkl created.'
    return 0

class Nsaba(object):
    aba = { 
        'exp_mat': None,
        'exp_df': None,
        'probe_mat': None,   
        'probe_df': None,
        'si_mat': None,
        'si_df': None,
    }

    mni_coords = None
    
    ns_db = None
    
    @classmethod
    def aba_load(cls, path='.',
            npy_names=['ABA_MicroExpression.npy',
                       'ABA_SampleAnnot.npy',
                       'ABA_Probes.npy']):
                           
        if len(npy_names) != 3:
            raise IndexError, "'npy_names' must a list of 3 'str' variables"
            
        print 'This may take a minute or two ...'
  
        npy_path = os.path.join(path, npy_names[0])
        cls.aba['exp_mat'] = np.load(npy_path)
        cls.aba['exp_df']  = pd.DataFrame(cls.aba['exp_mat'])
        print '%s loaded.' % npy_names[0]
        
        npy_path = os.path.join(path, npy_names[1])
        cls.aba['si_mat'] = np.load(npy_path)
        cls.aba['si_df']  = pd.DataFrame(cls.aba['si_mat'])
        print '%s loaded.' % npy_names[1]

        npy_path = os.path.join(path, npy_names[2])
        cls.aba['probe_mat'] = np.load(npy_path)
        cls.aba['probe_df'] = pd.DataFrame(cls.aba['probe_mat'])
        print '%s loaded.' % npy_names[2]
        
        cls.mni_coords = cls.aba['si_mat'][1:,10:].astype(float)
        print 'Nsaba.mni_coords intialized.'
        return 0
        
    @classmethod
    def ns_load(cls, file_path):
        cls.ns_db = Dataset.load(file_path)
        print '%s loaded' % os.path.basename(file_path)
        return 0
         
    def __init__(self):
        
        self.ge = {}
       # ...
       # ...
  
    def get_ge(self, entrez_ids):
        
        if self.__check_static_members() == 1:
            return 1
            
        try:
            iter(entrez_ids)
        except TypeError:
            print "Invalid parameter form; please contain entrez ids as 'str' types in iterable container"
            return 1
        else:
            if isinstance(entrez_ids, str):
                print "Invalid parameter form; please contain entrez ids as 'str' types in iterable container"
                return 1
                

        for entrez_id in entrez_ids:
            probe_ids = pd.DataFrame(self.aba['probe_df'].loc[self.aba['probe_df'][5] == \
                entrez_id]).index.tolist()
                
            if len(probe_ids) == 0:
                print 'Entrez ID: %s not registered with ABA database' % entrez_id
                continue
            
            ge_df = pd.DataFrame(self.aba['exp_df'].loc[probe_ids.pop(0)-1, :])
            for probe_id in probe_ids:
                ge_df = ge_df.join(self.aba['exp_df'].loc[probe_id-1, :])
            
            ge_mat = ge_df.as_matrix().astype(float)[1:,:]
            self.ge[entrez_id] = []
             
            self.ge[entrez_id].append([np.mean(row) for row in ge_mat])
            diff_arr = [np.amax(row) - np.amin(row) for row in ge_mat] 
            self.ge[entrez_id].append(diff_arr)
            self.ge[entrez_id].append(np.amax(diff_arr))
            self.ge[entrez_id].append(len(ge_mat[0]))

        return 0

              
    def make_ge_ns_mat(self):
        if self.__check_static_members() == 1:
            return 1 
            
    def get_ns_act(self, query, term_freq=.005):
         if self.__check_static_members() == 1:
            return 1
            
    def __check_static_members(self):
        for val in self.aba.itervalues():
            if val is None:
                print "Unassigned Nsaba 'aba' static variable: see Nsaba.aba_load(path)"
                return 1
        if self.ns_db is None:
            print "Unassigned Nsaba 'ns_db' static variable: see Nsaba.ns_load(path)"
            return 1
        return 0

        
# -----------------------------------------------------------

# Example Usage

'''
        
from nsaba import Nsaba

path 'blab' #path is folder location
Nsaba.load(path)    
A = Nsaba()

A.get_ge('7889')
A.get_ns_act('Alzhiemers')
mat = A.make_ge_ns_mat()

B = Nsaba()

'''

# Run science here

# -----------------------------------------------------------
