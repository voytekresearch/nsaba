# -*- coding: utf-8 -*-
"""
Nsaba Unit-Testing Framework
Created on Mon Jul 13 16:33:15 2015

Designed to be run with nosetests

@author: Simon Haxby
"""
import unittest
from nsaba import Nsaba

class NsabaLoaderTest(unittest.TestCase):

    def test_aba_load(self):
        aba_dir = 'data_dir'
        self.assertEquals(Nsaba.aba_load(aba_dir), 0)

    def test_static_si_intialize(self):
        self.assertIsNotNone(Nsaba.aba['si_mat'])
        self.assertIsNotNone(Nsaba.aba['si_df'])

    def test_static_exp_intialize(self):
        self.assertIsNotNone(Nsaba.aba['exp_mat'])
        self.assertIsNotNone(Nsaba.aba['exp_df'])

    def test_static_probe_intialize(self):
       self.assertIsNotNone(Nsaba.aba['probe_mat'])
       self.assertIsNotNone(Nsaba.aba['probe_df'])

    def test_mni_coords(self):
        self.assertIsNotNone(Nsaba.mni_coords)

class NsabaTestCase(unittest.TestCase):
      
    def test_get_ge_fail(self):
        nsaba_fail1 = Nsaba()
        nsaba_fail2 = Nsaba()

        self.assertEquals(nsaba_fail1.get_ge('2'), 1)
        self.assertEquals(nsaba_fail2.get_ge(2), 1)

    def test_get_ge_pass_single(self):
        nsaba_pass = Nsaba()
        self.assertEquals(nsaba_pass.get_ge(['2']), 0)
        self.assertEquals(len(nsaba_pass.ge), 1)
        self.assertEquals(len(nsaba_pass.ge['2']), 4)
        self.assertEquals(len(nsaba_pass.ge['2'][0]), 893)
        self.assertEquals(nsaba_pass.ge['2'][3], 4)

    def test_get_ge_pass_multi(self):
        nsaba_pass = Nsaba()
        self.assertEquals(nsaba_pass.get_ge(['2', '733', '-1', '1', 2]), 0)
        self.assertEquals(len(nsaba_pass.ge), 3)
        self.assertEquals(len(nsaba_pass.ge['1']), 4)
        self.assertEquals(len(nsaba_pass.ge['733'][0]), 893)
        self.assertEquals(nsaba_pass.ge['1'][3], 2)

if __name__ == '__main__':
    unittest.main()