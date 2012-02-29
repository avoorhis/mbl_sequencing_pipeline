#!/usr/local/www/vamps/software/python/bin/python
# -*- coding: utf-8 -*-

# Copyright (C) 2011, Marine Biological Laboratory
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.
#

import unittest
from pipeline.run import Run
from pipeline.trim_run import TrimRun
import pipeline.utils as utils

class unittestargs:
    def __init__(self):
        self.VERBOSE = False
        self.QUIET = False


       
class TestTrimFunctions(unittest.TestCase):

    def setUp(self):
        """

        set up data used in the tests.

        setUp is called before each test function execution.

        """
        
        args = unittestargs()
        run = Run('./unittesting/unittest10.ini')
        self.unittest_trimrun = TrimRun(run,'./unittesting',args)
        trim_code = self.unittest_trimrun.trimrun()
        self.id1 = 'GFADN4I02JYSYM' 
        self.seq1 = 'CGTGATGGGCGTAAAGTGGGTTTAAAGGGTGCGTAGGCGGATTTATAAGTCAGTGGTGAAAGCCTGCGGCTCAACCGTAG\
AATTGCCATTGATACTGTAGATCTTGAGTGCAATCGAGGTGGTTGGAATACGTAGTGTAGCTGTGAAATGCATAGATATTACGTAGAACACC\
AATTGCGAAGGCAGATCACTAGATTGTAACTGACGCTGAGGCACGAAAGCGTGGGGAGCGAACAGGATTAGATACCCTGGTAGTCCACGCCG\
TAAACGATGGATACTAGGTGTTCGGGATAATTGAGTCCTGAGTGCCCAAGCGAAAGCGATAAGTATCCCACCTGGGGAGTACGTCCGCAAGG\
ATGAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGATACGCGAGGAACCTTACCTGGGCTAAAT\
GTTAGTGACGTACCGTGAAAGCGGTATTTCTTCGGACACGAAACA'
        self.trim_info1 = self.unittest_trimrun.do_trim(self.id1,self.seq1)
        
        # proximal deleted:
        self.id2 = 'GFADN4I02HD8KF'
        self.seq2 = 'TGGGCCAACGCGAAAAACCTTACCTCGGAGGTTACCCCCTTCGGTGCCGCAGCTAACGCATTAAGTACTCCG\
CCTGGGGAGTACGCACGCAAGTGTGAAACTCAAAGGAATTGACGGGGACCCGCACAAGTAGCGGAGCATGTGGTTTAATTCGAAGCAACGCGAA\
GAACCTTACCTAGGCTTGACATCCTTCTGACCGAGGATTAATCTCCTCAGGTGGTGCATGGGCTGTCG'
        self.trim_info2 = self.unittest_trimrun.do_trim(self.id2,self.seq2)
        
        # deleted for key:
        self.id3 = 'GFADN4I02HD8KF'
        self.seq3 = 'TGGCCAACGCGAAAAACCTTACCTCGGAGGTTACCCCCTTCGGTGCCGCAGCTAACGCATTAAGTACTCCG\
CCTGGGGAGTACGCACGCAAGTGTGAAACTCAAAGGAATTGACGGGGACCCGCACAAGTAGCGGAGCATGTGGTTTAATTCGAAGCAACGCGAA\
GAACCTTACCTAGGCTTGACATCCTTCTGACCGAGGATTAATCTCCTCAGGTGGTGCATGGGCTGTCG'
        self.trim_info3 = self.unittest_trimrun.do_trim(self.id3,self.seq3)
        
        # deleted for N:
        self.id4 = 'GFADN4I02JYSYM'
        self.seq4 = 'CGTGATGGGCGTAAAGTGGGTTTAAAGGGTGCGTAGGCGGATTTATAAGTCAGTGGTGAAAGCCTGCGGCTCAACCGTAG\
AATTGCCATTGATACTGTAGATCTTGAGTGCAATCGAGGTGGTTGGAATACGTAGTGTAGCTGTGAAATGCATAGATATTACGTAGAACACC\
AATTGCGAAGGCAGATCACTAGATTGTAACTGACGCTGAGGCACGAAAGCGTGGGGAGCGAACAGGATTAGATACCCTGGTAGTCCACGCCG\
TAAACGATGGATACTAGGTGTTCGGGATAATTGAGTCCNTGAGTGCCCAAGCGAAAGCGATAAGTATCCCACCTGGGGAGTACGTCCGCAAGG\
ATGAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGATACGCGAGGAACCTTACCTGGGCTAAAT\
GTTAGTGACGTACCGTGAAAGCGGTATTTCTTCGGACACGAAACA'
        self.trim_info4 = self.unittest_trimrun.do_trim(self.id4,self.seq4)
        
        # deleted for distal:
        self.id5 = 'GFADN4I02I64C9'
        self.seq5 = 'TGGGCTGGGCGTAAAGCGAAGGGTGCGTAGGCGGATTTATAAGTCAGTGGTGAAAGCCTGCGGCTCAACCGTAGAATTG\
CCATTGATACTGTAGATCTTGAGTGCAATCGAGGTTTTGGTTGGAATACGTAGTGTAGCGGTGAAATGCATAGATATTACGTAGAACACCAATTGCGAAGGCA\
GATCACTAGATTGTAACTGACGCTGAGGCACGAAAGCGTGGGGAGCGAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGGATACTAGGTGT\
TCGGGATAATTGAGTCCTGAGTGCCCAAGCGAAAGCGATAAGTATCCCACCTGGGGAGTACGTCCGCAAGGATGAAACTCAAAGGAATTGACGGGGGCCC\
GCACAAGCGGTGGAGCATGTGGTTTAATTCGATGATACGCGAGGAACCTTACCTGGGCTAAATGTTAGT'
        self.trim_info5 = self.unittest_trimrun.do_trim(self.id5,self.seq5)
        
        
        
    def test_lane_key(self):        
        self.assertEqual(self.trim_info1['lane_key'],'1_CGTGA')
    
    def test_primer_name(self):        
        self.assertEqual(self.trim_info1['primer_name'],'565F')
        
    def test_eaxct_right(self):        
        self.assertEqual(self.trim_info1['exact_right'],'AGGTGGT')
        
    def test_exact_left(self):        
        self.assertEqual(self.trim_info1['exact_left'],'TGGGCGTAAAG')
        
    def test_deleted(self):        
        self.assertEqual(self.trim_info1['deleted'],False)  
        
    def test_orientation(self):        
        self.assertEqual(self.trim_info1['orientation'],'+')     
    
    ####################################################################
    def test_deleted_proximal(self):        
        self.assertEqual(self.trim_info2['deleted'],True)  
        self.assertEqual(self.trim_info2['delete_reason'],'proximal')
        
    def test_deleted_key(self):        
        self.assertEqual(self.trim_info3['deleted'],True)  
        self.assertEqual(self.trim_info3['delete_reason'],'runkey')   
        
    def test_deleted_N(self):        
        self.assertEqual(self.trim_info4['deleted'],True)  
        self.assertEqual(self.trim_info4['delete_reason'],'N')
        
    def test_deleted_distal(self):    
        print self.trim_info5
        self.assertEqual(self.trim_info5['deleted'],True)  
        self.assertEqual(self.trim_info5['delete_reason'],'distal') 
        
        
if __name__ == '__main__':
    #unittest.main()
    suite = unittest.TestLoader().loadTestsFromTestCase(TestTrimFunctions)
    unittest.TextTestRunner(verbosity=2).run(suite)