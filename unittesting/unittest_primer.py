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

import random
import unittest
from   suites.primer import PrimerSuite
import pipeline.utils as utils
import pipeline.primer_utils as putils
import constants as C

class TestPrimerFunctions(unittest.TestCase):

    def setUp(self):
        """

        set up data used in the tests.

        setUp is called before each test function execution.

        """
        self.psuite = PrimerSuite('Bacterial','v6v4')

    def test_domain(self):     
        self.assertEqual(self.psuite.domain, 'Bacterial')
        
    def test_region(self): 
        self.assertEqual(self.psuite.region, 'v6v4')
        
    def test_name(self):
        self.assertEqual(self.psuite.name, 'Bacterial:v6v4')
        
    def test_forward_primers(self):
        self.assertEqual(self.psuite.primer_seq_list['F'],['TGGGCGTAAAG'])
        
    def test_reverse_primers(self):
        self.assertEqual(sorted(self.psuite.primer_seq_list['R']),sorted(['AGGTG.TGCATGGCTGTCG','AGGTG.TGCATGGTTGTCG','AGGTG.TGCATGGCCGTCG','AGGTG.TGCATGGTCGTCG']))
    
    def test_revcomp(self):
        self.assertEqual(utils.revcomp('TGGGCGTAAAG'),'CTTTACGCCCA')
    
    
    # Testing the expand primers def
    def test_expand_primer(self):
        self.assertEqual(sorted(putils.expand('TG[GA]GNGTAA?AG')),sorted(['TGGGGGTAAG',  'TGGGAGTAAAG', 
                                                                          'TGGGCGTAAG',  'TGGGTGTAAAG', 
                                                                          'TGGGGGTAAAG', 'TGAGAGTAAAG', 
                                                                          'TGAGGGTAAG',  'TGAGCGTAAG', 
                                                                          'TGAGTGTAAAG', 'TGAGTGTAAG', 
                                                                          'TGGGCGTAAAG', 'TGGGTGTAAG', 
                                                                          'TGGGAGTAAG',  'TGAGAGTAAG', 
                                                                          'TGAGCGTAAAG', 'TGAGGGTAAAG']))
    def test_expand_dot(self):
        self.assertEqual( sorted(putils.expand('T.G')),     sorted( ['TAG','TTG','TGG','TCG']) )
    def test_expand_dot_qmark(self):
        self.assertEqual( sorted(putils.expand('T.?G')),    sorted( ['TG','TAG','TTG','TGG','TCG']) )
    def test_expand_N(self):
        self.assertEqual(sorted(putils.expand('TNG')),      sorted( ['TAG','TTG','TGG','TCG']) )  
    def test_expand_R(self):
        self.assertEqual( sorted(putils.expand('TRG')),     sorted( ['TAG','TGG']) )  
    def test_expand_Y(self):
        self.assertEqual( sorted(putils.expand('TYG')),     sorted( ['TCG','TTG']) )     
    def test_expand_W(self):
        self.assertEqual( sorted(putils.expand('TWG')),     sorted( ['TAG','TTG']) )     
    def test_expand_S(self):
        self.assertEqual( sorted(putils.expand('TSG')),     sorted( ['TGG','TCG']) )
    def test_expand_M(self):
        self.assertEqual( sorted(putils.expand('TMG')),     sorted( ['TAG','TCG']) )
    def test_expand_K(self):
        self.assertEqual( sorted(putils.expand('TKG')),     sorted( ['TTG','TGG']) )
    def test_expand_Lbracket(self):
        self.assertEqual( sorted(putils.expand('T[TG]G')),  sorted( ['TTG','TGG']) ) 
    def test_expand_qmark(self):
        self.assertEqual( sorted(putils.expand('TT?G')),    sorted( ['TTG','TG']) )
    def test_expand_plus(self):
        self.assertEqual( sorted(putils.expand('TT+G')),    sorted( ['TTG','TTTG']) )
    def test_expand_star(self):
        self.assertEqual( sorted(putils.expand('TT*G')),    sorted( ['TG','TTG','TTTG']) ) 
    def test_expand_Lbrace(self):
        self.assertEqual( sorted(putils.expand('TC{2,5}G')),sorted( ['TCCG','TCCCG','TCCCCG','TCCCCCG']) )
        
        
if __name__ == '__main__':
    #unittest.main()
    suite = unittest.TestLoader().loadTestsFromTestCase(TestPrimerFunctions)
    unittest.TextTestRunner(verbosity=2).run(suite)