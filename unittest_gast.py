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
from   chimera import Chimera
import pipeline.utils as utils
import constants as C

class TestChimeraFunctions(unittest.TestCase):

    def setUp(self):
        """

        set up data used in the tests.

        setUp is called before each test function execution.

        """
        self.chimera = Chimera('Bacterial','v6v4')

    def test_chimera(self):     
        pass
        
        
if __name__ == '__main__':
    #unittest.main()
    suite = unittest.TestLoader().loadTestsFromTestCase(TestChimeraFunctions)
    unittest.TextTestRunner(verbosity=2).run(suite)