# -*- coding: utf-8 -*-
#
# Copyright (C) 2011, Marine Biological Laboratory
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.
#

class Sample:
    """Doc string here.."""
    def __init__(self, run_key):
        self.run_key = run_key
        self.primers = {}
        self.anchor = None
        self.pool = None
        self.lane = None
        self.direction = None # enum('B','F','F RC','R')
        self.adaptor = None # enum('A','B')
        self.dna_region = None
        self.taxonomic_domain = None
        self.project = None
        self.dataset = None
        

