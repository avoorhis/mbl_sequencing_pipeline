# Copyright (C) 2011, Marine Biological Laboratory
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.
#

from pipeline.sample import Sample

class RunConfig:
    """Doc string here."""
    def __init__(self, config_file_path = None):
        self.config_file_path = config_file_path

        self.run_date   = None
        self.platform   = None # enum('454','illumina','ion_torrent','')
        self.input_dir  = None
        self.output_dir = None
        self.sff_files  = []

        self.run_keys = []
        self.samples = {}

        if config_file_path:
            self.configFromFile(config_file_path)

    def configFromFile(self, config_file_path):
        import ConfigParser
        
        user_config = ConfigParser.ConfigParser()
        user_config.read(config_file_path)

        # take care of the general section
        G = lambda v: user_config.get('general', v)
        self.run_date   = G('run_date')
        self.platform   = G('platform')
        self.input_dir  = G('input_dir')
        self.output_dir = G('output_dir')

        self.input_files  = [file.strip() for file in G('input_files').split(',')]
        self.input_file_type = G('input_file_type')
 
        # populate sample information for every run_key
        for run_key in [s for s in user_config.sections() if s != 'general']:
            #print run_key    # looks like:  1:ACACT
            S = lambda v: user_config.get(run_key, v)
            sample = Sample(run_key)
            
            # has defaults -not required
            try:
                sample.proximal_primers = S('forward_primers').strip("'").strip('"').split(',')
            except:
                sample.proximal_primers = []
            try:
                sample.distal_primers = S('reverse_primers').strip("'").strip('"').split(',')
            except:
                sample.distal_primers = []
            try:
                sample.stop_sequences = S('stop_sequences').strip("'").strip('"').split(',')
            except:
                sample.stop_sequences = []
            try:
                sample.anchor = S('anchor')
            except:
                sample.anchor = ''
            # required
            sample.direction = S('direction')
            sample.project = S('project_name')
            sample.dataset = S('dataset_name')
            sample.dna_region = S('dna_region')
            sample.taxonomic_domain = S('taxonomic_domain')
            
            # a list of run_keys
            # convert: change ':' to '_'
            key = run_key[:1]+'_'+run_key[2:]
            self.run_keys.append(key)
            # a dictionary of samples
            self.samples[key] = sample


