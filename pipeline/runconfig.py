# Copyright (C) 2011, Marine Biological Laboratory
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.
#

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
        self.sff_files  = [sff.strip() for sff in G('sff_files').split(',')]

        # populate sample information for every run_key
        for run_key in [s for s in user_config.sections() if s != 'general']:
            S = lambda v: user_config.get(run_key, v)
            sample = Sample(run_key)
            sample.forward_primer = S('forward_primer')
            sample.distal_primer = S('distal_primer')
            sample.direction = S('direction')
            # (...)

            self.run_keys.append(run_key)
            self.samples[run_key] = sample


