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
from pipeline.runconfig import RunConfig

class Run(RunConfig):
    def __init__(self, config_file_path = None):
        RunConfig.__init__(self, config_file_path)

