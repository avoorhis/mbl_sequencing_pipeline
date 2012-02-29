       
from subprocess import call
import sys, os
import time
import constants as C

class Vamps:
    """Doc string here.."""
    Name = "VAMPS"
    def __init__(self, run = None, outdir= None, args = None):       
        self.run = run
        self.outdir = outdir
        self.rundate = self.run.run_date
        self.outdir    = outdir
        self.QUIET     = args.QUIET
        self.VERBOSE   = args.VERBOSE
        if self.QUIET: self.VERBOSE = False
        
    def info(self):
        pass
    def projects(self):
        pass
    def taxonomy(self):
        pass
    def sequences(self):
        pass        
    def exports(self):
        pass