
import os
import sys

Map = {'A': 'T',
       'T': 'A',
       'C': 'G',
       'G': 'C',
       'N': 'N'}
       
class lineReader:
    pass





class FastaReader:
    def __init__(self,file_name=None):
        self.file_name = file_name
        self.h = open(self.file_name)
        self.seq = ''
        self.id = None
        self.revcomp_seq = None
        self.base_counts = None

    def next(self): 
        def read_id():
            return self.h.readline().strip()[1:]

        def read_seq():
            ret = ''
            while True:
                line = self.h.readline()
                
                while len(line) and not len(line.strip()):
                    # found empty line(s)
                    line = self.h.readline()
                
                if not len(line):
                    # EOF
                    break
                
                if line.startswith('>'):
                    # found new defline: move back to the start
                    self.h.seek(-len(line), os.SEEK_CUR)
                    break
                    
                else:
                    ret += line.strip()
                    
            return ret
        
        self.id = read_id()
        self.seq = read_seq()
        self.revcomp_seq = self.rev_comp()
        self.base_counts = self.counts()
        
        if self.id:
            return True        

    def counts(self):
        self.cnt = {}
        for b in self.seq:           
            self.cnt[b] = self.cnt.get(b,0) +1
        return self.cnt

  
    def rev_comp(self):
        return ''.join(reversed([Map[x] for x in self.seq]))

        
