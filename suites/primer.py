# -*- coding: utf-8 -*-
# name
# sequence
# direction
# domain
# dna_region
import pipeline.constants as C
import pipeline.primer_utils as putils


class Primer:
    """Doc string here.."""
    Name = "Primer"
    def __init__(self,name,dir,domain,region,seq):
        self.name = name
        self.direction = dir
        self.original_seq = seq
        self.domain = domain
        self.region = region
        self.expanded_seqs = []
        self.primer_names = {}
        self.expanded_seqs = putils.expand(self.original_seq)
        
	    
class PrimerSuite:
    """Doc string here.."""
    Name = "PrimerSuite"
    def __init__(self,domain,region):
 
        self.domain = domain
        self.region = region
        
        # list of included primer classes
        self.primer_list ={} 
        self.primer_list['F']=[]        
        self.primer_list['R']=[]
        
        # list of primer sequences
        self.primer_seq_list={}
        self.primer_seq_list['F'] = []
        self.primer_seq_list['R'] = [] 
        # list of expanded primer sequences
        self.primer_expanded_seq_list={}
        self.primer_expanded_seq_list['F'] = []
        self.primer_expanded_seq_list['R'] = []        
        self.primer_names ={}
        # changes Bacteria to Bacterial for the name
        if self.domain[-1:] != 'l': self.domain = self.domain + 'l'
        self.name = self.domain+':'+self.region
        
        from ConfigParser import SafeConfigParser

        parser = SafeConfigParser()
        
        #
        # Reads default primer FILE listed in constants.py
        #
#         parser.read(C.default_primers_file)
#         for section_name in parser.sections():
#             if section_name == self.name:
#                 for name, value in parser.items(section_name):
#                     # for each PrimerSuite instance
#                     # we want to make a list of names,dirs and seqs
#                     #
#                     direction = value[:1]
#                     sequence  = value[2:]
#                     #self.names.append(name)
#                     
#                     #self.dirs.append(direction)
#                     #self.orig_seqs.append(sequence)
#                     #self.expanded_seqs[name]={}
#                     #self.expanded_seqs[value[:1]] = expand(value[2:])
#                     p = Primer(name,direction,self.domain,self.region,sequence)
#                     self.primer_list[direction].append(p)
#                     self.primer_seq_list[direction].append(sequence)
#                     self.primer_expanded_seq_list[direction] = self.primer_expanded_seq_list[direction] + p.expanded_seqs
#                     #self.primer_names.append(p.primer_names)
#                     for eseq in p.expanded_seqs:
#                         self.primer_names[eseq] = p.name
                        
                        
        #
        # Reads default primers from dictionary in constants.py
        #
        for suite in C.mbl_primer_suites:
            if suite == self.name:
                
                for pname in C.mbl_primer_suites[suite]:
                    #print pname
                    direction = C.mbl_primer_suites[suite][pname]['direction']
                    sequence  = C.mbl_primer_suites[suite][pname]['sequence']
                    domain    = C.mbl_primer_suites[suite][pname]['domain']
                    region    = C.mbl_primer_suites[suite][pname]['region']
                    #self.names.append(name)
                    
                    p = Primer(pname,direction,domain,region,sequence)
                    self.primer_list[direction].append(p)
                    self.primer_seq_list[direction].append(sequence)
                    self.primer_expanded_seq_list[direction] = self.primer_expanded_seq_list[direction] + p.expanded_seqs
                    #self.primer_names.append(p.primer_names)
                    for eseq in p.expanded_seqs:
                        self.primer_names[eseq] = pname
                
            
            
            
            
            