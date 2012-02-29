#!/usr/local/www/vamps/software/python/bin/python

##!/usr/bin/env python
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
import os
import shutil

from suites.primer import PrimerSuite
from pipeline.primer_utils import *
from pipeline.utils import *
from pipeline.fastalib import *
from pipeline.Fasta import sfasta
from pipeline.anchortrimming_mbl import *

class TrimRun( object ):
    """ Define here"""
    def __init__(self, run = None, outdir= None, args = None):
    
        self.run       = run
        self.outdir    = outdir
        self.args       = args
        self.QUIET     = args.QUIET
        self.VERBOSE   = args.VERBOSE
        if self.QUIET: self.VERBOSE = False        
        
        # do something with 'run'.
        self.rundate = self.run.run_date
        if self.VERBOSE: print "Rundate:",self.rundate
        
        # sff or fasta or fastq ?
        self.input_file_type = self.run.input_file_type.strip("'").strip('"')
        if(self.input_file_type == 'sff'): self.input_file_type = 'sff-trim'
        self.input_files = self.run.input_files
        self.run_keys = [key[2:] for key in self.run.run_keys]
      
    
        #self.projects           = {}
        #self.datasets           = {}
        self.seqDirs            = {}
        self.dna_regions        = {}
        self.taxonomic_domain   = {}
        self.expanded_primers   = {}
        self.proximal_primers   = {}
        self.distal_primers     = {}
        
        self.anchor_name            = {}
        self.adtnl_anchors ={}
        
        self.anchors ={}
        
        
        self.psuite             = {}
        #self.id_list_all        = {}
        self.id_list_passed     = {}
        #self.unique_fasta_files = {}
        #self.names_files        = {}
        self.deleted_ids        = {}
        self.trimmed_ids        = {}
        self.uniques            = {}
        self.names              = {}
        self.uniquefa ={}
        self.abundfa ={}
        self.fa ={}
        self.statsFileName     = 'run_trim_stats'
        
                   
        self.runbin={}
        #if VERBOSE: print "Run Keys from .ini file: ", run_keys
        for lane_key in self.run.run_keys:
            # strip off surrounding single or double quotes
            #self.projects[lane_key]         = self.run.samples[lane_key].project.strip("'").strip('"')
            #self.datasets[lane_key]         = self.run.samples[lane_key].dataset.strip("'").strip('"')
            self.seqDirs[lane_key]          = self.run.samples[lane_key].direction.strip("'").strip('"')
            self.dna_regions[lane_key]      = self.run.samples[lane_key].dna_region.strip("'").strip('"')
            self.taxonomic_domain[lane_key] = self.run.samples[lane_key].taxonomic_domain.strip("'").strip('"')
            
            # this should be defaiult 'suite' of anchors
            # but also in ini file: anchor=XXXXX will be added to defaults
            self.anchor_name[lane_key]         = self.run.samples[lane_key].anchor.strip("'").strip('"')
            self.adtnl_anchors[lane_key]         = self.run.samples[lane_key].stop_sequences  #list
           
            self.anchors[lane_key]       = {}
            
            #self.id_list_all[lane_key]      = []
            self.id_list_passed[lane_key]   = []
            self.deleted_ids[lane_key]      = {}
            self.deleted_ids['nokey']       = {}
            self.trimmed_ids[lane_key]      = {}
            self.uniques[lane_key]          = {}
            self.names[lane_key]            = {}
            self.fa[lane_key]       = FastaOutput(self.outdir + '/' + lane_key + ".trimmed.fa")
            
            
            #####################
            #
            #  PrimerSuite Class
            #
            #####################
            self.psuite[lane_key] = PrimerSuite(self.taxonomic_domain[lane_key],self.dna_regions[lane_key])
            #self.runbin['psuite'][lane_key]= PrimerSuite(self.taxonomic_domain[lane_key],self.dna_regions[lane_key])
                
            if(self.seqDirs[lane_key] == 'F' or self.seqDirs[lane_key] == 'B'):
                self.proximal_primers[lane_key] = self.psuite[lane_key].primer_expanded_seq_list['F']
                self.distal_primers[lane_key]   = self.psuite[lane_key].primer_expanded_seq_list['R']
                if self.anchor_name[lane_key]:
                    self.anchors[lane_key] = get_anchor_list(self.anchor_name[lane_key], self.adtnl_anchors[lane_key])
 
                
            if(self.seqDirs[lane_key] == 'R' or self.seqDirs[lane_key] == 'B'):
                self.proximal_primers[lane_key] = revcomp( self.psuite[lane_key].primer_expanded_seq_list['F'] )
                self.distal_primers[lane_key]   = revcomp( self.psuite[lane_key].primer_expanded_seq_list['R'] ) 
                
                #self.proximal_primers[lane_key] = revcomp( self.psuite[lane_key].primer_expanded_seq_list['F'] )
                #self.distal_primers[lane_key]   = revcomp( self.psuite[lane_key].primer_expanded_seq_list['R'] )   
                if self.anchor_name[lane_key]:
                    self.anchors[lane_key] = revcomp( get_anchor_list(self.anchor_name[lane_key], self.adtnl_anchors[lane_key]) )
                    
            #print lane_key,revcomp(psuite[lane_key].primer_seq_list['F'])
            # this should be a list of primers, but it may be entered as a filename
            # the primer direction preceeds the sequence (colon sep)
            # ie:  F:CTAACCGA.GAACCT[CT]ACC
            
 
            
        # next we need to start reading sequences
        # 
        # END of Initialization
        #
        
    def trimrun(self, write_files = False):
        
        # use BioPython to read sff file or fasta file
        # depending on switch in config file: self.input_file_type
        from Bio import SeqIO
        #from Bio.Alphabet import DNAAlphabet
        #from Bio.Seq import Seq
        # 
        #    
    
        self.stats_fp = open( self.outdir+'/'+self.statsFileName,"w" )
        self.number_of_raw_sequences  = 0
        self.number_of_good_sequences = 0
        success_code = ()
        self.deleted_count = {}
        self.deleted_count['Runkey'] =0
        self.deleted_count['Proximal'] =0
        self.deleted_count['Distal'] =0
        self.deleted_count['N'] =0
        self.deleted_count['Quality'] =0
        self.deleted_count['No Insert'] =0
        self.deleted_count['Minimum Length'] =0
        # input type is defined in the config file and can be sff, fasta or fastq
        # need to figure out how if to get quality data with fasta    
        # save a list of read_ids: all_ids, passed_ids
        for file in self.input_files:  # usually just one file: a list of one?
            
            for record in SeqIO.parse(file, self.input_file_type): 
            
                self.number_of_raw_sequences += 1
                id = record.id
                seq = record.seq.tostring().upper()
                if(self.input_file_type == 'fasta'):
                    q_scores = ''
                else:
                    q_scores = record.letter_annotations["phred_quality"]  
                
    
                #==========================================
                #
                # here each sequence is trimmed
                #
                #==========================================
                trim_data = self.do_trim(id, seq, q_scores)
                #
                #
                #
                #
                
                
                
                
                deleted         = trim_data['deleted']
                lane_tag        = trim_data['lane_key']
                seq             = trim_data['trimmed_sequence']        
                exact_left      = trim_data['exact_left']   
                exact_right     = trim_data['exact_right']        
                delete_reason   = trim_data['delete_reason']
                primer_name     = trim_data['primer_name']
                
                ###################################################################
                #
                #   9-print out file of trimmed sequence
                #
                ###################################################################
                # output fasta file for each lane_tag: RAW or TRIMMED for chimera checking?
                # these will be input files for chimera checking
                # write to fasta file: need fasta 
     
                
                # print out a fasta file for each lane_tag
                
                if(lane_tag and not deleted):
                    
                    # for the names file of unique ids
                    try:
                        self.names[lane_tag][self.uniques[lane_tag][seq]].append(id)
                    except KeyError:
                        self.names[lane_tag][id] = [id]
                    
                    # hash for uniques -- this will use the first found as key
                    if seq not in  self.uniques[lane_tag]:
                        self.uniques[lane_tag][seq] = id
                    
    
                    if write_files:        
                        #fa = FastaOutput(read_id, trimmed_sequence, ('|').join(['project:'+projects[lane_tag],'dataset:'+datasets[lane_tag]]))
                        self.fa[lane_tag].write_id(id)
                        self.fa[lane_tag].write_seq(seq) 
                    
                   
                    
                    self.number_of_good_sequences += 1
                    if self.VERBOSE: print record.id,"Passed"
                elif(lane_tag and deleted):
                    self.deleted_ids[lane_tag][id] = delete_reason
                    if self.VERBOSE:
                        print id,"1-Deleted:",delete_reason
                else:
                    if self.VERBOSE:
                        print id,"2-Deleted:",delete_reason
                        
                        
        # after reading sff files                
                        
        count_uniques = 0
        good_lane_keys = []
        for lane_key in self.run.run_keys:
            count = count_keys(self.uniques[lane_key])
            if count > 0:
                good_lane_keys.append(lane_key)
            count_uniques = count_uniques + count               
        if not success_code:
            success_code = ('SUCCESS','',good_lane_keys)
        return success_code
 
    ###################################################################
    #
    #   do_trim - Trims each sequence
    #
    ################################################################### 
 
 
    def do_trim(self, read_id, raw_sequence, quality_scores=None):
    
        trim_collector = {}    
        
        
        if self.VERBOSE: print "\nTrimming read",   read_id
        
        
        
        trimmed_sequence = raw_sequence
        
        deleted = False
        delete_reason = ''
        
    

        if self.VERBOSE: print read_id,'RAW:',raw_sequence
        ###################################################################
        #
        #   1-check/count for Ns
        #
        ###################################################################
        # save count for later
        countNs = check_for_Ns(raw_sequence)
        

            
        ###################################################################
        #    
        #   3-find/remove runkey
        #
        ###################################################################
        
        tag      = ''
        lane     = ''
        lane_tag = ''
        tag, trimmed_sequence = remove_runkey(raw_sequence, self.run_keys)

        if( not tag ):
            deleted = True
            delete_reason = 'runkey'
            self.deleted_count['Runkey'] += 1
            self.deleted_ids['nokey'][read_id] = delete_reason
            if self.VERBOSE: print "deleted: tag"
        else:
            lane = self.run.run_keys[self.run_keys.index(tag)][:1]
            lane_tag = lane+'_'+tag
            if self.VERBOSE: print "lane and tag found:",lane_tag
            
        ###################################################################
        #
        #   4-determine read direction 
        #
        ###################################################################

        seq_direction = ''
        #print seqDirs[lane_tag]
        if( not deleted ):
            seq_direction = find_sequence_direction( self.seqDirs[lane_tag] )
            
        ###################################################################
        #
        #   5-find primers
        #     defs:
        #       proximal === closest to runkey
        #       distal   === furthest from run_key
        #       forward 'F' === sequences was read R->L  ie v3v5
        #       reverse 'R' === sequence was read L->R   ie v6v4
        ###################################################################
        exactRight = ''
        exactLeft=''
        orientation=''
        primer_name=''
        if seq_direction != 'F' and seq_direction != 'B' and seq_direction != 'R':
            print 'Unknown seq direction:',seq_direction
        else:           
            if seq_direction == 'F' or seq_direction == 'B':
                
                ###################################################################
                #
                #   5-find/remove proximal primer -exactLeft for forward reads
                #
                ###################################################################
                # proximal primer - exactLeft
                exactLeft, offset, trimmed_sequence \
                      = trim_proximal_primer( self.proximal_primers[lane_tag], trimmed_sequence, 'F' )
                primer_name = ''
                #primer_name = self.psuite[lane_tag].primer_names[exactLeft]
                if ( exactLeft and (seq_direction == 'F' or seq_direction == 'B' )):  
                    orientation = '+'
                    primer_name = self.psuite[lane_tag].primer_names[exactLeft] 
                    
                if( not exactLeft ):
                    deleted = True
                    if self.VERBOSE: print 'deleted: proximal'                        
                    if(not delete_reason):  
                        delete_reason = 'proximal'
                        self.deleted_count['Proximal'] += 1
                    #self.deleted_ids[lane_tag][self.read_id] = delete_reason
                else:
                    
                    if self.VERBOSE: print 'ExactLeft:',exactLeft, offset, orientation, primer_name
                    if offset:
                        primer_name = "offset " + primer_name
                ###################################################################
                #
                #   5-find/remove anchor (if present)
                #
                ###################################################################
                # anchor
                if self.anchor_name[lane_tag]:
                    print 'Have Anchor name: ',self.anchor_name[lane_tag]                    
                    exactRight, exactTrimmedOff, trimmed_sequence \
                                = trim_anchor( self.anchor_name[lane_tag],self.anchors[lane_tag], trimmed_sequence, 'F' )
                    if exactRight: 
                        anchor_found = True 
                        print 'found exactRight-anchor',exactRight
                ###################################################################
                #
                #   5-find/remove distal primer exactRight or fuzzy
                #
                ###################################################################
                # distal primer - exactRight or fuzzyRight
                if not exactRight:               
                    exactRight, exactTrimmedOff, trimmed_sequence \
                                = trim_distal_primer( self.distal_primers[lane_tag], trimmed_sequence, 'F' )
                    if exactRight: 
                        exact_found = True 
                        print 'found exactLeft-primer',exactLeft
                        
                if not exactRight:
                    deleted = True
                    if self.VERBOSE: print 'deleted: distal'                        
                    if(not delete_reason):  
                        delete_reason = 'distal'
                        self.deleted_count['Distal'] += 1
                        
            if seq_direction == 'R' or seq_direction == 'B':
                pass
                
                ###################################################################
                #
                #   5-find/remove distal primer -exactRight for reverse reads
                #
                ###################################################################
                # proximal primer - exactRight
                exactRight, offset, trimmed_sequence \
                      = trim_proximal_primer( self.proximal_primers[lane_tag], trimmed_sequence, 'R' ) 
                
                if ( exactRight and (seq_direction == 'R' or seq_direction == 'B' )):  
                    orientation = '+'
                    primer_name = self.psuite[lane_tag].primer_names[revcomp(exactRight)]
                    
                if( not exactRight ):
                    deleted = True
                    if self.VERBOSE: print 'deleted: proximal'                        
                    if(not delete_reason):  
                        delete_reason = 'proximal'
                        self.deleted_count['Proximal'] += 1
                    #self.deleted_ids[lane_tag][self.read_id] = delete_reason
                else:
                    
                    if self.VERBOSE: print 'ExactRight:',exactRight, offset, orientation, primer_name
                    if offset:
                        primer_name = "offset " + primer_name
                ###################################################################
                #
                #   5-find/remove anchor 
                #
                ###################################################################
                 # anchor   
                if self.anchor_name[lane_tag]:
                    print 'Have Anchor name: ',self.anchor_name[lane_tag]                    
                    exactLeft, exactTrimmedOff, trimmed_sequence \
                                = trim_anchor( self.anchor_name[lane_tag],self.anchors[lane_tag], trimmed_sequence, 'R' )
                    if exactLeft: 
                        anchor_found = True 
                        print 'found exactLeft-anchor',exactLeft
                ###################################################################
                #
                #   5-find/remove proximal primer exactLeft or fuzzy
                #
                ###################################################################
                # distal primer - exactLeft or fuzzyLeft
                if not exactLeft:               
                    exactLeft, exactTrimmedOff, trimmed_sequence \
                                = trim_distal_primer( self.distal_primers[lane_tag], trimmed_sequence, 'R' )
                    if exactLeft: 
                        exact_found = True
                        print 'found exactLeft-primer',exactLeft
                
                if not exactLeft:
                    deleted = True
                    if self.VERBOSE: print 'deleted: distal'                        
                    if(not delete_reason):  
                        delete_reason = 'distal'
                        self.deleted_count['Distal'] += 1
        ###################################################################
        #
        #   5-find/remove proximal primer
        #
        ###################################################################
        
#         offset = 0
#         exactLeft = 0
#         primer_name = ''
#         orientation = ''
#         
#         if( not deleted ):
#             #print forward_primers[lane_tag]
#             
#             exactLeft, offset, trimmed_sequence \
#                       = trim_proximal_primer( self.proximal_primers[lane_tag], trimmed_sequence,'F')
#                       
#             if ( exactLeft and (seq_direction == 'F' or seq_direction == 'B' )):  
#                 orientation = '+'
#                 primer_name = self.psuite[lane_tag].primer_names[exactLeft]
#                 
#             if ( exactLeft and (seq_direction == 'R' or seq_direction == 'B') ):
#                 orientation = '-'
#                 primer_name = self.psuite[lane_tag].primer_names[revcomp(exactLeft)]
#                 
#             if( not exactLeft ):
#                 deleted = True
#                 if self.VERBOSE: print 'deleted: proximal'                        
#                 if(not delete_reason):  
#                     delete_reason = 'proximal'
#                     self.deleted_count['Proximal'] += 1
#                 #self.deleted_ids[lane_tag][self.read_id] = delete_reason
#             else:
#                 
#                 if self.VERBOSE: print 'ExactLeft:',exactLeft, offset, orientation, primer_name
#                 if offset:
#                     primer_name = "offset " + primer_name
#             
#             print 'found prox primer',trimmed_sequence
#                 
#                 
#                
#         
#         ###################################################################
#         #
#         #   6-find  & remove distal primer or anchor
#         #
#         ###################################################################
#         
#         distal_primer = ''
#         exactRight   = ''
#         fuzzyRight = ''
#         loc           = 0
#         
#         if( delete_reason != 'runkey'):
#         
#             minLength = C.trim_lengths[self.dna_regions[lane_tag]]['length']
#             
#             if( len(trimmed_sequence) < minLength ):
#                 deleted = True
#                 if (not delete_reason): 
#                     delete_reason = 'minimum length'
#                     self.deleted_count['Minimum Length'] += 1
#                     
#             else:
#                 anchor_start = C.trim_lengths[self.dna_regions[lane_tag]]['start'];
#                 anchor_end   = C.trim_lengths[self.dna_regions[lane_tag]]['end'];
#                 trim_type = C.trim_lengths[self.dna_regions[lane_tag]]['trim_type']  # internal or distal
#                 exact_found  = False
#                 fuzzy_found  = False
#                 anchor_found = False
#                 
#  #                if(self.anchors[lane_tag]):
# #                     # here anchors are revcomed if read seq is reversed
# #                     print 'stop seqs:',self.anchors[lane_tag]
# #                     # for reverse reads anchor_start and anchor end are negative (in constants.py)
# #                     # but we have to reversese complement them
# #                     exactRight, exactTrimmedOff, trimmed_sequence \
# #                                 = trim_stop_seq( self.anchors[lane_tag], trimmed_sequence, trim_type, anchor_start, anchor_end )
# #                     if exactRight: 
# #                         anchor_found = True 
#                 if self.anchor_name[lane_tag]:
#                     print 'Have Anchor name: ',self.anchor_name[lane_tag]
#                     
#                     exactRight, exactTrimmedOff, trimmed_sequence \
#                                 = trim_anchor( self.anchor_name[lane_tag],self.anchors[lane_tag], trimmed_sequence )
#                     if exactRight: 
#                         anchor_found = True 
#                         
#                     
#                 if not exactRight:               
#                     exactRight, exactTrimmedOff, trimmed_sequence \
#                                 = trim_distal_primer( self.distal_primers[lane_tag], trimmed_sequence )
#                     if exactRight: 
#                         exact_found = True  
#                         
#                         
#                 if not exactRight:
#                     fuzzyRight, fuzzyDist, fuzzyTrimmed, fuzzy_match \
#                                 = trim_fuzzy_distal(  self.distal_primers[lane_tag], trimmed_sequence, trim_type, anchor_start, anchor_end)
#                     if exactRight: 
#                         fuzzy_found = True  
#                         
#                 # choose the longest trimmed off section: exactRight vs fuzzyRight
#                 if (len(exactTrimmedOff) < len(fuzzyRight) ):
#                     exactRight = fuzzyRight
#                     trimmed_sequence = fuzzyTrimmed
#                 
#                     
#                 if not exactRight:
#                     if not deleted:
#                         deleted = True
#                         if self.VERBOSE: print 'deleted distal'
#                         if (not delete_reason):  
#                             delete_reason = 'distal'
#                             self.deleted_count['Distal'] += 1
#                         #self.deleted_ids[lane_tag][self.read_id] = delete_reason
#                      
#                 else:                    
#                     if self.VERBOSE and anchor_found: print "ExactRight:", exactRight, "->Found Anchor<-"
#                     if self.VERBOSE and exact_found: print "ExactRight:", exactRight, "->Found Exact<-"
#                     if self.VERBOSE and fuzzy_found: print "ExactRight:", exactRight, "->Found Fuzzy<-"
#                     if self.VERBOSE: print "TrimmedOff:",exactTrimmedOff
#                     if self.VERBOSE: print "Trimmed:", trimmed_sequence
#         ###################################################################
        #
        #   7-check length and other sundry things
        #
        ###################################################################
        if self.VERBOSE: print 'length raw:',len(raw_sequence),'length trimmed',len(trimmed_sequence)
        if( delete_reason != 'runkey' and not trimmed_sequence ):
            deleted = True
            if not delete_reason:  
                delete_reason = 'no insert'
                self.deleted_count['No Insert'] += 1
            #self.deleted_ids[lane_tag][self.read_id] = delete_reason
        elif( delete_reason != 'runkey' and len(trimmed_sequence) < C.minimumLength):
            deleted = True
            if not delete_reason: 
                delete_reason = 'minimum length'
                self.deleted_count['Minimum Length'] += 1
            #self.deleted_ids[lane_tag][read_id] = delete_reason
        
        if ( (not deleted) and (countNs > C.maxN) ):     
            deleted = True
            if (not delete_reason):  
                delete_reason = 'N'
                self.deleted_count['N'] += 1
            #self.deleted_ids[lane_tag][read_id] = delete_reason
            if self.VERBOSE: print 'deleted N'
        
        # save passed ids
        if(not deleted):
            self.id_list_passed[lane_tag].append(read_id)
        
        ###################################################################
        #
        #   8-check quality
        #
        ###################################################################
        if(quality_scores and not deleted):
            average_score = check_for_quality(raw_sequence, trimmed_sequence, quality_scores)
            if average_score < C.minAvgQual:
                deleted = True
                if (not delete_reason):  
                    delete_reason = 'quality'
                    self.deleted_count['Quality'] += 1
                
                if self.VERBOSE: print 'deleted quality'
        
        trim_collector['raw_sequence']      = raw_sequence
        trim_collector['trimmed_sequence']  = trimmed_sequence
        trim_collector['exact_left']        = exactLeft
        trim_collector['exact_right']       = exactRight
        trim_collector['deleted']           = deleted
        trim_collector['delete_reason']     = delete_reason
        trim_collector['primer_name']       = primer_name
        trim_collector['lane_key']          = lane_tag
        trim_collector['orientation']       = orientation
       
        return trim_collector
####################################################################
#
#
#
####################################################################


                
    def write_data_files(self, lane_keys): 
    
        ###################################################################
        #
        #   10-print out files of unique trimmed sequence
        #           and deleted ids
        #    for each amplification:
        #       fasta file          is this needed?
        #       names file          (just like mothur)
        #       uniques fasta file  (just like mothur)
        #       abundance fasta file -formated for usearch (-uc): sorted most abund first
        #       deleted read_ids file
        #   Also there should be one file each for rawaseq, trimseq, primers and runkeys
        #       Not divided by amplification for rapidly puting in db.
        #   These are printed in out directory: './out'+rundate
        #
        ###################################################################  
        #rawseqFileName = self.outdir + '/rawseq_file.txt'
        #trimseqFileName = self.outdir + '/trimseq_file.txt'
        #f_rawseq  = open(rawseqFileName, "w")
        #f_trimseq = open(trimseqFileName,"w")
        for lane_key in self.run.run_keys:
            self.fa[lane_key].close()  
            uniquesFileName = self.outdir + '/' + lane_key + ".unique.fa"
            abundFileName   = self.outdir + '/' + lane_key + ".abund.fa"
            namesFileName   = self.outdir + '/' + lane_key + ".names"
            delFileName     = self.outdir + '/' + lane_key + ".deleted.txt"
            f_names = open(namesFileName,"w") 
            #  if we order the uniques by length of self.uniques[lane_key][seq] then we have abundance file
            
            # Write abund.fa file
            # mysort returns a list of tuples: (read_id, count, seq) sorted highest to lowest freq
            try:
                sorted_uniques = mysort( self.uniques[lane_key], self.names[lane_key] )
                for item in sorted_uniques:
                    read_id = item[0]
                    count = item[1]
                    seq = item[2]
                    
                    sfastaRead = read_id + ";size="+str(count)                
                    abundfa = sfasta(sfastaRead, seq)
                    abundfa.write(abundFileName,'a')

            except:
                success_code = ('FAIL','abund',lane_key)
               
            # Write uniques.fa file
            try:
                for seq in self.uniques[lane_key]:
                    read_id = self.uniques[lane_key][seq]
                    uniquefa = sfasta(read_id, seq)
                    uniquefa.write(uniquesFileName,'a')
                if self.VERBOSE: print "\nwrote uniques file",uniquesFileName
 
            except:
                success_code = ('FAIL','unique',lane_key)    
                
            # Write names file
            try:
                for id in self.names[lane_key]:
                    others = ','.join(self.names[lane_key][id])                
                    f_names.write(id+"\t"+others+"\n")
                f_names.close()
                if self.VERBOSE: print "wrote names file",namesFileName
            except:
                success_code = ('FAIL','names',lane_key) 
                
            # Write deleted.txt file   
            if lane_key in self.deleted_ids and self.deleted_ids[lane_key]:
                f_del   = open(delFileName,  "w") 
                for id in self.deleted_ids[lane_key]:
                    reason = self.deleted_ids[lane_key][id]                
                    f_del.write(id+"\t"+reason+"\n")
                f_del.close()
                if self.VERBOSE: print "wrote deleted file",delFileName
            
            
            
        # print out readids that failed the key test: one file only
        if 'nokey' in self.deleted_ids and self.deleted_ids['nokey']:
            delfileName = self.outdir + '/nokey.deleted.txt'
            f_del = open(delfileName,"w")
            for id in  self.deleted_ids['nokey']:                
                f_del.write(id+"\tnokey\n")
            f_del.close()
            
        if not self.QUIET: 
            print   
            print 'Output Directory:', './'+self.outdir
            print self.number_of_raw_sequences,  "raw sequences read"  
            pct = '%.1f' % ((float(self.number_of_good_sequences)/self.number_of_raw_sequences) *100)
            print self.number_of_good_sequences, "sequences passed" ,pct+'%'
            print "Unique Counts:"
            count_uniques = 0
            good_lane_keys = []
            for lane_key in self.run.run_keys:
                count = count_keys(self.uniques[lane_key])
                if count > 0:
                    good_lane_keys.append(lane_key)
                count_uniques = count_uniques + count
                print "   ",lane_key,self.dna_regions[lane_key],count
            print "   Total Uniques:",count_uniques
 
 
        #####
        #
        #  Write to stats file for this run
        #
        self.stats_fp.write("Run: "+self.rundate+"\n")
        self.stats_fp.write("Unique Counts:\n")
        #stats_fp.write("Run: "+self.rundate)
        count_uniques = 0
        for lane_key in self.run.run_keys:
            count = count_keys(self.uniques[lane_key])
            count_uniques = count_uniques + count
            self.stats_fp.write("   " + str(count)+"\t"+lane_key+"\n")
        
        self.stats_fp.write("Total Uniques: "+str(count_uniques)+"\n")
        self.stats_fp.write("\nDeleted Counts (before chimera check):\n")
        for reason in self.deleted_count:
            self.stats_fp.write("   " + str(self.deleted_count[reason]) +"\t"+reason+"\n")
        self.stats_fp.write("Total Deleted: "+str(self.number_of_raw_sequences-self.number_of_good_sequences)+"\n")    
        self.stats_fp.close()
        
        success_code=''
        if not success_code:
            success_code = ('SUCCESS','',good_lane_keys)
        return success_code
        

    
if __name__=='__main__':
    pass