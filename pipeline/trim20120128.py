
from pipeline.utils import *
from pipeline.primer_utils import *
from suites.primer import PrimerSuite

class Trim:
    def __init__(self, record = None, run = None, args = None):
        self.run = run
        input_file_type = self.run.input_file_type.strip("'").strip('"')

        self.VERBOSE = args.VERBOSE
        self.QUIET = args.QUIET
        self.read_id = record.id
        self.raw_sequence = record.seq.tostring().upper()
        
        if(input_file_type == 'fasta'):
            self.quality_scores = ''
        else:
            self.quality_scores = record.letter_annotations["phred_quality"]  
            
        self.run_keys = [key[2:] for key in self.run.run_keys]
        self.seqDirs            = {}
        self.psuite             = {}
        self.proximal_primers   = {}
        self.distal_primers     = {}
        self.seqDirs            = {}
        self.dna_regions        = {}
        self.taxonomic_domain   = {}
        self.id_list_passed     = {}
        for lane_key in self.run.run_keys:
            self.id_list_passed[lane_key]   = []
            self.seqDirs[lane_key]          = self.run.samples[lane_key].direction.strip("'").strip('"')
            self.dna_regions[lane_key]      = self.run.samples[lane_key].dna_region.strip("'").strip('"')
            self.taxonomic_domain[lane_key] = self.run.samples[lane_key].taxonomic_domain.strip("'").strip('"')
                    #####################
            #
            #  PrimerSuite Class
            #
            #####################
            self.psuite[lane_key] = PrimerSuite(self.taxonomic_domain[lane_key],self.dna_regions[lane_key])
            #self.runbin['psuite'][lane_key]= PrimerSuite(self.taxonomic_domain[lane_key],self.dna_regions[lane_key])
                
            if(self.seqDirs[lane_key] == 'F' or self.seqDirs[lane_key] == 'B'):
                self.proximal_primers[lane_key] = self.psuite[lane_key].primer_expanded_seq_list['F']
                self.distal_primers[lane_key] = self.psuite[lane_key].primer_expanded_seq_list['R']
                 
            if(self.seqDirs[lane_key] == 'R' or self.seqDirs[lane_key] == 'B'):
                self.proximal_primers[lane_key] = revcomp( self.psuite[lane_key].primer_expanded_seq_list['R'] )
                self.distal_primers[lane_key] = revcomp( self.psuite[lane_key].primer_expanded_seq_list['F'] )
                
                
    def do_trim(self):
    
        if self.VERBOSE: print "Starting do_trim()"        
        
        
        trimmed_sequence = self.raw_sequence
        deleted = False
        delete_reason = ''
        
    

        if self.VERBOSE: print self.read_id,'RAW:',self.raw_sequence
        ###################################################################
        #
        #   1-check/count for Ns
        #
        ###################################################################
        # save count for later
        countNs = check_for_Ns(self.raw_sequence)
        

            
        ###################################################################
        #    
        #   3-find/remove runkey
        #
        ###################################################################
        
        tag      = ''
        lane     = ''
        lane_tag = ''
        tag, trimmed_sequence = remove_runkey(self.raw_sequence, self.run_keys)

        if( not tag ):
            deleted = True
            delete_reason = 'runkey'
            self.deleted_count['Runkey'] += 1
            self.deleted_ids['nokey'][self.read_id] = delete_reason
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
            trim_type = C.trim_lengths[self.dna_regions[lane_tag]]['trim_type']  # internal or distal
        
        
        ###################################################################
        #
        #   5-find/remove proximal primer
        #
        ###################################################################
        
        offset = 0
        exactLeft = 0
        primer_name = ''
        orientation = ''
        
        if( not deleted ):
            #print forward_primers[lane_tag]
            
            exactLeft, offset, trimmed_sequence \
                      = trim_proximal_primer( self.proximal_primers[lane_tag], trimmed_sequence)
                      
            if ( exactLeft and (seq_direction == 'F' or seq_direction == 'B' )):  
                orientation = '+'
                primer_name = self.psuite[lane_tag].primer_names[exactLeft]
                
            if ( exactLeft and (seq_direction == 'R' or seq_direction == 'B') ):
                orientation = '-'
                primer_name = self.psuite[lane_tag].primer_names[revcomp(exactLeft)]
                
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
                #print 'found primer',self.read_id,exactLeft, 'Offset:',offset, tag,raw_sequence[:40],trimmed_sequence[:40]
                
                
               
        
        ###################################################################
        #
        #   6-find  & remove distal primer or anchor
        #
        ###################################################################
        
        distal_primer = ''
        exactRight   = ''
        fuzzyRight = ''
        loc           = 0
        if( delete_reason != 'runkey'):
        
            minLength = C.trim_lengths[self.dna_regions[lane_tag]]['length']
            
            if( len(trimmed_sequence) < minLength ):
                deleted = True
                if (not delete_reason): 
                    delete_reason = 'minimum length'
                    self.deleted_count['Minimum Length'] += 1
                #self.deleted_ids[lane_tag][self.read_id] = delete_reason
                    
            else:
                anchor_start = C.trim_lengths[self.dna_regions[lane_tag]]['start'];
                anchor_end   = C.trim_lengths[self.dna_regions[lane_tag]]['end'];
                exact_found = None
                fuzzy_found = None
                
                exactRight, exactTrimmedOff, trimmed_sequence \
                                = trim_distal_primer( self.distal_primers[lane_tag], trimmed_sequence, trim_type, anchor_start, anchor_end )
                
                if not exactRight:
                    fuzzyRight, fuzzyDist, fuzzyTrimmed, fuzzy_match \
                                = trim_fuzzy_distal(  self.distal_primers[lane_tag], trimmed_sequence, trim_type, anchor_start, anchor_end)
                    
                # choose the longest trimmed off section: exactRight vs fuzzyRight
                if (len(exactTrimmedOff) < len(fuzzyRight) ):
                    exactRight = fuzzyRight
                    trimmed_sequence = fuzzyTrimmed
                
                    
                if not exactRight:
                    if not deleted:
                        deleted = True
                        if self.VERBOSE: print 'deleted distal'
                        if (not delete_reason):  
                            delete_reason = 'distal'
                            self.deleted_count['Distal'] += 1
                        #self.deleted_ids[lane_tag][self.read_id] = delete_reason
                     
                else:            
                    if self.VERBOSE: print "ExactRight:", exactRight, "TrimmedOff:",exactTrimmedOff,"\nSeq:", trimmed_sequence
                    
        ###################################################################
        #
        #   7-check length and other sundry things
        #
        ###################################################################
        if self.VERBOSE: print 'length raw:',len(self.raw_sequence),'length trimmed',len(trimmed_sequence)
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
            self.id_list_passed[lane_tag].append(self.read_id)
        
        ###################################################################
        #
        #   8-check quality
        #
        ###################################################################
        if(self.quality_scores and not deleted):
            average_score = check_for_quality(self.raw_sequence, trimmed_sequence, self.quality_scores)
            if average_score < C.minAvgQual:
                deleted = True
                if (not delete_reason):  
                    delete_reason = 'quality'
                    self.deleted_count['Quality'] += 1
                
                if self.VERBOSE: print 'deleted quality'