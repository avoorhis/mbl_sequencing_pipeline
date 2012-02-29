import os,sys
from pipeline.fastalib import *
from pipeline.Fasta import sfasta
from pipeline.utils import *
from time import sleep
# class to receive fasta file
# and write out a unique fasta file, names file and abundance
# should subclass Meren's fastalib.py


class MBLPipelineFastaUtils:
    """
    will get either a list of keys and a directory to look for fasta files
    or a single fasta file
    
    """
    def __init__(self, lane_keys=None,outputdir=None,input_fasta_file_path=None):
        self.inputFileName = {}
        self.orphans = {}
        if input_fasta_file_path:
            self.inputFileName['none'] = input_fasta_file_path
            self.fasta = SequenceSource(self.inputFileName['none'], True, True)
        else:
            self.lane_keys = lane_keys
            self.outputdir = outputdir
            for lane_key in lane_keys:
                self.inputFileName[lane_key] = self.outputdir + "/" + lane_key + ".trimmed.fa"
                self.orphans[lane_key] = {}
        
    def write_clean_fasta_file(self):
        """
        def to write a new fasta from the original fasta file 
                using the deleted file
                
        The deleted file contains the trimming deleted as well
        as the chimera deleted
        Then write the uniques from Meren's fastalib
        """
        sleep(2)
        for lane_key in self.lane_keys:
            
            deleted_id_list = []
            original_trimmed_file   = self.outputdir + "/" + lane_key + ".trimmed.fa" 
            new_trimmed_file_name   = self.outputdir + "/" + lane_key + ".newtrimmed.fa"
            new_trimmed_file        = FastaOutput(new_trimmed_file_name)
            deleted_file            = self.outputdir + "/" + lane_key + ".deleted.txt" 
               

            if not (os.path.exists(deleted_file) and os.path.getsize(deleted_file) > 0):
                continue
            del_fh = open(deleted_file,"r")
            for line in del_fh.readlines():
                #print 'line:',line
                lst = line.strip().split()                
                deleted_id_list.append(lst[0])
            
            # open trimmed file and read a line             
            trimmedfasta = SequenceSource(original_trimmed_file)
 
            while trimmedfasta.next():
                if trimmedfasta.id not in deleted_id_list:
                    new_trimmed_file.store(trimmedfasta)
            new_trimmed_file.close()
            
            # rename to newtrimmed => trimmed
            os.rename(original_trimmed_file, self.outputdir + "/" + lane_key + ".trimmed_dirty.fa" )
            os.rename(new_trimmed_file_name, original_trimmed_file )
            
            
    def write_clean_names_file(self):
        """
        Writes a new names file from the deleted file.
        The deleted file contains the trimming delted as well
        as the chimera deleted
        """
        
        for lane_key in self.lane_keys:
            
            deleted_id_list = []
            original_names_file = self.outputdir + "/" + lane_key + ".names" 
            new_names_file      = self.outputdir + "/" + lane_key + ".newnames"
            deleted_file        = self.outputdir + "/" + lane_key + ".deleted.txt" 
            
            if not (os.path.exists(deleted_file) and os.path.getsize(deleted_file) > 0):
                continue
            del_fh = open(deleted_file,"r")
            for line in del_fh.readlines():
                lst = line.strip().split()                
                deleted_id_list.append(lst[0])
            
            newnames_fh = open( new_names_file,     "w" )
            names_fh    = open( original_names_file,"r" )
            for line in names_fh.readlines():
                #print lane_key,'line',line
                lst = line.strip().split()
                id = lst[0]                 # a read_id
                dupes = lst[1].split(',')   # a list of read_ids
                newdupes = []
                if id in deleted_id_list:
                    #print 'found id',id,dupes
                    if len(dupes) == 1:
                        continue
                    for subid in dupes:
                        if subid in deleted_id_list:
                            pass
                        else:
                            newdupes.append(subid)
                                            
                    self.orphans[lane_key][id] = newdupes
                    
                    id = newdupes[0] 
                    
                else:
                    # This id is not in deleted list
                    for subid in dupes:
                        #print 'looking through dupes',subid
                        if subid in deleted_id_list:
                            pass
                        else:
                            newdupes.append(subid)
                    
                newnames_fh.write(id+"\t"+','.join(newdupes)+"\n")
            
            newnames_fh.close()
            names_fh.close()
            
            # rename to newnames => names
            os.rename(original_names_file, self.outputdir + "/" + lane_key + ".names_dirty" )
            os.rename(new_names_file, original_names_file )            
            
                    
    def write_clean_uniques_file(self):
        """
        Write out a new unique file with all the deleted ids removed
           especially the chimeras which were detected after the original unique file 
           was created.
        """
        for lane_key in self.lane_keys:
       
            deleted_id_list = []
            new_unique_file_name    = self.outputdir + "/" + lane_key +".newunique.fa"
            new_unique_file         = FastaOutput(new_unique_file_name)            
            original_unique_file    = self.outputdir + "/" + lane_key + '.unique.fa'
            deleted_file            = self.outputdir + "/" + lane_key + ".deleted.txt" 

            if not (os.path.exists(deleted_file) and os.path.getsize(deleted_file) > 0):
                continue
            del_fh = open(deleted_file,"r")
            for line in del_fh.readlines():
                lst = line.strip().split()                
                deleted_id_list.append(lst[0])
            #print deleted_id_list
            # open unique file and read a line             
            uniquesfasta = SequenceSource(original_unique_file)
            while uniquesfasta.next():
                #print uniquesfasta.id,self.orphans[lane_key]
                
                if uniquesfasta.id in self.orphans[lane_key].keys():
                    #print "found orphan",uniquesfasta.id
                    uniquesfasta.id = self.orphans[lane_key][uniquesfasta.id][0]
                    #print "new id",uniquesfasta.id
                if uniquesfasta.id not in deleted_id_list:
                    new_unique_file.store(uniquesfasta)
            new_unique_file.close()
            
            # rename to newuniques => uniques
            os.rename(original_unique_file, self.outputdir + "/" + lane_key + ".unique_dirty.fa" )
            os.rename(new_unique_file_name, original_unique_file ) 
            
            
    def write_clean_abundance_file(self):
        """
        Writes the abundance file from the new names file and new unique file.
           Thes file have already had their ids checked from the deleted file
        """
        for lane_key in self.lane_keys:
            original_abundance_file = self.outputdir + "/" + lane_key + ".abund.fa"   
            new_abundance_file      = self.outputdir + "/" + lane_key + ".newabund.fa"            
            new_names_file          = self.outputdir + "/" + lane_key + ".names"
            new_unique_file         = self.outputdir + "/" + lane_key + ".unique.fa" 
            deleted_file            = self.outputdir + "/" + lane_key + ".deleted.txt" 
            names = {}
            uniques = {}
            if not (os.path.exists(deleted_file) and os.path.getsize(deleted_file) > 0):
                continue
                
            newnames_fh = open( new_names_file,     "r" )
            for line in newnames_fh.readlines():
                lst = line.strip().split()
                
                names[lst[0]] = lst[1].split(',')
            #print names  
            fasta = SequenceSource(new_unique_file)
            
            
            while fasta.next():
                fasta.id
                uniques[fasta.seq] = fasta.id
            #print uniques    
            sorted_uniques = mysort( uniques, names )
            
            for item in sorted_uniques:
                read_id = item[0]
                count = item[1]
                seq = item[2]
             
                sfastaRead = read_id + ";size="+str(count)                
                abundfa = sfasta(sfastaRead, seq)
                abundfa.write(new_abundance_file,'a')
            
            # rename to newuniques => uniques
            os.rename(original_abundance_file, self.outputdir + "/" + lane_key + ".abund_dirty.fa" )
            os.rename(new_abundance_file, original_abundance_file )            
     
            
    def write_clean_files_to_database(self):
        # Write new clean files to the database: env454
        # rawseq table not used
        # trimseq
        # runkeys
        # primers
        # run primers
        for lane_key in self.lane_keys:
            abundance_file      = self.outputdir + "/" + lane_key + ".abund.fa"            
            names_file          = self.outputdir + "/" + lane_key + ".names"
            unique_file         = self.outputdir + "/" + lane_key + ".unique.fa"     
            deleted_file        = self.outputdir + "/" + lane_key + ".deleted.txt"
            
    # primers table
    # primer_name:967F	direction:F	sequence:CAACGCGAAGAACCTTACC	dna_region:v6	originalseq:CAACGCGAAGAACCTTACC	domain:bacteria
    # run_primer_view
    # rn_primer
    # run_key
    
            
            
            