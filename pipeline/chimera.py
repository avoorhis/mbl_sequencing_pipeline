import subprocess
import sys, os
import time

import pipeline.constants as C

class Chimera:
    """ Define here """
    def __init__(self, run = None, outdir= None, args = None):
        self.run       = run
        self.run_keys  = self.run.run_keys
        self.rundate   = self.run.run_date
        self.outdir    = outdir
        self.QUIET     = args.QUIET
        self.VERBOSE   = args.VERBOSE
        if self.QUIET: self.VERBOSE = False
        
        self.usearch_cmd = 'usearch'
        self.abskew = '1.9'
        self.refdb = '/xraid2-2/g454/blastdbs/rRNA16S.gold.fasta'
        self.files = {}
        self.prefix = {}
        time.sleep(10)
        for lane_key in self.run.run_keys:
            #if file exists and if dna_region is v6v4 or v3v5:
            self.prefix[lane_key] = './' + outdir + '/' + lane_key
            self.files[lane_key] = {}
            self.files[lane_key]['names']       = self.prefix[lane_key] + '.names'
            self.files[lane_key]['unique']      = self.prefix[lane_key] + '.unique.fa'
            self.files[lane_key]['abund']       = self.prefix[lane_key] + '.abund.fa'
            self.files[lane_key]['deleted']     = self.prefix[lane_key]  + ".deleted.txt"
            #self.files[lane_key]['chimera_del'] = self.prefix[lane_key] + ".chimera.deleted.txt"
            self.files[lane_key]['chimera_db']  = self.prefix[lane_key] + ".chimeras.db"
            self.files[lane_key]['chimera_txt'] = self.prefix[lane_key] + ".chimeras.txt"
            
            
            
            
            if os.path.isfile(self.files[lane_key]['names']) and os.path.isfile(self.files[lane_key]['abund']):
                pass
            else:
                pass
        
       
    def chimera_denovo(self,lane_keys):
        
        chimera_region_found = False
        output = {}
        cluster_id_list = []
        for lane_key in lane_keys:
            
            dna_region  = self.run.samples[lane_key].dna_region.strip("'").strip('"')
            if dna_region in C.regions_to_chimera_check:
                chimera_region_found = True
            else:
                if self.VERBOSE: print 'region not checked', dna_region
                continue
            
            
            # file existance has already been checked in __init__             
  
            
            out_fileName = self.prefix[lane_key] + ".chimeras.txt"        
            #clusterize uchime454 -replace -r self.rundate -t chimeras_denovo

            
            uchime_cmd = ["./clusterize/clusterize"]
            uchime_cmd.append(self.usearch_cmd)
            uchime_cmd.append("--uchime")
            uchime_cmd.append(self.files[lane_key]['abund'])
            uchime_cmd.append("--uchimeout")
            uchime_cmd.append(out_fileName)
            uchime_cmd.append("--abskew")
            uchime_cmd.append(self.abskew)
            
            try:
                output[lane_key] = subprocess.check_output(uchime_cmd)
                #print output[lane_key]
                #print output[lane_key].split()[2]
                cluster_id_list.append(output[lane_key].split()[2])
                #print 'Have %d bytes in output' % len(output)
                #print 'denovo',lane_key,output,len(output)
                # len(output) is normally = 47
                if len(output[lane_key]) < 50 and len(output[lane_key]) > 40:
                    if self.VERBOSE: print  lane_key,"uchime denovo seems to have been submitted successfully"
                else:
                    print >>sys.stderr, "uchime denovo may have broken"                    

            except OSError, e:
                print >>sys.stderr, "Execution failed:", e               
        
        if not chimera_region_found:            
            return ('NOREGION','No regions found that need checking','')
        
        for lane_key in output:
            if len(output[lane_key]) > 50 or len(output[lane_key]) < 40:
                return ('ERROR','uchime ref may have broken or empty',lane_key)  
        
        # finally
        return ('SUCCESS','uchime ref seems to have been submitted successfully',cluster_id_list)
        
    def chimera_reference(self,lane_keys):
    
        chimera_region_found = False
        output = {}
        cluster_id_list = []
        for lane_key in lane_keys:
            
            dna_region  = self.run.samples[lane_key].dna_region.strip("'").strip('"')
            if dna_region in C.regions_to_chimera_check:
                chimera_region_found = True
            else:
                if self.VERBOSE: print 'region not checked', dna_region
                continue
            
            out_fileName = self.prefix[lane_key] + ".chimeras.db"      
            
            uchime_cmd = ["./clusterize/clusterize"]
            uchime_cmd.append(self.usearch_cmd)
            uchime_cmd.append("--uchime")
            uchime_cmd.append(self.files[lane_key]['abund'])
            uchime_cmd.append("--uchimeout")
            uchime_cmd.append(out_fileName)
            uchime_cmd.append("--db")
            uchime_cmd.append(self.refdb)
            
            
            try:
                output[lane_key] = subprocess.check_output(uchime_cmd)
                #print 'outsplit',output[lane_key].split()[2]
                cluster_id_list.append(output[lane_key].split()[2])
                #print 'Have %d bytes in output' % len(output)
                #print 'ref',lane_key,output,len(output)
                if len(output[lane_key]) < 50 and len(output[lane_key]) > 40:
                    if self.VERBOSE: print  lane_key,"uchime ref seems to have been submitted successfully"
                else:
                    print >>sys.stderr, "uchime ref may be broke"
               
            except OSError, e:
                print >>sys.stderr, "Execution failed:", e 
        
        if not chimera_region_found:            
            return ('NOREGION','No regions found that need checking','')
              
        for lane_key in output:
            if len(output[lane_key]) > 50 or len(output[lane_key]) < 40:
                return ('ERROR','uchime ref may have broken or empty',lane_key)  
        
        return ('SUCCESS','uchime ref seems to have been submitted successfully',cluster_id_list)
        
            
    def write_chimeras_to_deleted_file(self,lane_keys): 
    
        for lane_key in lane_keys:
            # open  deleted file and append chimera to it
            # open and read both chimeras files: chimeras.db and chimeras.txt
            
            # hash to remove dupes
            chimera_deleted = {}
            for file in [self.files[lane_key]['chimera_db'], self.files[lane_key]['chimera_txt']]:            
                fh = open(file,"r") 
                # make a list of chimera deleted read_ids            
                for line in fh.readlines():
                    lst = line.strip().split()
                    id = lst[1].split(';')[0]
                    chimera_yesno = lst[-1]
                    if(chimera_yesno) == 'Y':
                        chimera_deleted[id] = 'chimera'
        
            fh_del = open(self.files[lane_key]['deleted'],"a")
            for id in chimera_deleted:
                fh_del.write(id+"\tchimera\n") 
            
                      
                    
                            
                
            
            
 