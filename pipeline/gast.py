from subprocess import call, check_output
import sys, os
import time
import constants as C

class Gast:
    """Doc string here.."""
    Name = "GAST"
    def __init__(self, run = None, outdir= None, args = None):

        self.run = run
        self.outdir = outdir
        self.rundate = self.run.run_date
        self.outdir    = outdir
        self.QUIET     = args.QUIET
        self.VERBOSE   = args.VERBOSE
        if self.QUIET: self.VERBOSE = False
        
        
        # 1) clustergast
        # 2) gast cleanup
        # 3) gast2tax
        
        
    def clustergast(self,lane_keys):
        print 'in cg',lane_keys
        # use fasta.uniques file
        # split into smaller files
        # usearch --cluster each
        #######################################
        #
        # Split the uniques fasta and run UClust per node
        #
        #######################################
        for lane_key in lane_keys:
            nodes = 100
            qsub_prefix = 'clustergast_sub_'
            
            unique_file = self.outdir +'/'+lane_key+'.unique.fa'
            gast_dir = self.outdir +'/'+lane_key+'_gast/'
            grep_cmd = ['grep','-c','>',unique_file]
            print grep_cmd
            facount = check_output(grep_cmd)
            print lane_key,'count',facount
            calcnode_cmd = ['calcnodes','-t',str(facount),'-n',str(nodes),'-f','1']
            
            calcout = check_output(calcnode_cmd)
            print 'calcout',calcout
            #calcout:
            # node=1 start=1 end=1 rows=1
            # node=2 start=2 end=2 rows=1
            # node=3 start=3 end=3 rows=1           
            lines = clacout.split("\n")
            i = 1
            for line in lines:
                script_filename = gast_dir + qsub_prefix + i
                fh = open(script_filename,'w')
                log_file = gast_dir + 'clustergast.log_' + i
                qstat_name = "gast" + lane_key + '_' + self.rundate + "_" + i
                fh.write("#!/bin/csh\n")
                fh.write("#\$ -j y\n" )
                fh.write("#\$ -o " + log_file + "\n")
                fh.write("#\$ -N " + qstat_name + "\n\n")
                
                data = line.split()
                start = data[1].split('=')[1]
                end  = data[2].split('=')[1]

                fastasampler_cmd = ['fastasampler']
                fastasampler_cmd.append('-n')
                fastasampler_cmd.append(str(start)+','+str(end))
                fastasampler_cmd.append(unique_file)
                fastasampler_cmd.append(tmp_out_file)
                
                

                usearch_cmd = ['usearch']
                usearch_cmd.append('--global')
                usearch_cmd.append('--query')
                usearch_cmd.append(tmp_fasta_filename)
                usearch_cmd.append('--iddef')
                usearch_cmd.append('3')
                usearch_cmd.append('--gapopen')
                usearch_cmd.append('6I/1E')
                usearch_cmd.append('--db')
                usearch_cmd.append()
                usearch_cmd.append('--uc')
                usearch_cmd.append(tmp_usearch_filename)
                usearch_cmd.append('--maxaccepts')
                usearch_cmd.append(str(max_accepts))
                usearch_cmd.append('--maxrejects')
                usearch_cmd.append(str(max_rejects))
                usearch_cmd.append('--id')
                usearch_cmd.append(str(pctid_threshold))
  
                clean_usearch_cmd = []
                
                
                i += 1
                fh.close()
                
            print "exiting"
            sys.exit()
# 
#     # Count the number of sequences so the job can be split for nodes
#     $facount = `grep -c \">\" $fasta_uniqs_filename`;
#     chomp $facount;
#     my $calcs = `/bioware/seqinfo/bin/calcnodes -t $facount -n $nodes -f 1`;
#     my @lines = split(/\n/, $calcs);
#     warn "facount= $facount\ncalcnodes command= /bioware/seqinfo/bin/calcnodes -t $facount -n $nodes -f 1\n";
#     my $i = 1;
#     foreach my $l (@lines)
#     {
#         # calcnodes returns a series of lines describing each set, parse the set and use for the qsub
#         my @data = split (/\s+/, $l);
# 
#         my $start = $data[1];
#         $start =~ s/start=//;
#         my $end = $data[2];
#         $end =~ s/end=//;
# 
#         # Create qsub script
#         my $script_filename = $gastDir .'/'. $qsub_prefix . $file_prefix . $rand_file . $i . ".sh";
#         my $qstat_name = "gast" . $vRegion .'_'. $runcode . "_" . $i;
#         my $tmp_usearch_filename = $gastDir .'/'. $usearch_filename . $rand_file . "$i.txt";
#         my $usearch_filename = $usearch_filename . "_$i.txt";
#         my $tmp_fasta_filename = $gastDir .'/'. $fasta_uniqs_filename . $rand_file . $i;
#         
#         # here the prefix to $gast_filename should be the same as $gast_table name for mysqlImport to work
#         my $gast_filename = $gastDir."/gast_" . $file_prefix . ".$i";
#         open (QSUB, ">$script_filename") || die "Unable to write to qsub script: $script_filename.  Exiting\n";
#         my $log_file2 = $gastDir.'/clustergast.log_'.$i;
#         print QSUB
#         "#!/bin/csh\n" .
#         "#\$ -j y\n" .
#         "#\$ -o $log_file2\n" .
#         "#\$ -N $qstat_name\n" .
#         "\n";
# 
#         # Run UClust
#         my $script_text = "";
#         my $fastasampler_cmd = "/bioware/seqinfo/bin/fastasampler -n $start,$end ${gastDir}/${fasta_uniqs_filename} $tmp_fasta_filename";
#         $script_text .= "$fastasampler_cmd\n\n";
#         if ($verbose) {warn "$fastasampler_cmd\n";}
# 
# #       my $uclust_cmd = "/bioware/uclust/uclust --mbl_id --input $tmp_fasta_filename --lib $refhvr_fa --uc 
# #$tmp_uclust_filename --libonly --allhits --maxaccepts $max_accepts --maxrejects $max_rejects --id $pctid_threshold";
# #         my $uclust_cmd = "/bioware/uclust/usearch --cluster $tmp_fasta_filename --iddef 3 --gapopen 6I/1E --db 
# #$refhvr_fa --uc $tmp_uclust_filename --libonly --allhits --maxaccepts $max_accepts --maxrejects $max_rejects --id 
# #$pctid_threshold";
# 
#        # updated 9/13/2011 for usearch 5.0
#        #my $uclust_cmd =  "/bioware/uclust/usearch --cluster $tmp_fasta_filename --usersort --iddef 3 --gapopen 6I/1E --db $refhvr_fa --uc $tmp_usearch_filename --maxaccepts $max_accepts --maxrejects $max_rejects --id $pctid_threshold";
#        my $usearch_cmd = "/bioware/uclust/usearch --global --query $tmp_fasta_filename --iddef 3 --gapopen 6I/1E --db $refhvr_fa --uc $tmp_usearch_filename --maxaccepts $max_accepts --maxrejects $max_rejects --id $pctid_threshold";
# 
# print "$usearch_cmd\n";
#         $script_text .= "$usearch_cmd\n\n";
#         if ($verbose) {warn "$usearch_cmd\n";}
#     
#         # sort the results for valid hits, saving only the ids and pct identity
#         my $clean_usearch_cmd = "grep -P \"^H\\t\" $tmp_usearch_filename | sed -e 's/|.*\$//' | awk '{print \$9 \"\\t\" \$4 \"\\t\" \$10 \"\\t\" \$8}' | sort -k1,1b -k2,2gr | clustergast_tophit > $gast_filename";
#         $script_text .= "$clean_usearch_cmd\n\n";
#         if ($verbose) {warn "Ridiculous commandline uclust clean:\n\t$clean_usearch_cmd\n\n";}
# 
#         # Load the resulting file into the database
#         my $load_cmd = "$sqlImportCommand -u $db_user -p$db_password -C -v -L -h $db_hostname $dbName $gast_filename";      
#         #my $load_cmd = "mysql -h $db_hostname -u $db_user -p$db_password -D $dbName -e \"load data local infile '$gast_filename' into table $gast_table\"";
#         $script_text .= "$load_cmd\n\n";
#         if($verbose) {warn "Load to database: $load_cmd\n";}
# 
#         # Clean up the tmp files
#         my $clean_qsubfile_cmd = 
#         "/bin/rm -f $tmp_fasta_filename\n" .
#         "/bin/rm -f $script_filename\n" ;
# #        "/bin/rm -f $tmp_usearch_filename\n";
# #        $script_text .= "$clean_qsubfile_cmd\n\n";
#         if ($verbose) {warn "remove qsub file: $clean_qsubfile_cmd\n";}
# 
#         # Print the commands to the qsub script and then
#         # submit the script to the cluster!
#         print QSUB $script_text;
#         close(QSUB);
#         system("chmod +x $script_filename");
# 
#         my $qsub_error = system("/usr/local/sge/bin/lx24-amd64/qsub $script_filename");
#         #my $qsub_jobid = `qsub -terse $script_filename`;
#         if ($qsub_error) {my $err_txt = "Error submitting $script_filename to the cluster, $qsub_error.\n"; warn $err_txt; warn $err_txt;}
#         
#         $i++;
#     }
# }
     
    def gast_cleanup(self,lane_keys):
        pass
        
    def gast2tax(self,lane_keys):
        pass