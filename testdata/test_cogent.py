#!/usr/local/www/vamps/software/python/bin/python

#from cogent.parse import binary_sff

#fn = binary_sff('testdata/sff_reads_1050.sff')
#seqs = LoadSeqs(fn, moltype=DNA, aligned=False)
#print seqs
import qiime.split_libraries
from cogent import LoadSeqs, DNA
from cogent.parse.binary_sff import ( 
    seek_pad, parse_common_header, parse_read_header, parse_read_data,
    validate_common_header, parse_read, parse_binary_sff, UnsupportedSffError,
    write_pad, write_common_header, write_read_header, write_read_data,    
    parse_binary_sff, write_binary_sff 
    ) 
sff_in = open("testdata/sff_reads_1050.sff") 
#sff_out = open("filtered.sff", "wb") 
# Returns generator of reads 
header, reads = parse_binary_sff(sff_in, native_flowgram_values=True) 
aln = LoadSeqs(data=reads)
#header, reads = parse_read(sff_in, native_flowgram_values=True) 
for read in reads:
    print read["Name"], read["Bases"]
# Force evaluation of reads 
reads = [r for r in reads if r["number_of_bases"] > 504] 
# Adjust number of reads in SFF header 
header['number_of_reads'] = len(reads) 
# No index written by write_binary_sff 
header['index_offset'] = 0 
header['index_length'] = 0 
#write_binary_sff(sff_out, header, reads) 
sff_in.close() 
#sff_out.close() 
