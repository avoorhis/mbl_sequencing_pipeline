# constants.py
# for mbl sequncing pipeline

minimumLength    = 50
minAvgQual       = 30
maxN             = 0

# this is the maximum distance from the end of the sequence where script
# will accept a distal primer if found
distal_from_end  = 12

#default_primers_file = 'suites/mbl_primers.txt'

complementary_bases = {'A':'T',
                       'T':'A',
                       'C':'G',
                       'G':'C'}
                       
# ==================================================================================================
# cluster and chimera checking
regions_to_chimera_check = ['v6v4','v3v5','v4v5','v4v6','v3v1','v5v3']
cluster_max_wait                = 1*60*60  # 1 hour
cluster_check_interval          = 2
cluster_initial_check_interval  = 10
# ==================================================================================================

# 0: anchor begin -- where to start looking for an internal anchor (if distal, doesn't really matter)
# 1: anchor end -- where to stop looking for an internal anchor, set to -1 for distal trimming
# 2: minimum length -- delete if sequence is shorter than this length
# 3: trim type is it looking internal = inside for an anchor or distal = at the end
trim_lengths = {
    'v3'    : {'start' : 110, 'end' : -1,  'length' : 110,  'trim_type' : "distal"},
    'v4'    : {'start' : 112, 'end' : -1,  'length' : 112,  'trim_type' : "distal"},
    'v6'    : {'start' : 50,  'end' : -1,  'length' : 50,   'trim_type' : "distal"},
    'v9'    : {'start' : 70,  'end' : -1,  'length' : 70,   'trim_type' : "distal"},
    'v6v4'  : {'start' : -420, 'end' : -525, 'length' : 400,  'trim_type' : "internal"},
    'v6v4a' : {'start' : -325, 'end' : -425, 'length' : 325,  'trim_type' : "internal"},
    'v3v5'  : {'start' : 375, 'end' : 450, 'length' : 350,  'trim_type' : "internal"}
    }
    

# from Meren's anchortrimming.py script
mbl_anchors = {
              'v6v4-361': {'reversed' : True,  'start'   : 361, 'freedom' : 50, 'length'  : 13,
                           'sequence' : 'G[T,G]AG.[A,G]GT[A,G][A,G]AAT'},
                            # previously determiend anchor consensus: G[T,G]AG.[A,G]GT[A,G][A,G]AAT
              'v6v4-4xx': {'reversed' : True,  'start'   : 480, 'freedom' : 60, 'length'  : 11,
                           'sequence' : 'G[T,G]AG.[A,G]GT[A,G][A,G]AAT'},
                            # previously determiend anchor consensus: G[T,G]AG.[A,G]GT[A,G][A,G]AAT
              'v3v5-440': {'reversed' : False, 'start'   : 440, 'freedom' : 30, 'length'  : 13,
                           'sequence' : 'GGATTAGA[T,G]ACCC'},
                            # previously determined anchor consensus: GGATTAGA[T,G]ACCC
              'v3v5-370': {'reversed' : False, 'start'   : 370, 'freedom' : 50, 'length'  : 12,
                           'sequence' : '[A,T,C][A,T,G]GCGAA[A,G]GC[A,G][A,C,G]'},
                            # previously determined anchor consensus: [A,T,C][A,T,G]GCGAA[A,G]GC[A,G][A,C,G]
}

# for anchor trimming: lowenshtein distance
max_divergence = 0.9

mbl_primer_suites = {
        "Archaeal:v6" :  
            { 
                '958F'  : {'domain':'archaea', 'region':'v6', 'direction':'F', 'sequence':'AATTGGA.?TCAACGCC.G'},                        
                '1048R' : {'domain':'archaea', 'region':'v6', 'direction':'R', 'sequence':'GWGGTRCATGGCY?GY?CG'}
            },
        "Archaeal:v6v4" : 
            {
                '685F-a' : {'domain':'archaea', 'region':'v6v4', 'direction':'F', 'sequence':'G[TA]AG[GA][GA]GT[GA]AAAT' },
                '1048R'  : {'domain':'archaea', 'region':'v6v4', 'direction':'R', 'sequence':'GWGGTRCATGGCY?GY?CG' }
            },
        "Bacterial:v3" : 
            {
                '338F' : {'domain':'bacteria', 'region':'v3', 'direction':'F', 'sequence':'ACTCCTACGGGAGGCAGCAG' },
                '533R' : {'domain':'bacteria', 'region':'v3', 'direction':'R', 'sequence':'GTGCCAGCAGCCGCGGTAA' }
            },
        "Bacterial:v3v1" : 
            {
                '534R-1' : {'domain':'bacteria', 'region':'v3v1', 'direction':'R', 'sequence':'CCAGCAGCCGCGGTAAT' },
                '534R-2' : {'domain':'bacteria', 'region':'v3v1', 'direction':'R', 'sequence':'CCAGCAGCTGCGGTAA.' },
                '8F-2'   : {'domain':'bacteria', 'region':'v3v1', 'direction':'F', 'sequence':'AGAGTTTGATCCTGGCTCAG' }
            },
        "Bacterial:v3v5" : 
            {
                '341F-1' : {'domain':'bacteria', 'region':'v3v5', 'direction':'F', 'sequence':'CCTACGGGAGGCAGCAG' },
                '341F-2' : {'domain':'bacteria', 'region':'v3v5', 'direction':'F', 'sequence':'CCTACGGG.GGC[AT]GCAG' },
                '341F-3' : {'domain':'bacteria', 'region':'v3v5', 'direction':'F', 'sequence':'TCTACGGAAGGCTGCAG' },
                '785F-a' : {'domain':'bacteria', 'region':'v3v5', 'direction':'R', 'sequence':'GGATTAG.TACCC' }
            },
        "Bacterial:v4v5" : 
            {
                '518F-1' : {'domain':'bacteria', 'region':'v4v5', 'direction':'F', 'sequence':'CCAGCAGCCGCGGTAA.' },
                '518F-2' : {'domain':'bacteria', 'region':'v4v5', 'direction':'F', 'sequence':'CCAGCAGCTGCGGTAA.' },
                '926R-1' : {'domain':'bacteria', 'region':'v4v5', 'direction':'R', 'sequence':'ACT[CT]AAA.GAATTGACGG' },
                '926R-2' : {'domain':'bacteria', 'region':'v4v5', 'direction':'R', 'sequence':'ACTCAAAAGAATTGACGG' },
                '926R-3' : {'domain':'bacteria', 'region':'v4v5', 'direction':'R', 'sequence':'ACTCAAAGAAATTGACGG' }
            },
        "Bacterial:v4v6" : 
            {
                '1046F-1' : {'domain':'bacteria', 'region':'v4v6', 'direction':'F', 'sequence':'CGACAGCCATGCA.CACCT' },
                '1046F-2' : {'domain':'bacteria', 'region':'v4v6', 'direction':'F', 'sequence':'CGACAACCATGCA.CACCT' },
                '1046F-3' : {'domain':'bacteria', 'region':'v4v6', 'direction':'F', 'sequence':'CGACGGCCATGCA.CACCT' },
                '1046F-4' : {'domain':'bacteria', 'region':'v4v6', 'direction':'F', 'sequence':'CGACGACCATGCA.CACCT' },
                '534R-1'  : {'domain':'bacteria', 'region':'v4v6', 'direction':'R', 'sequence':'CCAGCAGCCGCGGTAAT' },
                '534R-2'  : {'domain':'bacteria', 'region':'v4v6', 'direction':'R', 'sequence':'CCAGCAGCTGCGGTAA.' }
            },
        "Bacterial:v5v3" : 
            {
                '909F-1' : {'domain':'bacteria', 'region':'v5v3', 'direction':'F', 'sequence':'CCGTCAATTC.?TTT.AGT' },
                '909F-2' : {'domain':'bacteria', 'region':'v5v3', 'direction':'F', 'sequence':'CCGTCAATTCTTTTGAGT' },
                '909F-3' : {'domain':'bacteria', 'region':'v5v3', 'direction':'F', 'sequence':'CCGTCAATTTCTTTGAGT' },
                '357R-1' : {'domain':'bacteria', 'region':'v5v3', 'direction':'R', 'sequence':'CTGCTGCCTCCCGTAGG' },
                '357R-2' : {'domain':'bacteria', 'region':'v5v3', 'direction':'R', 'sequence':'CTGC.GCC.CCCGTAGG' },
                '357R-3' : {'domain':'bacteria', 'region':'v5v3', 'direction':'R', 'sequence':'CTGCAGCCTTCCGTAGA' }
            },
        "Bacterial:v6" : 
            {
                '967F-AQ'   : {'domain':'bacteria', 'region':'v6', 'direction':'F', 'sequence':'CTAACCGA.GAACCT[CT]ACC' },
                '967F-PP'   : {'domain':'bacteria', 'region':'v6', 'direction':'F', 'sequence':'C.ACGCGAAGAACCTTA.C' },
                '967F-UC1'  : {'domain':'bacteria', 'region':'v6', 'direction':'F', 'sequence':'CAACGCGAAAA+CCTTACC' },
                '967F-UC2'  : {'domain':'bacteria', 'region':'v6', 'direction':'F', 'sequence':'CAACGCGCAGAACCTTACC' },
                '967F-UC3'  : {'domain':'bacteria', 'region':'v6', 'direction':'F', 'sequence':'ATACGCGA[AG]GAACCTTACC' },
                '1046R'     : {'domain':'bacteria', 'region':'v6', 'direction':'R', 'sequence':'AGGTG.?TGCATGG*CTGTCG' },
                '1046R-AQ1' : {'domain':'bacteria', 'region':'v6', 'direction':'R', 'sequence':'AGGTG.?TGCATGG*CCGTCG' },
                '1046R-AQ2' : {'domain':'bacteria', 'region':'v6', 'direction':'R', 'sequence':'AGGTG.?TGCATGG*TCGTCG' },
                '1046R-PP'  : {'domain':'bacteria', 'region':'v6', 'direction':'R', 'sequence':'AGGTG.?TGCATGG*TTGTCG' }
            },
        "Bacterial:v6v4" : 
            {
                '565F'    : {'domain':'bacteria', 'region':'v6v4', 'direction':'F', 'sequence':'TGGGCGTAAAG' },
                '1064R-1' : {'domain':'bacteria', 'region':'v6v4', 'direction':'R', 'sequence':'AGGTG.TGCATGGCTGTCG' },
                '1064R-2' : {'domain':'bacteria', 'region':'v6v4', 'direction':'R', 'sequence':'AGGTG.TGCATGGTTGTCG' },
                '1064R-3' : {'domain':'bacteria', 'region':'v6v4', 'direction':'R', 'sequence':'AGGTG.TGCATGGCCGTCG' },
                '1064R-4' : {'domain':'bacteria', 'region':'v6v4', 'direction':'R', 'sequence':'AGGTG.TGCATGGTCGTCG' }
            },
        "Eukaryal:v9" :
            {
                '1380F' : {'domain':'eukarya', 'region':'v9', 'direction':'F', 'sequence':'CCC+TGCC.TTT+GTACACAC..?CCC+' },
                '1389F' : {'domain':'eukarya', 'region':'v9', 'direction':'F', 'sequence':'TTGTACACACCGCCC+' },
                '1510R' : {'domain':'eukarya', 'region':'v9', 'direction':'R', 'sequence':'GTAGGTGAACCTGC.?GAAGG' }
            }
       }


  
    
# primer_hash = {
#         "1046F-1" : {"direction" : "F", "sequence" : "CGACAGCCATGCA.CACCT",  "regions" : "v6",  "location" : ""},
#         "1046F-2" : {"direction" : "F", "sequence" : "CGACAACCATGCA.CACCT",  "regions" : "v6",  "location" : ""},
#         "1046F-3" : {"direction" : "F", "sequence" : "CGACGGCCATGCA.CACCT",  "regions" : "v6",  "location" : ""},
#         "1046F-4" : {"direction" : "F", "sequence" : "CGACGACCATGCA.CACCT",  "regions" : "v6",  "location" : ""},
#         "1046R" : {"direction" : "R", "sequence" : "AGGTG.?TGCATGG*CTGTCG",  "regions" : "v6",  "location" : "bacteria"},
#         "1046R-AQ1" : {"direction" : "R", "sequence" : "AGGTG.?TGCATGG*CCGTCG",  "regions" : "v6",  "location" : "bacteria"},
#         "1046R-AQ2" : {"direction" : "R", "sequence" : "AGGTG.?TGCATGG*TCGTCG",  "regions" : "v6",  "location" : "bacteria"},
#         "1046R-PP" : {"direction" : "R", "sequence" : "AGGTG.?TGCATGG*TTGTCG",  "regions" : "v6",  "location" : "bacteria"},
#         "1046R-set" : {"direction" : "R", "sequence" : "AGGTG.TGCATGG[CT][CT]?GTCG",  "regions" : "v6",  "location" : "bacteria"},
#         "1048R" : {"direction" : "R", "sequence" : "GWGGTRCATGGCY?GY?CG",  "regions" : "v6",  "location" : "archaea"},
#         "1059R" : {"direction" : "R", "sequence" : "GTCGTCAG?CTCG?TG[TC]?CG?TGA",  "regions" : "v6_dutch_old",  "location" : "bacteria"},
#         "1061R" : {"direction" : "R", "sequence" : "[CG]TCGTCAGCTCGTG[TC]CGTGA?",  "regions" : "v6_dutch",  "location" : "bacteria"},
#         "1064R-1" : {"direction" : "R", "sequence" : "AGGTG.TGCATGGCTGTCG",  "regions" : "v6",  "location" : ""},
#         "1064R-2" : {"direction" : "R", "sequence" : "AGGTG.TGCATGGTTGTCG",  "regions" : "v6",  "location" : ""},
#         "1064R-3" : {"direction" : "R", "sequence" : "AGGTG.TGCATGGCCGTCG",  "regions" : "v6",  "location" : ""},
#         "1064R-4" : {"direction" : "R", "sequence" : "AGGTG.TGCATGGTCGTCG",  "regions" : "v6",  "location" : ""},
#         "10F-1" : {"direction" : "F", "sequence" : "AGTTTGATC.TGGCTCA",  "regions" : "v1",  "location" : ""},
#         "10F-2" : {"direction" : "F", "sequence" : "A[AG]TTTGATCTT[AG]GTTCA",  "regions" : "v1",  "location" : ""},
#         "10F-3" : {"direction" : "F", "sequence" : "AGTTTGATCCTGGCTTA",  "regions" : "v1",  "location" : ""},
#         "1194R" : {"direction" : "R", "sequence" : "GAGGAAGG.GGGGA[CT]GACGT",  "regions" : "v5v7",  "location" : "bacteria"},
#         "1380F" : {"direction" : "F", "sequence" : "CCC+TGCC.TTT+GTACACAC..?CCC+",  "regions" : "v9",  "location" : "eukarya"},
#         "1389F" : {"direction" : "F", "sequence" : "TTGTACACACCGCCC+",  "regions" : "v9",  "location" : "eukarya"},
#         "1510R" : {"direction" : "R", "sequence" : "GTAGGTGAACCTGC.?GAAGG",  "regions" : "v9",  "location" : "eukarya"},
#         "1513R" : {"direction" : "R", "sequence" : "AAGTC[AG]TAACAAGGTA[AG]CCGTA",  "regions" : "v9",  "location" : "bacteria"},
#         "26R-1" : {"direction" : "R", "sequence" : "TGAGCCA.GATCAAACT",  "regions" : "v1",  "location" : ""},
#         "26R-2" : {"direction" : "R", "sequence" : "TGAAC[CT]AAGATCAAA[CT]T",  "regions" : "v1",  "location" : ""},
#         "26R-3" : {"direction" : "R", "sequence" : "TAAGCCAGGATCAAACT",  "regions" : "v1",  "location" : ""},
#         "338F" : {"direction" : "F", "sequence" : "ACTCCTACGGGAGGCAGCAG",  "regions" : "v3",  "location" : "bacteria"},
#         "341F-1" : {"direction" : "F", "sequence" : "CCTACGGGAGGCAGCAG",  "regions" : "v3",  "location" : "bacteria"},
#         "341F-2" : {"direction" : "F", "sequence" : "CCTACGGG.GGC[AT]GCAG",  "regions" : "v3",  "location" : "bacteria"},
#         "341F-3" : {"direction" : "F", "sequence" : "TCTACGGAAGGCTGCAG",  "regions" : "v3",  "location" : "bacteria"},
#         "341F-A" : {"direction" : "F", "sequence" : "CCTACGGG[AG][GC]GCAGCAG",  "regions" : "v3",  "location" : "archaea"},
#         "357R-1" : {"direction" : "R", "sequence" : "CTGCTGCCTCCCGTAGG",  "regions" : "v3",  "location" : "bacteria"},
#         "357R-2" : {"direction" : "R", "sequence" : "CTGC.GCC.CCCGTAGG",  "regions" : "v3",  "location" : "bacteria"},
#         "357R-3" : {"direction" : "R", "sequence" : "CTGCAGCCTTCCGTAGA",  "regions" : "v3",  "location" : "bacteria"},
#         "518F-1" : {"direction" : "F", "sequence" : "CCAGCAGCCGCGGTAA.",  "regions" : "v4",  "location" : ""},
#         "518F-2" : {"direction" : "F", "sequence" : "CCAGCAGCTGCGGTAA.",  "regions" : "v4",  "location" : ""},
#         "520F" : {"direction" : "F", "sequence" : "A.TGGG..TAAAG.G",  "regions" : "v4",  "location" : "bacteria"},
#         "525F" : {"direction" : "F", "sequence" : "[CG]GT[CT][CT]GTAGC[CT][GT]G",  "regions" : "v4",  "location" : "archaea"},
#         "531R-A" : {"direction" : "R", "sequence" : "ACCGCGGC[GT]GCTGGC",  "regions" : "v3",  "location" : "archaea"},
#         "533R" : {"direction" : "R", "sequence" : "GTGCCAGCAGCCGCGGTAA",  "regions" : "v3",  "location" : "bacteria"},
#         "534R-1" : {"direction" : "R", "sequence" : "CCAGCAGCCGCGGTAAT",  "regions" : "v3",  "location" : ""},
#         "534R-2" : {"direction" : "R", "sequence" : "CCAGCAGCTGCGGTAA.",  "regions" : "v5",  "location" : ""},
#         "565F" : {"direction" : "F", "sequence" : "TGGGCGTAAAG",  "regions" : "v4",  "location" : "bacteria"},
#         "575R" : {"direction" : "R", "sequence" : "CTTTACGCCCA",  "regions" : "v4",  "location" : "bacteria"},
#         "680R-Vib" : {"direction" : "R", "sequence" : "CTGTAGAGGGGG+TAGAA",  "regions" : "v4",  "location" : ""},
#         "685F-a" : {"direction" : "F", "sequence" : "G[TA]AG[GA][GA]GT[GA]AAAT",  "regions" : "v4",  "location" : "archaea"},
#         "784F-DEG" : {"direction" : "F", "sequence" : "GGNTTAGATACCC",  "regions" : "v5",  "location" : "bacteria"},
#         "785F" : {"direction" : "F", "sequence" : "GGATTAGATACCC",  "regions" : "v4",  "location" : "bacteria"},
#         "785F-a" : {"direction" : "R", "sequence" : "GGATTAG.TACCC",  "regions" : "v3v5",  "location" : "bacteria"},
#         "785F-D" : {"direction" : "F", "sequence" : "GGATTAGATACCC.[AG]GTA[CG]TC",  "regions" : "v6_dutch",  "location" : "bacteria"},
#         "787F" : {"direction" : "F", "sequence" : "ATTAGATACCC.GGTAG",  "regions" : "quince_v6",  "location" : "bacteria"},
#         "797R" : {"direction" : "R", "sequence" : "GGGTATCTAATCC",  "regions" : "v4",  "location" : ""},
#         "802R" : {"direction" : "F", "sequence" : "GGATTAGATACCC..GTA",  "regions" : "v4",  "location" : "bacteria"},
#         "880R" : {"direction" : "R", "sequence" : "CNCCTGNGNAGTANG",  "regions" : "v5",  "location" : "bacteria"},
#         "8F" : {"direction" : "F", "sequence" : "AGAGTTTGATC[AC]TGGC",  "regions" : "v1",  "location" : ""},
#         "8F-1" : {"direction" : "F", "sequence" : "",  "regions" : "v3v1",  "location" : "bacteria"},
#         "8F-2" : {"direction" : "F", "sequence" : "AGAGTTTGATCCTGGCTCAG",  "regions" : "v1",  "location" : ""},
#         "8F-3" : {"direction" : "F", "sequence" : "",  "regions" : "v3v1",  "location" : "bacteria"},
#         "906F" : {"direction" : "F", "sequence" : "CCGTCAATT[CT][CT]TTT[GA]AGTTT",  "regions" : "v6_dutch_old",  "location" : "bacteria"},
#         "909F-1" : {"direction" : "F", "sequence" : "CCGTCAATTC.?TTT.AGT",  "regions" : "v5",  "location" : ""},
#         "909F-2" : {"direction" : "F", "sequence" : "CCGTCAATTCTTTTGAGT",  "regions" : "v5",  "location" : ""},
#         "909F-3" : {"direction" : "F", "sequence" : "CCGTCAATTTCTTTGAGT",  "regions" : "v5",  "location" : ""},
#         "909F-4" : {"direction" : "F", "sequence" : "",  "regions" : "v5v3",  "location" : "bacteria"},
#         "916F" : {"direction" : "F", "sequence" : "GAATTGACGGGG.CCCGCA",  "regions" : "v6_dutch_old",  "location" : "bacteria"},
#         "926R-1" : {"direction" : "R", "sequence" : "ACT[CT]AAA.GAATTGACGG",  "regions" : "v5",  "location" : ""},
#         "926R-2" : {"direction" : "R", "sequence" : "ACTCAAAAGAATTGACGG",  "regions" : "v5",  "location" : ""},
#         "926R-3" : {"direction" : "R", "sequence" : "ACTCAAAGAAATTGACGG",  "regions" : "v5",  "location" : ""},
#         "958F" : {"direction" : "F", "sequence" : "AATTGGA.?TCAACGCC.G",  "regions" : "v6",  "location" : "archaea"},
#         "967F" : {"direction" : "F", "sequence" : "CAACGCGAAGAACCTTACC",  "regions" : "v6",  "location" : "bacteria"},
#         "967F-AQ" : {"direction" : "F", "sequence" : "CTAACCGA.GAACCT[CT]ACC",  "regions" : "v6",  "location" : "bacteria"},
#         "967F-PP" : {"direction" : "F", "sequence" : "C.ACGCGAAGAACCTTA.C",  "regions" : "v6",  "location" : "bacteria"},
#         "967F-UC1" : {"direction" : "F", "sequence" : "CAACGCGAAAA+CCTTACC",  "regions" : "v6",  "location" : "bacteria"},
#         "967F-UC2" : {"direction" : "F", "sequence" : "CAACGCGCAGAACCTTACC",  "regions" : "v6",  "location" : "bacteria"},
#         "967F-UC3" : {"direction" : "F", "sequence" : "ATACGCGA[AG]GAACCTTACC",  "regions" : "v6",  "location" : "bacteria"},
#         "968F" : {"direction" : "F", "sequence" : "AACGCGAAGAACCTTACC",  "regions" : "v6",  "location" : "bacteria"},
#         "A-adaptor" : {"direction" : "F", "sequence" : "GCCTCCCTCGCGCCATCAG",  "regions" : "454",  "location" : ""},
#         "B-adaptor" : {"direction" : "R", "sequence" : "GCCTTGCCAGCCCGCTCAG",  "regions" : "454",  "location" : ""},
#         "BSL" : {"direction" : "F", "sequence" : "GGCTTATTACAACTTA.CAAG",  "regions" : "",  "location" : "eukarya"},
#         "CDSIII" : {"direction" : "R", "sequence" : "AAAAAAAAAAAAAA*TCGGCCCGCCTCGGCCTCTAG",  "regions" : "est",  "location" : "eukarya"},
#         "EST-GG" : {"direction" : "F", "sequence" : "GG?",  "regions" : "est",  "location" : "eukarya"},
#         "TopoF" : {"direction" : "F", "sequence" : "GTTTAAACGAATTCGCCCTT",  "regions" : "pCR4",  "location" : "eukarya"},
#         "TopoR" : {"direction" : "R", "sequence" : "AAGGGCGATTCGCGGCCGC",  "regions" : "pCR4",  "location" : "eukarya"}
#       }
