
import os, sys
import pipeline.anchortrimming_mbl as anchortrim

#from pipeline.utils import revcomp

import constants as C
try:
    import Levenshtein
except:
    print '''
    You need Levenshtein module installed to run this software.

    Here is a fast implementation of Levenshtein distance for Python:

        http://code.google.com/p/pylevenshtein/

'''
    sys.exit(-1)

def count_keys(hash):
    return len([x for x in hash])
        
def trim_stop_seq( stop_seqs, seq, trim_type, start, end ):
    
    for anchor in stop_seqs:
        anchor_length = len(anchor)
        print anchor,start,end,len(seq)
        for pos in range(start,end):
            seq_window = seq[pos:pos+anchor_length]

            
            dist = abs( Levenshtein.ratio( anchor,     seq_window ) )
            #dist2 = abs( Levenshtein.ratio( seq_window, anchor )     )
            if dist == 1.0:
                # perfect match
                # do I trim off before or after anchor?
                #before
                return anchor,seq[pos+anchor_length:],seq[:pos]
                #after (include anchor in trimmed seq):
                #return anchor,seq[pos:],seq[:pos+anchor_length]
            if dist >= C.max_divergence:
                pass
            #print pos,seq_window,dist1,dist2
    
    return '','',seq
    
def trim_anchor(anchor_name, anchor_list, sequence, readDir):
    print 'looking for anchor:',anchor_name
    exact = ''
    exactTrimmedOff = ''    
    
    #print sequence
    freedom         = C.mbl_anchors[anchor_name]['freedom']
    length          = C.mbl_anchors[anchor_name]['length']
    start           = C.mbl_anchors[anchor_name]['start']
    reversed        = C.mbl_anchors[anchor_name]['reversed']
    max_divergence  = C.max_divergence
    if (reversed == True and readDir == 'F') or (reversed == False and readDir == 'R'):
        print "anchor reversed vs readDir mismatch -Skipping",reversed,readDir
        return '','',sequence
    #print 'Reversed',reversed
    list_of_tuples = anchortrim.generate_tuples(start, freedom, length, reversed_read = reversed)
    #print list_of_tuples
    if readDir == 'R':
        import pipeline.utils as utils
        anchor_list = utils.revcomp(anchor_list)
    print anchor_list
    anchor, location = anchortrim.find_best_distance(sequence, anchor_list, max_divergence, list_of_tuples)
    
    if anchor and location:
        print 'anchor:',anchor,'loc:',location
              
        # busines time.
        if reversed:
            trimmed_sequence = sequence[location[0]:] # this includes the anchor in trimmed sequence..
        else:
            trimmed_sequence = sequence[:location[1]] # same thing here for the reversed == False
        #trimmed_sequence = anchortrim.trim_sequence(seq, location, s)
        
        #print 'trimmed',trimmed_sequence
        exact = anchor
        if reversed:
            exactTrimmedOff = sequence[:location[0]]
            #return sequence[location[0]:] # this includes the anchor in trimmed sequence..
        else:
            exactTrimmedOff = sequence[location[1]:]
            #return sequence[:location[1]] # same thing here for the reversed == False
    else:
        print 'no anchor location found'
        trimmed_sequence = sequence
    
    
    return exact, exactTrimmedOff, trimmed_sequence

def trim_distal_primer(primers_list, seq, readDir):
    """Doc string here"""
    d_primer_found = ""
    trimmedPortion = ""
    
    loc = 0
 
    for p in primers_list:
        #print 'dist primer:',p
        # find the furthest RIGHT p in seq
        if seq.rfind(p) != -1:
            pos = seq.rfind(p)
            #print 'pos',pos
                
            # found whole exact primer
            primer_found = p
            loc = seq.find(p)
            #print 'loc',loc
            trimmedPortion = seq[loc+len(p):] # 
            #print 'len1',trimmedPortion
            if len(trimmedPortion) < C.distal_from_end:
                #print 'len',trimmedPortion
                seq = seq[:loc]
            
                return primer_found, trimmedPortion, seq
            
        
        else:
            truncLength = len(p) + 1 #add one just in case!
            while truncLength >= 5:
                dSeq = p[:truncLength] # cut primer to the appropriate length
                # 123456789
                if dSeq in seq:
                    primer_found = dSeq
                    loc = seq.find(dSeq)
                    
                    trimmedPortion = seq[loc-len(dSeq):]
                    
                    if len(trimmedPortion) - len(dSeq) < C.distal_from_end:
                        seq = seq[:loc]
                        return primer_found, trimmedPortion, seq
                    
                truncLength = truncLength - 1
                
    return '', '', seq
    
def trim_fuzzy_distal(anchors_list, seq, trim_type, start, end):
    """Doc string here.."""
    max_distance = 3
    best_distance = max_distance + 1
    found_fuzzy = 0
    fuzzy_match = ""
    for anchor in anchors_list:
        anchor_length = len(anchor)
        for pos in range(start,end):
        
            seq_window = seq[pos:anchor_length]
            
            dist = 0
            
            #dist1 = abs( Levenshtein.ratio( anchor,     seq_window ) )
           # dist2 = abs( Levenshtein.ratio( seq_window, anchor )     )
            dist1 = abs( levenshtein( anchor,     seq_window ) )
            dist2 = abs( levenshtein( seq_window, anchor )     )
            if dist1 >= dist2:  dist = dist1
            else:               dist = dist2
            
            if (dist <= max_distance) and (dist < best_distance) and (seq_window[:2] == anchor[:2]):
                if seq_window[-3:] != anchor[-3:]:
                    
                    # check for deletion
                    if(seq_window[-4:][:3] == anchor[-3:]):                        
                        seq_window.strip()
                        print "Fuzzy with deletion",seq_window
                    # check for insertion
                    elif(seq_window[-3:] == anchor[-4:][:3]):                        
                        seq_window = seq_window + anchor[-1:]
                        print "fuzzy with insertion", seq_window
                        
                # Found a fuzzy match within tolerances, so store it
                found_fuzzy = 1;
                best_distance = dist;
                best_position = pos;
                fuzzy_match = seq_window;        
                if dist == 0: 
                    found_exact = 1
                    break
    fuzz_right = ''
    if found_fuzzy:
        fuzzy_right = seq
        if( trim_type == 'internal'):
            seq = seq[:best_position + len(fuzzy_match)]    
        else:
            seq = seq[:best_position]
            
        fuzzy_right = fuzz_right[len(seq):]
    return fuzz_right, best_distance, seq, fuzzy_match
        
def levenshtein(s1, s2):
    """ Use the Levenshtein.py module instead """
    if len(s1) < len(s2):
        return levenshtein(s2, s1)
    if not s1:
        return len(s2)
 
    previous_row = xrange(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1 # j+1 instead of j since previous_row and current_row are one character longer
            deletions = current_row[j] + 1       # than s2
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row
 
    return previous_row[-1]

def trim_proximal_primer(primer_list, seq, readDir):
    """Doc string here.."""
    p_primer_found = "";
    offset         = 0;
    #oppositDir     = "";
    #primerDir      = "F";
    #print primer
    seq_len = len(seq)
    #print seq
    for primer in primer_list:
        #print primer
        prim_len = len(primer)
        idx = seq.find(primer)
        if(primer in seq):
            if readDir == 'F':  # look on the left
                if(idx == 0):
                    p_primer_found = primer
                    seq = seq[prim_len:]
                    break
                elif(idx > 0 and idx < 10):
                    p_primer_found = primer
                    offset = idx
                    seq = seq[idx + prim_len:]
                    break
            else:             # look on the right
                hope = seq_len - prim_len
                if(idx == hope):
                    p_primer_found = primer
                    seq = seq[:-prim_len]
                    break
                elif(idx < hope and idx > hope - 10):
                    p_primer_found = primer
                    offset = seq_len - prim_len - idx 
                    seq = seq[:-(prim_len + offset)]
                    break
    return p_primer_found, offset, seq
   

def expand_primers(primer_list):
    """Takes a list of ambiguous primer sequences and expands them to remove ambiguities.
    Adds them to a hash (prevents dupes) and returns them.
    """
    # each primer has p_direction prepended   
    
    
    expandedPrims={}  # fully clean distals - hash
    workingPrimers=[]    # still cleaning distals
    bases = ['A','C','G','T']

    for dir_primer in primer_list:
        #print 'Primer:\n',dir_primer
        # change all '.' to N 
        primer = dir_primer.replace('.','N')
        #print '1',primer
        workingPrimers.append(primer.upper())
    
  
    while (len(workingPrimers) > 0):
        # this is the only place items are removed from working primers
        t = workingPrimers.pop()
        #print 'primer\n',d
        
        # remove prepended direction to expand 
        # then replace it again when re-appending to workingPrimers
        d   = t[2:]  
        dir = t[:1]
        

#     # For each N (blast doesn't like them) expand to 4 distals, one for each base
        if ('N' in d ):    
            for b in bases:      
                workingPrimers.append(dir+':'+d.replace('N',b,1))
                #print '1',d.replace('N',b,1)

#     # For R, Y, W, S, M, K expand to the pair of bases
        elif ('R' in d):
            workingPrimers.append(dir+':'+d.replace('R','A',1))
            workingPrimers.append(dir+':'+d.replace('R','G',1))
            #print '1',d.replace('R','G',1)

        elif ('Y' in d):
            workingPrimers.append(dir+':'+d.replace('Y','C',1))
            workingPrimers.append(dir+':'+d.replace('Y','T',1))        

        elif ('W' in d):
            workingPrimers.append(dir+':'+d.replace('W','A',1))
            workingPrimers.append(dir+':'+d.replace('W','T',1))     

        elif ('S' in d):
            workingPrimers.append(dir+':'+d.replace('S','G',1))
            workingPrimers.append(dir+':'+d.replace('S','C',1))     

        elif ('M' in d):
            workingPrimers.append(dir+':'+d.replace('M','A',1))
            workingPrimers.append(dir+':'+d.replace('M','C',1))     

        elif ('K' in d):
            workingPrimers.append(dir+':'+d.replace('K','G',1))
            workingPrimers.append(dir+':'+d.replace('K','T',1))     

#     # For each [CT] or [AT] ... expand to each base
        elif ('[' in d):
            #print d.find('[')
            #print d.find(']')
            if( d.find('[') + 3 == d.find(']') ):
                base1 = d[d.find('[')+1:d.find('[')+2]
                #print base1
                base2 = d[d.find('[')+2:d.find('[')+3]
                #print base2
                replace = '['+base1+base2+']'
                workingPrimers.append(dir+':'+d.replace(replace,base1,1))
                workingPrimers.append(dir+':'+d.replace(replace,base2,1)) 
                #print '1',d.replace(replace,base1,1)
                #print '2',d.replace(replace,base2,1)

#     # expand \?
        elif ('?' in d):
            if(d.find('?') == 0):
                workingPrimers.append(dir+':'+d[1:])
                #print '1',d[1:]
            else:
                preceder_plus = d[d.find('?')-1:d.find('?')+1] 
                #print preceder_plus
                # the preceing base exists: remove '?' only
                workingPrimers.append(dir+':'+d.replace('?','',1))
                # preceding base doesn't exist: 'remove '?' and preceding base
                workingPrimers.append(dir+':'+d.replace(preceder_plus,'',1))
                #print '2',d.replace('?','',1)
                #print '3',d.replace(preceder_plus,'',1)

#     # next expand + to 1 or 2
        elif ('+' in d):
            preceder = d[d.find('+')-1:d.find('+')]
            # the preceing base exists once: remove '+' only
            workingPrimers.append(dir+':'+d.replace('+','',1))
            # the preceing base exists twice: change '+' to preceding base
            workingPrimers.append(dir+':'+d.replace('+',preceder,1))

#     # expand * to 0,1,2
        elif ('*' in d):
            if(d.find('*') == 0):
                workingPrimers.append(dir+':'+d[1:])
            else:
                preceder = d[d.find('*')-1:d.find('*')]
                preceder_plus = d[d.find('*')-1:d.find('*')+1] 
                #print preceder,preceder_plus
                # the preceding base doesn't exist': remove '*' and preceding base
                workingPrimers.append(dir+':'+d.replace(preceder_plus,'',1))
                # the preceding base exists once: remove '*' only
                workingPrimers.append(dir+':'+d.replace('*','',1))
                # the preceding base exists twice: change '*' to preceding base
                workingPrimers.append(dir+':'+d.replace('*',preceder,1))

#     # For repeat bases, e.g., C{5,8} becomes 4 new primers: Cx5, Cx6, Cx7, Cx8
        elif('{' in d):
            if(d.find('}') == d.find('{') + 4 ):
                repeatBase = d[d.find('{')-1:d.find('{')]
                toReplace = d[d.find('{')-1:d.find('}')+1]
                #print toReplace
                #print 'repeatbase',repeatBase
                minCount = d[ d.find('{') + 1:d.find('{') + 2 ]
                maxCount = d[ d.find('{') + 3:d.find('{') + 4 ]
                #print minCount, maxCount
                for i in range (int(minCount),int(maxCount) + 1 ):
                    homopolymer = i * repeatBase
                    workingPrimers.append(dir+':'+d.replace(toReplace, homopolymer))
                
 
        # If it made it through everything else, move it to the final set
        # Use hash so can filter out duplicates
        else:
            expandedPrims[d] = dir
    
#   if ($verbose) {print join("\n", "Orig @$dist", keys %expandedDistals) . "\n";}
    return expandedPrims
    
def expand(seq):
    """Takes a single ambiguous primer sequence and expands it to remove ambiguities.
    Adds them to a hash (prevents dupes) and returns the hash.
    """
    expandedPrims={}  # fully clean distals - hash to prevent dupes
    workingPrimers=[]    # still cleaning holder
    bases = ['A','C','G','T']

    
    seq = seq.replace('.','N')
        #print '1',primer
    workingPrimers.append(seq.upper())
    
  
    while (len(workingPrimers) > 0):
        # this is the only place items are removed from working primers
        d = workingPrimers.pop() 
        
        # remove all commas
        d = d.replace(',','')
        
#     # For each N (blast doesn't like them) expand to 4 distals, one for each base
        if ('N' in d ):    
            for b in bases:      
                workingPrimers.append(d.replace('N',b,1))

#     # For R, Y, W, S, M, K expand to the pair of bases
        elif ('R' in d):
            workingPrimers.append(d.replace('R','A',1))
            workingPrimers.append(d.replace('R','G',1))
            #print '1',d.replace('R','G',1)

        elif ('Y' in d):
            workingPrimers.append(d.replace('Y','C',1))
            workingPrimers.append(d.replace('Y','T',1))        

        elif ('W' in d):
            workingPrimers.append(d.replace('W','A',1))
            workingPrimers.append(d.replace('W','T',1))     

        elif ('S' in d):
            workingPrimers.append(d.replace('S','G',1))
            workingPrimers.append(d.replace('S','C',1))     

        elif ('M' in d):
            workingPrimers.append(d.replace('M','A',1))
            workingPrimers.append(d.replace('M','C',1))     

        elif ('K' in d):
            workingPrimers.append(d.replace('K','G',1))
            workingPrimers.append(d.replace('K','T',1))     

#     # For each [CT] or [AGT] ... expand to each base
        elif ('[' in d):
            
            if( d.find('[') + 3 == d.find(']') ):
                # ie [AC]
                base1 = d[d.find('[')+1:d.find('[')+2]
                base2 = d[d.find('[')+2:d.find('[')+3]
                replace = '['+base1+base2+']'
                workingPrimers.append(d.replace(replace,base1,1))
                workingPrimers.append(d.replace(replace,base2,1)) 
            
            if( d.find('[') + 4 == d.find(']') ):
                # ie: [CTA]
                base1 = d[d.find('[')+1:d.find('[')+2]
                base2 = d[d.find('[')+2:d.find('[')+3]
                base3 = d[d.find('[')+3:d.find('[')+4]
                replace = '['+base1+base2+base3+']'
                workingPrimers.append(d.replace(replace,base1,1))
                workingPrimers.append(d.replace(replace,base2,1))       
                workingPrimers.append(d.replace(replace,base3,1))  

#     # expand \?
        elif ('?' in d):
            if(d.find('?') == 0):
                workingPrimers.append(d[1:])
            else:
                preceder_plus = d[d.find('?')-1:d.find('?')+1] 
                #print preceder_plus
                # the preceing base exists: remove '?' only
                workingPrimers.append(d.replace('?','',1))
                # preceding base doesn't exist: 'remove '?' and preceding base
                workingPrimers.append(d.replace(preceder_plus,'',1))

#     # next expand + to 1 or 2
        elif ('+' in d):
            preceder = d[d.find('+')-1:d.find('+')]
            # the preceing base exists once: remove '+' only
            workingPrimers.append(d.replace('+','',1))
            # the preceing base exists twice: change '+' to preceding base
            workingPrimers.append(d.replace('+',preceder,1))

#     # expand * to 0,1,2
        elif ('*' in d):
            if(d.find('*') == 0):
                workingPrimers.append(d[1:])
            else:
                preceder = d[d.find('*')-1:d.find('*')]
                preceder_plus = d[d.find('*')-1:d.find('*')+1] 
                #print preceder,preceder_plus
                # the preceding base doesn't exist': remove '*' and preceding base
                workingPrimers.append(d.replace(preceder_plus,'',1))
                # the preceding base exists once: remove '*' only
                workingPrimers.append(d.replace('*','',1))
                # the preceding base exists twice: change '*' to preceding base
                workingPrimers.append(d.replace('*',preceder,1))

#     # For repeat bases, e.g., C{5,8} becomes 4 new primers: Cx5, Cx6, Cx7, Cx8
        elif('{' in d):
            if(d.find('}') == d.find('{') + 4 ):
                repeatBase = d[d.find('{')-1:d.find('{')]
                toReplace = d[d.find('{')-1:d.find('}')+1]
                minCount = d[ d.find('{') + 1:d.find('{') + 2 ]
                maxCount = d[ d.find('{') + 3:d.find('{') + 4 ]
                #print minCount, maxCount
                for i in range (int(minCount),int(maxCount) + 1 ):
                    homopolymer = i * repeatBase
                    workingPrimers.append(d.replace(toReplace, homopolymer))
                
 
        # If it made it through everything else, move it to the final set
        # Use hash so can filter out duplicates
        else:
            if(d):
                expandedPrims[d] = 1

    return expandedPrims.keys()  
    
def get_expanded_primers(lane_key, sample_primers):
    
    orig_primer_seqs = {}  
    expanded_primers = {}
    if(os.path.exists(sample_primers) and os.path.isfile(sample_primers)):
        f = open(sample_primers, 'r')
        x = f.readlines()
        orig_primer_seqs[lane_key] = [p.strip() for p in x]
        
    else: # this has to be a single primer seq or comma separated list
        orig_primer_seqs[lane_key] = sample_primers.split(',')
        
    #for dir_raw_primer in proximal_primers[lane_key]:
    #p_direction, raw_primer = dir_raw_primer.split(':')
    # goal: expanded_primers[lane_key][p_direction] = [list of expanded primers]
    expanded_primers[lane_key] = {}
    #print lane_key
    expanded = expand_primers(orig_primer_seqs[lane_key])
    expanded_primers[lane_key]['F'] = []
    expanded_primers[lane_key]['R'] = []
    for p_list in expanded:
        if not p_list: continue
        if (expanded[p_list] == 'F'):
            expanded_primers[lane_key][expanded[p_list]].append(p_list) 
            # reverse complement here and add to this list also
            base_list = p_list[::-1]  # reverses primer
            xxx = [complementary_bases[base] for base in list(base_list)]
            expanded_primers[lane_key][expanded[p_list]].append(''.join(xxx) )
        elif(expanded[p_list] == 'R'):
            expanded_primers[lane_key][expanded[p_list]].append(p_list)
            # reverse complement here and add to this list also                
            base_list = p_list[::-1]  # reverses primer
            xxx = [complementary_bases[base] for base in list(base_list)]
            expanded_primers[lane_key][expanded[p_list]].append(''.join(xxx) )
            #print 'comp',item,''.join(xx)
            
        else:
            pass
            #print "ERROR expanding primers"
            print "no primer direction found for", p_list
    
    return expanded_primers[lane_key]
    
def get_anchor_list(anchor_name, adtnl_anchors_list):
    raw_sequence = C.mbl_anchors[anchor_name]['sequence']
    expanded_seqs = expand(raw_sequence)

    for a in adtnl_anchors_list:
        expanded_seqs = expanded_seqs + expand(a)

    return expanded_seqs