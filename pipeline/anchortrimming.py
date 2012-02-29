#!/usr/local/www/vamps/software/python/bin/python
# -*- coding: utf-8 -*-
#
# See README file for details. Send your questions / remarks to "A. Murat Eren", <meren / mbl.edu>


import sys

try:
    import Levenshtein
except:
    print '''
    You need Levenshtein module installed to run this software.

    Here is a fast implementation of Levenshtein distance for Python:

        http://code.google.com/p/pylevenshtein/

'''
    sys.exit(-1)


class Settings:
    """
    Settings class with templates for different regions of 16S rRNA gene (see __main__ for client side usage).

    * reversed: if the sequence is reverse sequenced.
    * start   : nucleotide position to start the search (if reversed is true, it is -(start),
                means not from the beginning but from the end of the sequence)
    * freedom : the length of the search space from both directions.
    * length  : expected length of the anchor sequence
    """

    def __init__(self, region = None):
        self.general_settings = {
              'v6v4-361': {'reversed': True,  'start'   : 361, 'freedom' : 50, 'length'  : 13},
                            # previously determiend anchor consensus: G[T,G]AG.[A,G]GT[A,G][A,G]AAT
              'v6v4-4xx': {'reversed': True,  'start'   : 480, 'freedom' : 60, 'length'  : 11},
                            # previously determiend anchor consensus: G[T,G]AG.[A,G]GT[A,G][A,G]AAT
              'v3v5-440': {'reversed': False, 'start'   : 440, 'freedom' : 30, 'length'  : 13},
                            # previously determined anchor consensus: GGATTAGA[T,G]ACCC
              'v3v5-370': {'reversed': False, 'start'   : 370, 'freedom' : 50, 'length'  : 12},
                            # previously determined anchor consensus: [A,T,C][A,T,G]GCGAA[A,G]GC[A,G][A,C,G]
        }

        if region:
            self.region_settings = DictDotNotationWrapper(self.general_settings[region])

    def available_regions(self):
       return self.general_settings.keys()


class DictDotNotationWrapper(object):
    """
    This is just a wrapper class for syntactic convenience (if your dictionary instance is from this class,
    basically you can reach items in the dictionary with dot notation) (so you can ignore this one if
    you're trying to understand what this program does).
    """
    
    def __init__(self, dictionary):
        self.dictionary = dictionary

    def __getattr__(self, key):
        return self.dictionary.get(key, None)

    def __setattr__(self, key, val):
        super.__setattr__(self, key, val)


def generate_tuples(start, freedom, length, direction = -1, step = 0, list_of_tuples = [], reversed_read = False):
    """
    This recursive function crates a list of tuples that are expanding from 'start' in both
    directions until they reach the borders of 'freedom' that sent as a parameter. 'length' is
    the distance from the start. if you call generate_tuples(10, 5, 3), you would get this:

    [(10, 14), (11, 15), (9, 13), (12, 16), (8, 12), (13, 17), (7, 11), (14, 18), (6, 10), (15, 19), (5, 9)]

    maybe this would give a better idea:

                                       (10, 14)
                                (9, 13)        (11, 15)
                         (8, 12)                       (12, 16)
                  (7, 11)                                      (13, 17)
           (6, 10)                                                     (14, 18)
     (5, 9)                                                                    (15, 19)

    so, if you know where it is likely to find the pattern, you can start from there, and expand search to both
    directions step by step in a somewhat optimized way.
    """

    if reversed_read:
        list_of_tuples.append((-(start), -(start - length)))
    else:
        list_of_tuples.append((start, start + length))
    
    if step == freedom * 2:
        return list_of_tuples

    direction *= -1
    step += 1

    return generate_tuples(start + direction * step, freedom, length, direction, step, list_of_tuples, reversed_read = reversed_read)


def trim_sequence(sequence, location, s):
    if s.reversed:
        return sequence[location[0]:] # this includes the anchor in trimmed sequence..
    else:
        return sequence[:location[1]] # same thing here for the reversed == False


def find_best_distance(sequence, valid_anchor_sequences, max_divergence, list_of_tuples):
    best_loc_for_every_anchor = []
    for anchor in valid_anchor_sequences:
        non_perfect_matches = []
        for tpl in list_of_tuples:
            r = Levenshtein.ratio(sequence[tpl[0]:tpl[1]], anchor)
            if r == 1:
                # r is 1, means that we found the perfect match to one of the valid anchors.
                # we'll return it imediately.
                return (anchor, tpl)
            elif r >= max_divergence:
                # r is not 1, but it is larger than the 'max_divergence' accepted.
                # we'll put it in a list and keep testing.
                non_perfect_matches.append((r, tpl))

        if len(non_perfect_matches):
            best_loc_for_every_anchor.append((anchor, sorted(non_perfect_matches, reverse=True)[0]))
    
    # potential FIXME: OK. at this point, 'best_loc_for_every_anchor' is such a list that every member of this 
    # list is an anchor from 'valid_anchor_sequences' and the location within the 'sequence' where the sequence
    # is most similar to this anchor. The data structure for every item in this list looks like this:
    #
    #    ('GGAGCGGTGGAAT', (0.7692307692307693, (-344, -331)))
    #
    # it reads "the best distance of GGAGCGGTGGAAT to every oligonucleotide in the given sequence was at sequence[-344:-331]"
    # 
    # but it is important to remember that this 'best' is coming from a sorted list. so, there may be equally good ones that left out
    # and maybe some of them were better candidates. this issue becomes even more of a challenge when we pick 'best_anchor' in the next
    # line by sorting this list. just by chance, from two equaly good candidates, the one that could be more prefferable in terms of the
    # trimming location in the sequence might be beaten by another one during the sorting. at some point we might want to plug
    # in a probabilistic logic here to pick competing locations (maybe equally distant options could be ranked based on a 
    # pre-computed mixture of gaussian curves for a given region to pick the most appealing location to use as an anchor):

    best_loc_for_every_anchor_sorted = sorted(best_loc_for_every_anchor, key = lambda k: k[1][0], reverse=True)
    if len(best_loc_for_every_anchor_sorted):
        best_anchor = best_loc_for_every_anchor_sorted[0]
        return (best_anchor[0], best_anchor[1][1])
    else:
        return (None, None)

def colorize(sequence, location):
    Green = lambda s: '\033[30m\033[42m' + s + '' + '\033[0m'
   
    return sequence[0:location[0]] + Green(sequence[location[0]:location[1]]) + sequence[location[1]:]

def main(s):
    # define the search space
    list_of_tuples = generate_tuples(s.start, s.freedom, s.length, reversed_read = s.reversed)

    input = s.input_fasta

    while input.next():
        anchor, location = find_best_distance(input.seq, s.valid_anchor_sequences, s.max_divergence, list_of_tuples)

        if anchor and location:
            # lets get that anchor at the beginning of the list..
            # this step will make the algorithm get faster as it carries frequently
            # used anchors at the beginning of the list
            if s.valid_anchor_sequences.index(anchor) != 0:
                s.valid_anchor_sequences.insert(0, s.valid_anchor_sequences.pop(s.valid_anchor_sequences.index(anchor)))
            
            # busines time.
            trimmed = trim_sequence(input.seq, location, s)

            if s.screen:
                s.output.write('>' + input.id + '\n')
                s.output.write(colorize(input.seq, location) + '\n')
            else:
                s.output.write('>' + input.id + '\n')
                s.output.write(trimmed + '\n')
                if input.pos % 100 == 0:
                    sys.stderr.write('\rTRIMMING: %.2d%% -- pos: %d' % (input.pos * 100 / input.total_seq, input.pos))
                    sys.stderr.flush()
        else:
            if s.screen:
                s.failed.write('>' + input.id + '\n')
                s.failed.write(input.seq + '\n')
                if input.pos % 100 == 0:
                    sys.stderr.write('\rTRIMMING: %.2d%% -- pos: %d' % (input.pos * 100 / input.total_seq, input.pos))
                    sys.stderr.flush()
            else:
                s.failed.write('>' + input.id + '\n')
                s.failed.write(input.seq + '\n')
    print



if __name__ == '__main__':
    import argparse
    import fastalib as u

    parser = argparse.ArgumentParser(description='Fuzzy anchor trimming for 454 Sequences')
    parser.add_argument('-i', '--input-fasta', required=True, metavar = 'FASTA_FILE',
                        help = 'Sequences file to be trimmed in FASTA format')
    parser.add_argument('-r', '--region', required=True, metavar = 'REGION',
                        help = 'Region in the 16S rRNA gene. Available options: %s' % 
                                    ', '.join(Settings().available_regions()), choices = Settings().available_regions())
    parser.add_argument('-a', '--anchor-sequences', required=True, metavar = 'ANCHORS_FILE',
                        help = 'Input file that contains the list of valid anchor sequences')
    parser.add_argument('-o', '--output', help = 'Where trimmed sequences will be written (default: standart output)')
    parser.add_argument('-d', '--max-divergence', type=float, default=0.9,
                        help = 'Maximum Levenshtine distance allowed candidate trimming site from one of the valid\
                                anchor sequence (default: 0.90). Please see README file in order to get more\
                                information on maximum divergence.')

    args = parser.parse_args()
   
    s = Settings(args.region).region_settings

    s.input_fasta = u.SequenceSource(args.input_fasta)
    s.valid_anchor_sequences = [sequence.strip() for sequence in open(args.anchor_sequences).readlines()]
    s.max_divergence = args.max_divergence

    if args.output:
        s.screen = False
        s.output = open(args.output, 'w')
        s.failed = open(args.output + '-FAILED', 'w')
    else:
        s.screen = True
        s.output = sys.stdout
        s.failed = sys.stderr

    sys.exit(main(s))
