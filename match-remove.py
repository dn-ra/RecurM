# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 15:05:15 2019

A script to remove repeatedly assembled sequences from bins. For use in determining coverage patterns between these repeat fragments and associated bins

Important for dev:
    the central feature of 

@author: dan_r
@email: daniel.rawlinson@uqconnect.edu.au
"""

'''bipython'''
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC

'''repeatm modules'''
import delta_parse
import single_linkage_cluster


'''other'''
import intervals

'''for testing'''
file = 'C:\\Users\\dan_r\\Documents\\Honours_Data\\2.delta'

bin_dir = '/srv/home/s4204666/abisko/aterrible_bins/12_assembly73_individuals_flat20160115/bins' #location of bin fasta files

exit_bin_dir = '' #location to save fasta file of edited bins

d = delta_parse.deltaread(file)

g_file = {}



seqs = {}
bin_matches = {} #dictionary of each repeatedly assembled sequence with it's associated match from the bins
#forge connection between contigs and their bins
for key, value in d.items():
    bin_ref = key.split("---")[1]
    for m in value:
        m.seqs[1] = bin_ref+"__"+m.seqs[1]
        #extract sequences for manipulation
       # seqs [m.seqs[1]] = retrieve_bin_seq(bin_dir, m) #seq_object dict to be passed into remove_span function
        try:
            bin_matches[m.seqs[0]].append(m)
        except KeyError:
            bin_matches[m.seqs[0]] = [m]

bin_counts = {} #dictionary of # of bins where each node (the key) is found at threshold value
for key, value in bin_matches.items():
    count = 0
    for m in value:
        if m.apply_threshold(threshold = 0.90) == True:
            count +=1
    bin_counts[key] = count

bin_found = [k for k,v in bin_counts.items() if v ==1]

for k in bin_found:
    all_stats = [m.gen_statistics() for m in bin_matches[k]]
    #find any more matches in which the whole of the fragment forms part of a larger assembled contig. That would be bad news.
    
    
'''functions in use here'''    
def retrieve_bin_seq(bin_dir, match_obj):
    
    bin_file, node = match_obj.seqs[1].split("__")
    seq = ''
    with open("/".join([bin_dir, bin_file])) as f:
        for i in f:
            if i.find(node)!=-1:
                wholeseq =False #flag to tell me if I've taken the whole sequence yet
                while wholeseq ==False:
                    line = next(f)
                    if line.startswith('>'):
                        wholeseq = True
                    else:
                        seq += line.strip()
                        
    return SeqIO.SeqRecord(Seq(seq, IUPAC.unambiguous_dna), id= node)

def remove_span(match_object, seq_object_dict, target = 2):
    #TODO - add options in here to completely remove a sequence if there is barely any of it left
    #target is 2 if bin_file sequences is the second seuqence of the match object (expected)
    #else, set as 1
    for i in len(match_object):
        span = intervals.closed(getattr(m, 'hitstarts_{}'.format(target))[i] -1, getattr(m, 'hitstops_{}'.format(target))[i] -1) #minus one for 0-based indexing
        if i ==0:
            c_span = span
        else:
            c_span = c_span | span
    total_length = sum([ivl.upper - ivl.lower for ivl in c_span])
    if total_length > (match_object.lengths[target-1] * 0.90): #if whole match area is >90% of the sequence
        spliced_seq =  Seq('', IUPAC.unambiguous_dna) #remove the whole sequence
    else:
        spliced_seq = seqs[match_object.seqs[ target -1 ]]
        for ivl in c_span:
            #TODO - here to set a boundary for how small a match to remove
            if ivl.upper - ivl.lower >0:
                print(ivl)
                spliced_seq = spliced_seq.seq[ 0:ivl.lower ] + spliced_seq.seq[ ivl.upper:: ]
    
    return spliced_seq             
        
def copy_edited_bin():
    '''export whole new bin file with edited contigs'''
    return

def coverm_genome_file():
    '''retain associaitons between bins and contigs to feed into coverm'''
    return