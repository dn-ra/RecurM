# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 15:05:15 2019

A script to remove repeatedly assembled sequences from bins. For use in determining coverage patterns between these repeat fragments and associated bins

Important for dev:
    the central feature of 

@author: dan_r
@email: daniel.rawlinson@uqconnect.edu.au
"""

'''biopython modules'''
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


'''repeatm modules'''
import delta_parse


'''other'''
#import intervals
import os
import sys

'''set locations and files'''
bin_dir = '/srv/home/s4204666/abisko/aterrible_bins/12_assembly73_individuals_flat20160115/bins_unique_contig_names.shorter' #location of bin fasta files
derep_bins_file = '/srv/home/s4204666/abisko/aterrible_bins/12_assembly73_individuals_flat20160115/dereplicated_genomes.txt'
exit_bin_file = 'for_coverm.fa'
repseqs_loc = '/srv/home/s4204666/abisko/dan/repeatm_tests/all_assemblies/cluster_test_unionfind/size_50+_all_repseqs.fa' #disk location of all_repseqs.fa


'''functions in use here'''
def remove_contig(match_object, target = 2):
    #target is 2 if bin_file sequences is the second seuqence of the match object (expected)
    #else, set as 1
    remove = False
    for i in len(match_object):
        span = intervals.closed(getattr(m, 'hitstarts_{}'.format(target))[i] -1, getattr(m, 'hitstops_{}'.format(target))[i] -1) #minus one for 0-based indexing
        if i ==0:
            c_span = span
        else:
            c_span = c_span | span
    total_length = sum([ivl.upper - ivl.lower for ivl in c_span])
    if total_length > (match_object.lengths[target-1] * 0.90): #if whole match area is >90% of the sequence
        remove = True
    
   
    #this stuff below is from when I was trying to slice out segments     
        #spliced_seq =  Seq('', IUPAC.unambiguous_dna) #remove the whole sequence
        
#    else:
#        spliced_seq = seqs[match_object.seqs[ target -1 ]]
#        for ivl in c_span:
#            #TODO - here to set a boundary for how small a match to remove
#            if ivl.upper - ivl.lower >0:
#                print(ivl)
#                spliced_seq = spliced_seq.seq[ 0:ivl.lower ] + spliced_seq.seq[ ivl.upper:: ]
    
    return remove


'''start process'''

print('building bin linkages')
sys.stdout.flush()

derep_bins = []
f = open(derep_bins_file, 'r')
for line in f:
    derep_bins.append(line.strip()+'.fna')
f.close()



alldeltas = []
for file in os.listdir():
    if file.endswith('.delta'):
        alldeltas.append(delta_parse.deltaread(file))

bin_contigs_remove = {}

bin_matches = {} #dictionary of each repeatedly assembled sequence with it's associated match from the bins
#forge connection between contigs and their bins
for d in alldeltas:
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

bin_finds = {} #dictionary of # of bins where each node (the key) is found at threshold value
for key, value in bin_matches.items():
    match_list = []
    for m in value:
        if m.apply_threshold(threshold = 0.90) == True:
            match_list.append(m)
    bin_finds[key] = match_list
        

'''for reporting on stats of matches to bins'''
single_bin_match = []
multiple_bin_match = []
no_bins = []

for k,v in bin_finds.items():
    if len(v) ==1:
        single_bin_match.append(k)
    elif len(v) >2:
        multiple_bin_match.append(k)
    elif len(v) == 0:
        no_bins.append(k)

'''construct dictionary of sequences to remove from bins'''
for matches in bin_finds.values():
    for m in matches:
        #if remove_contig(m): # a more robust measurement that accounts for overlapping regions. Probably not necessary at this point seeing that it has already passed the apply_threshold step. but this intervals method should really be replacing the apply_threshold method in future
        bin_ref, seq_name = m.seqs[1].split("__")
        try:
            bin_contigs_remove[bin_ref] += seq_name
        except KeyError:
            bin_contigs_remove[bin_ref] = [ seq_name ]
       


'''Go through all bin files. Copy and amend names. remove repeat elements.'''

f = open(exit_bin_file, 'w')
for file in derep_bins:
    print('processing {} in dereplicated bins'.format(file))
    sys.stdout.flush()
    for seq in SeqIO.parse(handle = os.path.join(bin_dir, file), format = 'fasta', alphabet=IUPAC.unambiguous_dna):
        if seq.id in bin_contigs_remove[file]:
            continue
        else:
            seq.id = file+"~"+seq.id
            SeqIO.write(seq, f, format='fasta')
    
                    
for num, seq in enumerate(SeqIO.parse(handle = repseqs_loc, format='fasta', alphabet = IUPAC.unambiguous_dna)):
    seq.id = 'cluster_{}~{}'.format(num, seq.id)
    SeqIO.write(seq, f, format = 'fasta')

f.close()







'''unused functions'''
        
def copy_edited_bin():
    '''export whole new bin file with edited contigs'''
    return

def find_larger(): #and find smaller? #is this better to have in the cluster object?
    return
