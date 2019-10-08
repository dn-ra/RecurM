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
from Bio.Alphabet import IUPAC


'''repeatm modules'''
import delta_parse


'''other'''
#import intervals
import os
import sys
import intervals
import csv

'''set locations and files'''
#directory of bin files
bin_dir = '/srv/home/s4204666/abisko/aterrible_bins/12_assembly73_individuals_flat20160115/bins_unique_contig_names.shorter' #location of bin fasta files
#names of derep bins
derep_bins_file = '/srv/home/s4204666/abisko/aterrible_bins/12_assembly73_individuals_flat20160115/dereplicated_genomes.txt'
#fasta file to save to
exit_bin_file = 'for_coverm.fa'
#disk location of all_repseqs.fa
repseq_locs = {'Perfect':'/srv/home/s4204666/abisko/dan/repeatm_tests/all_assemblies/cluster_complete/FINAL_perfectlinear_vsbins/ALL_PEAK_PERFECT.fa', 
               'Circular':'/srv/home/s4204666/abisko/dan/repeatm_tests/all_assemblies/set_scores/other_tools/circular_superior.fa'}
#directories where binvscluster nucmer results are stored
delta_dirs = ['/srv/home/s4204666/abisko/dan/repeatm_tests/all_assemblies/cluster_complete/FINAL_circular_vsbins', 
              '/srv/home/s4204666/abisko/dan/repeatm_tests/all_assemblies/cluster_complete/FINAL_perfectlinear_vsbins'] 


'''functions in use here'''
def remove_contig(match_object, target = 2):
    #target is 2 if bin_file sequences is the second seuqence of the match object (expected)
    #else, set as 1    
    
    remove = False
#    if match_object.get_ani() < 0.9:
#        return remove #early exit if ani fails threshold
    
    for i in range(len(match_object)):
        coords = [getattr(m, 'hitstarts_{}'.format(target))[i] -1, getattr(m, 'hitstops_{}'.format(target))[i] -1] #minus one for 0-based indexing
        span = intervals.closed(min(coords),max(coords))
        if i ==0:
            c_span = span
        else:
            c_span = c_span | span
    total_length = sum([ivl.upper - ivl.lower for ivl in c_span])
    if total_length > (match_object.lengths[target-1] * 0.90): #if whole match area is >90% of the sequence
        remove = True
    
    return remove

    

'''start process'''

print('building bin linkages')
sys.stdout.flush()

derep_bins = []
f = open(derep_bins_file, 'r')
for line in f:
    derep_bins.append(line.strip()+'.fna')
f.close()

print('{} bins submitted for processing'.format(len(derep_bins)))

print('processing bin delta files', flush=True)
alldeltas = []
for directory in delta_dirs:
    for file in os.listdir(directory):
        if file.endswith('.delta'):
            alldeltas.append(delta_parse.deltaread(os.path.join(directory, file)))


'''link match objects to each bin'''

print('building cluster-bin linkages',flush=True)
sys.stdout.flush()

bin_matches = {} #dictionary of each repeatedly assembled sequence with it's associated match from the bins
#forge connection between contigs and their bins
for d in alldeltas:
    for key, value in d.items():
        bin_ref = key.split("---")[1]
        for m in value:
            m.seqs[1] = bin_ref+"__"+m.seqs[1]
           
           #get bin location for each match object
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

f = open('bin_found_summary.txt', 'w')
for value in bin_finds.values():
    for m in value:
        f.write(m.seqs[0]+'\t'+m.seqs[1]+'\n')

f.close()



single_bin_match = []
multiple_bin_match = []
no_bins = []
bin_multi_seqs = {}

for k,v in bin_finds.items():
    if len(v) ==1:
        single_bin_match.append(k)
    elif len(v) >1:
        multiple_bin_match.append(k)
    elif len(v) == 0:
        no_bins.append(k)
    for m in v:
        try:
            bin_multi_seqs[m.seqs[1].split('__')[0]].append(m.seqs[0])
        except KeyError:
            bin_multi_seqs[m.seqs[1].split('__')[0]] = [m.seqs[0]]

f = open('bins_with_multiple_seqs', 'w')
w= csv.writer(f, delimiter="\t")
for key, value in bin_multi_seqs.items():
    if len(value) >1:
        w.writerow([key, value])
f.close()

print('Identified bin-cluster linkages')
print()
print('{} sequences map to bins'.format(len(single_bin_match)+len(multiple_bin_match)))

if multiple_bin_match:
    for key in multiple_bin_match:
        if len(set([m.seqs[1] for m in bin_finds[key]]))>1:
            print('Warning! {} was found in more than one bin'.format(key))

sys.stdout.flush()  

'''construct dictionary of sequences to remove from bins'''
bin_contigs_remove = {}

for matches in bin_finds.values():
    for m in matches:
        if remove_contig(m): # a more robust measurement that accounts for overlapping regions. Probably not necessary at this point seeing that it has already passed the apply_threshold step. but this intervals method should really be replacing the apply_threshold method in future
            bin_ref, seq_name = m.seqs[1].split("__")
            try:
                bin_contigs_remove[bin_ref].append(seq_name)
            except KeyError:
                bin_contigs_remove[bin_ref] = [ seq_name ]
       


'''Go through all bin files. Copy and amend names to contain linkages. remove repeat elements. Save total fasta file as for_coverm.fa'''


f = open(exit_bin_file, 'w')
cnt=0
for file in derep_bins:
    print('processing {} in dereplicated bins'.format(file), flush=True)
    sys.stdout.flush()
    
    if file in bin_contigs_remove.keys():
        print('Removing {} sequences from {}'.format(len(bin_contigs_remove[file]), file), flush = True)
        for seq in SeqIO.parse(handle = os.path.join(bin_dir, file), format = 'fasta', alphabet=IUPAC.unambiguous_dna):
                if seq.id in bin_contigs_remove[file]:
                    cnt+=1
                    continue #don't copy. It's in my repeated elements
                    
                else:
                    seq.id = file+"~"+seq.id
                    SeqIO.write(seq, f, format='fasta')
    else:
        for seq in SeqIO.parse(handle = os.path.join(bin_dir, file), format = 'fasta', alphabet=IUPAC.unambiguous_dna):
            seq.id = file+"~"+seq.id
            SeqIO.write(seq, f, format='fasta')    
            
print('{} sequences removed from bins'.format(cnt), flush=True)
#has to be a more succint way to write all this??^^

#write in all the cluster sequences that weren't thrown in with the bins

seen_set = set() #there are two sets of clusters - perfect and circular - so there will be duplicates. Save ids in here to skip if they are seen again

for typ, file in repseq_locs.items():
    for num, seq in enumerate(SeqIO.parse(handle = file, format='fasta', alphabet = IUPAC.unambiguous_dna)):
        if seq.id not in seen_set:
            seen_set.add(seq.id)
            seq.id = '{}~cluster_{}_{}'.format(seq.id,typ, num+1)
            seq.name = None
            SeqIO.write(seq, f, format = 'fasta')
    
 
f.close()
print('Bin linking Complete! File for processing in CoverM stored at {}'.format(exit_bin_file))
