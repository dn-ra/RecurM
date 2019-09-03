# -*- coding: utf-8 -*-
"""
Main pipeline for Nucmer_Matching from .delta files to cluster generation

author: Daniel Rawlinson, ACE
email: daniel.rawlinson@uqconnect.edu.au
"""
import single_linkage_cluster
import delta_parse
import os
import pickle
import union_find_cluster

#set locaiton of delta files
delta_dir = '/srv/home/s4204666/abisko/dan/repeatm_tests/all_assemblies/filtered_2000bp/nucmer_feed_out'

#set location of assemblies
assembly_dir = '/srv/home/s4204666/abisko/dan/repeatm_tests/all_assemblies/filtered_2000bp'

#set location of bam files
bam_dir = ''

cluster_method = 'union-find' #set to union-find or single-linkage

#set empty vector of all nucmer_matches that pass threshold (=0.90 by default)
collated_sig_matches = []
sigmatch_set = set()
#set empty vector to collect all nucmer_matches that are fragments of larger assemblies
fragments = []

#read in nucmer matches of each deltafile
for file in os.listdir(delta_dir):
#if file.endswith('.delta'):
    print('parsing {}'.format(file))
    delta = delta_parse.deltaread(delta_dir+'/'+file)
    #exit type is dictionary
    
    #threshold all matches and collate
    print('thresholding {}'.format(file))
    firstpass_fragments = []
    for m in next(iter(delta.values())):
        stats = m.gen_statistics()
        if m.apply_threshold(treshold = 0.90, stats = stats) == True:
            collated_sig_matches.append(m)
            sigmatch_set.add([m.seqs])
        elif m.is_fragment(upperthreshold = 0.90, lowerthreshold = 0.90, stats = stats):
            firstpass_fragments.append(m)

    #all match objects are read in. keep only those fragments that map to full matches
    for m in firstpass_fragments:
        if m.seqs[0] in sigmatch_set:
            fragments.append(m)
        elif m.seqs[1] in sigmatch_set:
            fragments.append(m)
                    
        #exit type is list of Nucmer_Match objects (sig_matches) + list of Nucmer_Match objects (fragments) that map to the sig_match objects

#save progress
f = open('pickled_sigmatches', 'wb')
pickle.dump(collated_sig_matches, f)
f.close()

f = open('pickled_fragmatches', 'wb')
pickle.dump(fragments, f)
f.close()

#cluster all significant matches & sort
if cluster_method == 'single-linkage': #now deprecated 
    clusters = single_linkage_cluster.cluster_nucmer_matches(collated_sig_matches)
    
elif cluster_method == 'union-find':
    clusters = union_find_cluster.union_find_pipe(collated_sig_matches)

single_linkage_cluster.sort_clusters(clusters)

f = open('pickled_cluster_objs', 'wb')
pickle.dump(clusters, f)
f.close()
#exit type is list of Contig_Cluster objects

#find fragments linked to clusters
#frag_seqs = [m.seqs.index(min(m.lengths)) for m in fragments]
#
#
#
#
##run analysis on clusters
#for c in clusters:
#    c.retrieve_seqs(assembly_dir = assembly_dir) #locations set at beginning of script. Does not return anything
#    
#f = open('pickled_clusters', 'wb')
#pickle.dump(clusters, f)
#f.close()


    
##determine orientation of matches in the cluster
