# -*- coding: utf-8 -*-
"""
For use in RepeatM. See: github.com/wwood/RepeatM
Main pipeline for Nucmer_Matching from .delta files to cluster generation


author: Daniel Rawlinson, ACE
email: daniel.rawlinson@uqconnect.edu.au
"""
'''imports'''

import single_linkage_cluster
import delta_parse
import os
import pickle
import union_find_cluster
import sys


'''constants'''

#set locaiton of delta files
delta_dir = '/srv/home/s4204666/abisko/dan/repeatm_tests/all_assemblies/filtered_2000bp/nucmer_feed_out'

#set location of assemblies
assembly_dir = '/srv/home/s4204666/abisko/dan/repeatm_tests/all_assemblies/filtered_2000bp'

#set location of bam files
bam_dir = ''


#set empty vector of all nucmer_matches that pass threshold (=0.90 by default)
collated_sig_matches = []
sigmatch_set = set()
#set empty vector to collect all nucmer_matches that are fragments of larger assemblies
fragments = []

'''-------------------------begin pipe-------------------------------------'''


#read in nucmer matches of each deltafile
for file in os.listdir(delta_dir):
    if file.endswith('.delta'):
        print('parsing {}'.format(file))
        delta = delta_parse.deltaread(delta_dir+'/'+file)
        #exit type is dictionary
        
        #threshold all matches and collate
        print('thresholding {}'.format(file))
        firstpass_fragments = []
        for m in next(iter(delta.values())):
            if m.seqs[0] ==m.seqs[1]:
                continue #it's a self match
            else:
                stats = m.gen_statistics()
                if m.apply_threshold(threshold = 0.90, stats = stats) == True:
                    collated_sig_matches.append(m)
                    sigmatch_set.update(set(m.seqs))
                elif m.is_fragment(upperthreshold = 0.90, lowerthreshold = 0.90, stats = stats):
                    firstpass_fragments.append(m)
    
        #all match objects are read in. keep only those fragments that map to full matches
        for m in firstpass_fragments:
            if m.seqs[0] in sigmatch_set:
                fragments.append(m)
            elif m.seqs[1] in sigmatch_set:
                fragments.append(m)
        sys.stdout.flush()
            #exit type is list of Nucmer_Match objects (sig_matches) + list of Nucmer_Match objects (fragments) that map to the sig_match objects

#save progress
f = open('pickled_sigmatches', 'wb')
pickle.dump(collated_sig_matches, f)
f.close()

f = open('pickled_fragmatches', 'wb')
pickle.dump(fragments, f)
f.close()

print('clustering...')
#cluster all significant matches & sort
clusters = single_linkage_cluster.cluster_nucmer_matches(collated_sig_matches)

#sort clusters
single_linkage_cluster.sort_clusters(clusters)

#further refine fragments to find those only that map to the trimmed clusters (ie. size >2)
cluster_nodes = []
for c in clusters:
    cluster_nodes +=c.nodes
cluster_node_set = set(cluster_nodes)

cluster_frags = []

for m in fragments:
    if m.seqs[0] in cluster_node_set:
        cluster_frags.append(m)
    elif m.seqs[1] in cluster_node_set:
        cluster_frags.append(m)

f = open('pickled_cluster_objs', 'wb')
pickle.dump(clusters, f)
f.close()
#exit type is list of Contig_Cluster objects

f = open('pickled_fragments', 'wb')
pickle.dump(cluster_frags, f)
f.close()
#exti type is list of match objects that map to clusters

##run analysis on clusters
#for c in clusters:
#    c.retrieve_seqs(assembly_dir = assembly_dir) #locations set at beginning of script. Does not return anything
#    
#f = open('pickled_clusters', 'wb')
#pickle.dump(clusters, f)
#f.close()


    
##determine orientation of matches in the cluster

'''-------------------------------------end pipe---------------------------'''