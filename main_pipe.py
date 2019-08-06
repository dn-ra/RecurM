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

#set locaiton of delta files
delta_dir = '/srv/home/s4204666/abisko/dan/repeatm_tests/all_assemblies/filtered_2000bp/nucmer_feed_out'

#set location of assemblies
assembly_dir = '/srv/home/s4204666/abisko/dan/repeatm_tests/all_assemblies/filtered_2000bp'

#set location of bam files
bam_dir = ''



#set empty vector of all nucmer_matches that pass threshold (=0.90 by default)
collated_sig_matches = []


#read in nucmer matches of each deltafile
for file in os.listdir(delta_dir):
    #if file.endswith('.delta'):
    print('parsing {}'.format(file))
    delta = delta_parse.deltaread(delta_dir+'/'+file)
    #exit type is dictionary
    
    #threshold all matches and collate
    print('thresholding {}'.format(file))
    collated_sig_matches += delta_parse.dict_threshold(delta, threshold = 0.90, collate = True)
    #exit type is list of Nucmer_Match objects

#save progress
f = open('pickled_sigmatches', 'wb')
pickle.dump(collated_sig_matches, f)
f.close()

#cluster all significant matches & sort
clusters = single_linkage_cluster.cluster_nucmer_matches(collated_sig_matches)
single_linkage_cluster.sort_clusters(clusters)
#exit type is list of Contig_Cluster objects

#run analysis on clusters
for c in clusters:
    assembly_node_dict = c.split_names()
    c.retrieve_seqs(assembly_node_dict, assembly_dir = assembly_dir) #locations set at beginning of script
    
f = open('pickled_clusters', 'wb')
pickle.dump(clusters, f)
f.close()


    
##determine orientation of matches in the cluster
