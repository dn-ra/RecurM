# -*- coding: utf-8 -*-
"""
Main pipeline for Nucmer_matching from .delta files to cluster generation

author: Daniel Rawlinson, ACE
email: daniel.rawlinson@uqconnect.edu.au
"""
import single_linkage_cluster
import delta_parse


#set locaiton of delta files
delta_dir = 'C:\Users\dan_r\Documents\Honours_Data\nucmertest\Nucmer_results'

#read in nucmer matches of each deltafile
for file in os.listdir(delta_dir):
    if file.endswith('.delta'):
        alldeltas.append(delta_parse.deltaread(file))
#exit type is dictionary

        
#set empty vector of all nucmer_matches that pass threshold (=0.90 by default)  
collated_sig_matches = []

#threshold all matches and collate
for delta_dict in alldeltas:
    collated_sig_matches += delta_parse.dict_threshold(delta_dict, threshold = 0.90, collate = True)
#exit type is list of Nucmer_Match objects

#cluster all significant matches & sort
clusters = single_linkage_cluster.cluster_nucmer_matches(collated_sig_matches)
single_linkage_cluster.sort_clusters(clusters)
#exit type if list of Contig_Cluster objects

#run analysis on clusters
for c in clusters:
    c.retrieve_seqs(get_bam_for_nodes = True)
