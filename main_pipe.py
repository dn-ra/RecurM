# -*- coding: utf-8 -*-
"""
Main pipeline for Nucmer_Matching from .delta files to cluster generation

author: Daniel Rawlinson, ACE
email: daniel.rawlinson@uqconnect.edu.au
"""
import single_linkage_cluster
import delta_parse
import os

#set locaiton of delta files
delta_dir = 'C:\\Users\\dan_r\\Documents\\Honours_Data\\nucmertest\\Nucmer_results'

#set location of bam files
bam_dir = ''

alldeltas = []
#read in nucmer matches of each deltafile
for file in os.listdir(delta_dir):
    if file.endswith('.delta'):
        alldeltas.append(delta_parse.deltaread(delta_dir+'/'+file))
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
#exit type is list of Contig_Cluster objects

#run analysis on clusters
for c in clusters:
    assembly_node_dict = c.retrieve_seqs(assembly_dir = ??, return_node_assembly_dict = True) #bam_location set at beginning of script
    c.gen_minibam(assembly_node_dict, bam_location = bam_dir)
    c.label_cluster()
    
##determine orientation of matches in the cluster
