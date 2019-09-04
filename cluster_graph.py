# -*- coding: utf-8 -*-
"""
A module for use in RepeatM. see github.com/wwood/repeatm

Uses partial match Nucmer_Match object to form connections between clusters and exogenous sequences.
Uses this to determine fidelity of repeated sequences.

@author: Daniel Rawlinson, Australian Centre for Ecogenomics (ACE)
@email: daniel.rawlinson@uqconnect.edu.au
"""

#clusters are the cluster objects
#cluster_frags are the framents (one align_len > 0.9 and ani > 0.9) that map to the clusters
       


'''imports'''
import delta_parse

'''constants'''


'''------------------------begin class definition---------------------------'''
class Cluster_Graph(object):
    
    def __init__(self, cluster_objs):
        self.vertices = []
        self.edges = []
        self.c_n_d = {} # dictionary with key = node, value = cluster it belongs to
        self.pointer_dict = {} #adjacency link list array
        self.degrees = {}
        for c in cluster_objs:
            self.pointer_dict[c] = [] #initate pointers for clusters
            for node in c.nodes:
                self.c_n_d[node] = c #set up node-cluster dict
        
        
    def add_edge(self, m): #where m is a frag_match object
        '''process match object and add into graph'''
        if not isinstance(m, delta_parse.Nucmer_Match):
            raise Exception('add_edge function requires a Nucmer_Match object as input')
        long_link = m.seqs[m.lengths.index(max(m.lengths))]
        short_link = m.seqs[m.lengths.index(min(m.lengths))]
        if long_link in self.c_n_d: # if that node belongs to a cluster get that cluster as the end being pointed to
            long_point = self.c_n_d[long_link]
        else:
            if long_link not in self.pointer_dict:
                self.pointer_dict[long_link] = [] #initiate with no pointers
            long_point = long_link
            
        if short_link in self.c_n_d:
            short_point = self.c_n_d[short_link]
        else:
            if short_link not in self.pointer_dict:
                self.pointer_dict[short_link] = [] #initiate with no pointers
            short_point = short_link
        
    
        self.pointer_dict[short_point].append(long_point)
    
        
    def split(self):
        '''split graph into smaller connected ones'''
        return
    def identify_larger(self, cluster_obj):
        return
       
#these come from fragments, not from cluster_objs. instantiate as they are read
#        n_p_d[node] = []
    
#have now generated directed graph. Need to separate into discrete graph structures

#some of these clusters point to themselves. That indicates that there are some 
#match objects that fail the first pass threshold (a consequence of single 
#linkage clsutering)

'''DFS'''

queue = []
visited = {k:False for k in pointer_dict.keys()}
results = []

while False in visited.values(): #driver
    a = next(vertex for vertex in visited.keys() if visited[vertex] == False)
    visited[a] = True
    connected_graph = set(str(a))
    for edge in pointer_dict[a]:
        visited[edge] = True
        connected_graph.add(edge)
    
    
    
