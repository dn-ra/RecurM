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
import nucmer_cluster


'''constants'''



'''constituent classes'''
#TODO - keep pointers within classes here?
#class cluster_vertex(object):
#    
#    def __init__(self, cluster_obj):
#        self.in_d = 0
#        self.out_d = 0
#        self.c = cluster_obj #python default behaviour is to copy memory pointer, not deep copy

#class single_vertex(object):
#    
#    def __init__(self, node_name):
#        self.name = node_name

'''------------------------begin class definition---------------------------'''
class Cluster_Graph(object):
    
    def __init__(self, cluster_objs): #initialise with all cluster objects first
        self.c_n_d = {} # dictionary with key = node (sequence), value = cluster it belongs to
        self.pointers = {} #adjacency link list array
        self.edges = {} #dictionary holding node objects, and a dictionary with edges in, out, and to self.
        
        #main initialisation of clusters as nodes.
        for c in cluster_objs:
            self.pointers[c] = [] #initate pointers for clusters
            self.edges[c] = {'in': 0, 'out': 0, 'self': 0}
            for node in c.nodes:
                self.c_n_d[node] = c #set up node-cluster dict. instantiation of vertex_objects happens 
                
    def add_vertex(self, node_name): #should only be adding single sequence vertices after graph has been initialised.
        if not isinstance(node_name, str):
            raise Exception('This is not a string. Clusters must be instantiated as vertices at the point of Graph creation')
        self.pointers[node_name] = []
        self.edges[node_name] = {'in':0,'out':0, 'self':None} #don't want self-matches. That's only for cluster objects

    def add_edge(self, m): #where m is a frag_match object
        '''process match object and add into graph'''
        long_link = m.seqs[m.lengths.index(max(m.lengths))]
        short_link = m.seqs[m.lengths.index(min(m.lengths))]
        if long_link in self.c_n_d: # if that node belongs to a cluster get that cluster as the end being pointed to
            long_point = self.c_n_d[long_link] # this will be a cluster object
        
        else: #if not match to a cluster
            if long_link not in self.pointers: #if it hasn't already been added to the graph 
                self.add_vertex(long_link)
            long_point = long_link #set pointer as the single sequence
            
        if short_link in self.c_n_d:
            short_point = self.c_n_d[short_link]
        else:
            if short_link not in self.pointers:
                self.add_vertex(short_link)
            short_point = short_link
            
        self.pointers[short_point].append(long_point)
        
        if short_point == long_point:
            self.edges[short_point]['self'] +=1
        else:
            self.edges[short_point]['out'] +=1
            self.edges[long_point]['in'] +=1

        return
    
        
    def quantify_subraphs(self):
        '''split graph into smaller connected ones'''
        #TODO - don't know how to achieve this. Do I need to add birectionality into the graph?
        subgraph_sets = []
        
        return
    
       
#these come from fragments, not from cluster_objs. instantiate as they are read
#        n_p_d[node] = []
    
#have now generated directed graph. Need to separate into discrete graph structures

#some of these clusters point to themselves. That indicates that there are some 
#match objects that fail the first pass threshold (a consequence of single 
#linkage clsutering)

    '''BFS'''

    def BFS(self, cluster_obj): #find all vertices that a partciular vertex points to (ie. is a fragment of)
        queue = [cluster_obj]
        visited = set()
        connected_graph = []
        
        
        while queue: #driver
            v = queue.pop(0)
            if v not in visited:
                points_to = [p for p in self.pointers[v] if p not in visited]
                queue += points_to
                if v != cluster_obj: #don't add query cluster into output
                    connected_graph.append(v)
                
                visited.add(v)
        return connected_graph


#'''---------try it with networkx data structure----------------------------'''
#Cluster_G = nx.MultiDiGraph(c_n_dict = None)
#
#for c in clusters:
#    Cluster_G.add_node(c, 'typ' = 'cluster','size' = c.size) #new node for each cluster
#    for node in c.nodes:
#        c_n_d[node] = c #set up node-cluster dict.
#
#Cluster_G.c_n_dict = c_n_d
#
#for m in frag_matches:
#    pass
    
