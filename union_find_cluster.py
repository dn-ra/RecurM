# -*- coding: utf-8 -*-
"""
Created on Fri Aug  9 15:43:35 2019

Alternative method for clustering match objects that come out of nucmer for RepeatM
Based on Disjoint Set / Union-Find algorithm.
See:
    https://www.geeksforgeeks.org/union-find/
    https://en.wikipedia.org/wiki/Disjoint-set_data_structure
    https://www.youtube.com/watch?v=wU6udHRIkcc
    
TODO:
    - make collapsed tree instead? will make extra_cluster step easier

@author: Daniel Rawlinson, Australian Centre for Ecogenomics (ACE)
@email: daniel.rawlinson@uqconnect.edu.au
"""

def find_parent(node, node_parent_array):
    #print('finding parent of', node)
    try:
        parent = node_parent_array[node]
        if isinstance(parent,int):
            return node
        else:
            return find_parent(parent, node_parent_array)
    except KeyError:
        #print('adding', node)
        node_parent_array[node] = -1
        return node


def union(node_1, node_2, node_parent_array):
    set_1 = find_parent(node_1, node_parent_array)
    set_2 = find_parent(node_2, node_parent_array)
    if set_1 == set_2:
        pass
    else:
        node_parent_array[set_2] = set_1
    

def cluster_conversion(node_parent_array):
    '''convert node_parent array into dictionary of node names:cluster_ids'''
    i = 0 #initial cluster counter
    for key, value in node_parent_array.items():
        if not isinstance(value, int): #if value is re-direct to another node
            #print('calling find_parent for', value)
            c_num = node_parent_array[find_parent(value, node_parent_array)] #find cluster number further up the graph
            #print('found', c_num, 'from', value, 'at cluster conversion')
            if c_num == -1:
                c_num = i
                i+=1
            cluster_chain(key, c_num, node_parent_array)
        elif value == -1:
            node_parent_array[key] = i #use this cluster number for this grouping
            i+=1 #iterate up the cluster counter
        else: #last option should be if the value is an already-decided cluster number. Can't trace back from this! so return to dict iterator
            continue

def extract_clusters(node_parent_array): #use max_c as the known number of clusters?
    '''pull out node names for each cluster id. Return list of id lists'''
    reverse = {} #dictionary of reversed name: cluster dict
    for name, number in node_parent_array.items():
        if number not in reverse:
            reverse[number] = [] #add that cluster number into dict if not there yet
        reverse[number].append(name) #add id into that cluster dict
    return list(reverse.values()) #output as list of lists I GET IT!
            
def cluster_chain(node, c_num, node_parent_array):
    '''follow path through and change values to cluster number'''
    while not isinstance(node, int): #when I know the algorithm is working, change this to: while node != c_num
        new_n = node_parent_array[node]
        #print('changing', node, 'to', c_num)
        node_parent_array[node] = c_num
        node = new_n
    
    
    
def find_rep():
    '''find the most represented node to pass on for analysis'''
    return None


def union_find_pipe(collated_sig_matches):
    disjoint_set_array = {}
    links = [m.seqs for m in collated_sig_matches]
    for link in links:
        union(link[0], link[1], disjoint_set_array)
    cluster_conversion(disjoint_set_array)
    clusters = extract_clusters(disjoint_set_array)
    
    return clusters