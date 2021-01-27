# -*- coding: utf-8 -*-

"""
Implementation of the Louvain method for community detection.
"""

import random, math
import networkx as nx

from asymmetree.datastructures.LinkedList import LinkedList
from asymmetree.datastructures.AVLTree import TreeSet


__author__ = 'David Schaller'


class Supernode:
    
    def __init__(self):
        
        self.nodes = []
    
    def extend(self, nodes):
        
        self.nodes.extend(nodes)
        
    # def remove(self, node):
        
    #     self.nodes.remove(node)


class Community:
    
    def __init__(self, nodes):
        
        self.nodes = {item for item in nodes}
        
        # sum of the weights of the links inside the community
        self._in = 0
        
        # sum of the weights of the links incident to nodes in the community
        self._tot = 0
    
    
    def __contains__(self, node):
        
        return node in self.nodes
    
    
    def add(self, node):
        
        self.nodes.add()
    
    
    def remove(self, node):
        
        self.nodes.remove(node)
        
        
class Level:
    
    def __init__(self, graph):
        
        self.graph = graph
        self.communities = set()
        self.node_to_com = {}
        self.nodes = [x for x in graph.nodes()]
        
        # sum of the weights of all edges incident to nodes
        self.k = {x: 0 for x in self.nodes}
        
        # sum of of all loop weights (times 2) for all nodes
        self.loops = {x: 0 for x in self.nodes}
        
        # sum of all edge weight
        self.m = 0
        
        for x in self.nodes:
            com = Community([nx])
            self.communities.add(com)
            self.node_to_com[x] = com
        
        for x, y, data in graph.edges(data=True):
            weight = data.get('weight', 1.0)
            self.k[x] += weight
            self.k[y] += weight
            self.node_to_com[x]._tot += weight
            self.node_to_com[y]._tot += weight
            if self.node_to_com[x] is self.node_to_com[y]:
                self.node_to_com[x]._in += 2 * weight
                self.loops[x]+= 2 * weight
            self.m += weight
        


class Louvain:
    
    def __init__(self, graph):
        
        if not isinstance(graph, nx.Graph) or graph.is_directed():
            raise TypeError("input graph must be an undirected NetworkX graph")
        
        # self.graph = graph
        
        levels = []
        
    
def compute_level(graph):
    
    communities = set()
    node_to_com = {}
    nodes = [x for x in graph.nodes()]
    
    # sum of the weights of all edges incident to nodes
    k = {}
    
    # sum of all edge weight
    m = 0
    
    for x in nodes:
        com = Community([nx])
        communities.add(com)
        node_to_com[x] = com
        k[x] = 0
    
    for x, y, data in graph.edges(data=True):
        weight = data.get('weight', 1.0)
        k[x] += weight
        k[y] += weight
        node_to_com[x].Sigma_tot += weight
        node_to_com[y].Sigma_tot += weight
        if node_to_com[x] is node_to_com[y]:
            node_to_com[x] += 2 * weight
    
    random.shuffle(nodes)
    
    while True:
        moved_node = False
        
        for x in nodes:
            
            neighbor_com_weights = _neighbor_communities(x, node_to_com, graph)
            
            for y in graph.neighbors(x):
                if node_to_com[x] is node_to_com[y]:
                    continue
                
                C_x = node_to_com[x]
                C_y = node_to_com[y]
        
        if not moved_node:
            break
    
    
def _move_into_C_gain(x, C, graph, m, k):
    
    # the sum of the weights of the links from i to nodes in C
    k_x_in = 0
    for y in graph.neighbors(x):
        if y in C:
            weight = graph.get_edge_data(x, y).get('weight', 1.0)
            k_x_in += weight
    
    delta_Q = ((C._in + 2 * k_x_in) / (2 * m) - \
               ((C._tot + k[x]) / (2 * m)) ** 2) - \
              (C._in / (2 * m) - \
               (C._tot / (2 * m)) ** 2 - \
               (k[x] / (2 * m)) ** 2 )
    
    return k_x_in, delta_Q
            

def _neighbor_communities(x, node_to_com, graph):
    
    # sum of the weight from x to every community C ( \{x} )
    neighbor_com_weights = {}
    
    for y in graph.neighbors(x):
        if x == y:
            continue
        C = node_to_com[y]
        weight = graph.get_edge_data(x, y).get('weight', 1.0)
        neighbor_com_weights[C] = neighbor_com_weights.get(C, 0.0) + weight
        

def initial_communities(self):
    
    communities = set()
    
    for x in self.graph.nodes:
        
        community = Community([x])
        
    pass
    
    
