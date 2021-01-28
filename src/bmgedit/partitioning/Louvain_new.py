# -*- coding: utf-8 -*-

"""
Implementation of the Louvain method for community detection.
"""

import random, math
import networkx as nx

from asymmetree.datastructures.LinkedList import LinkedList


__author__ = 'David Schaller'


class Louvain:
    
    def __init__(self, graph):
        
        if not isinstance(graph, nx.Graph) or graph.is_directed():
            raise TypeError("input graph must be an undirected NetworkX graph")
        
        # self.graph = graph
        
        levels = []


class Supernode:
    
    def __init__(self):
        
        self.nodes = []
    
    def extend(self, nodes):
        
        self.nodes.extend(nodes)
        
    # def remove(self, node):
        
    #     self.nodes.remove(node)


class Community:
    
    def __init__(self, initial_node):
        
        self.nodes = {initial_node}
        
        # sum of the weights of the links inside the community
        self._in = 0
        
        # sum of the weights of the links incident to nodes in the community
        self._tot = 0
        
        
    def __len__(self):
        
        return len(self.nodes)
    
    # def __contains__(self, node):
    #
    #     return node in self.nodes
    
    
    def insert(self, x, x_degree, x_loop, k_x_in):
        
        self.nodes.add(x)
        self._in += 2 * k_x_in + x_loop
        self._tot += x_degree
    
    
    def remove(self, x, x_degree, x_loop, k_x_in):
        
        self.nodes.remove(x)
        self._in -= 2 * k_x_in + x_loop
        self._tot -= x_degree
        
        
class Level:
    
    def __init__(self, graph):
        
        self.graph = graph
        self.communities = set()
        self.node_to_com = {}
        self.nodes = [x for x in graph.nodes()]
        
        self.moved_on_level = False
        
        # sum of the weights of all edges incident to nodes
        self.k = {x: 0 for x in self.nodes}
        
        # sum of of all loop weights (times 2) for all nodes
        self.loops = {x: 0 for x in self.nodes}
        
        # sum of all edge weight
        self.m = 0
        
        for x in self.nodes:
            com = Community(x)
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
            
    
    def _cluster(self):
        
        random.shuffle(self.nodes)
    
        while True:
            moved_node = False
            
            for x in self.nodes:
                
                C_x = self.node_to_com[x]
                
                # remove x from its community
                C_x.remove(x)
                
                # maps community C to the sum of weights of the links of x to
                # elements in C \ {x}
                k_x_in = self._weight_sums_to_communities(x)
                
                # C_x is preferred in case of ties, i.e. stays in its original
                # community
                best_gain = self._modularity_gain(x, C_x, k_x_in[C_x])
                best_C = C_x
                visited = {C_x}
                
                for y in self.graph.neighbors(x):
                    
                    C_y = self.node_to_com[y]
                    if C_y in visited:
                        continue
                    
                    new_gain = self._modularity_gain(x, C_y, k_x_in[C_y])
                    if new_gain > best_gain:
                        best_gain = new_gain
                        best_C = C_y
                        moved_node = True
                        self.moved_on_level = True
                
                best_C.insert(x)
                self.node_to_com[x] = best_C
                
                # remove the original community if it is now empty
                if len(C_x) == 0:
                    self.communities.remove(C_x)
            
            # exit the loop when all nodes stayed in their community
            if not moved_node:
                break
    
    
    def _weight_sums_to_communities(self, x):
    
        # sum of the weight from x to every community C ( \{x} )
        weight_sums = {}
        
        for y in self.graph.neighbors(x):
            if x == y:
                continue
            C = self.node_to_com[y]
            weight = self.graph.get_edge_data(x, y).get('weight', 1.0)
            weight_sums[C] = weight_sums.get(C, 0.0) + weight
        
        return weight_sums
    
    
    def _modularity_gain(self, x, C, k_x_in):
        
        delta_Q = ((C._in + 2 * k_x_in) / (2 * self.m) - \
                   ((C._tot + self.k[x]) / (2 * self.m)) ** 2) - \
                  (C._in / (2 * self.m) - \
                   (C._tot / (2 * self.m)) ** 2 - \
                   (self.k[x] / (2 * self.m)) ** 2 )
        
        return delta_Q
        
