# -*- coding: utf-8 -*-

"""
Implementation of the Louvain method for community detection.
"""

import random, itertools
import networkx as nx

# from asymmetree.datastructures.LinkedList import LinkedList


__author__ = 'David Schaller'


class Louvain:
    
    def __init__(self, graph, weight='weight'):
        
        if not isinstance(graph, nx.Graph) or graph.is_directed():
            raise TypeError("input graph must be an undirected NetworkX graph")
        
        self.orig_graph = graph
        self.weight = weight
        
        self.level_graphs = [self._level_zero_graph()]
               
        self.level_mods = [modularity(self.level_graphs[0],
                                      [[x] for x in self.level_graphs[0].nodes()],
                                      weight=self.weight)
                      ]
        print(self.level_mods[0])
        
    
    def run(self):
        
        while True:
            
            graph = self.level_graphs[-1]
            level = Level(graph, weight=self.weight)
            level.cluster()
            
            part = level.communities
            print([[x.nodes for x in com] for com in part])
            mod = modularity(graph, part, weight=self.weight)
            print('+', level.total_modularity_gain,
                  '=', self.level_mods[-1] + level.total_modularity_gain,
                  '=', mod)
            
            nl_graph = self._next_level_graph(graph, part)
            self.level_graphs.append(nl_graph)
            self.level_mods.append(mod)
            print(modularity(nl_graph, [[x] for x in nl_graph.nodes()],
                             weight=self.weight))
            for u,v, data in nl_graph.edges(data=True):
                print(u.nodes, v.nodes, data)
            
            
            if not level.moved_on_level:
                break
        
    
    def _level_zero_graph(self):
        
        partition = [[x] for x in self.orig_graph.nodes()]
        
        return self._next_level_graph(self.orig_graph, partition)
    
    
    def _next_level_graph(self, graph, partition):
        
        nl_graph = nx.Graph()
        
        # maps the old supernode to the supernode in the next level
        old_to_new = {}
        
        for i, part_set in enumerate(partition):
            new_node = Supernode()
            nl_graph.add_node(new_node)
            for old_node in part_set:
                if isinstance(old_node, Supernode):
                    new_node.extend(old_node)
                else:
                    new_node.append(old_node)
                old_to_new[old_node] = new_node
        
        for x, y, data in graph.edges(data=True):
            w = data.get(self.weight, 1.0)
            u = old_to_new[x]
            v = old_to_new[y]
            if nl_graph.has_edge(u, v):
                nl_graph[u][v][self.weight] += w
            else:
                nl_graph.add_edge(u, v)
                nl_graph[u][v][self.weight] = w
            
            # (new) loops count twice 
            if (u == v) and x != y:
                nl_graph[u][v][self.weight] += w
        
        return nl_graph
        


class Supernode:
    
    def __init__(self):
        
        self.nodes = []
    
    
    def __iter__(self):
        
        return iter(self.nodes)
    
    
    def __next__(self):
        
        pass
    
    
    def append(self, node):
        
        self.nodes.append(node)
    
    
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
    
    
    def __iter__(self):
        
        return iter(self.nodes)
    
    
    def __next__(self):
        
        pass
    
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
    
    def __init__(self, graph, weight='weight'):
        
        self.graph = graph
        self.weight = weight        
        
        self.communities = set()
        self.node_to_com = {}
        self.nodes = [x for x in graph.nodes()]
        
        self.moved_on_level = False
        self.total_modularity_gain = 0.0
        
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
            weight = data.get(self.weight, 1.0)
            self.k[x] += weight
            self.k[y] += weight
            self.node_to_com[x]._tot += weight
            self.node_to_com[y]._tot += weight
            if self.node_to_com[x] is self.node_to_com[y]:
                self.node_to_com[x]._in += 2 * weight
                self.loops[x] += 2 * weight
            self.m += weight
            
    
    def cluster(self):
        
        # for an edgeless graph, every node is in its own cluster
        if self.graph.size() == 0:
            return
        
        random.shuffle(self.nodes)
    
        while True:
            moved_node = False
            
            for x in self.nodes:
                
                C_x = self.node_to_com[x]
                
                # maps community C to the sum of weights of the links of x to
                # elements in C \ {x}
                k_x_in = self._weight_sums_to_communities(x)
                
                # remove x from its community
                C_x.remove(x, self.k[x], self.loops[x], k_x_in[C_x])
                
                # C_x is preferred in case of ties, i.e. stays in its original
                # community
                cost_removal = self._modularity_gain(x, C_x, k_x_in[C_x])
                best_gain = cost_removal
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
                
                best_C.insert(x, self.k[x], self.loops[x], k_x_in[best_C])
                self.node_to_com[x] = best_C
                self.total_modularity_gain += best_gain - cost_removal
                
                # remove the original community if it is now empty
                if len(C_x) == 0:
                    self.communities.remove(C_x)
            
            # exit the loop when all nodes stayed in their community
            if not moved_node:
                break
    
    
    def _weight_sums_to_communities(self, x):
    
        # sum of the weight from x to every community C ( \{x} )
        # the community of x must be present even if x is its only element
        weight_sums = {self.node_to_com[x]: 0.0}
        
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


def modularity(G, partition, weight='weight'):
    
    k = {x: 0.0 for x in G.nodes()}
    m = G.size(weight=weight)
    
    for x, y, data in G.edges(data=True):
        w = data.get(weight, 1.0)
        k[x] += w
        k[y] += w
        # if x == y:
        #     k[x] -= w
    
    nodes = [x for x in G.nodes()]
    
    if isinstance(partition, dict):
        node_to_part = partition
    else:
        node_to_part = {}
        for i, part_set in enumerate(partition):
            for x in part_set:
                node_to_part[x] = i
    
    total_sum = 0.0
    for x, y in itertools.product(nodes, nodes):
        if node_to_part[x] == node_to_part[y]:
            if G.has_edge(x, y):
                A_xy = G.get_edge_data(x, y).get(weight, 1.0)
            else:
                A_xy = 0.0
            total_sum += A_xy - k[x] * k[y] / (2 * m)
    
    return total_sum / (2 * m)


if __name__ == '__main__':
    
    # G = nx.erdos_renyi_graph(1000, 0.02)
    
    
    G = nx.Graph()
    G.add_nodes_from([x for x in range(1,11)])
    G.add_edges_from([(1,2), (2,3), (1,3),
                      (4,5), (5,6), (4,6),
                      (7,8), (8,9), (7,9),
                      (3,10), (6,10), (9,10),])
    
    part = [[x] for x in G.nodes()]
    # part = [[1,2,3], [4,5,6,10], [7,8,9]]
    
    
    print(modularity(G, part))
    
    louv = Louvain(G)
    louv.run()
