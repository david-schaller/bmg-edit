# -*- coding: utf-8 -*-

"""
Implementation of the Louvain method for community detection.

References
----------
.. [1] Blondel V.D., Guillaume J.-L., Lambiotte R., Lefebvre E. (2008)
   Fast unfolding of communities in large networks. J. Stat. Mech.
   doi:10.1088/1742-5468/2008/10/P10008
"""

import random, itertools
import networkx as nx
from networkx.algorithms.community import modularity


__author__ = 'David Schaller'


class Louvain:
    """
    Implementation of the Louvain method for community detection.
    
    References
    ----------
    .. [1] Blondel V.D., Guillaume J.-L., Lambiotte R., Lefebvre E. (2008)
       Fast unfolding of communities in large networks. J. Stat. Mech.
       doi:10.1088/1742-5468/2008/10/P10008
    """
    
    def __init__(self, graph, weight='weight', print_info=False):
        """Constructor for the Louvain class.
        
        Runs the identification of communities.
        
        Parameters
        ----------
        graph : networkx.Graph
            The graph in which communities shall be found.
        weight : str, optional
            The name of the attribute used for edge weighting. The edges of
            `graph` must have this attribute (the default is 'weight', in which
            case edges without this attribute have weight 1.0).
        print_info : bool, optional
            Print the modularity results after each level (the default is
            False).
        """
        
        if not isinstance(graph, nx.Graph) or graph.is_directed():
            raise TypeError("input graph must be an undirected NetworkX graph")
        
        self.orig_graph = graph
        self.weight = weight
        self.print_info = print_info
        
        self._run()
        
    
    def _run(self):
        
        self.partitions = [ [{x} for x in self.orig_graph.nodes()] ]
        self.modularities = [ modularity(self.orig_graph,
                                         self.partitions[0],
                                         weight=self.weight) ]
        
        # construct the graph for the first level
        graph = self._next_level_graph(self.orig_graph, self.partitions[0])
        
        while True:
            
            level = _Level(graph, weight=self.weight)
            
            if not level.moved_on_level:
                break
            
            # partition of the supernodes (= found communities)
            sn_part = list(level.communities.values())
            mod = modularity(graph, sn_part, weight=self.weight)
            
            self.partitions.append(level.get_partition())
            self.modularities.append(mod)
            
            graph = self._next_level_graph(graph, sn_part)
            
            if self.print_info:
                print('----- Level {} -----'.format(len(self.partitions)-1))
                print('computed gain:',
                      self.modularities[-2], '+', level.total_modularity_gain,
                      '=', self.modularities[-2] + level.total_modularity_gain)
                print('mod. based on supernodes:', mod)
                print('mod. based on original graph',
                      modularity(self.orig_graph, self.partitions[-1],
                                 weight=self.weight))
    
    
    def _next_level_graph(self, graph, partition):
        
        nl_graph = nx.Graph()
        
        # maps the old supernode to the supernode in the next level
        old_to_new = {}
        
        for i, part_set in enumerate(partition):
            new_node = _Supernode()
            nl_graph.add_node(new_node)
            for old_node in part_set:
                if isinstance(old_node, _Supernode):
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
        
        return nl_graph


class _Supernode:
    """Wrapper class for communities on the previous levels.
    
    The communities found in each level are the nodes in the next level. This
    class contains all nodes of the original graph that belong to such a 
    supernode.
    """
    
    __slots__ = ('nodes',)
    
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
        
        
class _Level:
    """Clustering on a single level in the Louvain method."""
    
    def __init__(self, graph, weight='weight'):
        
        self.graph = graph
        self.weight = weight  
        
        self._initialize()
        self._cluster()
    
    
    def _initialize(self):
        
        self.nodes = [x for x in self.graph.nodes()]
        self.node_to_com = {x: i for i, x in enumerate(self.nodes)}
        self.communities = {i: {x} for x, i in self.node_to_com.items()}
        
        # sum of the weights of the links incident to nodes in the community
        self.com_tot = {i: 0 for i in self.communities.keys()}
        
        # sum of the weights of all edges incident to nodes
        self.k = {x: 0 for x in self.nodes}
        
        # sum of all edge weight
        self.m = self.graph.size(weight=self.weight)
        
        self.moved_on_level = False
        self.total_modularity_gain = 0.0
        
        for x, y, data in self.graph.edges(data=True):
            weight = data.get(self.weight, 1.0)
            self.k[x] += weight
            self.k[y] += weight
            self.com_tot[self.node_to_com[x]] += weight
            self.com_tot[self.node_to_com[y]] += weight
            
    
    def _cluster(self):
        
        # for an edgeless graph, every node is in its own cluster
        if self.m == 0:
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
                self.communities[C_x].remove(x)
                self.com_tot[C_x] -= self.k[x]
                
                # C_x is preferred in case of ties, i.e. stays in its original
                # community
                
                # equation in Blondel et al. can be simplified to this
                cost_removal = (k_x_in[C_x] - 
                                self.com_tot[C_x] * self.k[x] / (2 * self.m)) \
                               / self.m
                best_gain = cost_removal
                best_C = C_x
                visited = {C_x}
                
                for y in self.graph.neighbors(x):
                    
                    C_y = self.node_to_com[y]
                    if C_y in visited:
                        continue
                    
                    new_gain = (k_x_in[C_y] - 
                                self.com_tot[C_y] * self.k[x] / (2 * self.m)) \
                               / self.m
                    
                    if new_gain > best_gain:
                        best_gain = new_gain
                        best_C = C_y
                        moved_node = True
                        self.moved_on_level = True
                
                self.communities[best_C].add(x)
                self.com_tot[best_C] += self.k[x]
                self.node_to_com[x] = best_C
                self.total_modularity_gain += best_gain - cost_removal
                
                # remove the original community if it is now empty
                if len(self.communities[C_x]) == 0:
                    del self.communities[C_x]
            
            # exit the loop when all nodes stayed in their community
            if not moved_node:
                break
            
    
    def get_partition(self):
        """Convert the communities of supernodes into a partition of the 
        orginal nodes."""
        
        return [set(itertools.chain.from_iterable(com))
                for com in self.communities.values()]
    
    
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


if __name__ == '__main__':
    
    # random graph
    G = nx.erdos_renyi_graph(100, 0.02)
    louv = Louvain(G, print_info=False)
    for mod, part in zip(louv.modularities, louv.partitions):
        print(mod)
    
    # example from Blondel et al. 2008
    G2 = nx.Graph()
    G2.add_nodes_from([x for x in range(16)])
    G2.add_edges_from([(0,2), (0,3), (0,4), (0,5),
                       (1,2), (1,4), (1,7), (2,4), (2,5), (2,6),
                       (3,7), (4, 10), (5,7), (5,11), (6,7), (6,11),
                       (8,9), (8,10), (8,11), (8,14), (8,15), (9,12), (9,14),
                       (10,11), (10,12), (10,13), (10,14), (11,13),
                       ])
    
    louv = Louvain(G2, print_info=True)
    for mod, part in zip(louv.modularities, louv.partitions):
        print(part, mod)
