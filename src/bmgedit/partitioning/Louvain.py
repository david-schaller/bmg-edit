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


class Louvain:
    """
    Implementation of the Louvain method for community detection.
    
    References
    ----------
    .. [1] Blondel V.D., Guillaume J.-L., Lambiotte R., Lefebvre E. (2008)
       Fast unfolding of communities in large networks. J. Stat. Mech.
       doi:10.1088/1742-5468/2008/10/P10008
    """
    
    def __init__(self, graph, weight='weight',
                 cost_function=None,
                 cost_function_args=None,
                 print_info=False):
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
                      self.modularities[-2], '+', level.total_gain,
                      '=', self.modularities[-2] + level.total_gain)
                print('mod. based on supernodes:', mod)
                print('mod. based on original graph',
                      modularity(self.orig_graph, self.partitions[-1],
                                 weight=self.weight))
    
    
    def _next_level_graph(self, graph, partition):
        
        nl_graph = nx.Graph()
        
        # maps the old supernode to the supernode in the next level
        old_to_new = {}
        
        for part_set in partition:
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


class LouvainCustomCost:
    """
    Louvain method for community detection with a customized cost function.
    
    References
    ----------
    .. [1] Blondel V.D., Guillaume J.-L., Lambiotte R., Lefebvre E. (2008)
       Fast unfolding of communities in large networks. J. Stat. Mech.
       doi:10.1088/1742-5468/2008/10/P10008
    """
    
    def __init__(self, graph,
                 cost_function=None,
                 cost_function_args=None,
                 print_info=False):
        """Constructor for the Louvain class.
        
        Runs the identification of communities.
        
        Parameters
        ----------
        graph : networkx.Graph
            The graph in which communities shall be found.
        cost_function : function object, optional
            Custom function that is used instead of modularity optimization
            (the default is None).
        cost_function_args : tuple, optional
            Additional argument relevant for a custom cost function.
        print_info : bool, optional
            Print the modularity results after each level (the default is
            False).
        """
        
        if not isinstance(graph, nx.Graph) or graph.is_directed():
            raise TypeError("input graph must be an undirected NetworkX graph")
        
        self.orig_graph = graph
        self.cost_function = cost_function
        self.cost_function_args = cost_function_args
        self.print_info = print_info
        
        self._run()
        
    
    def _run(self):
        
        self.partitions = [ [{x} for x in self.orig_graph.nodes()] ]
        self.costs = [ self.cost_function(self.partitions[0],
                                          *self.cost_function_args) ]
        
        # construct the graph for the first level
        graph = self._next_level_graph(self.orig_graph, self.partitions[0])
        
        while True:
            
            level = _Level(graph,
                           cost_function=self.cost_function,
                           cost_function_args=self.cost_function_args)
            
            if not level.moved_on_level:
                break
            
            # partition of the supernodes (= found communities)
            part = level.get_partition()
            
            self.partitions.append(part)
            self.costs.append(self.costs[-1] - level.total_gain)
            
            graph = self._next_level_graph(graph, part)
            
            if self.print_info:
                print('----- Level {} -----'.format(len(self.partitions)-1))
                print('computed gain:',
                      self.costs[-2], '-', level.total_gain,
                      '=', self.costs[-2] + level.total_gain)
                print('mod. based on original graph',
                      self.cost_function(self.partitions[-1],
                                         *self.cost_function_args))
    
    
    def _next_level_graph(self, graph, partition):
        
        nl_graph = nx.Graph()
        
        # maps the old supernode to the supernode in the next level
        old_to_new = {}
        
        for part_set in partition:
            new_node = _Supernode()
            nl_graph.add_node(new_node)
            for old_node in part_set:
                if isinstance(old_node, _Supernode):
                    new_node.extend(old_node)
                else:
                    new_node.append(old_node)
                old_to_new[old_node] = new_node
        
        for x, y, data in graph.edges(data=True):
            u = old_to_new[x]
            v = old_to_new[y]
            if not nl_graph.has_edge(u, v):
                nl_graph.add_edge(u, v)
        
        return nl_graph
        
        
class _Level:
    """Clustering on a single level in the Louvain method."""
    
    def __init__(self, graph, weight='weight',
                 cost_function=None, cost_function_args=None):
        
        self.graph = graph
        self.weight = weight
        self.cost_function = cost_function
        self.cost_function_args = cost_function_args
        
        self.nodes = [x for x in self.graph.nodes()]
        random.shuffle(self.nodes)
        
        self.node_to_com = {x: i for i, x in enumerate(self.nodes)}
        self.communities = {i: {x} for x, i in self.node_to_com.items()}
        
        self.moved_on_level = False
        self.total_gain = 0.0
        
        if self.cost_function is None:
            self._cluster_by_modularity()
        else:
            self._cluster_by_cost()
    
    
    def _cluster_by_modularity(self):
        
        # sum of all edge weight
        m = self.graph.size(weight=self.weight)
        
        # for an edgeless graph, every node is in its own cluster
        if m == 0:
            return
        
        # sum of the weights of the links incident to nodes in the community
        com_tot = {i: 0 for i in self.communities.keys()}
        
        # sum of the weights of all edges incident to nodes
        k = {x: 0 for x in self.nodes}
        
        for x, y, data in self.graph.edges(data=True):
            w = data.get(self.weight, 1.0)
            k[x] += w
            k[y] += w
            com_tot[self.node_to_com[x]] += w
            com_tot[self.node_to_com[y]] += w
    
        while True:
            moved_node = False
            
            for x in self.nodes:
                
                C_x = self.node_to_com[x]
                
                # maps community C to the sum of weights of the links of x to
                # elements in C \ {x}
                k_x_in = self._weight_sums_to_communities(x)
                
                # remove x from its community
                self.communities[C_x].remove(x)
                com_tot[C_x] -= k[x]
                
                # C_x is preferred in case of ties, i.e. stays in its original
                # community
                
                # equation in Blondel et al. can be simplified to this
                cost_removal = (k_x_in[C_x] - com_tot[C_x] * k[x] / (2 * m)) / m
                best_gain = cost_removal
                best_C = C_x
                visited = {C_x}
                
                for y in self.graph.neighbors(x):
                    
                    C_y = self.node_to_com[y]
                    if C_y in visited:
                        continue
                    
                    new_gain = (k_x_in[C_y] - com_tot[C_y] * k[x] / (2 * m)) / m
                    
                    if new_gain > best_gain:
                        best_gain = new_gain
                        best_C = C_y
                        moved_node = True
                        self.moved_on_level = True
                
                self.communities[best_C].add(x)
                com_tot[best_C] += k[x]
                self.node_to_com[x] = best_C
                self.total_gain += best_gain - cost_removal
                
                # remove the original community if it is now empty
                if len(self.communities[C_x]) == 0:
                    del self.communities[C_x]
            
            # exit the loop when all nodes stayed in their community
            if not moved_node:
                break
    
    
    def _cluster_by_cost(self):
        
        # for an edgeless graph, every node is in its own cluster
        if self.graph.size() == 0:
            return
        
        # avoid recomputation of the partition in each iteration
        part = {i: set(x.nodes) for x, i in self.node_to_com.items()}
        
        cost = self.cost_function(self.get_partition(),
                                  *self.cost_function_args)
    
        while True:
            moved_node = False
            
            for x in self.nodes:
                
                C_x = self.node_to_com[x]
                
                # remove x from its community
                self.communities[C_x].remove(x)
                part[C_x].difference_update(x.nodes)
                
                # C_x is preferred in case of ties, i.e. stays in its original
                # community
                best_gain = 0.0
                best_C = C_x
                visited = {C_x}
                
                for y in self.graph.neighbors(x):
                    
                    C_y = self.node_to_com[y]
                    if C_y in visited:
                        continue
                    
                    part[C_y].update(x.nodes)
                    new_cost = self.cost_function(self.get_partition(),
                                                  *self.cost_function_args)
                    part[C_y].difference_update(x.nodes)
                    
                    new_gain = cost - new_cost
                    
                    if new_gain > best_gain:
                        best_gain = new_gain
                        best_C = C_y
                        moved_node = True
                        self.moved_on_level = True
                
                self.communities[best_C].add(x)
                part[best_C].update(x.nodes)
                self.node_to_com[x] = best_C
                self.total_gain += best_gain
                cost -= best_gain
                
                # remove the original community if it is now empty
                if len(self.communities[C_x]) == 0:
                    del self.communities[C_x]
                    del part[C_x]
            
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
            w = self.graph.get_edge_data(x, y).get('weight', 1.0)
            weight_sums[C] = weight_sums.get(C, 0.0) + w
        
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
