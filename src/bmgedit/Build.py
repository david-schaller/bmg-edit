# -*- coding: utf-8 -*-

import networkx as nx

from asymmetree.datastructures.PhyloTree import PhyloTree, PhyloTreeNode
from asymmetree.datastructures.Partition import Partition

from bmgedit.partitioning.Karger import Karger
from bmgedit.partitioning.GreedyBipartition import (greedy_bipartition,
                                                    gradient_walk_bipartition)


__author__ = 'David Schaller'


def aho_graph(R, L, weighted=False, triple_weights=None):
    """Construct the auxiliary graph (Aho graph) for BUILD.
        
    Edges {a,b} are optionally weighted by the number of occurrences, resp.,
    sum of weights of triples of the form ab|x or ba|x.
    """
    
    G = nx.Graph()
    G.add_nodes_from(L)

    for a, b, c in R:
        if not G.has_edge(a, b):
            G.add_edge(a, b, weight=0)
        
        if weighted:
            if triple_weights:
                G[a][b]['weight'] += triple_weights[a, b, c]
            else:
                G[a][b]['weight'] += 1
    
    return G


def _triple_connect(G, t):
    
    G.add_edge(t[0], t[1])
    G.add_edge(t[0], t[2])
    

def mtt_partition(L, R, F):
    
    # auxiliary graph initialized as Aho graph
    G = aho_graph(R, L, weighted=False)
    
    # auxiliary partition
    P = Partition(nx.connected_components(G))
    
    if len(P) == 1:
        return P, G
    
    # aux. set of forbidden triples
    S = {t for t in F if P.separated_xy_z(*t)}
    
    # lookup of forb. triples to which u belongs
    L = {u: [] for u in L}
    for t in F:
        for u in t:
            L[u].append(t)
    
    while S:
        t = S.pop()
        _triple_connect(G, t)
        
        # merge returns the smaller of the two merged sets
        smaller_set = P.merge(t[0], t[2])
        
        # update S by traversing the L(u)
        for u in smaller_set:
            for t in L[u]:
                if t in S and not P.separated_xy_z(*t):
                    S.remove(t)
                    _triple_connect(G, t)
                elif t not in S and P.separated_xy_z(*t):
                    S.add(t)
    
    return P, G



        

class Build2:
    """BUILD / MTT algorithm with minimal cost bipartition."""
    
    def __init__(self, R, L, F=None,
                 allow_inconsistency=True,
                 bipart_method='mincut',
                 cost_function=None, cost_function_args=None,
                 weighted_mincut=False, triple_weights=None):
        
        self.R = R
        self.L = L
        
        # forbidden triples --> activates MTT if non-empty
        self.F = F
        
        # allow inconsistencies or return False?
        self.allow_inconsistency = allow_inconsistency
        
        if bipart_method in ('mincut', 'karger', 'greedy', 'gradient_walk'):
            self.bipart_method = bipart_method
        else:
            raise ValueError("unknown bipartition method "\
                             "'{}'".format(bipart_method))
        
        self.cost_function = cost_function
        self.cost_function_args = cost_function_args
            
        # parameters if bipartition method is mincut
        self.weighted_mincut = weighted_mincut
        self.triple_weights = triple_weights
    
    
    def build_tree(self, return_root=False):
        """Build a tree displaying all triples in R if possible.
        
        Keyword arguments:
            return_root - if True, return 'PhyloTreeNode' instead of
                'PhyloTree' instance
        """
        
        self.total_cost = 0
        
        if self.F:
            root = self._mtt(self.L, self.R, self.F)
        else:
            root = self._aho(self.L, self.R)
            
        return root if return_root else PhyloTree(root)
    
    
    def _trivial_case(self, L):
        
        if len(L) == 1:
            leaf = L.pop()
            return PhyloTreeNode(leaf, label=str(leaf))
        
        elif len(L) == 2:
            node = PhyloTreeNode(-1)
            for _ in range(2):
                leaf = L.pop()
                node.add_child(PhyloTreeNode(leaf, label=str(leaf)))
            return node
    
    
    def _aho(self, L, R):
        """Recursive Aho-algorithm."""
        
        # trivial case: one or two leaves left in L
        if len(L) <= 2:
            return self._trivial_case(L)
            
        aux_graph = aho_graph(R, L, weighted=self.weighted_mincut,
                                    triple_weights=self.triple_weights)
        partition = list(nx.connected_components(aux_graph))
        
        if len(partition) < 2:
            if not self.allow_inconsistency:
                return False
            else:
                partition = self._bipartition(L, aux_graph)
        
        node = PhyloTreeNode(-1)            # place new inner node
        for s in partition:
            Li, Ri = set(s), []
            for t in R:                     # construct triple subset
                if Li.issuperset(t):
                    Ri.append(t)
            Ti = self._aho(Li, Ri)          # recursive call
            if not Ti:
                return False                # raise False to previous call
            else:
                node.add_child(Ti)          # add roots of the subtrees
   
        return node
    
    
    def _mtt(self, L, R, F):
        """Recursive MTT algorithm."""
        
        # trivial case: one or two leaves left in L
        if len(L) <= 2:
            return self._trivial_case(L)
        
        partition, aux_graph = mtt_partition(L, R, F)
        
        if len(partition) < 2:
            if not self.allow_inconsistency:
                return False
            else:
                partition = self._bipartition(L, aux_graph)
        
        node = PhyloTreeNode(-1)            # place new inner node
        for s in partition:
            Li, Ri, Fi = set(s), [], []
            for Xi, X in ((Ri, R), (Fi, F)):
                for t in X:
                    if Li.issuperset(t):
                        Xi.append(t)
            Ti = self._mtt(Li, Ri, Fi)      # recursive call
            if not Ti:
                return False                # raise False to previous call
            else:
                node.add_child(Ti)          # add roots of the subtrees
   
        return node
    
    
    def _bipartition(self, L, aux_graph):
        
        best_cost, best_bp = float('inf'), None
        
        if self.bipart_method == 'mincut':
            # Stoer–Wagner algorithm
            best_cost, best_bp = nx.stoer_wagner(aux_graph)
        
        elif self.bipart_method == 'karger':
            karger = Karger(aux_graph)
            
            for _, bp in karger.generate():
                cost = self.cost_function(bp, *self.cost_function_args)
                
                if cost < best_cost:
                    best_cost, best_bp = cost, bp
        
        elif self.bipart_method == 'greedy':
            
            for _ in range(5):
                cost, bp = greedy_bipartition(L, self.cost_function,
                                              args=self.cost_function_args)
                if cost < best_cost:
                    best_cost, best_bp = cost, bp
        
        elif self.bipart_method == 'gradient_walk':
            
            for _ in range(5):
                cost, bp = gradient_walk_bipartition(L, self.cost_function,
                               args=self.cost_function_args)
                if cost < best_cost:
                    best_cost, best_bp = cost, bp
        
        self.total_cost += best_cost
        return best_bp