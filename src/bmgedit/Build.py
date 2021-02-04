# -*- coding: utf-8 -*-

import networkx as nx

from asymmetree.datastructures.PhyloTree import PhyloTree, PhyloTreeNode
from asymmetree.tools.Build import aho_graph, mtt_partition

from bmgedit.partitioning.Karger import Karger
from bmgedit.partitioning.GreedyBipartition import (greedy_bipartition,
                                                    gradient_walk_bipartition)
from bmgedit.partitioning.Louvain import Louvain, LouvainCustomCost


__author__ = 'David Schaller'


class Build2:
    """BUILD / MTT algorithm with minimal cost bipartition."""
    
    def __init__(self, R, L, F=None,
                 allow_inconsistency=True,
                 part_method='mincut',
                 cost_function=None, cost_function_args=None,
                 greedy_repeats=5,
                 weighted_mincut=False, triple_weights=None,):
        
        self.R = R
        self.L = L
        
        # forbidden triples --> activates MTT if non-empty
        self.F = F
        
        # allow inconsistencies or return False?
        self.allow_inconsistency = allow_inconsistency
        
        if part_method in ('mincut', 'karger', 'greedy', 'gradient_walk',
                             'louvain', 'louvain_cost'):
            self.part_method = part_method
        else:
            raise ValueError("unknown bipartition method "\
                             "'{}'".format(part_method))
        
        self.cost_function = cost_function
        self.cost_function_args = cost_function_args
        self.greedy_repeats = greedy_repeats
            
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
        part = list(nx.connected_components(aux_graph))
        
        if len(part) < 2:
            if not self.allow_inconsistency:
                return False
            else:
                cost, part = partition(L, self.part_method,
                                       cost_function=self.cost_function,
                                       args=self.cost_function_args,
                                       aux_graph=aux_graph,
                                       greedy_repeats=self.greedy_repeats)
                self.total_cost += cost
        
        node = PhyloTreeNode(-1)            # place new inner node
        for s in part:
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
        
        part, aux_graph = mtt_partition(L, R, F)
        
        if len(part) < 2:
            if not self.allow_inconsistency:
                return False
            else:
                cost, part = partition(L, self.part_method,
                                       cost_function=self.cost_function,
                                       args=self.cost_function_args,
                                       aux_graph=aux_graph,
                                       greedy_repeats=self.greedy_repeats)
                self.total_cost += cost
        
        node = PhyloTreeNode(-1)            # place new inner node
        for s in part:
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
    
    
def partition(L, method,
              cost_function=None, args=None,
              aux_graph=None,
              greedy_repeats=1):
    
    best_cost, best_bp = float('inf'), None
    
    if method == 'mincut':
        # Stoer–Wagner algorithm
        best_cost, best_bp = nx.stoer_wagner(aux_graph)
    
    elif method == 'karger':
        karger = Karger(aux_graph)
        
        for _, bp in karger.generate():
            cost = cost_function(bp, *args)
            
            if cost < best_cost:
                best_cost, best_bp = cost, bp
    
    elif method == 'greedy':
        
        for _ in range(greedy_repeats):
            cost, bp = greedy_bipartition(L, cost_function, args=args)
            if cost < best_cost:
                best_cost, best_bp = cost, bp
    
    elif method == 'gradient_walk':
        
        for _ in range(greedy_repeats):
            cost, bp = gradient_walk_bipartition(L, cost_function,
                                                 args=args)
            if cost < best_cost:
                best_cost, best_bp = cost, bp
    
    elif method == 'louvain':
        
        best_mod = float('-inf')
        best_cost = 0.0             # cost makes no sense here
        
        for _ in range(greedy_repeats):
            louv = Louvain(aux_graph, at_least_two=True)
            mod, bp = louv.modularities[-1], louv.partitions[-1]
            if mod > best_mod:
                best_mod, best_bp = mod, bp
    
    elif method == 'louvain_cost':
        
        for _ in range(greedy_repeats):
            louv = LouvainCustomCost(aux_graph, cost_function, args=args,
                                     at_least_two=True)
            cost, bp = louv.costs[-1], louv.partitions[-1]
            if cost < best_cost:
                best_cost, best_bp = cost, bp
    
    return best_cost, best_bp