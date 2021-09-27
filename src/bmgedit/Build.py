# -*- coding: utf-8 -*-

import networkx as nx

from tralda.datastructures.Tree import Tree, TreeNode
from tralda.supertree.Build import aho_graph, mtt_partition

from bmgedit.partitioning.Karger import Karger
from bmgedit.partitioning.GreedyBipartition import (greedy_bipartition,
                                                    gradient_walk_bipartition)
from bmgedit.partitioning.Louvain import Louvain, LouvainCustomObj
from bmgedit.partitioning.NumberPartition import balanced_coarse_graining


__author__ = 'David Schaller'


class Build2:
    """BUILD / MTT algorithm with optimal objective partition."""
    
    def __init__(self, R, L, F=None,
                 allow_inconsistency=True,
                 binarize=False,
                 part_method='mincut',
                 obj_function=None,
                 minimize=True,
                 obj_function_args=None,
                 greedy_repeats=5,
                 weighted_mincut=False, triple_weights=None,):
        
        self.R = R
        self.L = L
        
        # forbidden triples --> activates MTT if non-empty
        self.F = F
        
        # allow inconsistencies or return False?
        self.allow_inconsistency = allow_inconsistency
        
        # mode for tree binarization
        if not binarize:
            self.binarize = False
        elif isinstance(binarize, str) and \
            binarize.lower() in ('balanced', 'b'):
            self.binarize = 'b'
        elif isinstance(binarize, str) and \
            binarize.lower() in ('caterpillar', 'c'):
            self.binarize = 'c'
        else:
            raise ValueError('unknown binarization mode:' \
                             ' {}'.format(binarize))
        
        # partitioning method in case of inconsistencies
        if part_method in ('mincut', 'karger', 'greedy', 'gradient_walk',
                           'louvain', 'louvain_obj'):
            self.part_method = part_method
        else:
            raise ValueError("unknown bipartition method "\
                             "'{}'".format(part_method))
        self.obj_function =      obj_function
        self.minimize =          minimize
        self.obj_function_args = obj_function_args
        self.greedy_repeats =    greedy_repeats
            
        # parameters if bipartition method is mincut
        self.weighted_mincut = weighted_mincut
        self.triple_weights = triple_weights
    
    
    def build_tree(self, return_root=False):
        """Build a tree displaying all triples in R if possible.
        
        Keyword arguments:
            return_root - if True, return 'TreeNode' instead of
                'Tree' instance
        """
        
        self.total_obj = 0
        
        if self.F:
            root = self._mtt(self.L, self.R, self.F)
        else:
            root = self._aho(self.L, self.R)
            
        return root if return_root else Tree(root)
    
    
    def _trivial_case(self, L):
        
        if len(L) == 1:
            leaf = L.pop()
            return TreeNode(label=leaf)
        
        elif len(L) == 2:
            node = TreeNode()
            for _ in range(2):
                leaf = L.pop()
                node.add_child(TreeNode(label=leaf))
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
                obj, part = partition(L, self.part_method,
                                      obj_function=self.obj_function,
                                      minimize=self.minimize,
                                      args=self.obj_function_args,
                                      aux_graph=aux_graph,
                                      greedy_repeats=self.greedy_repeats)
                self.total_obj += obj
        
        if self.binarize == 'b' and len(part) > 2:
            part = balanced_coarse_graining(part)
        
        root = TreeNode()                   # place new inner node
        node = root
        for i, s in enumerate(part):
            Li, Ri = set(s), []
            for t in R:                     # construct triple subset
                if Li.issuperset(t):
                    Ri.append(t)
            Ti = self._aho(Li, Ri)          # recursive call
            if not Ti:
                return False                # raise False to previous call
            else:
                node.add_child(Ti)          # add roots of the subtrees
            
            # resolve to binary (caterpillar)
            if self.binarize == 'c' and i < len(part)-2:
                new_node = TreeNode()
                node.add_child(new_node)
                node = new_node
   
        return root
    
    
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
                obj, part = partition(L, self.part_method,
                                      obj_function=self.obj_function,
                                      minimize=self.minimize,
                                      args=self.obj_function_args,
                                      aux_graph=aux_graph,
                                      greedy_repeats=self.greedy_repeats)
                self.total_obj += obj
        
        if self.binarize == 'b' and len(part) > 2:
            part = balanced_coarse_graining(part)
        
        root = TreeNode()                   # place new inner node
        node = root
        for i, s in enumerate(part):
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
            
            # resolve to binary (caterpillar)
            if self.binarize == 'c' and i < len(part)-2:
                new_node = TreeNode()
                node.add_child(new_node)
                node = new_node
   
        return root
    
    
def partition(L, method,
              obj_function=None, minimize=True, args=None,
              aux_graph=None,
              greedy_repeats=1):
    
    best_obj = float('inf') if minimize else float('-inf')
    best_part = None
    
    if method == 'mincut':
        
        # Stoer–Wagner algorithm
        best_obj, best_part = nx.stoer_wagner(aux_graph)
    
    elif method == 'karger':
        karger = Karger(aux_graph)
        
        for _, part in karger.generate(runs=greedy_repeats):
            obj = obj_function(part, *args)
            
            if ((minimize and obj < best_obj) or
                (not minimize and obj > best_obj)):
                best_obj, best_part = obj, part
    
    elif method == 'greedy':
        
        for _ in range(greedy_repeats):
            obj, part = greedy_bipartition(L, obj_function, 
                                           minimize=minimize,
                                           args=args)
            if ((minimize and obj < best_obj) or
                (not minimize and obj > best_obj)):
                best_obj, best_part = obj, part
    
    elif method == 'gradient_walk':
        
        for _ in range(greedy_repeats):
            obj, part = gradient_walk_bipartition(L, obj_function,
                                                  minimize=minimize,
                                                  args=args)
            if ((minimize and obj < best_obj) or
                (not minimize and obj > best_obj)):
                best_obj, best_part = obj, part
    
    elif method == 'louvain':
        
        best_obj = float('-inf')    # modularity is always maximized
        
        for _ in range(greedy_repeats):
            louv = Louvain(aux_graph, at_least_two=True)
            obj, part = louv.modularities[-1], louv.partitions[-1]
            if obj > best_obj:
                best_obj, best_part = obj, part
    
    elif method == 'louvain_obj':
        
        for _ in range(greedy_repeats):
            louv = LouvainCustomObj(aux_graph, obj_function,
                                    minimize=minimize,
                                    args=args,
                                    at_least_two=True)
            obj, part = louv.objectives[-1], louv.partitions[-1]
            if ((minimize and obj < best_obj) or
                (not minimize and obj > best_obj)):
                best_obj, best_part = obj, part
    
    return best_obj, best_part
