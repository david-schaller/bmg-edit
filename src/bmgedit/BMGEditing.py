# -*- coding: utf-8 -*-

"""Heuristics for BMG editing."""

import itertools

# import networkx as nx

from asymmetree.datastructures.PhyloTree import PhyloTree#, PhyloTreeNode

import asymmetree.best_matches.LeastResolvedTree as LRT
from asymmetree.best_matches.TrueBMG import bmg_from_tree
from asymmetree.tools.GraphTools import sort_by_colors
from asymmetree.datastructures.Tree import LCA
from asymmetree.tools.Build import best_pair_merge_first

from bmgedit.Build import Build2


__author__ = 'David Schaller'


class BMGEditor:
    """Wrapper for triple-based BMG editing heuristics."""
    
    def __init__(self, G):
        
        self.G = G
        self.color_dict = sort_by_colors(G)
        
        self.L = [v for v in G.nodes()]
        
        # informative and forbidden triples
        self.R, self.F = LRT.informative_forbidden_triples(G,
                             color_dict=self.color_dict)
        
        # current tree built by one of the heuristics
        self._tree = None
        
    
    def extract_consistent_triples(self):
        
        if not self._tree:
            raise RuntimeError('no tree has been built yet')
        
        lca = LCA(self._tree)
        
        return lca.consistent_triples(self.R)
    
    
    def build(self, mode):
        
        mode = mode.lower()
        
        if mode == 'bpmf':
            self._tree = best_pair_merge_first(self.R, self.L,
                                           triple_weights=None)
        elif mode in ('mincut', 'karger', 'greedy', 'gradient_walk',
                      'louvain', 'louvain_cost'):
            build = Build2(self.R, self.L,
                           allow_inconsistency=True,
                           part_method=mode,
                           cost_function=unsatisfiability_cost,
                           cost_function_args=(self.G,),
                           weighted_mincut=True)
            self._tree = build.build_tree()
    
    
    def get_bmg(self, extract_triples_first=False):
        
        if not extract_triples_first:
            self._tree.reconstruct_info_from_graph(self.G)
            return bmg_from_tree(self._tree)
        
        else:
            R_consistent = self.extract_consistent_triples()
            build = Build2(R_consistent, self.L, allow_inconsistency=False)
            tree = build.build_tree()
            tree.reconstruct_info_from_graph(self.G)
            return bmg_from_tree(tree)
        
        
def unsatisfiability_cost(partition, G):
    """Return the unsatisfiability cost of a partition.
    
    This cost is the number of unsatisfiable (non-)arcs induced by the
    bipartition.
    
    Parameters
    ----------
    partition : iterable of iterables of nodes
        A partition for (a subset of) the nodes in the graph.
    G : networkx.DiGraph
        A directed graph.
    
    Returns
    -------
    int
        The number of unsatisfiable best match arcs and non-arcs between
        elements in the partition, if this partition correspond to a split
        in a tree.
    """
    
    partition = list(partition)
    color_sets = [{} for V in partition]
    
    for V, colors in zip(partition, color_sets):
        for v in V:
            c = G.nodes[v]['color']
            if c not in colors:
                colors[c] = 1
            else:
                colors[c] += 1
    
    cost = 0
    
    for V, colors in zip(partition, color_sets):
        
        if not isinstance(V, set):
            V = set(V)
        
        for x, y in itertools.product(V, G.nodes()):
            
            y_color = G.nodes[y]['color']
            
            # skip pairs with the same color
            if G.nodes[x]['color'] == y_color:
                continue
            
            if y not in V:
                if G.has_edge(x, y):
                    if colors.get(y_color):
                        cost += 1
                elif not colors.get(y_color):
                    cost += 1
            elif not G.has_edge(x, y) and colors.get(y_color) == 1:
                cost += 1
                    
    return cost


def get_U1_U2_U3(partition, G):
    """Get the unsatisfiable relations of each type for a partition.
    
    U1: (x, y) in E,
        x in V_i, y in V \ V_i,
        color of y is present in V_i
    U2: (x, y) not in E,
        x in V_i, y in V \ V_i,
        color of y is not present in V_i
    U3: (x, y) not in E,
        distinct x, y in V_i
        y is the only vertex of its color in V_i.
    
    Parameters
    ----------
    partition : iterable of iterables of nodes
        A partition for (a subset of) the nodes in the graph.
    G : networkx.DiGraph
        A directed graph.
    
    Returns
    -------
    tuple of three lists of arcs
        The unsatisfiable arcs and non-arcs if the partition corresponds to a
        split in a tree sorted by the three types.
    """
    
    partition = list(partition)
    color_sets = [{} for V in partition]
    
    for V, colors in zip(partition, color_sets):
        for v in V:
            c = G.nodes[v]['color']
            if c not in colors:
                colors[c] = 1
            else:
                colors[c] += 1
    
    U1, U2, U3 = [], [], []
    
    for V, colors in zip(partition, color_sets):
        
        if not isinstance(V, set):
            V = set(V)
        
        for x, y in itertools.product(V, G.nodes()):
            
            y_color = G.nodes[y]['color']
            
            # skip pairs with the same color
            if G.nodes[x]['color'] == y_color:
                continue
            
            if y not in V:
                if G.has_edge(x, y):
                    if colors.get(y_color):
                        U1.append( (x, y) )
                elif not colors.get(y_color):
                    U2.append( (x, y) )
            elif not G.has_edge(x, y) and colors.get(y_color) == 1:
                U3.append( (x, y) )
                    
    return U1, U2, U3


if __name__ == '__main__':
    
    from time import time
    
    from asymmetree.tools.GraphTools import disturb_graph, symmetric_diff
    
    random_tree = PhyloTree.random_colored_tree(20, 10,
                                                force_all_colors=True)
    bmg = bmg_from_tree(random_tree)
    disturbed = disturb_graph(bmg, 0.1, 0.1, preserve_properly_colored=True)
    
    print('-------------')
    print([(v, bmg.nodes[v]['color']) for v in bmg.nodes()])
    print('orginal', bmg.edges(), bmg.size())
    print('disturbed', disturbed.edges(), disturbed.size())
    print('-------------')
    print('original BMG is BMG:   ', bool(LRT.is_bmg(bmg)  ) )
    print('disturbed graph is BMG:', bool(LRT.is_bmg(disturbed)) )
    print('original vs disturbed:', symmetric_diff(bmg, disturbed))
    print('-------------')
    
    editor = BMGEditor(disturbed)
    
    for mode in ('Mincut', 'BPMF', #'Karger',
                 'Greedy', 'Gradient_Walk',
                 'Louvain', 'Louvain_Cost'):
        start_time = time()
        editor.build(mode)
        end_time = time() - start_time
        bmg1 = editor.get_bmg(extract_triples_first=False)
        bmg2 = editor.get_bmg(extract_triples_first=True)
        print('-----', mode, '-----')
        print('time:', end_time)
        print('orig./dist. vs BMG1 (no extr):',
              symmetric_diff(bmg, bmg1), symmetric_diff(disturbed, bmg1))
        print('orig./dist. vs BMG2 (extr):',
              symmetric_diff(bmg, bmg2), symmetric_diff(disturbed, bmg2))
