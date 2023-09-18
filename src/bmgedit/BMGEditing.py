# -*- coding: utf-8 -*-

"""Heuristics for BMG editing."""

import itertools

from tralda.datastructures import LCA
from tralda.supertree.Build import best_pair_merge_first
from tralda.tools.GraphTools import sort_by_colors

from asymmetree.analysis.BestMatches import (informative_triples,
                                             binary_explainable_triples,
                                             bmg_from_tree)
from asymmetree.tools.PhyloTreeTools import (reconstruct_reconc_from_graph,)


from bmgedit.Build import Build2


__author__ = 'David Schaller'


class BMGEditor:
    """Wrapper for triple-based BMG editing heuristics."""
    
    def __init__(self, G, binary=False, binarization_mode='balanced', use_binary_triples=True):
        
        self.G = G
        self.color_dict = sort_by_colors(G)
        
        self.L = [v for v in G.nodes()]
        
        # informative triples
        if not binary:
            self.binarize = False
            if use_binary_triples:
                self.R = binary_explainable_triples(G, color_dict=self.color_dict)
            else:
                self.R = informative_triples(G, color_dict=self.color_dict)
        else:
            self.binarize = binarization_mode
            self.R = binary_explainable_triples(G, color_dict=self.color_dict)
            
        # current tree built by one of the heuristics
        self._tree = None


    def extract_consistent_triples(self):
        
        if not self._tree:
            raise RuntimeError('no tree has been built yet')
        
        lca = LCA(self._tree)
        
        return lca.consistent_triples(self.R)
    
    
    def build(self, method, objective='cost'):
        
        if objective == 'cost':
            minimize = True
            f_obj = unsatisfiability_cost
        elif objective == 'gain':
            minimize = False
            f_obj = satisfied_relations
        else:
            raise ValueError('unknown mode for objective ' \
                             'function: {}'.format(objective))

        method = method.lower()
        
        if method == 'bpmf':
            self._tree = best_pair_merge_first(self.R,self.L,
                                               triple_weights=None) 
        elif method in ('mincut', 'karger', 'greedy', 'gradient_walk',
                        'louvain', 'louvain_obj'):
            build = Build2(self.R, self.L,
                           allow_inconsistency=True,
                           binarize=self.binarize,
                           part_method=method,
                           obj_function=f_obj,
                           minimize=minimize,
                           obj_function_args=(self.G,),
                           weighted_mincut=True,
                          )
            self._tree = build.build_tree()
        else:
            raise ValueError('unknown partition method: {}'.format(method))
    
    
    def get_bmg(self, extract_triples_first=False,
                supply_inner_vertex_count=False):
        
        if not extract_triples_first:
            reconstruct_reconc_from_graph(self._tree, self.G)
            tree = self._tree
        else:
            R_consistent = self.extract_consistent_triples()
            build = Build2(R_consistent, self.L,
                           allow_inconsistency=False,
                           binarize=self.binarize,
                          )
            tree = build.build_tree()
            reconstruct_reconc_from_graph(tree, self.G)
        
        if supply_inner_vertex_count:
            return bmg_from_tree(tree), sum(1 for _ in tree.inner_nodes())
        else:
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


def satisfied_relations(partition, G):
    """Return the no. of arcs and non-arcs that are satisfied by a partition.
    
    Only considers pairs of vertices of different color.
    
    Parameters
    ----------
    partition : iterable of iterables of nodes
        A partition for (a subset of) the nodes in the graph.
    G : networkx.DiGraph
        A directed graph.
    
    Returns
    -------
    int
        The number of best match arcs and non-arcs between vertices of
        different colors that are satisfied if the partition corresponds to
        a split in a leaf-colored tree.
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
    
    count = 0
    
    for V, colors in zip(partition, color_sets):
        
        if not isinstance(V, set):
            V = set(V)
        
        for x, y in itertools.product(V, G.nodes()):
            
            y_color = G.nodes[y]['color']
            
            # skip pairs with the same color
            if G.nodes[x]['color'] == y_color:
                continue
            
            if y not in V:
                if not G.has_edge(x, y):
                    if colors.get(y_color):
                        count += 1
                elif not colors.get(y_color):
                    count += 1
            elif G.has_edge(x, y) and colors.get(y_color) == 1:
                count += 1
                    
    return count


def get_S1_S2_S3(partition, G):
    """Get the arcs and non-arcs that are satisfied by a partition of each type.
    
    S1: (x, y) notin E,
        x in V_i, y in V \ V_i,
        color of y is present in V_i
    S2: (x, y) in E,
        x in V_i, y in V \ V_i,
        color of y is not present in V_i
    S3: (x, y) in E,
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
        The arcs and non-arcs that are satisfied by a partition if this
        partition corresponds to a split in a tree sorted by the three types.
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
    
    S1, S2, S3 = [], [], []
    
    for V, colors in zip(partition, color_sets):
        
        if not isinstance(V, set):
            V = set(V)
        
        for x, y in itertools.product(V, G.nodes()):
            
            y_color = G.nodes[y]['color']
            
            # skip pairs with the same color
            if G.nodes[x]['color'] == y_color:
                continue
            
            if y not in V:
                if not G.has_edge(x, y):
                    if colors.get(y_color):
                        S1.append( (x, y) )
                elif not colors.get(y_color):
                    S2.append( (x, y) )
            elif G.has_edge(x, y) and colors.get(y_color) == 1:
                S3.append( (x, y) )
                    
    return S1, S2, S3


if __name__ == '__main__':
    
    from time import time
    
    from tralda.tools.GraphTools import disturb_graph, symmetric_diff
    
    from asymmetree.tools.PhyloTreeTools import random_colored_tree
    from asymmetree.analysis.BestMatches import is_bmg
    
    random_tree = random_colored_tree(30, 10, force_all_colors=True)
    bmg = bmg_from_tree(random_tree)
    disturbed = disturb_graph(bmg, 0.1, 0.1, preserve_properly_colored=True)
    
    print('-------------')
    print([(v, bmg.nodes[v]['color']) for v in bmg.nodes()])
    print('orginal', bmg.edges(), bmg.size())
    print('disturbed', disturbed.edges(), disturbed.size())
    print('-------------')
    print('original BMG is BMG:   ', bool(is_bmg(bmg)  ) )
    print('disturbed graph is BMG:', bool(is_bmg(disturbed)) )
    print('original vs disturbed:', symmetric_diff(bmg, disturbed))
    print('-------------')
    
    editor = BMGEditor(disturbed, binary=True)
    for method in ('Mincut', 'BPMF', 'Karger',
                   'Greedy', 'Gradient_Walk',
                   'Louvain', 'Louvain_Obj'):
        for objective in ('cost', 'gain'):
            start_time = time()
            editor.build(method, objective=objective)
            end_time = time() - start_time
            bmg1 = editor.get_bmg(extract_triples_first=False)
            bmg2 = editor.get_bmg(extract_triples_first=True)
            print('-----', method, objective, '-----')
            print('time:', end_time)
            print('orig./dist. vs BMG1 (no extr):',
                  symmetric_diff(bmg, bmg1), symmetric_diff(disturbed, bmg1))
            print('orig./dist. vs BMG2 (extr):',
                  symmetric_diff(bmg, bmg2), symmetric_diff(disturbed, bmg2))
