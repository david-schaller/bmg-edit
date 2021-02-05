# -*- coding: utf-8 -*-

import random, itertools
import numpy as np


__author__ = 'David Schaller'


# ----------------------------------------------------------------------------
#                           Greedy bipartition
# ----------------------------------------------------------------------------

def partition_cut_value(partition, G):
    
    cut_value = 0
    
    for Vi, Vj in itertools.combinations(partition, 2):
        for x, y in itertools.product(Vi, Vj):
            if G.has_edge(x, y):
                cut_value += 1
    
    return cut_value
    

def _move(from_set, to_set, x):
    
    from_set.remove(x)
    to_set.add(x)


def greedy_bipartition(V, f_obj, minimize=True, args=()):
    """Randomized greedy bipartitioning with custom objective function."""
    
    if len(V) < 2:
        raise ValueError('iterable must have >=2 elements')
    
    V1, V2 = set(), set(V)
    remaining = set(V)
    
    best_obj = float('inf') if minimize else float('-inf')
    best_bp = None
    
    for _ in range(len(V)-1):
        
        best_obj_local = float('inf') if minimize else float('-inf')
        best_x = None
        
        for x in remaining:
            _move(V2, V1, x)
            obj = f_obj([V1, V2], *args)
            if ((minimize and obj < best_obj_local) or
                (not minimize and obj > best_obj_local)):
                best_obj_local = obj
                best_x = [x]
            elif obj == best_obj_local:
                best_x.append(x) 
            _move(V1, V2, x)
            
        x = random.choice(best_x)
        _move(V2, V1, x)
        remaining.remove(x)
        
        if ((minimize and best_obj_local < best_obj) or
            (not minimize and best_obj_local > best_obj)):
            best_obj = best_obj_local
            best_bp = [(list(V1), list(V2))]
        elif best_obj_local == best_obj:
            best_bp.append( (list(V1), list(V2)) )
            
    return best_obj, random.choice(best_bp)


# ----------------------------------------------------------------------------
#                        Adaptive walk bipartition
# ----------------------------------------------------------------------------

def gradient_walk_bipartition(V, f_obj, 
                              minimize=True, args=(),
                              initial_bipartition=None):
    
    if len(V) < 2:
        raise ValueError('iterable must have >=2 elements')
    
    # partition P
    P = [set(), set()]
    lookup = {}
    
    if not isinstance(V, (list, tuple)):
        V = list(V)
    
    if not initial_bipartition:
        for k, m in enumerate(np.random.permutation(len(V))):
            x = V[m]
            i = 0 if k <= len(V) / 2 else 1
            P[i].add(x)
            lookup[x] = i
    elif len(initial_bipartition) == 2:
        for i in range(2):
            for x in initial_bipartition[i]:
                P[i].add(x)
                lookup[x] = i
    else:
        raise ValueError('not a bipartition')
        
    best_obj = f_obj([*P], *args)
    best_bp = [list(s) for s in P]
    
    while True:
        
        best_obj_local = float('inf') if minimize else float('-inf')
        best_x = None
        
        for x in V:
            V1 = P[lookup[x]]
            V2 = P[(lookup[x] + 1) % 2]
            
            # avoid that one set becomes empty
            if len(V1) == 1:
                continue
            
            _move(V1, V2, x)
            obj = f_obj([V1, V2], *args)
            if ((minimize and obj < best_obj_local) or
                (not minimize and obj > best_obj_local)):
                best_obj_local = obj
                best_x = [x]
            elif obj == best_obj_local:
                best_x.append(x) 
            _move(V2, V1, x)
        
        if ((minimize and best_obj_local < best_obj) or
            (not minimize and best_obj_local > best_obj)):
            x = random.choice(best_x)
            V1 = P[lookup[x]]
            V2 = P[(lookup[x] + 1) % 2]
            _move(V1, V2, x)
            lookup[x] = (lookup[x] + 1) % 2
            
            best_obj = best_obj_local
            best_bp = [list(s) for s in P]
        else:
            break
    
    return best_obj, best_bp
        
    

if __name__ == '__main__':
    
    import networkx as nx
    
    G = nx.Graph()
    G.add_edges_from([('a', 'b'), ('a', 'c'), ('a', 'd'),
                      ('b', 'd'), ('c', 'd'), ('c', 'e'), ('d', 'e')])
        
    V = list(G.nodes())
    
    for i in range(3):
        print( *greedy_bipartition(V, partition_cut_value, args=(G,)) )
    
    print('---------')
    for i in range(3):
        print( *gradient_walk_bipartition(V, partition_cut_value, args=(G,)) )

