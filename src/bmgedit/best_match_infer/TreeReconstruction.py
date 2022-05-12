# -*- coding: utf-8 -*-

import os, subprocess, itertools, time

from tralda.datastructures.Tree import TreeNode

from asymmetree.analysis.BestMatches import bmg_from_tree

from asymmetree.tools.PhyloTreeTools import (delete_and_reconnect,
                                             distances_from_root,
                                             parse_newick,)

from bmgedit.best_match_infer.ScenarioFileIO import matrix_to_phylip


__author__ = 'David Schaller'


# --------------------------------------------------------------------------
#                            TREE RECONSTRUCTION
#
# --------------------------------------------------------------------------

def bmg_by_tree_reconstruction(scenario, matrix_filename,
                               supply_rbmg=True, return_calltime=False):
    
    if return_calltime:
        start_time = time.time()
    
    nj_tree = neighbor_joining(scenario.genes, scenario.gene_index, matrix_filename)
    midpoint_rooting(nj_tree)
    bmg, rbmg = bmg_from_tree(nj_tree, supply_rbmg=supply_rbmg)
    
    if return_calltime:
        return bmg, rbmg, time.time() - start_time
    else:
        return bmg, rbmg
        
        

def neighbor_joining(leaves, leaf_index, matrix_filename,
                     return_calltime=False,
                     binary_path=None):
    """
    Keyword argument:
        return_calltime -- return the time needed for RapidNJ call
        binary_path -- path to 'qinfer' binary (if not available
                       within path)
    """
    
    if not binary_path:
        nj_command = "rapidnj"
    elif os.path.exists(binary_path):
         nj_command = binary_path
    else:
        raise FileNotFoundError("path to RapidNJ binary file '{}' does not "\
                                "exist".format(binary_path))
            
    start_time = time.time()
    
    try:
        output = subprocess.run([nj_command, matrix_filename, "-i", "pd"],
                                stdout=subprocess.PIPE)
    except:
        raise FileNotFoundError("calling RapidNJ failed")
        
    calltime = time.time() - start_time
    
    newick = output.stdout.decode()
    
    tree = parse_newick(newick)
    
    # restore (leaf) indeces and reconciliation
    leaf_dict = {leaf.label: leaf for leaf in leaves}
    index_counter = 0
    for v in tree.preorder():
        if not v.children:
            v.reconc = leaf_dict[v.label].reconc
        else:
            while index_counter in leaf_dict:
                index_counter += 1
            v.label = index_counter
            index_counter += 1
    
    if return_calltime:
        return tree, calltime
    else:
        return tree
    

def reroot(tree, node):
    
    edges = []                              # edges (u, v, distance) to be changed
    pos = node
    while pos.parent:
        edges.append( (pos.parent, pos, pos.dist) )    
        pos = pos.parent
    old_root = pos
    
    for u, v, dist in reversed(edges):      # change direction of edges
        u.remove_child(v)
        v.add_child(u)
        u.dist = dist
    
    node.detach()
    node.dist = 0.0
    tree.root = node
    
    if len(old_root.children) <= 1:         # delete old root if its out-degree is 1
        delete_and_reconnect(tree, old_root)


def midpoint_rooting(tree):
    
    # identify the two most distant leaves and their l.c.a.
    distance_dict = distances_from_root(tree)
    leaf_dict = tree.leaf_dict()
    max_dist, leaf1, leaf2, lca = float('-inf'), None, None, None
    max_label = 0
    
    for v in tree.preorder():
        if v.children:
            for c1, c2 in itertools.combinations(v.children, 2):
                for x in leaf_dict[c1]:
                    x_dist = distance_dict[x] - distance_dict[v]
                    for y in leaf_dict[c2]:
                        y_dist = distance_dict[y] - distance_dict[v]
                        if x_dist + y_dist > max_dist:
                            max_dist, leaf1, leaf2, lca = x_dist + y_dist, x, y, v
        max_label = max([max_label, v.label])
    
    # identify the edge for the new root
    edge, cut_at = None, 0
    
    pos, remaining = leaf1, max_dist/2
    while (not edge) and (pos is not lca):
        if remaining < pos.dist:
            edge = (pos.parent, pos)
            cut_at = pos.dist - remaining
        remaining -= pos.dist
        pos = pos.parent
    
    pos, remaining = leaf2, max_dist/2
    while (not edge) and (pos is not lca):
        if remaining < pos.dist:
            edge = (pos.parent, pos)
            cut_at = pos.dist - remaining
        remaining -= pos.dist
        pos = pos.parent
    
    # rerooting
    if edge:
        u, v = edge
        u.remove_child(v)
        new_root = TreeNode(label=max_label+1, dist=cut_at)
        u.add_child(new_root)
        new_root.add_child(v)
        v.dist -= new_root.dist
        reroot(tree, new_root)
    else:
        reroot(tree, lca)
        

def nj_from_numpy_matrix(leaves, leaf_index, matrix,
                         filename='temp.phylip'):

    matrix_to_phylip(filename, leaves, matrix)
    tree = neighbor_joining(leaves, leaf_index, filename)
    os.remove(filename)
    
    return tree
