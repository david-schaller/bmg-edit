# -*- coding: utf-8 -*-

"""
Scenario.

Wrapper class for species and gene tree. Compute statistics, BMG/RBMG.
"""

import networkx as nx

from asymmetree.treeevolve.GeneTree import observable_tree
from asymmetree.analysis.BestMatches import orthology_from_tree, bmg_from_tree
from asymmetree.tools.PhyloTreeTools import (reconc_sorted_leaves,
                                             distance_matrix,)


__author__ = 'David Schaller'


class Scenario:
    
    def __init__(self, S, TGT, dupl_rate, loss_rate, hgt_rate, OGT=None):
        
        self.S = S
        self.TGT = TGT
        self.DLH_rates = (dupl_rate, loss_rate, hgt_rate)
        
        self.OGT = OGT if OGT else observable_tree(TGT)
        
        self.reconc_dict, self.genes = reconc_sorted_leaves(self.OGT, return_list=True)
        self.gene_index = {gene: i for i, gene in enumerate(self.genes)}
        
        self._count_events()
        self._sort_species_to_subtrees()
        
        self.TOG = orthology_from_tree(self.OGT)
        self.bmg, self.rbmg = bmg_from_tree(self.OGT, supply_rbmg=True)
        
        self.bmg_subtrees, self.rbmg_subtrees = self.reduce_to_subtrees(self.bmg, self.rbmg)
    
    
    def _count_events(self):
        
        self.event_counts = [self.S.number_of_species, 0, 0, 0, 0]    # count S / D / L / H /
                                                                      # ancient duplications                
        for v in self.TGT.preorder():
            if v.event == 'D':
                self.event_counts[1] += 1
                if self.S.root.label == v.reconc[0]:
                    self.event_counts[4] += 1
            elif v.event == 'L':
                self.event_counts[2] += 1
            elif v.event == 'H':
                self.event_counts[3] += 1
    
    
    def _sort_species_to_subtrees(self):
        
        # sort the species into the subtrees of the first speciation event of S
        self.subtree_list, self.subtree_index = [], {}
        self.S_leaves = self.S.leaf_dict()
        for u in self.S.preorder():                          # exclude root of planted tree
            if len(u.children) > 1:
                for i in range(len(u.children)):
                    self.subtree_list.append([l.label for l in self.S_leaves[u.children[i]]])
                    self.subtree_index.update({item: i for item in self.subtree_list[i]})
                break
        
        # max. number of outgroup genes (in T) for each subtree of S
        self.outgroup_dict = {i: [] for i in range(len(self.subtree_list))}
        for gene in self.genes:
            for i in self.outgroup_dict.keys():
                if self.subtree_index[gene.reconc] != i:
                    self.outgroup_dict[i].append(gene)
    
    
    def get_data(self):
        
        return (self.S, self.genes, self.gene_index,
                self.subtree_list, self.subtree_index,
                self.outgroup_dict, self.reconc_dict)
    
    
    def rates_and_counts(self):
        """Return event rates and counts."""
        
        return list(self.DLH_rates) + self.event_counts
    
    
    def reduce_to_subtrees(self, full_bmg, full_rbmg):
        """Return a subgraph of the true RBMG with edges {u,v} for which the
        corresponding species are in the same subtree of root(S)."""
        
        bmg_subtrees = nx.DiGraph()
        rbmg_subtrees = nx.Graph()
        
        for G_subtrees, true_G in [(bmg_subtrees, full_bmg),
                                   (rbmg_subtrees, full_rbmg)]:
            G_subtrees.add_nodes_from(true_G.nodes(data=True))
            for u, v in true_G.edges:
                u_col = true_G.nodes[u]['color']
                v_col = true_G.nodes[v]['color']
                
                if self.subtree_index[u_col] == self.subtree_index[v_col]:
                    G_subtrees.add_edge(u, v)
        
        return bmg_subtrees, rbmg_subtrees
    
    
    def get_distance_matrix(self):
        
        _, D = distance_matrix(self.OGT, leaf_order=self.genes)
        
        return D
    
    
    def possible_edges_bmg(self):
        """Return the number of possible edges in the BMG i.e. for gene pairs for
        which an outgroup is available in S."""
        
        counts = [0 for i in range(len(self.subtree_list))]
        
        for g in self.genes:
            counts[self.subtree_index[g.reconc]] += 1
        
        possible_edges = 0
        for count in counts:
            possible_edges += (count * (count-1))
        
        return possible_edges
