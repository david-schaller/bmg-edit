# -*- coding: utf-8 -*-

import itertools

from docplex.mp.model import Model

import networkx as nx
import tralda.tools.GraphTools as gt


__author__ = 'David Schaller'


class BMGEditor:
    
    
    def __init__(self, graph):
        
        self.graph = graph
        self.n = self.graph.order()
        self.color_dict = gt.sort_by_colors(graph)
        
    
    def build_model(self):
        
        # constant values repr. the edges in the original graph
        self.E = {(x, y): int(self.graph.has_edge(x, y))
                     for x in self.graph.nodes()
                     for y in self.graph.nodes()
                     if x != y}
        
        # build model with variables and objective function
        self.model = Model('BMG_Editing')
        self._variables()
        
        # subject to constraints
        self._well_colored()
        self._inf_forb_triples()
        self._proper_hierarchy()
        
        self._objective()
        
        self.model.print_information()
        
    
    def _variables(self):
        
        # binary varibles pairs of vertices
        self.e = self.model.binary_var_dict(self.E, name='e')
        
        # binary variable for each triple ab|c (= ba|c)
        triples = ((a, b, c)
                   for a, b, c in itertools.permutations(self.graph.nodes(), 3)
                   if a < b)
        self.t = self.model.binary_var_dict(triples, name='t')
        
        # matrix for hierarchy repr. the tree
        self.M = self.model.binary_var_matrix(self.graph.nodes(),
                                              range(self.n-2),
                                              name='M')
        
        # a, b in custer p, but not c ?
        self.m = self.model.binary_var_dict(((a, b, c, p)
                                             for a, b, c in self.t
                                             for p in range(self.n-2)),
                                            name='m')
        
        # varibles for three-gamete condition
        self.C = self.model.binary_var_dict(
            ((p, q, gam) for p, q in itertools.combinations(range(self.n-2), 2)
                         for gam in range(3)),
             name='C')
        
    
    def _objective(self):
        """Sets the objective function."""
        
        self.obj = self.model.linear_expr()
        
        # add (1 - e_xy) * E_xy + (1 - E_xy) * e_xy for all x, y
        for x, y in self.E:
            if self.E[x, y]:
                self.obj += 1 - self.e[x, y]
            else:
                self.obj += self.e[x, y]
        
        self.model.minimize(self.obj)
        
    
    def _well_colored(self):
        """Ensures that each vertex has out-arcs to every other color."""
        
        for x in self.graph.nodes():
            for col in self.color_dict:
                if self.graph.nodes[x]['color'] != col:
                    self.model.add_constraint(
                        self.model.sum((self.e[x, y]
                                       for y in self.color_dict[col])) >= 1,
                        ctname='c[{},{}]'.format(x, col)
                        )
    
    
    def _inf_forb_triples(self):
        """Ensures that all informative triples of the graph are displayed
        by the tree but none of the forbidden triples."""
        
        for x in self.graph.nodes():
            for col in self.color_dict:
                if self.graph.nodes[x]['color'] != col:
                    for y1, y2 in itertools.permutations(self.color_dict[col],
                                                         2):
                        if x < y1:
                            # informative triple x y1 | y2
                            self.model.add_constraint(
                                (self.e[x, y1] - self.e[x, y2] -
                                 self.t[x, y1, y2] <= 0),
                                ctname='it[{},{},{}]'.format(x, y1, y2)
                                )
                            # forbidden triple x y1 | y2
                            self.model.add_constraint(
                                (self.e[x, y1] + self.e[x, y2] +
                                 self.t[x, y1, y2] <= 2),
                                ctname='ft[{},{},{}]'.format(x, y1, y2)
                                )
                        else:
                            # informative triple x y1 | y2
                            self.model.add_constraint(
                                (self.e[x, y1] - self.e[x, y2] -
                                 self.t[y1, x, y2] <= 0),
                                ctname='it[{},{},{}]'.format(x, y1, y2)
                                )
                            # forbidden triple x y1 | y2
                            self.model.add_constraint(
                                (self.e[x, y1] + self.e[x, y2] +
                                 self.t[y1, x, y2] <= 2),
                                ctname='ft[{},{},{}]'.format(x, y1, y2)
                                )
        
        # set m((ab|c), p) correctly
        for a, b, c, p in self.m:
            
            expr = (self.M[a, p] + self.M[b, p] + (1 - self.M[c, p]) - 
                    3 * self.m[a, b, c, p])
            self.model.add_constraint(expr >= 0)
            self.model.add_constraint(expr <= 2)
        
        # informative and forbidden triples
        for a, b, c in self.t:
            expr = self.model.sum((self.m[a, b, c, p] for p in range(self.n-2)))
            self.model.add_constraint(self.t[a, b, c] <= expr)
            self.model.add_constraint(expr <= (self.n - 2) * self.t[a, b, c])
        
        
    def _proper_hierarchy(self):
        """Ensures that a proper tree is constructed by the three-gamete
        condition."""
        
        self.model.add_constraints(
            (self.C[p, q, 0] >= self.M[a, q] - self.M[a, p]
             for p, q in itertools.combinations(range(self.n-2), 2)
             for a in self.graph.nodes()))
        self.model.add_constraints(
            (self.C[p, q, 1] >= self.M[a, p] - self.M[a, q]
             for p, q in itertools.combinations(range(self.n-2), 2)
             for a in self.graph.nodes()))
        self.model.add_constraints(
            (self.C[p, q, 2] >= self.M[a, q] + self.M[a, p] - 1
             for p, q in itertools.combinations(range(self.n-2), 2)
             for a in self.graph.nodes()))
                               
        self.model.add_constraints(
            (self.C[p, q, 0] + self.C[p, q, 1] + self.C[p, q, 2] <= 2
             for p, q in itertools.combinations(range(self.n-2), 2)))        
                                
                                    
    def optimize(self, time_limit=False):
        """Solves the editing problem."""
        
        # set new time limit temporarily
        if time_limit:
            original_time_limit = self.model.time_limit
            self.model.set_time_limit(time_limit)
        
        self.solution = self.model.solve(url=None, key=None)
        
        # reset time limit
        if time_limit:
            self.model.time_limit = original_time_limit
        
        print('Solve details:')
        print(self.model.get_solve_details())
        
        
    def get_solution(self):
        """Objective value and resulting digraph of the (current) solution."""
        
        sol_graph = nx.DiGraph()
        sol_graph.add_nodes_from(self.graph.nodes.data())
        
        for edge, variable in self.e.items():
            
            # handle potential numerical inaccuracies
            if self.solution.get_value(variable) > 0.5:
                sol_graph.add_edge(*edge)
        
        
        return self.solution.get_objective_value(), sol_graph
    
    
    def get_status(self):
        
        return self.solution.solve_details.status
    
    
    def get_solve_time(self):
        
        return self.solution.solve_details.time
        

if __name__ == '__main__':
    
    from asymmetree.tools.PhyloTreeTools import random_colored_tree
    from asymmetree.analysis import (bmg_from_tree,
                                     is_bmg,)
    
    random_tree = random_colored_tree(7, 4, force_all_colors=True)
    bmg = bmg_from_tree(random_tree)
    graph = gt.disturb_graph(bmg, 0.2, 0.2)
    
    solver = BMGEditor(graph)
    solver.build_model()
    solver.optimize(time_limit=3)
    solver.get_solution()
    opt_val, sol_graph = solver.get_solution()
    
    print('Status:', solver.get_status())
    print('Runtime:', solver.get_solve_time(), 's')
    
    print('-------------')
    print([(x, bmg.nodes[x]['color']) for x in bmg.nodes])
    print(bmg.edges)
    print('-------------')
    print([(x, sol_graph.nodes[x]['color']) for x in sol_graph.nodes])
    print(sol_graph.edges)
    print('-------------')
    
    print('original vs disturbed:', gt.symmetric_diff(bmg, graph))
    print('original vs new:', gt.symmetric_diff(sol_graph, bmg))
    print('disturbed vs new', gt.symmetric_diff(sol_graph, graph),
          'should be', opt_val,
          'and <=', gt.symmetric_diff(bmg, graph))
    print(sol_graph.order(), sol_graph.size())
    print('is BMG orig:', True if is_bmg(bmg) else False)
    print('is BMG dist:', True if is_bmg(graph) else False)
    success = True if is_bmg(sol_graph) else False
    print('is BMG new:', success)
    