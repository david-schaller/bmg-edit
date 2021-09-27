# -*- coding: utf-8 -*-

import itertools

import gurobipy as gp
from gurobipy import GRB

import networkx as nx
import tralda.tools.GraphTools as gt


__author__ = 'David Schaller'


class EditorSuperClass:
    
    status = {1: 'LOADED',
              2: 'OPTIMAL',
              3: 'INFEASIBLE',
              4: 'INF_OR_UNBD',
              5: 'UNBOUNDED',
              6: 'CUTOFF',
              7: 'ITERATION_LIMIT',
              8: 'NODE_LIMIT',
              9: 'TIME_LIMIT',
              10: 'SOLUTION_LIMIT',
              11: 'INTERRUPTED',
              12: 'NUMERIC',
              13: 'SUBOPTIMAL',
              14: 'INPROGRESS',
              15: 'USER_OBJ_LIMIT',}
    
    
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
        self.model = gp.Model()
        
        # binary varibles pairs of vertices
        self.e = self.model.addVars(self.E, vtype=GRB.BINARY, name='e')
        
        # finally sink all variables
        self.model.update()
        
        # set objective
        self._objective()
        
    
    def _objective(self):
        """Sets the objective function."""
        
        self.obj = gp.LinExpr()
        
        # add (1 - e_xy) * E_xy + (1 - E_xy) * e_xy for all x, y
        for x, y in self.E:
            if self.E[x, y]:
                self.obj += 1 - self.e[x, y]
            else:
                self.obj += self.e[x, y]
        
        self.model.setObjective(self.obj, GRB.MINIMIZE)
                                
                                    
    def optimize(self, time_limit=False):
        """Solves the editing problem."""
        
        # set new time limit temporarily
        if time_limit:
            original_time_limit = self.model.Params.timeLimit
            self.model.Params.timeLimit = time_limit
        
        self.model.optimize()
        
        # reset time limit
        if time_limit:
            self.model.Params.timeLimit = original_time_limit
    
    
    def _sink_free_colored(self):
        """Ensures that each vertex has out-arcs to every other color."""
        
        for x in self.graph.nodes():
            for col in self.color_dict:
                if self.graph.nodes[x]['color'] != col:
                    self.model.addConstr(
                        gp.quicksum((self.e[x, y]
                                     for y in self.color_dict[col])) >= 1,
                        name='c[{},{}]'.format(x, col)
                        )
        
        
    def get_solution(self):
        """Objective value and resulting digraph of the (current) solution."""
        
        # if self.model.status != GRB.OPTIMAL:
        #     raise RuntimeError('could not find optimal solution')
        
        sol_graph = nx.DiGraph()
        sol_graph.add_nodes_from(self.graph.nodes.data())
        
        for edge, value in self.model.getAttr('x', self.e).items():
            if value > 0.5:     # handle potential numerical inaccuracies
                sol_graph.add_edge(*edge)
        
        
        return self.model.ObjVal, sol_graph
    
    
    def get_status(self):
        
        return BMGEditor.status[self.model.status]
    
    
    def get_solve_time(self):
        
        return self.model.Runtime
    
    


class BMGEditor(EditorSuperClass):
    
    
    def __init__(self, graph):
        
        super().__init__(graph)
        
    
    def build_model(self):
        
        # arc variables and objective
        super().build_model()
        
        self._additional_variables()
        
        # subject to constraints
        self._sink_free_colored()
        self._inf_forb_triples()
        self._proper_hierarchy()
        
    
    def _additional_variables(self):
        
        # binary variable for each triple ab|c (= ba|c)
        triples = ((a, b, c)
                   for a, b, c in itertools.permutations(self.graph.nodes(), 3)
                   if a < b)
        self.t = self.model.addVars(triples, vtype=GRB.BINARY, name='t')
        
        # matrix for hierarchy repr. the tree
        self.M = self.model.addVars(((a, p) for a in self.graph.nodes()
                                            for p in range(self.n-2)),
                                vtype=GRB.BINARY, name='M')
        
        # a, b in custer p, but not c ?
        self.m = self.model.addVars(((a, b, c, p) for a, b, c in self.t
                                                  for p in range(self.n-2)),
                                    vtype=GRB.BINARY, name='m')
        
        # varibles for three-gamete condition
        self.C = self.model.addVars(
            ((p, q, gam) for p, q in itertools.combinations(range(self.n-2), 2)
                         for gam in range(3)),
            vtype=GRB.BINARY, name='C')
        
        # finally sink all variables
        self.model.update()
        
#        for v in self.model.getVars():
#            print(v)
    
    
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
                            self.model.addConstr(
                                (self.e[x, y1] - self.e[x, y2] -
                                 self.t[x, y1, y2] <= 0),
                                name='it[{},{},{}]'.format(x, y1, y2)
                                )
                            # forbidden triple x y1 | y2
                            self.model.addConstr(
                                (self.e[x, y1] + self.e[x, y2] +
                                 self.t[x, y1, y2] <= 2),
                                name='ft[{},{},{}]'.format(x, y1, y2)
                                )
                        else:
                            # informative triple x y1 | y2
                            self.model.addConstr(
                                (self.e[x, y1] - self.e[x, y2] -
                                 self.t[y1, x, y2] <= 0),
                                name='it[{},{},{}]'.format(x, y1, y2)
                                )
                            # forbidden triple x y1 | y2
                            self.model.addConstr(
                                (self.e[x, y1] + self.e[x, y2] +
                                 self.t[y1, x, y2] <= 2),
                                name='ft[{},{},{}]'.format(x, y1, y2)
                                )
        
        # set m((ab|c), p) correctly
        for a, b, c, p in self.m:
            
            expr = (self.M[a, p] + self.M[b, p] + (1 - self.M[c, p]) - 
                    3 * self.m[a, b, c, p])
            self.model.addConstr(expr >= 0)
            self.model.addConstr(expr <= 2)
        
        # informative and forbidden triples
        for a, b, c in self.t:
            expr = gp.quicksum((self.m[a, b, c, p] for p in range(self.n-2)))
            self.model.addConstr(self.t[a, b, c] <= expr)
            self.model.addConstr(expr <= (self.n - 2) * self.t[a, b, c])
        
        
    def _proper_hierarchy(self):
        """Ensures that a proper tree is constructed by the three-gamete
        condition."""
        
        self.model.addConstrs(
            (self.C[p, q, 0] >= self.M[a, q] - self.M[a, p]
             for p, q in itertools.combinations(range(self.n-2), 2)
             for a in self.graph.nodes()))
        self.model.addConstrs(
            (self.C[p, q, 1] >= self.M[a, p] - self.M[a, q]
             for p, q in itertools.combinations(range(self.n-2), 2)
             for a in self.graph.nodes()))
        self.model.addConstrs(
            (self.C[p, q, 2] >= self.M[a, q] + self.M[a, p] - 1
             for p, q in itertools.combinations(range(self.n-2), 2)
             for a in self.graph.nodes()))
                               
        self.model.addConstrs(
            (self.C[p, q, 0] + self.C[p, q, 1] + self.C[p, q, 2] <= 2
             for p, q in itertools.combinations(range(self.n-2), 2)))
    
    
class TwoBMGEditor(EditorSuperClass):

    def __init__(self, graph):

        super().__init__(graph)
        
        if len(self.color_dict) != 2:
            raise RuntimeError('not a 2-colored digraph')


    def build_model(self):
        
        # arc variables and objective
        super().build_model()
        
        # subject to constraints
        self._sink_free_colored()
        self._forbidden_subgraphs()


    def _forbidden_subgraphs(self):
        """Ensures that the graph does not contain a (2-colored) forbidden
        induced subgraph."""

        for col1, col2 in itertools.permutations(self.color_dict, 2):
            for x1, x2 in itertools.permutations(self.color_dict[col1], 2):

                for y1, y2 in itertools.permutations(self.color_dict[col2], 2):

                    # type (F1)
                    self.model.addConstr(
                        (self.e[x1, y1] + self.e[y1, x2] + self.e[y2, x2] - 
                         self.e[x1, y2] - self.e[y2, x1] <= 2),
                        name='F1[{},{},{},{}]'.format(x1, x2, y1, y2)
                        )

                    # type (F2)
                    self.model.addConstr(
                        (self.e[x1, y1] + self.e[y1, x2] + self.e[x2, y2] - 
                         self.e[x1, y2] <= 2),
                        name='F2[{},{},{},{}]'.format(x1, x2, y1, y2)
                        )

                for y1, y2, y3 in itertools.permutations(self.color_dict[col2],
                                                         3):
                    # since type (F3) subgraphs are symmetrical
                    if y1 > y2:
                        continue

                    # type (F3)
                    self.model.addConstr(
                        (self.e[x1, y1] + self.e[x1, y3] +
                         self.e[x2, y2] + self.e[x2, y3] - 
                         self.e[x1, y2] - self.e[x2, y1] <= 3),
                        name='F3[{},{},{},{},{}]'.format(x1, x2, y1, y2, y3)
                        )


class BinaryBMGEditor(EditorSuperClass):
    
    
    def __init__(self, graph):
        
        super().__init__(graph)
        
    
    def build_model(self):
        
        # arc variables and objective
        super().build_model()
        
        self._additional_variables()
        
        # subject to constraints
        self._sink_free_colored()
        self._ext_inf_triples()
        self._triple_consistency()
        
    
    def _additional_variables(self):
        
        # binary variable for each triple ab|c (= ba|c)
        triples = ((a, b, c)
                   for a, b, c in itertools.permutations(self.graph.nodes(), 3)
                   if a < b)
        self.t = self.model.addVars(triples, vtype=GRB.BINARY, name='t')
        
        # finally sink all variables
        self.model.update()
    
    
    def _ext_inf_triples(self):
        """Ensures that all informative triples of the graph are displayed
        by the tree but none of the forbidden triples."""
        
        for x in self.graph.nodes():
            for col in self.color_dict:
                if self.graph.nodes[x]['color'] != col:
                    for y1, y2 in itertools.permutations(self.color_dict[col],
                                                         2):
                        if x < y1:
                            # informative triple x y1 | y2
                            self.model.addConstr(
                                (self.e[x, y1] - self.e[x, y2] -
                                 self.t[x, y1, y2] <= 0),
                                name='it[{},{},{}]'.format(x, y1, y2)
                                )
                        else:
                            # informative triple x y1 | y2
                            self.model.addConstr(
                                (self.e[x, y1] - self.e[x, y2] -
                                 self.t[y1, x, y2] <= 0),
                                name='it[{},{},{}]'.format(x, y1, y2)
                                )
                        
                        # one direction is sufficient
                        if y1 < y2:
                            self.model.addConstr(
                                (self.e[x, y1] + self.e[x, y2] -
                                 self.t[y1, y2, x] <= 1),
                                name='it[{},{},{}]'.format(x, y1, y2)
                                )
    
    
    def _triple_consistency(self):
        """Ensures concistency of all triples."""
        
        # strictly dense triple set
        for triple in itertools.combinations(self.graph.nodes(), 3):
            a, b, c = sorted(triple)
            
            self.model.addConstr((self.t[a, b, c] + self.t[a, c, b] +
                                  self.t[b, c, a] == 1),
                                 name='d[{},{},{}]'.format(a, b, c))
                                    
        # compatibility based on 2-order inference rules
        for a, b, c, d in itertools.permutations(self.graph.nodes(), 4):
            #a, b, c, d = sorted(triple)
            
            ab = sorted([a, b])
            ad = sorted([a, d])
            bd = sorted([b, d])
                    
            self.model.addConstr((2 * self.t[(*ab, c)] + 2 * self.t[(*ad, b)] -
                                  self.t[(*bd, c)] - self.t[(*ad, c)] <= 2),
                                 name='cp[{},{},{},{}]'.format(a, b, c, d))
        

if __name__ == '__main__':
    
    from asymmetree.tools.PhyloTreeTools import random_colored_tree
    from asymmetree.analysis import (bmg_from_tree,
                                     is_bmg,
                                     binary_refinable_tree)
    
    # i = 0
    # while True:
    #     i += 1
    
    random_tree = random_colored_tree(12, 4,
                                      # binary=True,
                                      force_all_colors=True,)
    bmg = bmg_from_tree(random_tree)
    graph = gt.disturb_graph(bmg, 0.2, 0.2)
    
    # solver = BMGEditor(graph)
    solver = BinaryBMGEditor(graph)
    solver.build_model()
    solver.optimize(time_limit=None)
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
    be_success = True if binary_refinable_tree(sol_graph, mincut=False) else False
    print('is beBMG new:', be_success)
    
    colors = ['lightblue', 'red', 'green', 'yellow', 'black']
    color_map = []
    for v in sol_graph:
        color_map.append(colors[sol_graph.nodes[v]['color']-1])
    nx.draw(sol_graph,
            node_color=color_map,
            with_labels=True)
    
        # print('iteration', i)
        # if not success:
        #     break
    