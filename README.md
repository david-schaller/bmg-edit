# Best match graph editing

[![license: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

Implemention of various heuristics and ILP formulations for best match graph (BMG) editing.

Best match graphs (BMGs) are a class of colored digraphs that naturally appear in mathematical phylogenetics as a representation of the pairwise most closely related genes among multiple species. An arc connects a gene x with a gene y from another species (encoded as the color of the nodes) Y whenever it is one of the phylogenetically closest relatives of x.
This package contains various methods to edit an arbitrary vertex-colored digraph to a valid BMG, i.e., a graph that has a certain representation as (leaf-colored) tree.


## Installation and Dependencies

`bmg-edt` requires Python 3.7 or higher. It has the following dependencies:

* [NetworkX](https://networkx.github.io/)
* [Scipy and Numpy](http://www.scipy.org/install.html)
* [Matplotlib](https://matplotlib.org/)
* [AsymmeTree](https://github.com/david-schaller/AsymmeTree)
* [tralda](https://github.com/david-schaller/tralda)

In order to use the ILP versions for BMG editing, an installation of Gurobi Optimizer (9.0 or higher) or IBM ILOG CPLEX Optimization Studio (12.10 or higher) is required.
Moreover, the corresponding Python packages `gurobipy` or `docplex`, respectively, must be installed.

## Usage

The functions in `bmg-edit` require a NetworkX `DiGraph` as input.
Moreover, all nodes must have an attribute 'color'.

### ILP

The following classes for optimal BMG editing are available in the module `ilp.GurobiBMG` (requires an installation of Gurobi Optimzizer):

- **BMGEditor** edits the input graph with an arbitrary number of colors to the closest BMG.
- **BinaryBMGEditor** edits the input graph with an arbitrary number of colors to the closest BMG that can be explained by a binary tree.
- **TwoBMGEditor** edits the input graph with at most two distinct colors to the closest (2-)BMG.

<details>
<summary>Example usage: (Click to expand)</summary>

    solver = BMGEditor(input_graph)
    solver.build_model()
    
    # run the optimization with an optional time limit in seconds
    solver.optimize(time_limit=None)
    
    optimal_editing_cost, solution_graph = solver.get_solution()
    
</details>

The following classes for optimal BMG editing are available in the module `ilp.CplexBMG` (requires an installation of IBM ILOG CPLEX Optimization Studio):

- **BMGEditor** edits the input graph with an arbitrary number of colors to the closest BMG.

<details>
<summary>Example usage: (Click to expand)</summary>

    solver = BMGEditor(input_graph)
    solver.build_model()
    
    # run the optimization with an optional time limit in seconds
    solver.optimize(time_limit=3)
    solver.get_solution()
    
    optimal_editing_cost, solution_graph = solver.get_solution()
    
</details>

### Heuristics Algorithms

The package also implements various heuristic approaches for BMG editing.
Some of these methods are based on the unsatisfiable relation (UR) which are insertions or deletions of arcs that are associated with a certain inner node of the tree that explains the editing results.
More precisely, the heuristics construct this tree in a top-down manner (i.e., starting with the root) and attempt, in each step, to minimize the UR (see refrenced paper below for details).

The class `BMGEditor` in the module `BMGEditor` manages the editing:

    editor = BMGEditor(disturbed, binary=True)
    editor.build('Mincut', objective='cost')
    solution_graph = editor.get_bmg(extract_triples_first=False)

<details>
<summary>The following methods are available (first parameter of the `build` method): (Click to expand)</summary>

- 'Mincut'
- 'BPMF'
- 'Karger'
- 'Greedy'
- 'Gradient_Walk'
- 'Louvain'
- 'Louvain_Obj'

See the paper for an explanation of these methods.
</details>

## Citation and References

If you use `bmg-edit` in your project or code from it, please consider citing:

* **Schaller, D., Geiß, M., Hellmuth, M., Stadler, P. F. (2021) Heuristic Algorithms for Best Match Graph Editing.
Algorithms for Molecular Biology (in press).
arXiv:2103.07280 [math.CO].**
