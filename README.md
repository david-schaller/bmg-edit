# Best match graph editing

[![license: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

Implemention of various heuristics and ILP formulations for best match graph (BMG) editing.


## Installation and Dependencies

`bmg-edt` requires Python 3.7 or higher. It has the following dependencies:

* [NetworkX](https://networkx.github.io/)
* [Scipy and Numpy](http://www.scipy.org/install.html)
* [Matplotlib](https://matplotlib.org/)
* [AsymmeTree](https://github.com/david-schaller/AsymmeTree)

In order to use the ILP versions for BMG editing, an installation of Gurobi Optimizer (9.0 or higher) or IBM ILOG CPLEX Optimization Studio (12.10 or higher) is required.
Moreover, the corresponding Python packages `gurobipy` or `docplex` must be installed.