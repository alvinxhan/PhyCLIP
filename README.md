# PhyCLIP (_Phylogenetic Clustering by Linear Integer Programming_)

PhyCLIP is a integer linear programming (ILP) approach that assigns statistically-principled cluster membership to as many taxa as possible for a given **_rooted_** phylogenetic tree based on its pairwise patristic distance distribution, subject to the following statistical constraints: 
1. Minimum number of taxa in a cluster (_cs_)
2. Multiple (_gamma_) of deviations from the grand median of the mean pairwise patristic distance that defines the within-cluster limit.
3. False discovery rate (_fdr_) for rejecting the null hypotheses that the pairwise patristic distance distributions of every combinatorial pair of clusters are empirically equivalent to that should they form a single cluster.

A manuscript describing PhyCLIP is available here:  
MANUSCRIPT LINK

## Installation
PhyCLIP is written in Python 2.7 and is currently distributed in two versions depending on the user's accessibility to the following supported ILP solvers: 
1. **PuLP** is a free, open-source LP modeler written in python (https://github.com/coin-or/pulp). 
2. **Gurobi** optimizer (http://www.gurobi.com/) is a commercial linear and quadratic programming solver with free licenses available for academic users.
<!-- The GNU Linear Programmiong Kit (**GLPK**; https://www.gnu.org/software/glpk/) is a free and open-source package intended for solving large-scale linear programming, mixed integer programming, and other related problems. -->

### Dependencies 
Several dependencies are required in all versions, including: 
* numpy, scipy, statsmodels  (mathematical/statistical operations)
* ete3 (parsing phylogenetic trees in NEWICK format) 

To install the dependencies: 
```python
PyPi
$pip install numpy scipy ete3 statsmodels

Conda
$conda install numpy scipy ete3 statsmodels
```

### Gurobi version 
