# PhyCLIP (_Phylogenetic Clustering by Linear Integer Programming_)

PhyCLIP is a integer linear programming (ILP) approach that assigns statistically-principled cluster membership to as many taxa as possible for a given **_rooted_** phylogenetic tree based on its pairwise patristic distance distribution, subject to the following statistical constraints: 
1. Minimum number of taxa in a cluster (_cs_)
2. the false discovery rate (_fdr_) for rejecting the null hypotheses that the pairwise distance distributions associated with  sibling monophyletic clusters are empirically equivalent.
3. the multiple (_gamma_) of positive-skew median absolute deviations of pairwise distance admissible for all clusters 

PhyCLIP is written in C++/python and is currently distributed in two versions depending on the user's accessibility to ILP solvers: 
1. The **Gurobi** optimizer (http://www.gurobi.com/) is a commercial linear and quadratic programming solver with a free academic licenses available for university users.  
2. The GNU Linear Programmiong Kit (**GLPK**; https://www.gnu.org/software/glpk/) is a free and open-source package intended for solving large-scale linear programming, mixed integer programming, and other related problems.
