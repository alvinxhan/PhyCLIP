# PhyCLIP (_Phylogenetic Clustering by Linear Integer Programming_)

## Overview

PhyCLIP is an integer linear programming (ILP) approach that optimally delineates a tree into statistically-principled clusters. Other than a **_rooted_** phylogeny, 3 additional inputs are required from the user: 
1. Minimum number of sequences (_S_) that can be quantified as a cluster.
2. Multiple of deviations (_gamma_) from the grand median of the mean pairwise sequence patristic distance that defines the within-cluster divergence limit. 
3. False discovery rate (_FDR_) to infer that the diversity observed for every combinatorial pair of output clusters is significantly distinct from one another.

A manuscript describing PhyCLIP is currently available on bioRxiv:  
Phylogenetic Clustering by Linear Integer Programming (PhyCLIP)
Alvin Xiaochuan Han, Edyth Parker, Frits Scholer, Sebastian Maurer-Stroh, Colin Russell
bioRxiv 446716; doi: [https://doi.org/10.1101/446716](https://doi.org/10.1101/446716)

Full documentation: 
[https://github.com/alvinxhan/PhyCLIP/wiki](https://github.com/alvinxhan/PhyCLIP/wiki)

**We highly encourage that you go through the ENTIRE DOCUMENTATION before starting your analysis.**

If you are based in an academic institution, you may go to the [Quickstart](https://github.com/alvinxhan/PhyCLIP/wiki/I.-Quickstart-(Beginners)) guide to promptly install and run PhyCLIP using the recommended procedure.  


**Clustering results of PhyCLIP maximises inclusivity and thus should not be not be interpreted as taxa linked by rapid transmission events. We reccomend you to use [Phydelity](https://github.com/alvinxhan/Phydelity), a re-purposing of PhyCLIP for the identification of putative transmission clusters.**
