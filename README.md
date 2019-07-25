# PhyCLIP (_Phylogenetic Clustering by Linear Integer Programming_)

## Updates  

* 25-Jul-2019: v2.1 - Fixed int64 declaration for inter-platform compatibility; tested on Windows and MAC 
* 16-May-2019: PhyCLIP version 2.0 release  
  * PhyCLIP is now fully Cythonised. If you have been using PhyCLIP previously, please make sure that Cython is installed (if you are using Anaconda, type ```conda install -c anaconda cython``` into your terminal/command prompt)
  * **~270x improvement in computational speed** for a phylogenetic tree with ~1,200 taxa when compared to version 1.0.  
  * Decreased peak memory usage by >2x  
  * ```--force``` option available to heuristically force putative subtrees with low divergence (low information) to cluster under the stipulated statistical framework (More information to come in wiki documentation).  
  * Bug fixes.  

## Overview  

PhyCLIP is an integer linear programming (ILP) approach that optimally delineates a tree into statistically-principled clusters. Other than a **_rooted_** phylogeny, 3 additional inputs are required from the user: 
1. Minimum number of sequences (_S_) that can be quantified as a cluster.
2. Multiple of deviations (_gamma_) from the grand median of the mean pairwise sequence patristic distance that defines the within-cluster divergence limit. 
3. False discovery rate (_FDR_) to infer that the diversity observed for every combinatorial pair of output clusters is significantly distinct from one another.

A manuscript describing PhyCLIP (version 1.0) is published in *Molecular Biology and Evolution*:  
Phylogenetic Clustering by Linear Integer Programming (PhyCLIP)
Alvin Xiaochuan Han, Edyth Parker, Frits Scholer, Sebastian Maurer-Stroh, Colin Russell
*Molecular Biology and Evolution*, msz053; doi: [https://doi.org/10.1093/molbev/msz053](https://doi.org/10.1093/molbev/msz053)

Full documentation: 
[https://github.com/alvinxhan/PhyCLIP/wiki](https://github.com/alvinxhan/PhyCLIP/wiki)

**We highly encourage that you go through the ENTIRE DOCUMENTATION before starting your analysis.**

If you are based in an academic institution, you may go to the [Quickstart](https://github.com/alvinxhan/PhyCLIP/wiki/I.-Quickstart-(Beginners)) guide to promptly install and run PhyCLIP using the recommended procedure.  

**Clustering results of PhyCLIP maximises inclusivity and thus should not be interpreted as taxa linked by rapid transmission events. We reccomend you to use [Phydelity](https://github.com/alvinxhan/Phydelity), a re-purposing of PhyCLIP for the identification of putative transmission clusters.**  
