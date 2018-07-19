# PhyCLIP (_Phylogenetic Clustering by Linear Integer Programming_)

## Overview

PhyCLIP is an integer linear programming (ILP) approach that optimally delineates a tree into statistically-principled clusters. Other than a **_rooted_** phylogeny, 3 additional inputs are required from the user: 
1. Minimum number of sequences (_S_) that can be quantified as a cluster.
2. Multiple of deviations (_gamma_) from the grand median of the mean pairwise sequence patristic distance that defines the within-cluster divergence limit. 
3. False discovery rate (_FDR_) to infer that the diversity observed for every combinatorial pair of output clusters is significantly distinct from one another.

A manuscript describing PhyCLIP is available here:  
MANUSCRIPT LINK

Full documentation: 
http://github.com/alvinxhan/PhyCLIP/wiki

## Installation
PhyCLIP is written in Python 2.7 (no support for Python 3 currently) and depends on several python libraries and at least one ILP solver. 

To simplify the installation process, we highly reccomend that you use Anaconda, a free and open-source distribution of Python and package management system. Visit http://www.anaconda.com/download/ to download and install the **Python 2.7 version** distribution for your preferred OS. 

### Prerequisite: Python libraries    

PhyCLIP depends on several Python libraries: 
* numpy, scipy, statsmodels  (mathematical/statistical operations)
* pathos (multiprocessing)
* ete3 (parsing phylogenetic trees) 

To install the dependencies, go to Terminal (Mac/Linux) or Command/Anaconda Prompt (Windows): 
```
$ conda install -c etetoolkit ete3
$ conda install -c conda-forge pathos
$ conda install numpy scipy statsmodels

```

Alternatively, if you are using PyPi:
```
$ pip install ete3
$ pip install numpy scipy statsmodels pathos
```

### Prerequisite: ILP solver 
PhyCLIP currently supports two ILP solvers. You can choose either **_ONE_** to install depending on your access to these solvers: 

1. **Gurobi** optimizer (http://www.gurobi.com/) is a commercial linear and quadratic programming solver with FREE licenses available for academic users.
2. **GLPK** (GNU Linear Programming Kit, http://www.gnu.org/software/glpk/) is a free and open-source package intended for solving large-scale linear programming, mixed integer programming, and other related problems.

If you are a university user (i.e. you have internet access from a recognized academic domain, e.g. '.edu' addresss), we highly reccomend running PhyCLIP with the Gurobi solver. GLPK performs poorly in terms of both speed and solvability (GLPK version 4.65 solved only 2 of the 87 standard test-set mixed-integer programming models whereas Gurobi is the fastest solver for all 87 benchmark problems, see http://plato.asu.edu/ftp/milpc.html). 

Furthermore, as with any other linear programming problems, it is possible to obtain multiple optimal solutions. Currently, GLPK can only return ONE solution that is guaranteed to be the global optimal if and only if the feasible region is convex and bounded. However, this may not always be the case. Gurobi, on the other hand, generates a solution pool which may include > 1 optimal solution.


#### Gurobi
**IMPORTANT: Take note of the version of Gurobi you are using (printed on summary-stats_*.txt file, see Output files)**. Gurobi is updated periodically to enhance solver performance. Correspondingly, we do find minor changes to PhyCLIP's clustering results in some cases as a result of Gurobi updates.

The easiest way to install Gurobi is via the Anaconda platform:

1. Make sure you have Anaconda for Python 2.7 installed (see above). 

2. Install the Gurobi package via conda:

```
$ conda config --add channels http://conda.anaconda.org/gurobi
$ conda install gurobi
```

3. You need to install a Gurobi licence next. Visit http://www.gurobi.com/registration/academic-license-reg to register for a free Gurobi account. Follow the instructions in the verification email from Gurobi to set your password and login to your Gurobi account via http://www.gurobi.com/login. 

4. You can now access http://user.gurobi.com/download/licenses/free-academic to request for a free academic license. Once requested, you will be brought to the License Detail webpage.

5. To install the license, go to Terminal/Command Prompt:  ```$ grbgetkey XXXXXXXX-XXXX-XXXX-XXXX-XXXXXXXXXXXX``` where ```XXXXXXXX-XXXX-XXXX-XXXX-XXXXXXXXXXXX``` is your unique license key shown in the License Detail webpage. Note that an active internet connection from a recognized academic domain (e.g. '.edu' addresss) is required for this step. 

#### GLPK

You can easily install GLPK via Anaconda as well: 
```
$ conda install -c conda-forge glpk
```

Alternatively, you can also install GLPK from source, go to http://ftp.gnu.org/gnu/glpk/ and download the latest distribution of glpk as a tarball. You can find installation information in the documentation provided.


### Install PhyCLIP 

Finally, install phyclip.py by: 
```
$ cd PhyCLIP-master/ 
$ python setup.py install
```
You may need sudo privileges for system-wide installation. Otherwise, it is also possible to use phyclip.py locally by adding the phyclip_modules folder to your $PYTHONPATH.

## Usage 

### Input file format
Prior to running phyclip.py, prepare the input text file indicating the path of the **rooted** input phylogenetic tree in **NEWICK** format and list of parameter sets to run in the following format: 
```
/path/to/input_newick_tree.nwk
cs1,fdr1,gam1
cs2,fdr2,gam2
...

E.g. (input_example.txt in examples/ folder): 
examples/example.nwk
3,0.2,2
5,0.2,1
```

You can use programs such as FigTree (http://tree.bio.ed.ac.uk/software/figtree/) for tree rooting and/or conversion to the NEWICK format.

### Running phyclip.py

```
usage: phyclip.py [-h] -i INPUT_FILE [--treeinfo TREEINFO] [--no_treeinfo]
                  [--optimise {intermediate,high}] [--prior PRIOR]
                  [--prior_weights PRIOR_WEIGHTS] [--pdf_tree]
                  [--collapse_zero_branch_lengths {0,1}]
                  [--equivalent_zero_length EQUIVALENT_ZERO_LENGTH]
                  [--gam_method {mad,qn}] [--hypo_test {kuiper,kolsmi}]
                  [--preQ {0,1}]
                  [--subsume_sensitivity_induced_clusters {0,1}]
                  [--sensitivity_percentile SENSITIVITY_PERCENTILE]
                  [--subsume_subclusters {0,1}] [--solver {glpk,gurobi}]
                  [--solver_verbose {0,1}] [--solver_check]
                  [--threads THREADS]

Phylogenetic Clustering by Linear Integer Programming (PhyCLIP) v0.1

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input_file INPUT_FILE
                        Input file. See manual for format details.
  --treeinfo TREEINFO   *_treeinfo.txt file generated from PREVIOUS PhyCLIP
                        run of the SAME phylogenetic tree.
  --no_treeinfo         Stop PhyCLIP from generating *_treeinfo.txt file.
  --optimise {intermediate,high}
                        PhyCLIP automatically searches for the clustering
                        result from the OPTIMAL input parameter set. See
                        documentation for details. Must have >1 input
                        parameter set in input file.
  --prior PRIOR         Prior information on clustered taxa. File format same
                        as PhyCLIP cluster text file output (i.e. CLUSTER-
                        ID{tab}SEQUENCE-NAME). Only works with gurobi solver,
                        no support for glpk.
  --prior_weights PRIOR_WEIGHTS
                        Weights on importance/confidence between prior
                        clusters to remain clustered together. File format:
                        CLUSTER-ID{tab}INT/FLOAT_WEIGHT{\newline}. Equal
                        weights will be assumed if none is given.
  --pdf_tree            PDF tree output annotated with cluster results.
  --collapse_zero_branch_lengths {0,1}
                        Collapse internal nodes with zero branch lengths of
                        tree before running PhyCLIP (default = 0).
  --equivalent_zero_length EQUIVALENT_ZERO_LENGTH
                        Maximum branch length to be rounded to zero if the
                        --collapse_zero_branch_lengths flag is passed
                        (advanced option, default = 1e-06).
  --gam_method {mad,qn}
                        Method to estimate robust dispersion measure (default
                        = qn).
  --hypo_test {kuiper,kolsmi}
                        Hypothesis test to use for statistical differentiation
                        of distance distributions (default = kuiper).
  --preQ {0,1}          Perform Benjamini-Hochberg corrections of p-values
                        BEFORE filtering nodes that are < minimum cluster size
                        (advanced option, default = 0).
  --subsume_sensitivity_induced_clusters {0,1}
                        Subsume cluster-size sensitivity-induced clusters into
                        parent cluster (default = 1).
  --sensitivity_percentile SENSITIVITY_PERCENTILE
                        Percentile of cluster size distribution under which a
                        cluster is considered to be sensitivity-induced
                        (advanced option, default = 25%).
  --subsume_subclusters {0,1}
                        Subsume sub-clusters into their respective parent
                        clusters (default = 0).
  --solver {glpk,gurobi}
                        Preferred ILP solver IF more than one solvers are
                        available (default: gurobi).
  --solver_verbose {0,1}
                        ILP solver verbose (default: 0)
  --solver_check        Check available ILP solver(s) installed.
  --threads THREADS     Number of threads (default = all).
```

For example, running phyclip.py by the default settings: 
```
$phyclip.py --input_file input_example.txt
```

Several output files will be generated by phylcip.py: 
* {tree_filename}\_treeinfo.txt - As statistical evaluation of the pariwise sequence   patristic distance distributions associated with the phylogeny can take some time, you can quicken future analyses of the _same_ phylogenetic tree by passing this file via the --treeinfo flag. 
* summary-stats_{input_filename}.txt - Tab-delimited file summarizing the clustering output (e.g. % clustered, mean pairwise patristic distance of clusters and its dispersion, etc.) of all parameter sets.
* cluster\_{gam\_method}\_{hypo\_test}\_{sol\_index}_{_S_}\_{_FDR_}\_{_gamma_}\_{tree_filename}.txt - Tab-delimited file detailing the cluster-ID of every clustered taxon  for the input (_S_, _FDR_, _gamma_) parameter set denoted in the filename.
* tree\_{gam\_method}\_{hypo\_test}\_{sol\_index}\_{_S_}\_{_FDR_}\_{_gamma_}\_{tree_filename}.tre - NEXUS tree with FigTree annotations of clusters for the input (_S_, _FDR_, _gamma_) parameter set denoted in the filename.
