# PhyCLIP (_Phylogenetic Clustering by Linear Integer Programming_)

## Overview

PhyCLIP is an integer linear programming (ILP) approach that assigns statistically-principled cluster membership to as many taxa as possible for a given **_rooted_** phylogenetic tree based on its pairwise patristic distance distribution. Other than the phylogeny, 3 additional inputs are required from the user: 
1. Minimum number of taxa in a cluster (_cs_)
2. Multiple (_gamma_) of deviations from the grand median of the mean pairwise patristic distance that defines the within-cluster limit.
3. False discovery rate (_fdr_) for rejecting the null hypotheses that the pairwise patristic distance distributions of every combinatorial pair of clusters are empirically equivalent to that should they form a single cluster.

A manuscript describing PhyCLIP is available here:  
MANUSCRIPT LINK

## Installation
PhyCLIP is written in Python 2.7 (no support for Python 3 currently) and depends on several python libraries and at least one ILP solver. 

To simplify the installation process, we highly reccomend that you use Anaconda, a free and open-source distribution of Python and package management system. Visit https://www.anaconda.com/download/ to download and install the **Python 2.7 version** distribution for your preferred OS. 

### Prerequisite: Python libraries    

PhyCLIP depends on several Python libraries: 
* numpy, scipy, statsmodels  (mathematical/statistical operations)
* ete3 (parsing phylogenetic trees) 

To install through Anaconda: 
```
$ conda install numpy scipy ete3 statsmodels
```

Alternatively, if you are using PyPi:
```
$ pip install numpy scipy ete3 statsmodels
```

### Prerequisite: ILP solver 
PhyCLIP currently supports two ILP solvers. You can choose either **_ONE_** to install depending on your access to these solvers: 

1. **Gurobi** optimizer (http://www.gurobi.com/) is a commercial linear and quadratic programming solver with free licenses available for academic users.
2. **GLPK** (GNU Linear Programming Kit, https://www.gnu.org/software/glpk/) is a free and open-source package intended for solving large-scale linear programming, mixed integer programming, and other related problems.

#### Gurobi
If you are a university user (i.e. you have internet access from a recognized academic domain, e.g. '.edu' addresss), we reccomend running PhyCLIP with the Gurobi optimizer. The easiest way to install Gurobi is via the Anaconda platform:  

1. Make sure you have Anaconda for Python 2.7 installed (see above). 

2. Install the Gurobi package via conda: 
```$ conda install gurobi```

3. You need to install a Gurobi licence next. Visit http://www.gurobi.com/registration/academic-license-reg to register for a free Gurobi account. Follow the instructions in the verification email from Gurobi to set your password and login to your Gurobi account via https://www.gurobi.com/login. 

4. You can now access https://user.gurobi.com/download/licenses/free-academic to request for a free academic license. To install the license: 
```$ grbgetkey XXX-XXX-XXX-XXX``` 
where command along with the license key stipulated in the License Detail page in your command/terminal prompt. Note that an active internet connection from a recognized academic domain (e.g. '.edu' addresss) is required. 

#### GLPK

```
$ conda install -c conda-forge pyomo glpk
```


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

You can use progarms such as FigTree (http://tree.bio.ed.ac.uk/software/figtree/) for tree rooting and/or conversion to the NEWICK format.

### Running phyclip.py

```
usage: phyclip.py [-h] -i INPUT_FILE [--treeinfo TREEINFO]
                  [--collapse_zero_branch_lengths {0,1}]
                  [--equivalent_zero_length EQUIVALENT_ZERO_LENGTH]
                  [--gam_method {MAD,Qn}] [--hypo_test {Kuiper,KS}]
                  [--preQ {0,1}]
                  [--subsume_sensitivity_induced_clusters {0,1}]
                  [--sensitivity_percentile SENSITIVITY_PERCENTILE]
                  [--ilp_verbose {0,1}]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input_file INPUT_FILE
                        Input file.
  --treeinfo TREEINFO   Master tree information file.
  --collapse_zero_branch_lengths {0,1}
                        Collapse nodes with zero branch lengths of tree prior
                        to running PhyCLIP (default = 0).
  --equivalent_zero_length EQUIVALENT_ZERO_LENGTH
                        Maximum branch length to be rounded to zero if the
                        --collapse_zero_branch_lengths flag is passed
                        (advanced option, default = 1e-07).
  --gam_method {MAD,Qn}
                        Method to estimate robust dispersion measure (default
                        = Qn).
  --hypo_test {Kuiper,KS}
                        Hypothesis test to use for statistical differentiation
                        of distance distributions (default = Kuiper).
  --preQ {0,1}          Perform Benjamini-Hochberg corrections of p-values
                        BEFORE filtering nodes that are < minimum cluster size
                        (advanced option, default = 0).
  --subsume_sensitivity_induced_clusters {0,1}
                        Subsume cluster-size sensitivity-induced clusters into
                        parent cluster (advanced option, default = 1).
  --sensitivity_percentile SENSITIVITY_PERCENTILE
                        Percentile of cluster size distribution under which a
                        cluster is considered to be sensitivity-induced
                        (advanced option, default = 25 percent).
  --ilp_verbose {0,1}   ILP solver verbose (default: 0)

```

For example, running phyclip.py by the default settings: 
```
$phyclip.py --input_file input_example.txt
```

Several output files will be generated by phylcip.py: 
* {tree_filename}\_treeinfo.txt - As statistical evaluation of the pariwise patristic distance distributions associated with the phylogeny can take some time, you can quicken future analyses of the _same_ phylogenetic tree by passing this file via the --treeinfo flag. 
* summary-stats_{input_filename}.txt - Tab-delimited file summarizing the clustering output (e.g. % clustered, mean pairwise patristic distance of clusters and its dispersion, etc.) of all parameter sets.
* cluster\_{gam\_method}\_{hypo\_test}\_cs\_fdr\_gam\_{input_filename}.txt - Tab-delimited file detailing the cluster-ID of every clustered taxon based on the input (cs, fdr, gam) parameter set denoted in the file name.
* tree\_{gam\_method}\_{hypo\_test}\_cs\_fdr\_gam\_{input_filename}.tre - NEXUS tree with FigTree annotations of clusters based on the input (cs, fdr, gam) parameter set denoted in the file name.
