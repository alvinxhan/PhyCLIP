# PhyCLIP (_Phylogenetic Clustering by Linear Integer Programming_)

## Overview

PhyCLIP is an integer linear programming (ILP) approach that assigns statistically-principled clade membership to as many taxa as possible for a given **_rooted_** phylogenetic tree based on its pairwise patristic distance distribution, subject to the following statistical constraints: 
1. Minimum number of taxa in a cluster (_cs_)
2. Multiple (_gamma_) of deviations from the grand median of the mean pairwise patristic distance that defines the within-cluster limit.
3. False discovery rate (_fdr_) for rejecting the null hypotheses that the pairwise patristic distance distributions of every combinatorial pair of clusters are empirically equivalent to that should they form a single cluster.

A manuscript describing PhyCLIP is available here:  
MANUSCRIPT LINK

## Installation
PhyCLIP is written in Python 2.7 (no support for Python 3 currently) and depends on several python libraries and at least one ILP solver.  

### Prerequisite: Python libraries    

* numpy, scipy, statsmodels  (mathematical/statistical operations)
* ete3 (parsing phylogenetic trees) 

To install the dependencies: 
```
PyPi
$pip install numpy scipy ete3 statsmodels

Conda/Anaconda
$conda install numpy scipy ete3 statsmodels
```
### Prerequisite: ILP solver 
PhyCLIP currently supports two ILP solvers. You can choose either **_ONE_** to install depending on your access to these solvers: 

1. **GLPK** (GNU Linear Programming Kit, https://www.gnu.org/software/glpk/) is a free and open-source package intended for solving large-scale linear programming, mixed integer programming, and other related problems. PyGLPK (http://tfinley.net/software/pyglpk/) is the Python module used to encapsulates the functionality of GLPK. 
2. **Gurobi** optimizer (http://www.gurobi.com/) is a commercial linear and quadratic programming solver with free licenses available for academic users.

#### Gurobi
If you are a university user (i.e. you have internet access from a recognized academic domain, e.g. '.edu' addresss), we reccomend running PhyCLIP with the Gurobi optimizer. The easiest way to install Gurobi is via the Anaconda platform:  

1. If you are not using Anaconda for Python 2.7 or have **not** installed Anaconda in your operating system, visit http://www.gurobi.com/downloads/get-anaconda to download the appropriate installer (choose support for Python 2.7) and install the Anaconda platform to your system. 

2. Once Anaconda is installed, enter ```conda install gurobi``` in your command/terminal prompt to install the Gurobi package. 

3. You need to install a Gurobi licence next. Visit http://www.gurobi.com/registration/academic-license-reg to register for a free Gurobi account. Follow the instructions in the verification email from Gurobi to set your password and login to your Gurobi account via https://www.gurobi.com/login. 

4. You can now access https://user.gurobi.com/download/licenses/free-academic to request for a free academic license. To install the license, enter the ```grbgetkey``` command along with the license key stipulated in the License Detail page in your command/terminal prompt. Note that an active internet connection from a recognized academic domain (e.g. '.edu' addresss) is required. 

#### GLPK

### Install PhyCLIP 

Finally, install phyclip.py by: 
```
$cd phyclip-master/ 
$python setup.py install
```
You may need sudo privileges for system-wide installation. Otherwise, it is also possible to use phyclip.py locally by adding the phyclip_modules folder to your $PYTHONPATH.

## Usage 

### Input file

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
                        Collapse zero branch lengths of tree prior to running
                        PhyCLIP (default = 0).
  --equivalent_zero_length EQUIVALENT_ZERO_LENGTH
                        Maximum branch length to be rounded to zero if
                        --collapse_zero_branch_lengths is called (advanced
                        option, default = 1e-07).
  --gam_method {MAD,Qn}
                        Robust estimation of dispersion parameter (default =
                        Qn).
  --hypo_test {Kuiper,KS}
                        Hypothesis test to use for statistical differentiation
                        of distance distributions (default = Kuiper).
  --preQ {0,1}          Perform Benjamini-Hochberg corrections of p-values
                        BEFORE filtering nodes < minimum cluster size
                        (advanced option, default = 0).
  --subsume_sensitivity_induced_clusters {0,1}
                        Subsume cluster-size sensitivity-induced clades into
                        parent clade (advanced option, default = 1).
  --sensitivity_percentile SENSITIVITY_PERCENTILE
                        Percentile threshold of clade size distribution under
                        which a clade is considered to be sensitivity-induced
                        (advanced option, default = 25 percent).
  --ilp_verbose {0,1}   ILP solver verbose (default: 0)
```
