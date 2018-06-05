# PhyCLIP (_Phylogenetic Clustering by Linear Integer Programming_)

PhyCLIP is a integer linear programming (ILP) approach that assigns statistically-principled cluster membership to as many taxa as possible for a given **_rooted_** phylogenetic tree based on its pairwise patristic distance distribution, subject to the following statistical constraints: 
1. Minimum number of taxa in a cluster (_cs_)
2. Multiple (_gamma_) of deviations from the grand median of the mean pairwise patristic distance that defines the within-cluster limit.
3. False discovery rate (_fdr_) for rejecting the null hypotheses that the pairwise patristic distance distributions of every combinatorial pair of clusters are empirically equivalent to that should they form a single cluster.

A manuscript describing PhyCLIP is available here:  
MANUSCRIPT LINK

## Installation
PhyCLIP is written in Python 2.7 and is currently distributed in two versions (phyclip-glpk.py and phyclip-gurobi.py) depending on the user's accessibility to the following supported ILP solvers: 
1. **GLPK** (GNU Linear Programming Kit, https://www.gnu.org/software/glpk/) is a free and open-source package intended for solving large-scale linear programming, mixed integer programming, and other related problems. PyGLPK (http://tfinley.net/software/pyglpk/) is the Python module used to encapsulates the functionality of GLPK. 
2. **Gurobi** optimizer (http://www.gurobi.com/) is a commercial linear and quadratic programming solver with free licenses available for academic users.
<!-- The GNU Linear Programmiong Kit (**GLPK**; ) is a . -->

### Prerequisites  
Several python libraries are required for **all** versions, including: 
* numpy, scipy, statsmodels  (mathematical/statistical operations)
* ete3 (parsing phylogenetic trees in NEWICK format) 

To install the dependencies: 
```
PyPi
$pip install numpy scipy ete3 statsmodels

Conda/Anaconda
$conda install numpy scipy ete3 statsmodels
```

### Gurobi version (phyclip-gurobi.py)
If you are a university user (i.e. you have internet access from a recognized academic domain, e.g. '.edu' addresss), we reccomend running PhyCLIP with the Gurobi optimizer. The easiest way to install Gurobi is via the Anaconda platform:  

1. If you are not using Anaconda for Python 2.7 or have **not** installed Anaconda in your operating system, visit http://www.gurobi.com/downloads/get-anaconda to download the appropriate installer (choose support for Python 2.7) and install the Anaconda platform to your system. 

2. Once Anaconda is installed, enter ```conda install gurobi``` in your command/terminal prompt to install the Gurobi package. 

3. You need to install a Gurobi licence next. Visit http://www.gurobi.com/registration/academic-license-reg to register for a free Gurobi account. Follow the instructions in the verification email from Gurobi to set your password and login to your Gurobi account via https://www.gurobi.com/login. 

4. You can now access https://user.gurobi.com/download/licenses/free-academic to request for a free academic license. To install the license, enter the ```grbgetkey``` command along with the license key stipulated in the License Detail page in your command/terminal prompt. Note that an active internet connection from a recognized academic domain (e.g. '.edu' addresss) is required. 

Finally, install phyclip-gurobi.py by: 
```
$cd phyclip-master/ 
$python setup.py install
```

### GLPK version (phyclip-glpk.py)
