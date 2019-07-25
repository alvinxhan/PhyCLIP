from __future__ import division
import re
import numpy as np
import itertools
import cython
import random
import multiprocessing as mp
from scipy.misc import comb as nCr
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
cimport numpy as np

from ete3 import Tree
import time
import sys

from copy import deepcopy as dc

from phyclip_modulex.pyutilx import collapse_zero_branch_length
from phyilpx_stats import qn, hypotest

IF UNAME_SYSNAME == "Windows":
    ctypedef long long _int64
ELSE:
    ctypedef long _int64

cdef struct Node:
    _int64 parent
    _int64* children
    _int64 children_length
    float edge_distance

@cython.no_gc_clear
cdef class phyilpx_treeinfo:

    cdef Node* data
    cdef _int64 depth, total_nr_nodes
    cdef _int64 hytest_method
    cdef object leafname_to_leaf_nodeid, leaf_nodeid_to_leafname, internalnodes, np_buffer, node_to_leaves, node_to_pwdist
    #cdef object nodepair_w_ancrelation
    cdef str original_tree_string
    cdef str treefname
    #cdef np.ndarray node_to_pwdist
    cdef np.ndarray leaf_dist_to_node, nodepair_to_distance, nodepair_to_pval

    #cdef object ete_nodeid_to_node

    def __init__(self, str newick_tree, str treefname, outgroup, int collapse_zero_branch_length_binary, float eq_zero_branch_length):

        cdef _int64 node_id

        self.np_buffer = None # numpy memory buffer

        # read newick tree as ete tree
        tree = Tree(newick_tree)

        # changing tree root
        if outgroup != False:

            if outgroup == 'midpoint': # midpoint rooting
                root_node = tree.get_midpoint_outgroup()
            else:
                # search for outgroup by name
                try:
                    root_node = tree&"'{}'".format(outgroup)
                except:
                    raise Exception('\nGiven outgroup {} did not match any taxon node in tree.\n'.format(outgroup))

            tree.set_outgroup(root_node)
            print ('\nWARNING: Tree outgroup changed to {}.'.format(outgroup))

        tree.ladderize()  # ladderize tree

        # collapse zero branch length
        if collapse_zero_branch_length_binary == 1:
            tree, self.treefname = collapse_zero_branch_length(tree, eq_zero_branch_length, treefname)

            if tree == False:
                raise Exception('No branches were collapsed. Check the upper limit of zero-length branches on your tree and adjust accordingly using --equivalent_zero_length')
        else:
            self.treefname = treefname

        treesize = len(tree.get_descendants()) + 1 # size of tree (all nodes + leaves)
        self.data = <Node*> PyMem_Malloc(treesize*sizeof(Node)) # create data array and allocate memory
        if not self.data:
            raise MemoryError()

        self.leaf_nodeid_to_leafname = {}
        self.leafname_to_leaf_nodeid = {}
        self.internalnodes = []

        #self.ete_nodeid_to_node = {}
        cdef float dist
        for node_id, node in enumerate(tree.traverse(strategy='levelorder')):

            #self.ete_nodeid_to_node[node_id] = node
            node.add_feature('node_id', node_id)

            if node.is_leaf() :
                self.leaf_nodeid_to_leafname[node_id] = node.name
                self.leafname_to_leaf_nodeid[node.name] = node_id
            else:
                node.add_feature('name', node_id)
                self.internalnodes.append(node_id)

        self.total_nr_nodes = node_id + 1 # total number of nodes

        # get original tree string
        self.original_tree_string = tree.write(format=3)

        for node_id, node in enumerate(tree.traverse(strategy='levelorder')):
            try:
                self.data[node_id].parent = node.up.node_id
            except:
                self.data[node_id].parent = -1 # root

            dist = node.dist
            if dist <= eq_zero_branch_length:
                self.data[node_id].edge_distance = 0.
            else:
                self.data[node_id].edge_distance = dist

            if node.is_leaf():
                children_list = []
            else:
                children_list = [child_node.node_id for child_node in node.get_children()]

            self.data[node_id].children_length = len(children_list)

            self.data[node_id].children = <_int64*> PyMem_Malloc(len(children_list) * sizeof(_int64))
            if not self.data[node_id].children:
                raise MemoryError()

            for c, child in enumerate(children_list):
                self.data[node_id].children[c] = child

        # determine the maximum depth of the tree
        for node_id in self.leaf_nodeid_to_leafname.keys():
            n = 1
            while True :
                if self.data[node_id].parent == -1: # break when parent is root
                    break
                node_id = self.data[node_id].parent
                n += 1
            if n > self.depth :
                self.depth = n

    def get_children(self, node_id):
        # return node ids of child nodes for a given node
        return [self.data[node_id].children[c] for c in xrange(self.data[node_id].children_length)]

    def get_leaves(self, node_id):
        # return array of leaves subtended by a given node
        cdef _int64 i, n = 0

        if self.np_buffer is None:
            self.np_buffer = np.ndarray(len(self.leaf_nodeid_to_leafname), dtype=np.int64)

        to_visit = [node_id]
        for i in to_visit:
            children_list = self.get_children(i)
            if len(children_list) == 0: # is leaf
                self.np_buffer[n] = i
                n += 1
            else:
                to_visit += children_list
        return np.array(self.np_buffer[:n])

    def mrca(self, a, b) :
        # return mrca of 2 nodes
        visited = np.zeros(self.depth, dtype=np.int64)

        return self._mrca(visited, a, b)

    @cython.boundscheck(False)
    cdef _int64 _mrca(self, _int64[:] visited, _int64 a, _int64 b) nogil :
        cdef _int64 i, a_depth
        cdef _int64 n, mrca_node = -1

        n = a
        i = 0
        while True :
            visited[i] = n
            n = self.data[n].parent
            i += 1
            if n == -1 : break
        a_depth = i

        n = b
        while True :
            i = 0
            while True :
                if i >= a_depth : break
                if visited[i] == n :
                    mrca_node = visited[i]
                    break
                i += 1
            if mrca_node != -1 : break
            n = self.data[n].parent
            if n == -1 :
                mrca_node = n
                break

        for i in xrange(self.depth) :
            visited[i] = -1

        return mrca_node

    @cython.boundscheck(False)
    cdef float get_distance(self, _int64 a, _int64 b):
        cdef _int64 mrca_node
        cdef float d = 0.
        cdef _int64 n

        mrca_node = self.mrca(a, b)

        n = a
        while n != mrca_node :
            d += self.data[n].edge_distance
            n =  self.data[n].parent
        n = b
        while n != mrca_node :
            d += self.data[n].edge_distance
            n =  self.data[n].parent
        return d

    def is_ancestor(self, a, b):
        # returns 1 if a is ancestor of b, -1 if b is ancestor of a, 0 if otherwise 
        return self._is_ancestor( a, b )

    @cython.boundscheck(False)
    cdef _int64 _is_ancestor(self, _int64 a, _int64 b) nogil :
        cdef _int64 i
        cdef _int64 n

        # is a an ancestor of b?
        i = b
        while True :
            n = self.data[i].parent
            if n == -1 :
                break
            if n == a :
                return 1
            i = n

        # is b an ancestor of a?
        i = a
        while True :
            n = self.data[i].parent
            if n == -1 :
                break
            if n == b :
                return -1
            i = n

        # or neither?
        return 0

    def properties(self):
        return self.leaf_nodeid_to_leafname, self.internalnodes, self.original_tree_string, self.treefname

    def get_nodepair_distance(self):
        cdef _int64 i, j
        cdef _int64 N = self.total_nr_nodes # length of all nodes
        cdef float dist

        self.nodepair_to_distance = np.zeros((N,N), dtype=np.float32)

        cdef _int64 N_leafpairs = 2*nCr(len(self.leaf_nodeid_to_leafname), 2)

        # pairwise patristic distance array of all nodes present
        for (i, j) in itertools.combinations(range(N), 2):
            dist = self.get_distance(i, j)
            self.nodepair_to_distance[(i,j)] = dist
            self.nodepair_to_distance[(j,i)] = dist

        return self.nodepair_to_distance

    def get_leaf_dist_to_node(self):
        cdef _int64 i, node_id
        cdef _int64 k = 0
        #cdef _int64 N = 0
        cdef _int64 N = self.total_nr_nodes
        cdef float dist

        cdef object leaf_to_ancestors

        # get leaf to internal node distances
        leaf_to_ancestors = {}
        for i, node_id in itertools.product(self.leaf_nodeid_to_leafname.keys(), self.internalnodes):
            if self.is_ancestor(node_id, i) == 1:
                try:
                    leaf_to_ancestors[i].append(node_id)
                except:
                    leaf_to_ancestors[i] = [node_id]
                #N += 1

        # structured array of node-leaf distance
        #self.leaf_dist_to_node = np.zeros(N, dtype={'names':('leaf', 'node', 'dist'), 'formats':(np.int64, np.int64, np.float32)})
        self.leaf_dist_to_node = np.full((N, N), np.nan, dtype = np.float32) # leaf v node
        for i in leaf_to_ancestors.keys():
            for node_id in leaf_to_ancestors[i]:
                dist = self.nodepair_to_distance[(i, node_id)]
                #self.leaf_dist_to_node[k] = (i, node_id, dist)
                self.leaf_dist_to_node[(i, node_id)] = dist
                k += 1

        # dictionary of leaf arrays reverse sorted by distance to node
        self.node_to_leaves = {}
        cdef np.ndarray curr_struc_array, leafdist_to_node_arr, leaves_of_node
        for node_id in self.internalnodes:
            #curr_struc_array = self.leaf_dist_to_node[self.leaf_dist_to_node['node'] == node_id][['leaf', 'dist']]
            #self.node_to_leaves[node_id] = np.sort(curr_struc_array, order='dist')['leaf'][::-1]

            curr_struc_array = self.leaf_dist_to_node[:,node_id]
            leafdist_to_node_arr = curr_struc_array[np.isnan(curr_struc_array) == False]
            leaves_of_node = np.where(np.isnan(curr_struc_array) == False)[0]
            self.node_to_leaves[node_id] = leaves_of_node[leafdist_to_node_arr.argsort()][::-1]

        return self.leaf_dist_to_node, self.node_to_leaves

    def get_node_to_parent_node(self):
        cdef _int64 node_id, parent, i
        cdef _int64 N = self.total_nr_nodes - 1
        #cdef np.ndarray node_to_parent_node = np.zeros(N, dtype={'names':('node', 'parent'), 'formats':(np.int64, np.int64)})
        cdef object node_to_parent_node = {}

        for i, node_id in enumerate(range(1, self.total_nr_nodes)):
            parent = self.data[node_id].parent
            #node_to_parent_node[i] = (node_id, parent)
            node_to_parent_node[node_id] = parent

        return node_to_parent_node

    def get_node_to_mean_child_dist2root(self):

        cdef _int64 node_id, child_node_id, i, k
        cdef float dist

        cdef object children_of_node

        cdef _int64 N = len(self.internalnodes)
        cdef np.ndarray child_dist_to_node
        cdef np.ndarray node_to_mean_child_dist2root = np.zeros(N, dtype={'names':('node', 'dist'), 'formats':(np.int64, np.float32)})

        for k, node_id in enumerate(self.internalnodes):
            children_of_node = self.get_children(node_id)

            child_dist_to_node = np.zeros(len(children_of_node), dtype=np.float32)
            for i, child_node_id in enumerate(children_of_node):

                dist = self.nodepair_to_distance[(child_node_id, 0)]
                child_dist_to_node[i] = dist

            node_to_mean_child_dist2root[k] = (node_id, np.mean(child_dist_to_node))

        return node_to_mean_child_dist2root

    def get_ancestral_relations(self):

        cdef object entry
        cdef _int64 i, j, k, node_id
        cdef _int64 N = len(self.internalnodes) + len(self.leaf_nodeid_to_leafname)
        cdef float mean_dist_to_root

        cdef object leaf_to_ancestors, node_to_descendant_nodes, ancestors_to_n, leafpairs_to_node, node_to_mean_pwdist

        cdef np.ndarray leaves_to_node
        cdef np.ndarray node_to_mean_child_dist2anc = np.zeros((N,N), dtype=np.float32)

        node_to_mean_child_dist2root = self.get_node_to_mean_child_dist2root() # get mean child node distances to root for every internal node

        leaf_to_ancestors = {}
        self.node_to_pwdist = {}
        node_to_mean_pwdist = {}
        node_to_descendant_nodes = {}

        for entry in np.sort(node_to_mean_child_dist2root, order='dist'): # starting node with the lowest mean children distance to root
            node_id, mean_dist_to_root = entry

            leaves_to_node = self.node_to_leaves[node_id]
            for i in leaves_to_node:
                try:
                    leaf_to_ancestors[i].append(node_id)
                except:
                    leaf_to_ancestors[i] = [node_id]

            ancestors_to_n = [i for i in self.internalnodes[::-1] if self.is_ancestor(i, node_id) == 1]

            for i in ancestors_to_n:
                try:
                    node_to_descendant_nodes[i].append(node_id)
                except:
                    node_to_descendant_nodes[i] = [node_id]

                node_to_mean_child_dist2anc[(node_id, i)] = mean_dist_to_root - self.nodepair_to_distance[i][0] # indexed by node to anc

            leafpairs_to_node = list(itertools.combinations(leaves_to_node, 2))
            pwdist_to_node = np.zeros(len(leafpairs_to_node), dtype=np.float32)
            for k, (i,j) in enumerate(leafpairs_to_node):
                pwdist_to_node[k] = self.nodepair_to_distance[(i,j)]

            self.node_to_pwdist[node_id] = np.sort(pwdist_to_node)
            node_to_mean_pwdist[node_id] = np.mean(pwdist_to_node)

        return node_to_descendant_nodes, leaf_to_ancestors, node_to_mean_child_dist2anc, node_to_mean_pwdist

    def _get_pval(self, main_node_index):

        cdef _int64 _, j, x, y, k
        cdef _int64 main_node = self.internalnodes[main_node_index]
        cdef _int64 n = len(self.internalnodes)
        cdef float pval0, pval1
        cdef np.ndarray ij_product_pwdist, ij_pwdist
        cdef np.ndarray pval_for_main_node = np.full(n, np.nan, dtype=np.float32)
        cdef object leaves_product

        for _ in xrange(main_node_index+1, n): # self.internalnodes should be sorted already
            j = self.internalnodes[_]

            if self.is_ancestor(main_node, j) == 0: # neither i or j are ancestors of each other's

                leaves_product = list(itertools.product(self.node_to_leaves[main_node], self.node_to_leaves[j]))

                ij_product_pwdist = np.zeros(len(leaves_product), dtype=np.float32)
                for k, (x,y) in enumerate(leaves_product):
                    ij_product_pwdist[k] = self.nodepair_to_distance[(x,y)]
                ij_pwdist = np.concatenate((self.node_to_pwdist[main_node], self.node_to_pwdist[j], ij_product_pwdist))
                ij_pwdist = np.sort(ij_pwdist)

                pval0 = hypotest(self.node_to_pwdist[main_node], ij_pwdist, self.hytest_method)
                pval1 = hypotest(self.node_to_pwdist[j], ij_pwdist, self.hytest_method)
                if pval0 > pval1:
                    pval_for_main_node[_] = pval0
                else:
                    pval_for_main_node[_] = pval1
            else:
                pval_for_main_node[_] = hypotest(self.node_to_pwdist[main_node], self.node_to_pwdist[j], self.hytest_method)

        return pval_for_main_node

    def get_global_pval(self, hytest_method, treeinfo_file, no_treeinfo_binary):

        cdef _int64 N = len(self.internalnodes) # number of indices to analyse
        cdef _int64 N_intnodes = N # number of internal nodes
        cdef _int64 i
        cdef object indices = []

        cdef _int64 errstatus = 2
        cdef _int64 progressN

        try:
            with np.load(treeinfo_file) as npzfhandle:
                self.nodepair_to_pval = npzfhandle['x']

            errstatus = 1 # treeinfo file parsed

            # check that the correct treeinfo file is given
            if self.nodepair_to_pval.shape[0] == self.nodepair_to_pval.shape[1] == N:
                errstatus = 0 # treeinfo file parsed and correct
            else:
                raise Exception() # invalid treeinfo file

            # check that all rows are NOT filled with N (except for last row)
            for i in range(N-1):
                if len(set(np.isnan(self.nodepair_to_pval[i,:]))) == 1: # all elements are nan
                    indices.append(i)
        except:
            pass

        if errstatus == 0:  # treeinfo file given
            N = len(indices)
            if N == 0:
                print ("Treeinfo file given...OK.")
                return self.nodepair_to_pval # treeinfo file parsed and all p-values found
            else:
                print ("Treeinfo file given...")
                progressN = sum([N_intnodes-i for i in indices])

        elif errstatus == 1: # treeinfo file given
            raise Exception('\nGiven treeinfo file does not match with input tree.\n') # but incorrectly matched with tree

        elif errstatus == 2: # treeinfo file NOT given
            indices = range(N)
            progressN = sum(range(N+1))
            self.nodepair_to_pval = np.full((N, N), np.nan, dtype=np.float32)

        cdef _int64 _i
        cdef _int64 writebreak = np.ceil(N/20) # 20 updates to treeinfo file

        cdef _int64 toolbar_width = 8
        cdef _int64 progress = 0

        self.hytest_method = hytest_method

        # setup toolbar
        if errstatus == 2:
            sys.stdout.write("Calculating p-values... %s" % (" " * toolbar_width))
        else:
            sys.stdout.write("Resuming calculations of p-values... %s" % (" " * toolbar_width))
        sys.stdout.flush()
        sys.stdout.write("\b" * (toolbar_width+1)) # return to start of line,

        random.shuffle(indices)

        for _i, i in enumerate(indices):
            self.nodepair_to_pval[i,:] = self._get_pval(i)

            # update the bar
            progress += (N_intnodes-i)
            sys.stdout.write("{:<8.3%}".format((progress)/progressN))
            sys.stdout.flush()
            sys.stdout.write("\b" * (toolbar_width)) # return to start of line,

            if no_treeinfo_binary == True: # do not print treeinfo
                continue
            else:
                if (_i+1)%writebreak == 0 or _i == N-1:
                    np.savez(treeinfo_file, x=self.nodepair_to_pval) # save over every write break

        sys.stdout.write("\n")

        return self.nodepair_to_pval

    def multiple_testing_correction(self, *args):

        import statsmodels.api as sm

        cdef np.ndarray list_of_ancestral_nodes, indices_of_curr_nodes

        if len(args) > 0:
            list_of_ancestral_nodes = args[0]
            indices_of_curr_nodes = return_indices_of_y_in_x(np.array(self.internalnodes, dtype = np.int64), list_of_ancestral_nodes)
        else:
            indices_of_curr_nodes = np.arange(len(self.internalnodes), dtype = np.int64)

        cdef _int64 i, _i, j
        cdef np.ndarray pval_array = np.array([self.nodepair_to_pval[(i,j)] for _i, i in enumerate(indices_of_curr_nodes) for j in indices_of_curr_nodes[_i+1:]], dtype = np.float32)
        cdef np.ndarray qval_values = sm.stats.multipletests(pval_array, method='fdr_bh')[1]

        cdef k = 0
        cdef float qval
        cdef object qval_dict = {}
        for _i, i in enumerate(indices_of_curr_nodes):
            for j in indices_of_curr_nodes[_i+1:]:
                qval = qval_values[k]
                qval_dict[(self.internalnodes[i], self.internalnodes[j])] = qval
                k += 1

        return qval_dict

    def __dealloc__(self):
        PyMem_Free(self.data)  # no-op if self.data is NULL

# get indices of elements in y array in x array
@cython.boundscheck(False)
cdef np.ndarray return_indices_of_y_in_x(np.ndarray x, np.ndarray y):

    cdef np.ndarray index = np.argsort(x)
    cdef np.ndarray sorted_x = x[index]
    cdef np.ndarray sorted_index = np.searchsorted(sorted_x, y)

    cdef np.ndarray yindex = np.take(index, sorted_index, mode="clip")
    cdef np.ndarray mask = x[yindex] != y

    return np.ma.array(yindex, mask=mask)

@cython.no_gc_clear
cdef class distal_dissociation:
    cdef _int64 min_cluster_size
    cdef float within_cluster_limit
    cdef str gam_method
    cdef np.ndarray nodes_list, nodepair_to_distance, leaf_dist_to_node
    cdef object node_to_parent_node, node_to_leaves, node_to_descendant_nodes, node_to_mean_pwdist, node_to_mean_child_dist2anc, leaf_to_ancestors, leaf_node_id_to_leafname, node_to_children_node

    def __init__(self, _int64 min_cluster_size, float within_cluster_limit, str gam_method, np.ndarray nodes_list, object node_to_leaves, object node_to_descendant_nodes,
                 object node_to_mean_pwdist, object node_to_mean_child_dist2anc, object node_to_parent_node, np.ndarray nodepair_to_distance, np.ndarray leaf_dist_to_node,
                 object leaf_to_ancestors, object leaf_node_id_to_leafname):

        self.min_cluster_size = min_cluster_size
        self.within_cluster_limit = within_cluster_limit
        self.gam_method = gam_method
        self.node_to_mean_child_dist2anc = node_to_mean_child_dist2anc
        self.node_to_parent_node = node_to_parent_node
        self.nodepair_to_distance = nodepair_to_distance
        self.leaf_dist_to_node = leaf_dist_to_node
        self.leaf_to_ancestors = leaf_to_ancestors
        self.leaf_node_id_to_leafname = leaf_node_id_to_leafname
        self.node_to_leaves = node_to_leaves
        self.node_to_descendant_nodes = node_to_descendant_nodes
        self.node_to_mean_pwdist = node_to_mean_pwdist
        self.nodes_list = nodes_list

        cdef _int64 node, parent
        self.node_to_children_node = {}
        for node, parent in self.node_to_parent_node.items():
            try:
                self.node_to_children_node[parent].append(node)
            except:
                self.node_to_children_node[parent] = [node]

    cdef float get_pwdist_from_leaf_distances_to_node(self, object leaves_dist_to_node, object desc_node_to_leaves):
        cdef _int64 n_i = len(leaves_dist_to_node)

        cdef float term_a = (n_i-1)*sum(leaves_dist_to_node)
        cdef float term_b = 0.

        cdef _int64 desc_node, n_desc, parent_node
        cdef object desc_node_leaves

        for desc_node, desc_node_leaves in desc_node_to_leaves.items():
            n_desc = len(desc_node_leaves)
            parent_node = self.node_to_parent_node[desc_node]
            term_b += (n_desc)*(n_desc-1)*self.nodepair_to_distance[(parent_node, desc_node)]

        return 2*(term_a-term_b)/(n_i*(n_i-1))

    def remove_outlier_leaves(self, _int64 main_node, *args):
        # main_node, rol_leaves, rol_node_to_descendant_nodes[main_node], rol_node_to_leaves
        # node, self.node_to_leaves[node], self.node_to_descendant_nodes, self.node_to_leaves

        cdef _int64 leaf
        cdef np.ndarray rol_leaves
        cdef object rol_node_to_leaves, descendant_nodes_to_main_node

        if len(args) > 0:
            rol_leaves = args[0]
            descendant_nodes_to_main_node = args[1]
            rol_node_to_leaves = args[2]
        else:
            rol_leaves = self.node_to_leaves[main_node]
            try:
                descendant_nodes_to_main_node = self.node_to_descendant_nodes[main_node]
            except:
                descendant_nodes_to_main_node = []
            rol_node_to_leaves = self.node_to_leaves

        cdef object rol_leaf_to_dist = {}
        cdef float dist, reduced_mean_pwdist
        cdef _int64 rdn, x

        for leaf in rol_leaves:
            #dist = self.leaf_dist_to_node[(self.leaf_dist_to_node['leaf'] == leaf) & (self.leaf_dist_to_node['node'] == main_node)]['dist'].sum()
            dist = self.leaf_dist_to_node[(leaf, main_node)]
            rol_leaf_to_dist[leaf] = dist

        cdef float med_leaf_dist_to_node = np.median(rol_leaf_to_dist.values())
        cdef float mad_leaf_dist_to_node = qn(rol_leaf_to_dist.values())

        cdef np.ndarray rol_remaining_leaves = np.array([leaf for leaf, dist in rol_leaf_to_dist.items() if dist <= (med_leaf_dist_to_node + (3*mad_leaf_dist_to_node))], dtype = np.int64)

        """
        cdef _int64 n_leaves = len(rol_leaves)
        cdef object leafpair_list = list(itertools.combinations(range(n_leaves), 2))
        cdef np.ndarray _pwdist = np.full((n_leaves, n_leaves), np.nan, dtype=np.float32)
        cdef np.ndarray cluster_pwdist = np.zeros(len(leafpair_list), dtype = np.float32)
        cdef _int64 _, _i, _j, i, j

        for _, (_i, _j) in enumerate(leafpair_list):
            i = rol_leaves[_i]
            j = rol_leaves[_j]
            _pwdist[(_i, _j)] = _pwdist[(_j, _i)] = self.nodepair_to_distance[(i, j)]
            cluster_pwdist[_] = self.nodepair_to_distance[(i, j)]

        cdef np.ndarray mean_d = np.zeros(n_leaves, dtype = np.float32)
        cdef np.ndarray mean_d_i
        for i in xrange(n_leaves):
            mean_d_i = _pwdist[(i,)]
            mean_d[i] = np.mean(mean_d_i[~np.isnan(mean_d_i)])

            leaf = rol_leaves[i]
            dist = self.leaf_dist_to_node[(leaf, main_node)]
            rol_leaf_to_dist[leaf] = dist

        cdef float med_leaf_dist_to_node = np.median(mean_d)
        cdef float mad_leaf_dist_to_node = qn(mean_d)
        cdef np.ndarray rol_remaining_leaves = np.array([rol_leaves[i] for i in range(n_leaves) if mean_d[i] <= (med_leaf_dist_to_node + (3*mad_leaf_dist_to_node))], dtype = np.int64)
        """

        cdef object rol_leaves_to_dissociate = list(set(rol_node_to_leaves[main_node]) - set(rol_remaining_leaves))
        cdef object remaining_descendant_nodes_to_leaves, desc_nodes_subtending_leaf, nodes_of_dissociated_leaves, rol_descendant_nodes_to_dissociate, y

        if len(rol_leaves_to_dissociate) == 0 or len(rol_remaining_leaves) == 1:
            return False
        else:
            remaining_descendant_nodes_to_leaves = {}

            for leaf in rol_remaining_leaves:
                desc_nodes_subtending_leaf = list(set(descendant_nodes_to_main_node)&set(self.leaf_to_ancestors[leaf]))

                for rdn in desc_nodes_subtending_leaf:
                    try:
                        remaining_descendant_nodes_to_leaves[rdn].append(leaf)
                    except:
                        remaining_descendant_nodes_to_leaves[rdn] = [leaf]

            reduced_mean_pwdist = self.get_pwdist_from_leaf_distances_to_node([rol_leaf_to_dist[leaf] for leaf in rol_remaining_leaves], remaining_descendant_nodes_to_leaves)

            nodes_of_dissociated_leaves = list(set(descendant_nodes_to_main_node)&set([x for y in [self.leaf_to_ancestors[leaf] for leaf in rol_leaves_to_dissociate] for x in y]))
            rol_descendant_nodes_to_dissociate = [x for x in nodes_of_dissociated_leaves if set(rol_node_to_leaves[x]) <= set(rol_leaves_to_dissociate)]

        return rol_remaining_leaves, rol_descendant_nodes_to_dissociate, reduced_mean_pwdist

    def leave_one_out_leaf_reduction(self, _int64 main_node, *args):
        # input arg = node, sorted_leaves
        cdef np.ndarray sorted_leaves # note that sorted_leaves is already reverse sorted by distance to main_node

        if len(args) > 0:
            sorted_leaves = args[0]
        else:
            sorted_leaves = self.node_to_leaves[main_node]

        # immediate return if len(sorted_leaves) < self.min_cluster_size
        if len(sorted_leaves) < self.min_cluster_size:
            return False

        cdef _int64 leaf, rdn, x
        cdef _int64 l_index
        cdef float dist, reduced_mean_pwdist
        cdef object remaining_leaves_to_node_dist, remaining_descendant_nodes_to_leaves, desc_nodes_subtending_leaf, leaves_dissociated, old_descendant_nodes_to_dissociate, y, nodes_of_dissociated_leaves

        for l_index in range(-1, len(sorted_leaves), 1):
            # return if number of leaves left < self.min_cluster_size
            if l_index + 1 == len(sorted_leaves) - self.min_cluster_size:
                return False

            remaining_leaves_to_node_dist = {}
            remaining_descendant_nodes_to_leaves = {}

            # start from not having any leaves removed
            for leaf in sorted_leaves[l_index+1:]:
                #dist = self.leaf_dist_to_node[(self.leaf_dist_to_node['leaf'] == leaf) & (self.leaf_dist_to_node['node'] == main_node)]['dist'].sum()
                dist = self.leaf_dist_to_node[(leaf, main_node)]
                remaining_leaves_to_node_dist[leaf] = dist

                try:
                    desc_nodes_subtending_leaf = list(set(self.node_to_descendant_nodes[main_node])&set(self.leaf_to_ancestors[leaf]))
                except:
                    desc_nodes_subtending_leaf = []

                for rdn in desc_nodes_subtending_leaf:
                    try:
                        remaining_descendant_nodes_to_leaves[rdn].append(leaf)
                    except:
                        remaining_descendant_nodes_to_leaves[rdn] = [leaf]

            reduced_mean_pwdist = self.get_pwdist_from_leaf_distances_to_node(remaining_leaves_to_node_dist.values(), remaining_descendant_nodes_to_leaves)

            # break if <= self.within_cluster_limit and dissociate
            if reduced_mean_pwdist <= self.within_cluster_limit:
                # find nodes of dissociated leaves
                leaves_dissociated = sorted_leaves[:l_index+1]

                if len(leaves_dissociated) == 0:
                    old_descendant_nodes_to_dissociate = []
                else:
                    try:
                        nodes_of_dissociated_leaves = list(set([x for y in [list(set(self.node_to_descendant_nodes[main_node])&set(self.leaf_to_ancestors[leaf])) for leaf in leaves_dissociated] for x in y]))
                        old_descendant_nodes_to_dissociate = [desc_node for desc_node in list(set(nodes_of_dissociated_leaves) - set(remaining_descendant_nodes_to_leaves.keys())) if set(self.node_to_leaves[desc_node]) <= set(leaves_dissociated)]
                    except:
                        old_descendant_nodes_to_dissociate = []

                # return leaves to keep for node AND descendant nodes (which all of its subtended leaves are to be removed) to dissociate
                return np.array(sorted_leaves[l_index+1:], dtype = np.int64), old_descendant_nodes_to_dissociate, reduced_mean_pwdist

        return False

    # distal dissociation main
    def dd_main(self):

        print ('Distal dissociation...')

        cdef object sd_node_to_leaves, sd_node_to_descendant_nodes, sd_node_to_mean_pwdist, rol_output, descendant_nodes_to_dissociate, loo_output, leaves_to_remove, y, loo_descendant_nodes_to_dissociate, desc_nodes_of_node
        cdef _int64 node, desc_node, x, leaf
        cdef float reduced_mean_pwdist, mean_pwdist
        cdef np.ndarray remaining_leaves, leaves_to_keep

        # create deep copies to be edited
        sd_node_to_leaves = dc(self.node_to_leaves)
        sd_node_to_descendant_nodes = dc(self.node_to_descendant_nodes)
        sd_node_to_mean_pwdist = {} # mean pairwise distance dictionary for current set of parameters

        # reassociate subtrees - starting from the root by level order
        for node in self.nodes_list: # nodes_list already sorted

            # --- current node <= self.within_cluster_limit --- #
            if self.node_to_mean_pwdist[node] <= self.within_cluster_limit:
                # check there are no outlying leaves
                rol_output = self.remove_outlier_leaves(node)

                if rol_output == False:
                    sd_node_to_mean_pwdist[node] = self.node_to_mean_pwdist[node]
                else:
                    leaves_to_keep, descendant_nodes_to_dissociate, mean_pwdist = rol_output
                    # update node to remaining leaves
                    sd_node_to_leaves[node] = leaves_to_keep[:]
                    try:
                        # dissociate descendant nodes from node if any
                        sd_node_to_descendant_nodes[node] = list(set(sd_node_to_descendant_nodes[node])-set(descendant_nodes_to_dissociate))
                    except:
                        pass
                    # update mean pwd dist
                    sd_node_to_mean_pwdist[node] = mean_pwdist

                continue
            # --- current node <= self.within_cluster_limit --- #

            # --- current node > self.within_cluster_limit --- #
            descendant_nodes_to_dissociate = []  # list to save nodes for dissociation
            try:
                ## --- node has descendant internal nodes --- ##
                # descendants nodes are already sorted by mean child-nodes' distance to node
                for desc_node in self.node_to_descendant_nodes[node]:
                    if desc_node in descendant_nodes_to_dissociate:
                        continue

                    elif self.node_to_mean_child_dist2anc[(desc_node, node)] > self.within_cluster_limit:
                        # append not only desc_node but all descendant nodes of desc_node itself
                        descendant_nodes_to_dissociate.append(desc_node)
                        try:
                            descendant_nodes_to_dissociate = list(set(descendant_nodes_to_dissociate)|set(self.node_to_descendant_nodes[desc_node]))
                        except:
                            pass

                ## --- node has descendant internal nodes --- ##
            except:
                ## --- dead-end node with no descendants --- ##
                loo_output = self.leave_one_out_leaf_reduction(node)

                # no leaves to remove by leave-one-out-wcl approach
                if loo_output == False:
                    # check there are no outlying leaves by distance to node
                    rol_output = self.remove_outlier_leaves(node)
                    # no outlying leaves
                    if rol_output == False:
                        sd_node_to_mean_pwdist[node] = self.node_to_mean_pwdist[node]
                    # outlying leaves found
                    else:
                        leaves_to_keep, descendant_nodes_to_dissociate, mean_pwdist = rol_output
                        # update node to remaining leaves
                        sd_node_to_leaves[node] = leaves_to_keep[:]
                        # update mean pwd dist
                        sd_node_to_mean_pwdist[node] = mean_pwdist

                # leave-one-out-wcl approach removes some leaves from dead-end node
                else:
                    leaves_to_keep, descendant_nodes_to_dissociate, mean_pwdist = loo_output
                    # update node_to_leaves as per output of leave-one-out-wcl approach first
                    sd_node_to_leaves[node] = leaves_to_keep[:]

                    # check that there are no outlying leaves by distance to node
                    try:
                        desc_nodes_of_node = sd_node_to_descendant_nodes[node]
                    except:
                        desc_nodes_of_node = []
                    rol_output = self.remove_outlier_leaves(node, sd_node_to_leaves[node], desc_nodes_of_node, sd_node_to_leaves)

                    # no outlying leaves
                    if rol_output == False:
                        # update mean_pwdist based on leave-one-out-wcl approach
                        sd_node_to_mean_pwdist[node] = mean_pwdist
                    # outlying leaves to node found
                    else:
                        leaves_to_keep, descendant_nodes_to_dissociate, mean_pwdist = rol_output
                        # update node to remaining leaves
                        sd_node_to_leaves[node] = leaves_to_keep[:]
                        # update mean pwd dist
                        sd_node_to_mean_pwdist[node] = mean_pwdist

                continue
                ## --- dead-end node with no descendants --- ##

            ## --- node has descendant internal nodes --- ##
            leaves_to_remove = list(set([x for y in [self.node_to_leaves[desc_node] for desc_node in descendant_nodes_to_dissociate] for x in y]))
            # remove all leaves from nodes that could be potentially dissociated (leaves in self.node_to_leaves[node] already reverse-sorted by distance to node)
            remaining_leaves = np.array([leaf for leaf in self.node_to_leaves[node] if leaf not in leaves_to_remove], dtype = np.int64)

            loo_output = self.leave_one_out_leaf_reduction(node, remaining_leaves) # leave-one-out-wcl approach removes leaves

            # no outlying leaves
            if loo_output != False:
                # update mean_pwdist from leave-one-out-wcl output
                leaves_to_keep, loo_descendant_nodes_to_dissociate, mean_pwdist = loo_output

                # update as per output of leave-one-out-wcl approach
                # update node to remaining leaves
                sd_node_to_leaves[node] = leaves_to_keep[:]
                # dissociate descendant nodes from node
                descendant_nodes_to_dissociate = list(set(descendant_nodes_to_dissociate)|set(loo_descendant_nodes_to_dissociate))
                sd_node_to_descendant_nodes[node] = list(set(sd_node_to_descendant_nodes[node])-set(descendant_nodes_to_dissociate))

                # now check that there are no outlying leaves by distance to node
                try:
                    desc_nodes_of_node = sd_node_to_descendant_nodes[node]
                except:
                    desc_nodes_of_node = []
                rol_output = self.remove_outlier_leaves(node, sd_node_to_leaves[node], desc_nodes_of_node, sd_node_to_leaves)

                # no outlying leaves
                if rol_output == False:
                    # update mean_pwdist from leave-one-out-wcl output
                    sd_node_to_mean_pwdist[node] = mean_pwdist
                # outlying leaves by distance to node found
                else:
                    leaves_to_keep, descendant_nodes_to_dissociate, mean_pwdist = rol_output
                    # update node to remaining leaves
                    sd_node_to_leaves[node] = leaves_to_keep[:]
                    # dissociate descendant nodes from node
                    sd_node_to_descendant_nodes[node] = list(set(sd_node_to_descendant_nodes[node])-set(descendant_nodes_to_dissociate))
                    # update mean pwd dist
                    sd_node_to_mean_pwdist[node] = mean_pwdist

            # leave-one-out-wcl approach did not remove any leaves
            else:
                # perform leave-one-out-wcl approach again, now based on self.node_to_leaves[node]
                loo_output = self.leave_one_out_leaf_reduction(node)

                # leave-one-out-wcl approach still did not remove any leaves
                if loo_output == False:
                    # check for outlying leaves by distance to node
                    rol_output = self.remove_outlier_leaves(node)
                    # no outlying leaves
                    if rol_output == False:
                        sd_node_to_mean_pwdist[node] = self.node_to_mean_pwdist[node]
                    # outlying leaves found
                    else:
                        leaves_to_keep, descendant_nodes_to_dissociate, mean_pwdist = rol_output
                        # update node to remaining leaves
                        sd_node_to_leaves[node] = leaves_to_keep[:]
                        # dissociate descendant nodes from node
                        sd_node_to_descendant_nodes[node] = list(set(sd_node_to_descendant_nodes[node])-set(descendant_nodes_to_dissociate))
                        # update mean pwd dist
                        sd_node_to_mean_pwdist[node] = mean_pwdist

                # leave-one-out-wcl approach removes leaves
                else:
                    leaves_to_keep, descendant_nodes_to_dissociate, mean_pwdist = loo_output
                    # update as per output of leave-one-out-wcl approach
                    # update node to remaining leaves
                    sd_node_to_leaves[node] = leaves_to_keep[:]
                    # dissociate descendant nodes from node
                    sd_node_to_descendant_nodes[node] = list(set(sd_node_to_descendant_nodes[node])-set(descendant_nodes_to_dissociate))

                    # check that there are no outlying leaves by distance to node
                    try:
                        desc_nodes_of_node = sd_node_to_descendant_nodes[node]
                    except:
                        desc_nodes_of_node = []
                    rol_output = self.remove_outlier_leaves(node, sd_node_to_leaves[node], desc_nodes_of_node, sd_node_to_leaves)

                    # no outlying leaves
                    if rol_output == False:
                        # this mean_pwdist is from loo_output
                        sd_node_to_mean_pwdist[node] = mean_pwdist
                    # outlying leaves found
                    else:
                        leaves_to_keep, descendant_nodes_to_dissociate, mean_pwdist = rol_output
                        # update node to remaining leaves
                        sd_node_to_leaves[node] = leaves_to_keep[:]
                        # dissociate descendant nodes from node
                        sd_node_to_descendant_nodes[node] = list(set(sd_node_to_descendant_nodes[node])-set(descendant_nodes_to_dissociate))
                        # update mean pwd dist
                        sd_node_to_mean_pwdist[node] = mean_pwdist

        """
        # check that no two nodes are subtending the same set of leaves
        cdef _int64 i, j
        cdef object leaf_tuple, list_of_nodes_to_del, leaf_tuple_to_nodes_to_del = {}
        for i, j in itertools.combinations(sd_node_to_mean_pwdist.keys(), 2):
            if set(sd_node_to_leaves[i]) == set(sd_node_to_leaves[j]):
                leaf_tuple = tuple(np.sort(sd_node_to_leaves[i]))

                if i < j:
                    try:
                        leaf_tuple_to_nodes_to_del[leaf_tuple].append(i)
                    except:
                        leaf_tuple_to_nodes_to_del[leaf_tuple] = [i]
                else:
                    try:
                        leaf_tuple_to_nodes_to_del[leaf_tuple].append(j)
                    except:
                        leaf_tuple_to_nodes_to_del[leaf_tuple] = [j]

        for leaf_tuple, list_of_nodes_to_del in leaf_tuple_to_nodes_to_del.items():
            for i in set(list_of_nodes_to_del):
                del sd_node_to_leaves[i]
                try:
                    del sd_node_to_descendant_nodes[i]
                except:
                    pass
                del sd_node_to_mean_pwdist[i]
        """
        """
        # ! --- addition --- #
        # get nodes <= curr_wcl
        cdef np.ndarray nodes_below_wcl = np.sort(np.array([node for node in sd_node_to_mean_pwdist.keys() if sd_node_to_mean_pwdist[node] <= self.within_cluster_limit], dtype = np.int64))
        cdef float med_pwdist, mad_pwdist, max_pwdist, dist_i, dist_j, upper_tol = -1.

        cdef object leaf_pairs, curr_nodes_and_leaves, children_of_node, nodes_to_separate
        cdef _int64 _, i, j, child
        cdef np.ndarray _pwdist, leaves


        for node in nodes_below_wcl:
            leaf_pairs = list(itertools.combinations(sd_node_to_leaves[node], 2))
            _pwdist = np.zeros(len(leaf_pairs), dtype = np.float32)

            for _, (i,j) in enumerate(leaf_pairs):
                _pwdist[_] = self.nodepair_to_distance[(i,j)]

            med_pwdist = np.median(_pwdist)
            mad_pwdist = qn(_pwdist)

            max_pwdist = med_pwdist + 3*mad_pwdist
            if max_pwdist > upper_tol:
                upper_tol = max_pwdist

        #print upper_tol, self.within_cluster_limit
        # homogeneity check
        for node, leaves in sd_node_to_leaves.items():
            try:
                curr_nodes_and_leaves = set(sd_node_to_descendant_nodes[node])|set(leaves)
            except:
                curr_nodes_and_leaves = set(leaves)

            children_of_node = list(set(self.node_to_children_node[node])&(curr_nodes_and_leaves))
            if len(children_of_node) < 2:
                continue

            nodes_to_separate = []
            for i, j in itertools.combinations(children_of_node, 2):
                if i in self.leaf_to_ancestors:
                    dist_i = self.nodepair_to_distance[(i, node)]
                else:
                    dist_i = self.node_to_mean_child_dist2anc[(i, node)]

                if j in self.leaf_to_ancestors:
                    dist_j = self.nodepair_to_distance[(j, node)]
                else:
                    dist_j = self.node_to_mean_child_dist2anc[(j, node)]

                if dist_i + dist_j > upper_tol:
                    #print node, i, j
                    nodes_to_separate.append(i)
                    nodes_to_separate.append(j)

            if len(nodes_to_separate) > 0:
                for child in set(nodes_to_separate):
                    if nodes_to_separate.count(child) == len(children_of_node) - 1:

                        if child in self.leaf_to_ancestors: # child is leaf
                            sd_node_to_leaves[node] = np.array(list(set(sd_node_to_leaves[node]) - set([child])), dtype = np.int64)
                        else:
                            try:
                                sd_node_to_leaves[node] = np.array(list(set(sd_node_to_leaves[node]) - set(sd_node_to_leaves[child])), dtype = np.int64)
                            except:
                                pass

                            try:
                                sd_node_to_descendant_nodes[node] = list(set(sd_node_to_descendant_nodes[node]) - set(sd_node_to_descendant_nodes[child]))
                            except:
                                pass

                if len(sd_node_to_leaves[node]) < self.min_cluster_size:
                    del sd_node_to_leaves[node]
                else:
                    sd_node_to_mean_pwdist[node] = np.mean([self.nodepair_to_distance[(i,j)] for i, j in itertools.combinations(sd_node_to_leaves[node], 2)])
            # ! --- addition --- #
            """

        return sd_node_to_leaves, sd_node_to_descendant_nodes, sd_node_to_mean_pwdist