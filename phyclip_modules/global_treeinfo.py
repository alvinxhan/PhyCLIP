from __future__ import division
import re
import itertools
import traceback
import math
import json
import numpy as np
from pathos.helpers import mp
from random import shuffle
import time

class inter_cluster_hytest(object):
    '''
    Hypothesis tests used to determine if the pairwise patristic distance distributions between two clusters are statistically distinct.

    Attributes:
         data1, data2 - SORTED lists/arrays of two distributions
    '''

    def __init__(self, data1=None, data2=None):
        self.data1 = data1
        self.data2 = data2

    def hytest(self, method='kuiper'):
        '''
        Performs hypothesis testing based on method called

        *Parameters*
        method (str): hypothesis test to use, either 'kuiper' or 'ks'.

        *Returns*
        prob (float): p-value
        '''

        if method == 'kuiper':
            prob = self._kuiper_2samp(self.data1, self.data2)
        else:
            prob = self._ks_2samp(self.data1, self.data2)

        return prob

    def _ks_2samp(self, data1, data2):
        from numpy import asarray
        from scipy.stats import kstwobign

        data1, data2 = map(asarray, (data1, data2)) # change data list to array

        n1 = len(data1) # number of elements in data1
        n2 = len(data2) # number of elements in data2

        #data1 = np.sort(data1) # sort data1 by ascending values
        #data2 = np.sort(data2) # sort data2 by ascending values

        data_all = np.concatenate([data1,data2]) # concatenate both data arrays and sort - this is basically the array of all possible values

        cdf1 = np.searchsorted(data1,data_all,side='right')/n1
        cdf2 = np.searchsorted(data2,data_all,side='right')/n2

        d = np.max(np.absolute(cdf1-cdf2))
        # Note: d absolute not signed distance
        en = np.sqrt(n1*n2/float(n1+n2))

        try:
            prob = kstwobign.sf((en + 0.12 + 0.11 / en) * d)
        except:
            prob = 1.0
        return prob

    def _kuiper_2samp(self, data1, data2):
        j1 = 0
        j2 = 0

        n1 = len(data1)  # number of elements in data1
        n2 = len(data2)  # number of elements in data2
        eff_n = np.sqrt((n1*n2)/(n1+n2))

        fn1 = 0
        fn2 = 0
        d = 0

        while j1 < n1 and j2 < n2:
            d1 = data1[j1]
            d2 = data2[j2]
            if d1 <= d2:
                j1 += 1
                fn1 = j1/n1

            if d2 <= d1:
                j2 += 1
                fn2 = j2/n2

            dt = abs(fn2-fn1)
            if dt > d:
                d = dt

        def kuiper_dist(q):
            EPS1 = 0.001
            EPS2 = 0.00000001
            a2 = -2.0*q*q
            sum = 0
            termbf = 0
            dist = 1

            for j in xrange(1, 100000, 1):
                term = 2.0*((4.0*j*j*q*q)-1)*np.exp(a2*j*j)
                sum += term

                if abs(term) <= EPS1*termbf or abs(term) <= EPS2*sum:
                    dist = sum
                    break

                termbf = abs(term)
            return dist

        lam = (eff_n + 0.155 + 0.24/eff_n)*d
        prob = kuiper_dist(lam)

        return prob

class get_global_tree_info(object):
    '''
    Get global tree information
    '''
    def __init__(self, tree_object=None, leaf_dist_to_node=None, leafpair_to_distance=None, nodepair_to_pval=None, treeinfo_fname=None, hypo_test_method=None, no_treeinfo=False, cores=False):
        self.tree_object = tree_object
        self.leaf_dist_to_node = leaf_dist_to_node
        self.leafpair_to_distance = leafpair_to_distance
        self.nodepair_to_pval = nodepair_to_pval
        self.treeinfo_fname = treeinfo_fname
        self.hypo_test_method = hypo_test_method
        self.no_treeinfo = no_treeinfo
        self.cores = cores

    #!-- multi-processing: get leaf distance to node for all leaves subtended  --#
    def get_leaf_distance_to_node(self, n_list, inode_list, ldn_queue, ntl_queue):

        sorted_leaves_to_node = {}
        currn_leaf_to_dist = {}

        for _, n_index in enumerate(n_list):
            inode = inode_list[_]

            # distance of leaf to each of its ancestral node
            for leaf_node in inode.get_leaves():
                leaf = leaf_node.name

                if self.treeinfo_file_given > 0:
                    try:
                        # check that we have correctly recovered leaf dist to node from treeinfo file (this will be the check that the correct treeinfo file is being used)
                        dist = self.leaf_dist_to_node[leaf][n_index]
                    except:
                        raise Exception('Mismatched leaf and node index found in TREEINFO file. Re-run without TREEINFO file.')
                else:
                    dist = leaf_node.get_distance(inode)

                try:
                    currn_leaf_to_dist[leaf].append((n_index, dist))
                except:
                    currn_leaf_to_dist[leaf] = [(n_index, dist)]

            # sort leaves by distance to node in reverse-order
            sorted_leaves_to_node[n_index] = sorted(inode.get_leaf_names(), key=lambda leaf: currn_leaf_to_dist[leaf][-1], reverse=True)

        # put result in result queue
        ldn_queue.put(currn_leaf_to_dist)
        ntl_queue.put(sorted_leaves_to_node)

    # !-- multi-processing: get leaf distance to node --#

    def node_indexing(self):
        '''
        1) Index tree nodes by level-order.
        2) Annotate node id to tree string.
        3) Get leaf to node distances.
        4) Calculate pairwise inter-node distances using leaf to node distances
        5) Calculate mean distance of child-nodes of each node to root
        '''

        tree_string = self.tree_object.write(format=5)  # append node id annotation

        node_to_leaves = {}
        nindex_to_node = {}
        node_to_nindex = {}
        node_to_parent_node = {}

        # binary indicating that treeinfo file was parsed
        self.treeinfo_file_given = 0
        if len(self.leaf_dist_to_node) > 0:
            self.treeinfo_file_given = 1
        else:
            if self.no_treeinfo:
                print ('\nWARNING: NO TREEINFO FILE WILL BE GENERATED.\n')

        # level-order traversal
        print ('\nIndexing internal nodes...')
        for n, node in enumerate(self.tree_object.traverse()):
            if node.is_leaf():
                continue

            nindex_to_node[n] = node
            node_to_nindex[node] = n

            # get parent node (except for root)
            try:
                node_to_parent_node[n] = node_to_nindex[node.up]
            except:
                pass

            # node annotation for final tree output
            node_string = re.sub('[^\)]+$', '', node.write(format=5))
            tree_string = tree_string.replace(node_string, '{}[&NODE_ID={}]'.format(node_string, n))

        # multi-proc setup
        manager = mp.Manager()
        # shared memory
        leaf_dist_to_node_queue = manager.Queue()
        node_to_leaves_queue = manager.Queue()
        # generate processes
        processes = []

        nindex_list = nindex_to_node.keys()[:]
        shuffle(nindex_list) # shuffle to make multi-processes more equitable
        increment = int(len(nindex_list)/self.cores)

        for p in xrange(self.cores):
            if p == self.cores-1:
                curr_nindex_list = nindex_list[p*increment:]
            else:
                curr_nindex_list = nindex_list[p*increment:(p*increment)+increment]

            #for n, node in nindex_to_node.items():
            proc = mp.Process(target=self.get_leaf_distance_to_node, args=(curr_nindex_list, [nindex_to_node[n] for n in curr_nindex_list], leaf_dist_to_node_queue, node_to_leaves_queue))
            processes.append(proc)
            proc.start()

        # collect results to dictionary
        for p in xrange(len(processes)):
            node_to_leaves.update(node_to_leaves_queue.get())

            for leaf_key, list_value in leaf_dist_to_node_queue.get().items():
                for (n, distance) in list_value:
                    try:
                        self.leaf_dist_to_node[leaf_key][n] = distance
                    except:
                        self.leaf_dist_to_node[leaf_key] = {n:distance}

        # wait for all processes to end
        for proc in processes:
            proc.join()

        # write to treeinfo file
        if self.treeinfo_file_given < 1 and self.no_treeinfo == False:
            print ('Writing to treeinfo file...')
            output = open(self.treeinfo_fname, 'w')
            json.dump(self.leaf_dist_to_node, output)
            output.write('\n')
            output.close()

        """
        # legacy single-thread code
        node_to_leaves_single = {}
        leaf_dist_to_node_single = {}

        # level-order traversal
        for n, node in enumerate(self.tree_object.traverse()):
            if node.is_leaf():
                continue

            # distance of leaf to each of its ancestral node
            for leaf_node in node.get_leaves():
                leaf = leaf_node.name
                dist = leaf_node.get_distance(node)

                try:
                    leaf_dist_to_node_single[leaf][n] = dist
                except:
                    leaf_dist_to_node_single[leaf] = {n: dist}

            # sort leaves by distance to node in reverse-order
            node_to_leaves_single[n] = sorted(node.get_leaf_names(), key=lambda leaf: self.leaf_dist_to_node[leaf][n], reverse=True)

        
        # check single vs multi-proc
        print ('Keys for leaf_dist_to_node: {}'.format(set(leaf_dist_to_node_single.keys()) == set(self.leaf_dist_to_node.keys())))
        for leaf, node_to_dist in self.leaf_dist_to_node.items():
            if set(node_to_dist.keys()) != set(leaf_dist_to_node_single[leaf].keys()):
                print (leaf, set(node_to_dist.keys())^set(leaf_dist_to_node_single[leaf].keys()))

            for node, dist in node_to_dist.items():
                if dist != leaf_dist_to_node_single[leaf][node]:
                    print (leaf, dist, leaf_dist_to_node_single[leaf][node])

        print ('Keys for node_to_leaves: {}'.format(set(node_to_leaves.keys()) == set(node_to_leaves_single.keys())))
        for node, leaves in node_to_leaves.items():
            if leaves != node_to_leaves_single[node]:
                print (node)
                print (leaves)
                print (node_to_leaves_single[node])
                print ('\n')
        """

        # get ancestral nodepair to dist
        # legacy single thread code (faster)
        ancestral_nodepair_to_dist = {}
        for leaf, node_to_dist in self.leaf_dist_to_node.items():
            ancestors_of_leaf = node_to_dist.keys()[:]
            for (i,j) in itertools.combinations(ancestors_of_leaf, 2):

                if (i in ancestral_nodepair_to_dist and j in ancestral_nodepair_to_dist[i]) or (j in ancestral_nodepair_to_dist and i in ancestral_nodepair_to_dist[j]):
                    continue
                else:
                    ij_dist = abs(node_to_dist[i] - node_to_dist[j])

                    try:
                        ancestral_nodepair_to_dist[i][j] = ij_dist
                    except:
                        ancestral_nodepair_to_dist[i] = {j:ij_dist}

                    try:
                        ancestral_nodepair_to_dist[j][i] = ij_dist
                    except:
                        ancestral_nodepair_to_dist[j] = {i:ij_dist}

        # get sibling nodepair to dist
        # legacy single thread code here (faster)
        sibling_nodepair_to_dist = {}
        for (i,j) in itertools.combinations(ancestral_nodepair_to_dist.keys(), 2):
            if (i in ancestral_nodepair_to_dist and j in ancestral_nodepair_to_dist[i]) or (i in sibling_nodepair_to_dist and j in sibling_nodepair_to_dist[i]) or (j in sibling_nodepair_to_dist and i in sibling_nodepair_to_dist[j]):
                continue
            else:
                ancestors_to_i = [node for node in ancestral_nodepair_to_dist[i].keys() if node < i]
                ancestors_to_j = [node for node in ancestral_nodepair_to_dist[j].keys() if node < j]
                common_ancestors = sorted(set(ancestors_to_i)&set(ancestors_to_j))
                common_ancestor = common_ancestors[-1]
                ij_dist = ancestral_nodepair_to_dist[i][common_ancestor] + ancestral_nodepair_to_dist[j][common_ancestor]

                try:
                    sibling_nodepair_to_dist[i][j] = ij_dist
                except:
                    sibling_nodepair_to_dist[i] = {j:ij_dist}

                try:
                    sibling_nodepair_to_dist[j][i] = ij_dist
                except:
                    sibling_nodepair_to_dist[j] = {i:ij_dist}

        nodepair_to_dist = ancestral_nodepair_to_dist.copy()

        for i in nodepair_to_dist.keys():
            nodepair_to_dist[i][i] = 0
            try:
                nodepair_to_dist[i].update(sibling_nodepair_to_dist[i])
            except:
                continue

        # get mean distance of children nodes of each node to root
        node_to_mean_child_dist2root = {n:np.mean([self.leaf_dist_to_node[child.name][0] if child.is_leaf() else nodepair_to_dist[node_to_nindex[child]][0] for child in nindex_to_node[n].get_children()]) for n in node_to_leaves.keys()}

        return tree_string, node_to_leaves, nindex_to_node, node_to_nindex, self.leaf_dist_to_node, nodepair_to_dist, node_to_parent_node, node_to_mean_child_dist2root

    # !-- remap non-string keys for json dump -- #
    def remap_keys(self, mapping):
        return [{'key':k, 'value': v} for k, v in mapping.iteritems()]

    def pwdist_dist_and_ancestral_trace(self, node_to_leaves, nindex_to_node, node_to_nindex, node_to_mean_child_dist2root, nodepair_to_dist):
        '''
        1) Get pairwise distances of all leaves
        2) Get ancestral/descendant traces
        3) Get pairwise distance distributions of nodes
        4) Get leaf to ancestor trace
        5) Get mean child-nodes distance to ancestral trace
        '''

        #! -- multiprocessing: calculate pairwise leaf distance -- #
        def get_pw_leaf_dist(lp_list, queue):
            lp_to_dist = {}
            for (x, y) in lp_list:
                lp_to_dist[(x.name, y.name)] = lp_to_dist[(y.name, x.name)] = x.get_distance(y)

            queue.put(lp_to_dist)
        # ! -- multiprocessing: calculate pairwise leaf distance -- #

        # get pairwise sequence patristic distance
        if self.treeinfo_file_given < 1:
            print ('\nParsing all pairwise distances between leaves...')

            """
            # multi-proc setup (pool)
            pool = mp.Pool(processes=self.cores)
            result = pool.map(get_pw_leaf_dist, list(itertools.combinations(self.tree_object.get_leaves(), 2)))

            for (leaf_x, leaf_y, dist) in result:
                self.leafpair_to_distance[(leaf_x, leaf_y)] = self.leafpair_to_distance[(leaf_y, leaf_x)] = dist
            """

            # multi-proc setup
            manager = mp.Manager()
            # shared memory
            leafpair_to_distance_queue = manager.Queue()
            # generate processes
            processes = []

            leafpair_list = list(itertools.combinations(self.tree_object.get_leaves(), 2))
            increment = int(len(leafpair_list)/self.cores)

            for p in xrange(self.cores):
                if p == self.cores-1:
                    curr_leafpair_list = leafpair_list[p*increment:]
                else:
                    curr_leafpair_list = leafpair_list[p*increment:(p*increment)+increment]

                #for n, node in nindex_to_node.items():
                proc = mp.Process(target=get_pw_leaf_dist, args=(curr_leafpair_list, leafpair_to_distance_queue))
                processes.append(proc)
                proc.start()

            # collect results to dictionary
            for p in xrange(len(processes)):
                self.leafpair_to_distance.update(leafpair_to_distance_queue.get())

            # wait for all processes to end
            for proc in processes:
                proc.join()

            if self.no_treeinfo == False:
                print ('Writing to treeinfo file...')
                output = open(self.treeinfo_fname, 'a')
                json.dump(self.remap_keys(self.leafpair_to_distance), output)
                output.write('\n')
                output.close()

            """# single thread legacy code
            leafpair_to_distance_single = {}
            for x, y in itertools.combinations(self.tree_object.get_leaves(), 2):
                leaf_x = x.name
                leaf_y = y.name
                dist = x.get_distance(y)
                leafpair_to_distance_single[(leaf_x, leaf_y)] = leafpair_to_distance_single[(leaf_y, leaf_x)] = dist

            print ('Keys to leafpair_to_distance: {}'.format(set(self.leafpair_to_distance.keys()) == set(leafpair_to_distance_single.keys())))
            for (leaf_x, leaf_y), dist in self.leafpair_to_distance.items():
                if dist != leafpair_to_distance_single[(leaf_x, leaf_y)]:
                    print (leaf_x, leaf_y, dist, leafpair_to_distance_single[(leaf_x, leaf_y)])"""

        node_to_ancestral_nodes = {}
        node_to_descendant_nodes = {}
        node_to_pwdist = {}
        node_to_mean_pwdist = {}
        leaf_to_ancestors = {}
        node_to_mean_child_dist2anc = {}

        # get ancestry and pairwise sequence distance distribution (single thread code faster)
        print ('\nSorting lineages and PWD distributions...')
        for n in sorted(node_to_mean_child_dist2root, key=node_to_mean_child_dist2root.get):
            leaves = node_to_leaves[n]
            mean_dist2root = node_to_mean_child_dist2root[n]

            # get leaf to ancestor nodes subtending it
            for leaf in leaves:
                try:
                    leaf_to_ancestors[leaf].append(n)
                except:
                    leaf_to_ancestors[leaf] = [n]

            ancestors_to_n = [node_to_nindex[anc] for anc in nindex_to_node[n].iter_ancestors()]

            node_to_ancestral_nodes[n] = ancestors_to_n
            for anc in ancestors_to_n:
                try:
                    node_to_descendant_nodes[anc].append(n)
                except:
                    node_to_descendant_nodes[anc] = [n]

                try:
                    node_to_mean_child_dist2anc[n][anc] = mean_dist2root - nodepair_to_dist[anc][0]
                except:
                    node_to_mean_child_dist2anc[n] = {anc:mean_dist2root - nodepair_to_dist[anc][0]}

            pwdist = sorted([self.leafpair_to_distance[(x,y)] for (x,y) in itertools.combinations(leaves, 2)])
            node_to_pwdist[n] = pwdist
            node_to_mean_pwdist[n] = np.mean(pwdist)

        return self.leafpair_to_distance, node_to_pwdist, node_to_mean_pwdist, node_to_ancestral_nodes, node_to_descendant_nodes, leaf_to_ancestors, node_to_mean_child_dist2anc

    def get_global_pval(self, hytest_method, node_to_leaves, node_to_ancestral_nodes, node_to_pwdist): #, leafpair_to_distance):
        '''
        Perform all inter-clusters' hypotheses tests
        '''
        if self.treeinfo_file_given < 1:
            from ctypes import c_char_p
            #import multiprocessing as mppy
            import multiprocessing.sharedctypes as mpsc

            # global
            lpd = self.leafpair_to_distance

            # worker
            def get_interclus_pval(np_list, method, ntl_dict, ntan_dict, ntpwd_dict, q):
                currp_np_to_pval = {}
                for (i,j) in np_list:
                    if (ntan_dict[j] != False and i in list(ntan_dict[j])) or (ntan_dict[i] != False and j in list(ntan_dict[i])):
                        pval = inter_cluster_hytest(list(ntpwd_dict[i]), list(ntpwd_dict[j])).hytest(method)
                    else:
                        ij_pwdist = sorted([lpd[(x, y)] for x, y in itertools.combinations(list(set(ntl_dict[i])|set(ntl_dict[j])), 2)])
                        # take the conservative (max) p-value comparing node i/j individually to i+j
                        pval = max([inter_cluster_hytest(list(ntpwd_dict[i]), ij_pwdist).hytest(method), inter_cluster_hytest(list(ntpwd_dict[j]), ij_pwdist).hytest(method)])
                    currp_np_to_pval[(i,j)] = pval
                q.put(currp_np_to_pval)

            print ('\nPerforming {} tests...'.format(hytest_method))

            # multi-proc setup
            manager = mp.Manager()

            # shared memory
            queue = manager.Queue()

            max_node = max(node_to_leaves.keys())
            node_to_leaves_shared = [mpsc.Array(c_char_p, node_to_leaves[n]) if n in node_to_leaves.keys() else False for n in xrange(max_node+1)]
            node_to_ancestral_nodes_shared = [mp.Array('i', node_to_ancestral_nodes[n]) if n in node_to_ancestral_nodes else False for n in xrange(max_node+1)]
            node_to_pwdist_shared = [mp.Array('d', node_to_pwdist[n]) if n in node_to_leaves.keys() else False for n in xrange(max_node+1)]

            # generate processes
            processes = []

            # split nodepair list into ncpu sets
            nodepair_list = list(itertools.combinations(node_to_leaves.keys(), 2))
            shuffle(nodepair_list) # shuffle to make more equitable

            increment = int(len(nodepair_list)/self.cores)
            for p in xrange(self.cores):
                if p == self.cores-1:
                    curr_nodepair_list = nodepair_list[p*increment:]
                else:
                    curr_nodepair_list = nodepair_list[p*increment:(p*increment)+increment]

                proc = mp.Process(target=get_interclus_pval, args=(curr_nodepair_list, hytest_method, node_to_leaves_shared, node_to_ancestral_nodes_shared, node_to_pwdist_shared, queue))
                processes.append(proc)
                proc.start()

            # collect results to dictionary
            for p in xrange(len(processes)):
                self.nodepair_to_pval.update(queue.get())

            # wait for all processes to end
            for proc in processes:
                proc.join()

            if self.no_treeinfo == False:
                print ('Writing to treeinfo file...')
                output = open(self.treeinfo_fname, 'a')
                json.dump(self.remap_keys(self.nodepair_to_pval), output)
                output.write('\n')
                output.close()

            """
            # single thread legacy code
            nodepair_to_pval_single = {}
            for i,j in itertools.combinations(node_to_leaves.keys(), 2):
                if (j in node_to_ancestral_nodes and i in node_to_ancestral_nodes[j]) or (i in node_to_ancestral_nodes and j in node_to_ancestral_nodes[i]):
                    pval = inter_cluster_hytest(node_to_pwdist[i], node_to_pwdist[j]).hytest(hytest_method)
                else:
                    ij_pwdist = sorted([leafpair_to_distance[(x,y)] for x,y in itertools.combinations(list(set(node_to_leaves[i])|set(node_to_leaves[j])), 2)])
                    # take the conservative (max) p-value comparing node i/j individually to i+j
                    pval = max([inter_cluster_hytest(node_to_pwdist[i], ij_pwdist).hytest(hytest_method), inter_cluster_hytest(node_to_pwdist[j], ij_pwdist).hytest(hytest_method)])

                nodepair_to_pval_single[(i,j)] = pval

            # check
            print ('Sets of nodepair_to_pval: {}'.format(set(self.nodepair_to_pval.keys()) == set(nodepair_to_pval_single.keys())))
            for (i, j), pval in self.nodepair_to_pval.items():
                if pval != nodepair_to_pval_single[(i,j)]:
                    print (i, j, pval, nodepair_to_pval_single[(i, j)])
            """

        return self.nodepair_to_pval
