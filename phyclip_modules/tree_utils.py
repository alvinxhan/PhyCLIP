from copy import deepcopy as dc
import numpy as np
import ete3
import re
import itertools

def collapse_zero_branch_lengths(tree_object, retain_length_bound):
    '''
    Collapse nodes with zero branch lengths and reorder tree
    '''
    for node in tree_object.traverse():
        if node.is_leaf():
            continue
        elif node.dist < retain_length_bound:
            node.delete()

    tree_object.ladderize()
    return tree_object

def parse_newick_tree(tree_file):

    fhandle = filter(None, [_.strip() for _ in open(tree_file, 'rU')])
    fhandle = fhandle.pop(0)
    new_fhandle = []
    prev_end = 0
    for expr in re.finditer('(\(|,)([^(,:]+):', fhandle):
        new_fhandle.append(fhandle[prev_end:expr.start()+1])
        new_fhandle.append("'{}'".format(re.sub("(^'|'$)", "", expr.group(2))))
        prev_end = expr.end()-1
    new_fhandle.append(fhandle[prev_end:])

    return ''.join(new_fhandle)

def parse_treeinfo_file(treeinfo_file, hytest_method):
    '''
    treeinfo file parser
    '''
    leaf_dist_to_node = {}
    leafpair_to_distance = {}
    nodepair_to_pval = {}

    fhandle = open(treeinfo_file, 'rU').readlines()
    for line in fhandle:
        try:
            leaf, n, dist = re.search('L(.+),N(\d+),D([\d\.eE-]+)', line).group(1,2,3)
            try:
                leaf_dist_to_node[leaf][int(n)] = float(dist)
            except:
                leaf_dist_to_node[leaf] = {int(n):float(dist)}
        except:
            try:
                leaf_x, leaf_y, dist = re.search('I(.+),J(.+),D([\d\.eE-]+)', line).group(1, 2, 3)
                leafpair_to_distance[(leaf_x, leaf_y)] = leafpair_to_distance[(leaf_y, leaf_x)] = float(dist)
            except:
                try:
                    i, j, ks_pval, kp_pval = re.search('I(\d+),J(\d+),KS([\d\.eE-]+),KP([\d\.eE-]+)', line).group(1, 2, 3, 4)
                    if hytest_method == 'Kuiper':
                        nodepair_to_pval[(int(i), int(j))] = float(kp_pval)
                    else:
                        nodepair_to_pval[(int(i), int(j))] = float(ks_pval)
                except:
                    i, j, test_type, pval = re.search('I(\d+),J(\d+),(KS|KP)([\d\.eE-]+)', line).group(1, 2, 3, 4)
                    if (test_type == 'KS' and hytest_method == 'Kuiper') or (test_type == 'KP' and hytest_method == 'KS'):
                        continue
                    else:
                        nodepair_to_pval[(int(i), int(j))] = float(pval)

    return leaf_dist_to_node, leafpair_to_distance, nodepair_to_pval

# class of modules to reassociate node and leaves by current min cluster size and within-cluster limit
class node_leaves_reassociation(object):

    def __init__(self, min_cluster_size=None, within_cluster_limit=None, gam_method=None, nodes_list=None, node_to_leaves=None, node_to_descendant_nodes=None, node_to_mean_pwdist=None, node_to_mean_child_dist2anc=None, node_to_parent_node=None, nodepair_to_dist=None, leaf_dist_to_node=None, leaf_to_ancestors=None):

        self.min_cluster_size = min_cluster_size
        self.within_cluster_limit = within_cluster_limit
        self.gam_method = gam_method
        self.nodes_list = nodes_list
        self.node_to_leaves = node_to_leaves
        self.node_to_descendant_nodes = node_to_descendant_nodes
        self.node_to_mean_pwdist = node_to_mean_pwdist
        self.node_to_mean_child_dist2anc = node_to_mean_child_dist2anc
        self.node_to_parent_node = node_to_parent_node
        self.nodepair_to_dist = nodepair_to_dist
        self.leaf_dist_to_node = leaf_dist_to_node
        self.leaf_to_ancestors = leaf_to_ancestors

    def get_pwdist_from_leaf_distances_to_node(self, leaves_dist_to_node, desc_node_to_leaves):
        n_i = len(leaves_dist_to_node)

        term_a = (n_i-1)*sum(leaves_dist_to_node)

        term_b = 0
        for desc_node, desc_node_leaves in desc_node_to_leaves.items():
            n_desc = len(desc_node_leaves)
            parent_node = self.node_to_parent_node[desc_node]
            term_b += (n_desc)*(n_desc-1)*self.nodepair_to_dist[parent_node][desc_node]

        return 2*(term_a-term_b)/(n_i*(n_i-1))

    def remove_outlier_leaves(self, rol_leaves, main_node, rol_node_to_descendant_nodes, rol_node_to_leaves):
        rol_leaf_to_dist = {leaf:self.leaf_dist_to_node[leaf][main_node] for leaf in rol_leaves}

        med_leaf_dist_to_node = np.median(rol_leaf_to_dist.values())
        if self.gam_method == 'Qn':
            from phyclip_modules.stats_utils import qn
            mad_leaf_dist_to_node = qn(rol_leaf_to_dist.values())
        else:
            mad_leaf_dist_to_node = np.median([abs(x-med_leaf_dist_to_node) for x in rol_leaf_to_dist.values() if x >= med_leaf_dist_to_node])

        rol_remaining_leaves = [leaf for leaf, dist in rol_leaf_to_dist.items() if dist <= (med_leaf_dist_to_node + (3*mad_leaf_dist_to_node))]
        rol_leaves_to_dissociate = list(set(rol_node_to_leaves[main_node]) - set(rol_remaining_leaves))

        if len(rol_leaves_to_dissociate) == 0 or len(rol_remaining_leaves) == 1:
            return False
        else:
            remaining_descendant_nodes_to_leaves = {}
            for leaf in rol_remaining_leaves:
                try:
                    desc_nodes_subtending_leaf = list(set(rol_node_to_descendant_nodes[main_node])&set(self.leaf_to_ancestors[leaf]))
                except:
                    desc_nodes_subtending_leaf = []

                for rdn in desc_nodes_subtending_leaf:
                    try:
                        remaining_descendant_nodes_to_leaves[rdn].append(leaf)
                    except:
                        remaining_descendant_nodes_to_leaves[rdn] = [leaf]

            reduced_mean_pwdist = self.get_pwdist_from_leaf_distances_to_node([rol_leaf_to_dist[leaf] for leaf in rol_remaining_leaves], remaining_descendant_nodes_to_leaves)

            try:
                # find descendant nodes which fully subtend leaves to be removed
                nodes_of_dissociated_leaves = list(set(rol_node_to_descendant_nodes[main_node])&set([x for y in [self.leaf_to_ancestors[leaf] for leaf in rol_leaves_to_dissociate] for x in y]))
                rol_descendant_nodes_to_dissociate = [desc_node for desc_node in nodes_of_dissociated_leaves if set(rol_node_to_leaves[desc_node]) <= set(rol_leaves_to_dissociate)]

            except:
                # no descendant nodes to main_node
                rol_descendant_nodes_to_dissociate = []

            return rol_remaining_leaves, rol_descendant_nodes_to_dissociate, reduced_mean_pwdist

    def leave_one_out_leaf_reduction(self, sorted_leaves, main_node):
        # note that sorted_leaves is already reverse sorted by distance to main_node
        # immediate return if len(sorted_leaves) < self.min_cluster_size
        if len(sorted_leaves) < self.min_cluster_size:
            return False

        for l_index in xrange(-1, len(sorted_leaves), 1):
            if l_index + 1 == len(sorted_leaves) - self.min_cluster_size:
                return False
                break

            remaining_leaves_to_node_dist = {}
            remaining_descendant_nodes_to_leaves = {}

            # start from not having any leaves removed
            for leaf in sorted_leaves[l_index+1:]:
                remaining_leaves_to_node_dist[leaf] = self.leaf_dist_to_node[leaf][main_node]

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
                return sorted_leaves[l_index+1:], old_descendant_nodes_to_dissociate, reduced_mean_pwdist
                break

        return False

    def nla_main(self):

        # create deep copies to be edited
        sd_node_to_leaves = dc(self.node_to_leaves)
        sd_node_to_descendant_nodes = dc(self.node_to_descendant_nodes)
        sd_node_to_mean_pwdist = {} # mean pairwise distance dictionary for current set of parameters

        # reassociate subtrees - starting from the root by level order
        for node in self.nodes_list: # nodes_list already sorted

            # current node <= self.within_cluster_limit
            if self.node_to_mean_pwdist[node] <= self.within_cluster_limit:
                # check there are no outlying leaves
                rol_output = self.remove_outlier_leaves(self.node_to_leaves[node], node, self.node_to_descendant_nodes, self.node_to_leaves)
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

            # current > self.within_cluster_limit
            descendant_nodes_to_dissociate = []  # list to save nodes for dissociation
            try:
                # descendants nodes are already sorted by mean child-nodes' distance to node
                for desc_node in self.node_to_descendant_nodes[node]:
                    if desc_node in descendant_nodes_to_dissociate:
                        continue

                    elif self.node_to_mean_child_dist2anc[desc_node][node] > self.within_cluster_limit:
                        # append not only desc_node but all descendant nodes of desc_node itself
                        descendant_nodes_to_dissociate.append(desc_node)
                        try:
                            descendant_nodes_to_dissociate = list(set(descendant_nodes_to_dissociate)|set(self.node_to_descendant_nodes[desc_node]))
                        except:
                            pass
                # code con't below...
            except:
                # dead-end node with no descendants but > self.within_cluster_limit
                loo_output = self.leave_one_out_leaf_reduction(self.node_to_leaves[node], node)

                # no leaves to remove by leave-one-out-wcl approach
                if loo_output == False:
                    # check there are no outlying leaves by distance to node
                    rol_output = self.remove_outlier_leaves(self.node_to_leaves[node], node, self.node_to_descendant_nodes, self.node_to_leaves)
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
                    rol_output = self.remove_outlier_leaves(sd_node_to_leaves[node], node, sd_node_to_descendant_nodes[node], sd_node_to_leaves)

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

            # ...resumed code from above
            leaves_to_remove = list(set([x for y in [self.node_to_leaves[desc_node] for desc_node in descendant_nodes_to_dissociate] for x in y]))
            # remove all leaves from nodes that could be potentially dissociated (leaves in self.node_to_leaves[node] already reverse-sorted by distance to node)
            remaining_leaves = [leaf for leaf in self.node_to_leaves[node] if leaf not in leaves_to_remove]

            loo_output = self.leave_one_out_leaf_reduction(remaining_leaves, node)

            # leave-one-out-wcl approach removes leaves
            if loo_output != False:

                leaves_to_keep, loo_descendant_nodes_to_dissociate, mean_pwdist = loo_output
                # update as per output of leave-one-out-wcl approach
                # update node to remaining leaves
                sd_node_to_leaves[node] = leaves_to_keep[:]
                # dissociate descendant nodes from node
                descendant_nodes_to_dissociate = list(set(descendant_nodes_to_dissociate)|set(loo_descendant_nodes_to_dissociate))
                sd_node_to_descendant_nodes[node] = list(set(sd_node_to_descendant_nodes[node])-set(descendant_nodes_to_dissociate))

                # now check that there are no outlying leaves by distance to node
                rol_output = self.remove_outlier_leaves(sd_node_to_leaves[node], node, sd_node_to_descendant_nodes[node], sd_node_to_leaves)
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

            # leave-one-out-wcl approach did not remove any leaves -- THIS SECTION OF THE CODE MAY BE REDUNDANT --
            else:
                # perform leave-one-out-wcl approach again, now based on self.node_to_leaves[node]
                loo_output = self.leave_one_out_leaf_reduction(self.node_to_leaves[node], node)
                # leave-one-out-wcl approach still did not remove any leaves
                if loo_output == False:
                    # check for outlying leaves by distance to node
                    rol_output = self.remove_outlier_leaves(self.node_to_leaves[node], node, self.node_to_descendant_nodes, self.node_to_leaves)
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
                    rol_output = self.remove_outlier_leaves(sd_node_to_leaves[node], node, sd_node_to_descendant_nodes[node], sd_node_to_leaves)
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

        return sd_node_to_leaves, sd_node_to_descendant_nodes, sd_node_to_mean_pwdist

# clean-up modules
class clean_up_modules(object):
    def __init__(self, current_node_to_descendant_nodes=None, node_to_parent_node=None, node_to_leaves=None, leafpair_to_distance=None, current_node_to_leaves=None, within_cluster_limit=None, min_cluster_size=None):
        self.current_node_to_descendant_nodes = current_node_to_descendant_nodes
        self.node_to_leaves = node_to_leaves
        self.leafpair_to_distance = leafpair_to_distance
        self.current_node_to_leaves = current_node_to_leaves
        self.within_cluster_limit = within_cluster_limit
        self.min_cluster_size = min_cluster_size
        self.node_to_parent_node = node_to_parent_node

    def most_desc_nodeid_for_cluster(self, clusterid_to_taxa, taxon_to_clusterid):
        # clean up cluster id - must follow most descendant node
        for cluster in sorted(clusterid_to_taxa.keys()): # start from most ancestral node, proceed by leve-order
            taxa = clusterid_to_taxa[cluster]
            try:
                # find the most descendant node which subtends all of clustered taxa
                for node in sorted([cluster] + self.current_node_to_descendant_nodes[cluster], reverse=True):
                    if set(taxa) <= set(self.current_node_to_leaves[node]):
                        break
            except:
                continue

            # if the most descendant node subtending all taxa != cluster-id
            if node != cluster:
                clusterid_to_taxa[node] = taxa[:]
                del clusterid_to_taxa[cluster]

                for taxon in taxa:
                    taxon_to_clusterid[taxon] = node

        return clusterid_to_taxa, taxon_to_clusterid

    def get_cluster_size_distribution(self, clusterid_to_taxa):
        clusterlen_to_frequency = {}
        for clusterid in clusterid_to_taxa.keys():
            try:
                clusterlen_to_frequency[len(clusterid_to_taxa[clusterid])] += 1
            except:
                clusterlen_to_frequency[len(clusterid_to_taxa[clusterid])] = 1

        return [i for j in [[clusterlen] * frequency for clusterlen, frequency in clusterlen_to_frequency.items()] for i in j]

    def subsume_subclusters_under_x_percentile(self, clusterid_to_taxa, taxon_to_clusterid, clusterlen_distribution, percentile):

        subsumed_taxa_to_clusterid = {}

        # regard any clusters of length <= 25 (default) percentile lacking evidence to be a potential trajectory
        clusterlen_cutoff = np.percentile(clusterlen_distribution, percentile)
        if percentile < 100:
            print ('Subsuming clusters <= {} taxa...'.format(int(clusterlen_cutoff)))

        # determine ancestry of clusters
        cluster_to_desc_clusters = {}
        for cluster in sorted(clusterid_to_taxa.keys()):
            try:
                desc_clusters = list(set(self.current_node_to_descendant_nodes[cluster])&set(clusterid_to_taxa.keys()))
                if len(desc_clusters) > 0:
                    cluster_to_desc_clusters[cluster] = desc_clusters
            except:
                continue

        """
        cluster_to_anc_clusters  = {}
        for cluster, desc_clusters in cluster_to_desc_clusters.items():
            for desc in desc_clusters:
                try:
                    cluster_to_anc_clusters[desc].append(cluster)
                except:
                    cluster_to_anc_clusters[desc] = [cluster]
        """

        # check nodes which are <= x-percentile and do not have any descending clusters
        for cluster in sorted(clusterid_to_taxa.keys()):
            taxa = clusterid_to_taxa[cluster]
            if len(taxa) <= clusterlen_cutoff and cluster not in cluster_to_desc_clusters:
                try:
                    #parent_cluster = sorted(cluster_to_anc_clusters[cluster])[-1]
                    #parent_taxa = clusterid_to_taxa[parent_cluster][:]
                    parent_node = self.node_to_parent_node[cluster]
                    parent_taxa = clusterid_to_taxa[parent_node][:]
                except:
                    continue

                # check that taxa set <= self.current_node_to_leaves[parent_cluster] (in other words, the parent cluster, a statistically significant cluster, is only so with the inclusion of the taxa set as well)
                if set(taxa) < set(self.current_node_to_leaves[parent_cluster]):
                    # check that if we subsume cluster into immediate parent cluster, it would still fulfill self.within_cluster_limit
                    combined_mean_pwd = np.mean([self.leafpair_to_distance[(x, y)] for x, y in itertools.combinations(parent_taxa + taxa, 2)])
                    if combined_mean_pwd <= self.within_cluster_limit:
                        for taxon in taxa:
                            taxon_to_clusterid[taxon] = parent_cluster
                            subsumed_taxa_to_clusterid[taxon] = cluster

                        clusterid_to_taxa[parent_cluster] = list(set(parent_taxa)|set(taxa))
                        del clusterid_to_taxa[cluster]

        return clusterid_to_taxa, taxon_to_clusterid, subsumed_taxa_to_clusterid

    def remove_clusters_below_cs(self, clusterid_to_taxa, taxon_to_clusterid):
        # remove any clusters < min cluster size
        for clusterid, taxa in clusterid_to_taxa.items():
            if len(taxa) < self.min_cluster_size:
                for taxon in taxa:
                    del taxon_to_clusterid[taxon]
                del clusterid_to_taxa[clusterid]
        return clusterid_to_taxa, taxon_to_clusterid

    def decluster_outlying_taxa_clustered_to_anc_clusters(self, clusterid_to_taxa, taxon_to_clusterid):
        clusterid_to_skip = []
        # start from the most ancestral to the most descendant clusters
        for clusterid in sorted(clusterid_to_taxa.keys()):
            taxa = clusterid_to_taxa[clusterid]
            # skip if all leaves subtended by cluster node is of the same set as those parsed in the tree
            if clusterid in clusterid_to_skip or set(taxa) == set(self.node_to_leaves[clusterid]):
                continue

            # get descendant cluster nodes of current cluster
            try:
                descendant_clusterids = list(set(self.current_node_to_descendant_nodes[clusterid])&set(clusterid_to_taxa.keys()))
            except:
                continue

            if len(descendant_clusterids) > 0:
                for desc_clusterid in descendant_clusterids:
                    # skip if all leaves subtended by cluster node is of the same set as those parsed in the tree
                    if set(self.node_to_leaves[desc_clusterid]) == set(clusterid_to_taxa[desc_clusterid]):
                        clusterid_to_skip.append(desc_clusterid)
                        continue

                    # leaf outliers of nodes appears in global node_to_leaves dict but not in current
                    outliers_of_desc_clusterid = list(set(self.node_to_leaves[desc_clusterid])-set(self.current_node_to_leaves[desc_clusterid]))
                    # decluster any outlying leaves of descending nodes that are clustered in current node
                    outliers_of_desc_clusterid_clustered_in_clusterid = list(set(taxa)&set(outliers_of_desc_clusterid))

                    if len(outliers_of_desc_clusterid_clustered_in_clusterid) > 0:
                        for outlier in outliers_of_desc_clusterid_clustered_in_clusterid:
                            del taxon_to_clusterid[outlier]
                            clusterid_to_taxa[clusterid].remove(outlier)
                            # delete any empty clusters
                            if len(clusterid_to_taxa[clusterid]) == 0:
                                del clusterid_to_taxa[clusterid]

        return clusterid_to_taxa, taxon_to_clusterid