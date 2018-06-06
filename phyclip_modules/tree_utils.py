import ete3
import re

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

def dissociate_subtrees(nodes_list, node_to_leaves, node_to_mean_pwdist, node_to_ancestral_nodes, node_to_descendant_nodes, nodepair_to_qval, within_cluster_limit, master_node_to_ancestral_nodes, master_node_to_descendant_nodes, master_node_to_leaves):
    '''
    # dissociate potential clades - starting from the root by level order
    '''
    nodes_to_skip = []

    for node in sorted(nodes_list):
        if (node in nodes_to_skip) or (node_to_mean_pwdist[node] > within_cluster_limit):
            continue

        ancestors_to_node = list(set(master_node_to_ancestral_nodes[node])&set(nodes_list))
        try:
            descendants_to_node_and_node = list(set(master_node_to_descendant_nodes[node])&set(nodes_list)) + [node]
        except:
            descendants_to_node_and_node = [node]

        # for all ancestors to node
        ancestors_to_dissociate = []
        for anc in ancestors_to_node:
            # dissociate descendant leaves
            node_to_leaves[anc] = list(set(node_to_leaves[anc]) - set(master_node_to_leaves[node]))
            # dissociate descendant nodes
            node_to_descendant_nodes[anc] = list(set(node_to_descendant_nodes[anc])-set(descendants_to_node_and_node))
            # append anc to dissociate
            ancestors_to_dissociate.append(anc)

        # dissociate ancestors from descendants
        if len(ancestors_to_dissociate) > 0:
            for desc in descendants_to_node_and_node:
                node_to_ancestral_nodes[desc] = list(set(node_to_ancestral_nodes[desc])-set(ancestors_to_dissociate))

            nodes_to_skip = list(set(nodes_to_skip)|set(descendants_to_node_and_node))

    return node_to_leaves, node_to_ancestral_nodes, node_to_descendant_nodes, sorted(nodes_to_skip)


def get_pwdist_from_leaf_distances_to_node(leaves_dist_to_node, desc_node_to_leaves, node_to_ancestral_nodes, nodepair_to_dist):
    n_i = len(leaves_dist_to_node)

    term_a = (n_i-1)*sum(leaves_dist_to_node)

    term_b = 0
    for desc_node, desc_node_leaves in desc_node_to_leaves.items():
        n_desc = len(desc_node_leaves)
        parent_node = sorted(node_to_ancestral_nodes[desc_node][:])[-1]
        term_b += (n_desc)*(n_desc-1)*nodepair_to_dist[parent_node][desc_node]

    return 2*(term_a-term_b)/(n_i*(n_i-1))

def dissociate_outlier_leaves(node_to_leaves, node_to_mean_pwdist, within_cluster_limit, master_leaf_dist_to_node, node_to_descendant_nodes, node_to_ancestral_nodes, master_nodepair_to_dist, master_node_to_descendant_nodes, master_node_to_ancestral_nodes, list_of_ancestral_node, leaf_to_ancestors):
    '''
    dissociate outlier leaves for clusters that are still > wcl
    '''
    nodes_with_leaf_membership_changes = []
    node_to_leaves_removed = {}

    for node, leaves in node_to_leaves.items():
        # only for those nodes with mean pairwise distance > wcl
        if node_to_mean_pwdist[node] > within_cluster_limit:
            # get distances of leaves to node
            leaf_to_node_dist = {leaf:master_leaf_dist_to_node[leaf][node] for leaf in leaves}
            # sort leaves starting from the longest distance to node
            sorted_leaves = sorted(leaf_to_node_dist.keys(), key=leaf_to_node_dist.get, reverse=True)

            if node in node_to_descendant_nodes:
                # descendant nodes of node
                descendant_nodes_of_node = node_to_descendant_nodes[node][:]

                leaves_to_remove = []
                # leave-one-out strategy
                for l in xrange(1, len(sorted_leaves), 1):
                    reduced_leaf_to_dist = {leaf:leaf_to_node_dist[leaf] for leaf in sorted_leaves[l:]}
                    # break if only 1 leaf left
                    if len(reduced_leaf_to_dist) == 1:
                        break

                    reduced_descendant_nodes_to_leaves = {}
                    for leaf in reduced_leaf_to_dist.keys():
                        for desc_node in list(set(leaf_to_ancestors[leaf])&set(descendant_nodes_of_node)):
                            try:
                                reduced_descendant_nodes_to_leaves[desc_node].append(leaf)
                            except:
                                reduced_descendant_nodes_to_leaves[desc_node] = [leaf]

                    # calculate mean pwdist of reduced set
                    reduced_mean_pwdist = get_pwdist_from_leaf_distances_to_node(reduced_leaf_to_dist.values(), reduced_descendant_nodes_to_leaves, node_to_ancestral_nodes, master_nodepair_to_dist)

                    # break if current mean pairwise distane <= wcl
                    if reduced_mean_pwdist <= within_cluster_limit:
                        leaves_to_remove = list(set(leaves)-set(reduced_leaf_to_dist.keys()))
                        break

                # save leaves to remove
                if len(leaves_to_remove) > 0:
                    nodes_with_leaf_membership_changes.append(node)
                    node_to_leaves_removed[node] = leaves_to_remove

    for node, leaves_removed in node_to_leaves_removed.items():
        # remove leaves from node
        node_to_leaves[node] = list(set(node_to_leaves[node])-set(leaves_removed))
        if node in node_to_descendant_nodes:
            for desc in sorted(node_to_descendant_nodes[node]):
                # if all leaves of descending node are to be removed in ancestor node
                if set(node_to_leaves[desc]) <= set(leaves_removed):
                    # dissociate descending node and future descending nodes from the ancestor
                    try:
                        descendants_to_remove = list(set(master_node_to_descendant_nodes[desc] + [desc])&set(list_of_ancestral_node))
                    except:
                        descendants_to_remove = [desc]
                    node_to_descendant_nodes[node] = list(set(node_to_descendant_nodes[node])-set(descendants_to_remove))
                    # dissociate all ancestor nodes for the descendants
                    try:
                        ancestors_to_remove = list(set(master_node_to_ancestral_nodes[node] + [node])&set(list_of_ancestral_node))
                    except:
                        ancestors_to_remove = [node]
                    node_to_ancestral_nodes[desc] = list(set(node_to_ancestral_nodes[desc])-set(ancestors_to_remove))

    return node_to_leaves, node_to_ancestral_nodes, node_to_descendant_nodes, nodes_with_leaf_membership_changes

def post_solver_cleanup(clusterid_to_taxa, taxon_to_clusterid, master_node_to_leaves, master_leaf_dist_to_node, node_to_leaves, master_node_to_descendant_nodes):
    '''
    Cluster potential cluster-able descendant clades to their parent clade which was clustered (supercluster)
    '''

    supercluster_to_taxa = {}
    taxon_to_supercluster = {}
    # start from most descendant cluster working towards the most ancestral one
    for clusterid in sorted(clusterid_to_taxa.keys(), reverse=True):
        taxa = clusterid_to_taxa[clusterid]
        # if set of clustered leaves != original as parsed from tree
        if set(taxa) != set(master_node_to_leaves[clusterid]):
            # get the distance of the farthest away leaf of node (i.e. minus outliers)
            max_dist = max([master_leaf_dist_to_node[leaf][clusterid] for leaf in node_to_leaves[clusterid] if leaf in clusterid_to_taxa[clusterid]])
            for leaf in master_node_to_leaves[clusterid]:
                # cluster to current cluster if leaf <= max_dist and not clustered in another other clusters
                if (master_leaf_dist_to_node[leaf][clusterid] <= max_dist) and (leaf not in taxon_to_clusterid):
                    taxon_to_clusterid[leaf] = clusterid

                    taxon_to_supercluster[leaf] = clusterid
                    try:
                        supercluster_to_taxa[clusterid].append(leaf)
                    except:
                        supercluster_to_taxa[clusterid] = [leaf]

    for node, taxa in supercluster_to_taxa.items():
        clusterid_to_taxa[node] = list(set(clusterid_to_taxa[node])|set(taxa))

    # de-cluster outliers of descendent cluster that got clustered in ancestor cluster
    clusterid_to_skip = []
    # start from the most ancestral to the most descendant clusters
    for clusterid in sorted(clusterid_to_taxa.keys()):
        taxa = clusterid_to_taxa[clusterid]
        # skip if all leaves subtended by cluster node is of the same set as those parsed in the tree
        if clusterid in clusterid_to_skip or set(taxa) == set(master_node_to_leaves[clusterid]):
            continue

        # get descendant cluster nodes of current cluster
        descendant_clusterids = list(set(master_node_to_descendant_nodes[clusterid])&set(clusterid_to_taxa.keys()))
        if len(descendant_clusterids) > 0:
            for desc_clusterid in descendant_clusterids:
                # skip if all leaves subtended by cluster node is of the same set as those parsed in the tree
                if set(master_node_to_leaves[desc_clusterid]) == set(clusterid_to_taxa[desc_clusterid]):
                    clusterid_to_skip.append(desc_clusterid)
                    continue
                # de-cluster outliers of descendent cluster that got clustered in ancestor cluster
                outliers_of_desc_clusterid = list(set(master_node_to_leaves[desc_clusterid])-set(node_to_leaves[desc_clusterid]))
                outliers_of_desc_clusterid_clustered_in_clusterid = list(set(taxa)&set(outliers_of_desc_clusterid))

                if len(outliers_of_desc_clusterid_clustered_in_clusterid) > 0:
                    for outlier in outliers_of_desc_clusterid_clustered_in_clusterid:
                        del taxon_to_clusterid[outlier]
                        clusterid_to_taxa[clusterid].remove(outlier)
                        if len(clusterid_to_taxa[clusterid]) == 0:
                            del clusterid_to_taxa[clusterid]

    # clean up cluster id - must follow most descendant node
    for cluster in sorted(clusterid_to_taxa.keys()):
        taxa = clusterid_to_taxa[cluster]
        try:
            for node in sorted([cluster] + master_node_to_descendant_nodes[cluster], reverse=True):
                if set(taxa) <= set(master_node_to_leaves[node]):
                    break
        except:
            continue

        if node != cluster:
            clusterid_to_taxa[node] = taxa[:]
            del clusterid_to_taxa[cluster]

    return clusterid_to_taxa, taxon_to_clusterid, taxon_to_supercluster

def clean_up_sensitivity_induced_clusters(clusterid_to_taxa, taxon_to_clusterid, sensitivity_percentile, master_node_to_descendant_nodes, master_leafpair_to_distance, within_cluster_limit):
    """
    Clean up sensitivity-induced clusters
    """
    import numpy as np
    import itertools

    clusterlen_to_frequency = {}
    sensitivity_subsumed_taxa_to_clusterid = {}

    # determine distribution of clusters
    for clusterid in clusterid_to_taxa.keys():
        try:
            clusterlen_to_frequency[len(clusterid_to_taxa[clusterid])] += 1
        except:
            clusterlen_to_frequency[len(clusterid_to_taxa[clusterid])] = 1

    clusterlen_distribution = [i for j in [[clusterlen]*frequency for clusterlen, frequency in clusterlen_to_frequency.items()] for i in j]

    # regard any clusters of length <= 25 percentile lacking evidence to be a potential trajectory
    clusterlen_cutoff = np.percentile(clusterlen_distribution, sensitivity_percentile)

    # determine ancestry of clusters
    cluster_to_desc_clusters = {}
    for cluster in sorted(clusterid_to_taxa.keys()):
        try:
            desc_clusters = list(set(master_node_to_descendant_nodes[cluster])&set(clusterid_to_taxa.keys()))
            if len(desc_clusters) > 0:
                cluster_to_desc_clusters[cluster] = desc_clusters
        except:
            continue

    cluster_to_anc_clusters  = {}
    for cluster, desc_clusters in cluster_to_desc_clusters.items():
        for desc in desc_clusters:
            try:
                cluster_to_anc_clusters[desc].append(cluster)
            except:
                cluster_to_anc_clusters[desc] = [cluster]

    # check nodes that are <= stipulated percentile do not have any descending clusters
    for cluster in sorted(clusterid_to_taxa.keys()):
        taxa = clusterid_to_taxa[cluster]
        if len(taxa) <= clusterlen_cutoff and cluster not in cluster_to_desc_clusters:
            try:
                parent_cluster = sorted(cluster_to_anc_clusters[cluster])[-1]
                parent_taxa = clusterid_to_taxa[parent_cluster][:]
            except:
                continue

            # check that if we subsume cluster into immediate parent cluster, it will still be <= wcl
            combined_mean_pwd = np.mean([master_leafpair_to_distance[(x, y)] for x, y in itertools.combinations(parent_taxa + taxa, 2)])
            if combined_mean_pwd <= within_cluster_limit:
                for taxon in taxa:
                    taxon_to_clusterid[taxon] = parent_cluster
                    sensitivity_subsumed_taxa_to_clusterid[taxon] = parent_cluster

                clusterid_to_taxa[parent_cluster] = list(set(parent_taxa)|set(taxa))
                del clusterid_to_taxa[cluster]

    return clusterid_to_taxa, taxon_to_clusterid, sensitivity_subsumed_taxa_to_clusterid