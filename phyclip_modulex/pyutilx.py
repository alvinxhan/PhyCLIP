from __future__ import division
import re
import itertools
import numpy as np
from copy import deepcopy as dc
import ete3
import time

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

def collapse_zero_branch_length(tree_object, retain_length_bound, treefname):
    '''
    Collapse nodes with zero branch lengths and reorder tree
    '''
    no_node_collapse = 1
    diff = max([retain_length_bound, 1e-6])
    for node in tree_object.traverse():
        if node.is_leaf():
            continue
        elif abs(np.float32(node.dist) - retain_length_bound) <= diff:
            no_node_collapse = 0
            node.delete()

    if no_node_collapse == 1:
        return False # no branches collapsed
    else:
        tree_object.ladderize()

        # write zero branch length collapsed tree as a newick file
        treefname = 'zero-branch-length-collapsed_{}'.format(treefname)
        tree_object.write(format=5, outfile='{}.nwk'.format(treefname))
        print ('Collapsed newick tree file...{}.nwk'.format(treefname))

        return tree_object, treefname

def parse_input_parameters(infhandle):
    # subsequent lines are parameters (cs, fdr, gam)
    parameters = []
    max_gam = -1

    int = np.int64
    float = np.float32

    for _, line in enumerate(infhandle):
        try:
            cs, fdr, gam = line.strip().split(',')
            try:
                parameters.append((int(cs), float(fdr), float(gam))) # single parameter sets
                if gam > max_gam:
                    max_gam = gam
                continue
            except:
                pass
        except:
            raise Exception('\nInvalid parameter set format in line {} of input file.\n'.format(_ + 2))

        # check for range
        cs_range = []
        try:
            lb, rb, increment = re.search('(\d+)-(\d+)\((\d+)\)', cs).group(1, 2, 3)
            if int(rb) <= int(lb):
                raise Exception('\nInvalid min cluster size format in line {} of input file.\n'.format(_ + 2))
            for i in range(int(lb), int(rb)+1, int(increment)):
                cs_range.append(i)
        except:
            try:
                cs_range.append(int(cs))
            except:
                raise Exception('\nInvalid min cluster size format in line {} of input file.\n'.format(_ + 2))

        fdr_range = []
        try:
            lb, rb, increment = re.search('(\d+\.\d+)-(\d+\.\d+)\((\d+\.\d+)\)', fdr).group(1, 2, 3)
            if float(rb) <= float(lb):
                raise Exception('\nInvalid fdr format in line {} of input file.\n'.format(_ + 2))
            currip = float(lb)
            while currip <= float(rb):
                fdr_range.append(currip)
                currip += float(increment)
        except:
            try:
                fdr_range.append(float(fdr))
            except:
                raise Exception('\nInvalid fdr format in line {} of input file.\n'.format(_ + 2))

        gam_range = []
        try:
            lb, rb, increment = re.search('(\d+(\.\d+)*)-(\d+(\.\d+)*)\((\d+(\.\d+)*)\)', gam).group(1, 3, 5)
            if float(rb) <= float(lb):
                raise Exception('\nInvalid gamma format in line {} of input file.\n'.format(_ + 2))
            currip = float(lb)
            while currip <= float(rb):
                gam_range.append(currip)
                currip += float(increment)
        except:
            try:
                gam_range.append(float(gam))
            except:
                raise Exception('\nInvalid gamma format in line {} of input file.\n'.format(_ + 2))

        max_gam = np.amax(gam_range) # maximum gamma

        for x in cs_range:
            for y in fdr_range:
                for z in gam_range:
                    parameters.append((x, y, z))

    print('Parameter sets...OK')

    return parameters, max_gam

def get_cluster_size_distribution(clusterid_to_taxa):
    clusterlen_to_frequency = {}
    for clusterid, taxa in clusterid_to_taxa.items():
        try:
            clusterlen_to_frequency[len(taxa)] += 1
        except:
            clusterlen_to_frequency[len(taxa)] = 1

    return [i for j in [[clusterlen] * frequency for clusterlen, frequency in clusterlen_to_frequency.items()] for i in j]

# clean-up modules
class clean_up_modules():

    def __init__(self, current_node_to_descendant_nodes=None, node_to_leaves=None, current_node_to_leaves=None, within_cluster_limit=None, min_cluster_size=None, leaf_dist_to_node=None, leaf_to_ancestors=None, node_to_parent_node=None, nodepair_to_dist=None, leaf_node_id_to_leaf_name=None):

        self.current_node_to_descendant_nodes = current_node_to_descendant_nodes
        self.node_to_leaves = node_to_leaves
        self.current_node_to_leaves = current_node_to_leaves
        self.within_cluster_limit = within_cluster_limit
        self.min_cluster_size = min_cluster_size
        self.leaf_dist_to_node = leaf_dist_to_node
        self.leaf_to_ancestors = leaf_to_ancestors
        self.node_to_parent_node = node_to_parent_node
        self.nodepair_to_dist = nodepair_to_dist
        self.leaf_node_id_to_leaf_name = leaf_node_id_to_leaf_name

    def get_pwdist_from_leaf_distances_to_node_cleanup(self, leaves_dist_to_node, desc_node_to_leaves):
        n_i = len(leaves_dist_to_node)

        term_a = (n_i-1)*sum(leaves_dist_to_node)
        term_b = 0
        for desc_node, desc_node_leaves in desc_node_to_leaves.items():
            n_desc = len(desc_node_leaves)
            parent_node = self.node_to_parent_node[desc_node]
            term_b += (n_desc)*(n_desc-1)*self.nodepair_to_dist[(parent_node, desc_node)]

        return np.float32(2*(term_a-term_b)/(n_i*(n_i-1)))

    def leave_one_out_leaf_reduction_cleanup(self, sorted_leaves, main_node):
        # note that sorted_leaves is already reverse sorted by distance to main_node
        # immediate return if len(sorted_leaves) < self.min_cluster_size
        if len(sorted_leaves) < self.min_cluster_size:
            return False, None

        for l_index in range(-1, len(sorted_leaves), 1):
            if len(sorted_leaves[l_index+1:]) < self.min_cluster_size:
                return False, None

            remaining_leaves_to_node_dist = {}
            remaining_descendant_nodes_to_leaves = {}

            # start from not having any leaves removed
            for leaf in sorted_leaves[l_index + 1:]:
                remaining_leaves_to_node_dist[leaf] = self.leaf_dist_to_node[(leaf, main_node)]
                try:
                    desc_nodes_subtending_leaf = list(set(self.current_node_to_descendant_nodes[main_node])&set(self.leaf_to_ancestors[leaf]))
                except:
                    desc_nodes_subtending_leaf = []

                for rdn in desc_nodes_subtending_leaf:
                    try:
                        remaining_descendant_nodes_to_leaves[rdn].append(leaf)
                    except:
                        remaining_descendant_nodes_to_leaves[rdn] = [leaf]

            reduced_mean_pwdist = self.get_pwdist_from_leaf_distances_to_node_cleanup(remaining_leaves_to_node_dist.values(), remaining_descendant_nodes_to_leaves)

            # break if <= self.within_cluster_limit and dissociate
            if reduced_mean_pwdist <= self.within_cluster_limit:
                # return leaves to keep for node
                return True, sorted_leaves[l_index+1:]

        return False, None

    def loo_wcl_violation(self, clusterid_to_taxa, taxon_to_clusterid, leaf_node_id_to_leafname):
        from phyilpx_stats import qn

        for clusterid, taxa in clusterid_to_taxa.items():

            taxa_copy = taxa[:]

            N = len(taxa)
            leafpair_list = list(itertools.combinations(range(N), 2))
            _pwdist = np.full((N, N), np.nan, dtype=np.float32)
            cluster_pwdist = np.zeros(len(leafpair_list), dtype = np.float32)
            for _, (_i, _j) in enumerate(leafpair_list):
                i = taxa[_i]
                j = taxa[_j]
                _pwdist[(_i, _j)] = _pwdist[(_j, _i)] = self.nodepair_to_dist[(i, j)]
                cluster_pwdist[_] = self.nodepair_to_dist[(i, j)]

            mean_d = np.zeros(N, dtype = np.float32)
            for i in xrange(N):
                mean_d_i = _pwdist[(i,)]
                mean_d[i] = np.mean(mean_d_i[~np.isnan(mean_d_i)])

            # violated wcl
            if np.mean(cluster_pwdist) > self.within_cluster_limit:

                iid = mean_d.argsort()[::-1]
                rsorted_taxa = list(np.array(taxa, dtype=np.int64)[iid])

                loo_output_binary, loo_output = self.leave_one_out_leaf_reduction_cleanup(rsorted_taxa, clusterid)
                if loo_output_binary == False:
                    # remvoe entire cluster since it entirely violates the within cluster limit
                    del clusterid_to_taxa[clusterid]
                    for taxon in taxa:
                        del taxon_to_clusterid[taxon]
                else:
                    # if we could still have a cluster after removing "outlying" taxa
                    for taxon in list(set(taxa) - set(loo_output)):
                        del taxon_to_clusterid[taxon]
                        clusterid_to_taxa[clusterid].remove(taxon)

            # check that there are no weirdly far away sequence
            median_mean_d = np.median(mean_d)
            mad_mean_d = np.median([_ for _ in mean_d if _ >= median_mean_d])
            tol_mean_d = median_mean_d + 3*mad_mean_d

            for i in xrange(N):
                if mean_d[i] > tol_mean_d:
                    taxon = taxa_copy[i]
                    try:
                        del taxon_to_clusterid[taxon]
                        clusterid_to_taxa[clusterid].remove(taxon)
                    except:
                        continue

            """
            # check mean pairwise distance
            leafpairs_list = list(itertools.combinations(taxa, 2))
            new_pwdist = np.zeros(len(leafpairs_list), dtype='f4')
            for _, (i,j) in enumerate(leafpairs_list):
                new_pwdist[_] = self.nodepair_to_dist[(i, j)]

            if np.mean(new_pwdist) > self.within_cluster_limit:
                # reverse sort clustered taxa by distance to node
                leaf_dist_of_cluster = {}
                #for leaf in self.node_to_leaves[clusterid]:
                for leaf in taxa:
                    leaf_dist_of_cluster[leaf] = self.leaf_dist_to_node[(leaf, clusterid)]
                rsorted_taxa = sorted(leaf_dist_of_cluster.keys(), key=leaf_dist_of_cluster.get, reverse=True)

                loo_output_binary, loo_output = self.leave_one_out_leaf_reduction_cleanup(rsorted_taxa, clusterid)
                if loo_output_binary == False:
                    # remvoe entire cluster since it entirely violates the within cluster limit
                    del clusterid_to_taxa[clusterid]
                    for taxon in taxa:
                        del taxon_to_clusterid[taxon]
                else:
                    # if we could still have a cluster after removing "outlying" taxa
                    for taxon in list(set(taxa) - set(loo_output)):
                        del taxon_to_clusterid[taxon]
                        clusterid_to_taxa[clusterid].remove(taxon)
            """
        return clusterid_to_taxa, taxon_to_clusterid

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

    def subsume_subclusters_under_x_percentile(self, clusterid_to_taxa, taxon_to_clusterid, clusterlen_distribution, percentile):

        subsumed_taxa_to_clusterid = {}

        # regard any clusters of length <= 25 (default) percentile lacking evidence to be a potential trajectory
        clusterlen_cutoff = np.percentile(clusterlen_distribution, percentile)
        if percentile < 100:
            print ('Subsuming cluster-size sensitivity-induced subclusters <= {} taxa (if any)...'.format(int(clusterlen_cutoff)))

        # determine ancestry of clusters
        cluster_to_desc_clusters = {}
        for cluster in sorted(clusterid_to_taxa.keys()):
            try:
                desc_clusters = list(set(self.current_node_to_descendant_nodes[cluster])&set(clusterid_to_taxa.keys()))
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

        # check nodes which are <= x-percentile and do not have any descending clusters
        for cluster in sorted(clusterid_to_taxa.keys()):
            taxa = clusterid_to_taxa[cluster]
            if len(taxa) <= clusterlen_cutoff and cluster not in cluster_to_desc_clusters:
                try:
                    parent_cluster = sorted(cluster_to_anc_clusters[cluster])[-1]
                    parent_taxa = clusterid_to_taxa[parent_cluster][:]
                except:
                    continue

                # check that taxa set <= self.current_node_to_leaves[parent_cluster] (in other words, the parent cluster, a statistically significant cluster, is only so with the inclusion of the taxa set as well)
                if set(taxa) < set(self.current_node_to_leaves[parent_cluster]):
                    # check that if we subsume cluster into immediate parent cluster, it would still fulfill self.within_cluster_limit
                    combined_mean_pwd = np.mean([self.nodepair_to_dist[(x, y)] for x, y in itertools.combinations(parent_taxa + taxa, 2)])
                    if combined_mean_pwd <= self.within_cluster_limit:
                        for taxon in taxa:
                            taxon_to_clusterid[taxon] = parent_cluster
                            subsumed_taxa_to_clusterid[taxon] = cluster

                        clusterid_to_taxa[parent_cluster] = list(set(parent_taxa)|set(taxa))
                        #print parent_cluster, cluster, [self.leaf_node_id_to_leaf_name[taxon] for taxon in taxa], '*'
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

                    #print clusterid, desc_clusterid, [self.leaf_node_id_to_leaf_name[taxon] for taxon in outliers_of_desc_clusterid_clustered_in_clusterid]

                    if len(outliers_of_desc_clusterid_clustered_in_clusterid) > 0:
                        #print clusterid, desc_clusterid, [self.leaf_node_id_to_leaf_name[outlier] for outlier in outliers_of_desc_clusterid_clustered_in_clusterid]

                        for outlier in outliers_of_desc_clusterid_clustered_in_clusterid:
                            del taxon_to_clusterid[outlier]
                            clusterid_to_taxa[clusterid].remove(outlier)
                            # delete any empty clusters
                            if len(clusterid_to_taxa[clusterid]) == 0:
                                del clusterid_to_taxa[clusterid]

        return clusterid_to_taxa, taxon_to_clusterid

    def force_lowinfo_cluster(self, clusterid_to_taxa, taxon_to_clusterid, node_to_mean_pwdist):
        from phyilpx_stats import qn

        # get inter-cluster distance distributions
        inter_cluster_dist = [self.nodepair_to_dist[(i, j)] for (i, j) in itertools.combinations(clusterid_to_taxa.keys(), 2)]
        med_inter_cluster_dist = np.float32(np.median(inter_cluster_dist))
        mad_inter_cluster_dist = np.float32(qn(inter_cluster_dist))
        inter_cluster_tol = np.float32(max([min(inter_cluster_dist), med_inter_cluster_dist - 3.*mad_inter_cluster_dist]))
        #print inter_cluster_tol, min(inter_cluster_dist)

        # determine what is left unclustered
        unclustered_taxa = list(set(self.leaf_to_ancestors.keys()) - set(taxon_to_clusterid.keys()))

        # get descendants of clusters
        clustered_nodes_and_desc_nodes = clusterid_to_taxa.keys()[:]
        for clusterid in clusterid_to_taxa.keys():
            try:
                clustered_nodes_and_desc_nodes = list(set(clustered_nodes_and_desc_nodes)|set(self.current_node_to_descendant_nodes[clusterid]))
            except:
                continue

        # putative forced clusters
        putative_forced_cluster_to_taxa = {}

        # nodes subtending unclustered taxa
        for node in sorted(self.current_node_to_leaves.keys()):
            leaves = self.current_node_to_leaves[node]
            if (len(leaves) < self.min_cluster_size) or (node_to_mean_pwdist[node] > self.within_cluster_limit) or (node in clustered_nodes_and_desc_nodes):
                continue

            overlapping_taxa_set = set(leaves)&set(unclustered_taxa)
            if overlapping_taxa_set and set(leaves) <= set(unclustered_taxa):
                # find nearest cluster to node
                nearest_cluster_to_dist = {clusterid:self.nodepair_to_dist[(clusterid, node)] for clusterid in clusterid_to_taxa.keys()}
                dist_of_nearest_cluster_to_node = np.min(nearest_cluster_to_dist.values())
                nearest_clusters = [clusterid for clusterid, dist in nearest_cluster_to_dist.items() if dist == dist_of_nearest_cluster_to_node]

                # only accept if the distnace to nearest clusters >= inter_cluster_tol
                if dist_of_nearest_cluster_to_node >= inter_cluster_tol:
                    #print node
                    putative_forced_cluster_to_taxa[node] = list(leaves)
                    clusterid_to_taxa[node] = list(leaves)
                    unclustered_taxa = list(set(unclustered_taxa) - set(leaves))

                else:
                    cluster_to_remove = list(set(nearest_clusters)&set(putative_forced_cluster_to_taxa.keys()))
                    if len(cluster_to_remove) > 0:
                        # remove cluster
                        #print cluster_to_remove
                        unclustered_taxa = list(set(unclustered_taxa)|set(clusterid_to_taxa[clusterid]))
                        for clusterid in cluster_to_remove:
                            del putative_forced_cluster_to_taxa[clusterid]
                            del clusterid_to_taxa[clusterid]

        for clusterid, taxa in putative_forced_cluster_to_taxa.items():
            for taxon in taxa:
                taxon_to_clusterid[taxon] = clusterid
        """
        # check if we could include any other "outlying sequences"
        putative_forced_taxon_to_cluster = {}
        for clusterid in sorted(clusterid_to_taxa.keys()):
            unclustered_taxa = list((set(self.node_to_leaves[clusterid]) - set(taxon_to_clusterid.keys()))&set(self.current_node_to_leaves[clusterid]))

            if len(unclustered_taxa) > 0:
                for taxon in unclustered_taxa:
                    putative_forced_taxon_to_cluster[taxon] = clusterid

        putative_forced_cluster_to_taxa = {}
        for taxon, clusterid in putative_forced_taxon_to_cluster.items():
            try:
                putative_forced_cluster_to_taxa[clusterid].append(taxon)
            except:
                putative_forced_cluster_to_taxa[clusterid] = [taxon]

        for clusterid, taxa in putative_forced_cluster_to_taxa.items():
            if clusterid in clusterid_to_taxa:
                clusterid_to_taxa[clusterid] = list(set(taxa)|set(clusterid_to_taxa[clusterid]))
                for taxon in taxa:
                    taxon_to_clusterid[taxon] = clusterid
            else:
                if len(taxa) >= self.min_cluster_size:
                    clusterid_to_taxa[clusterid] = taxa[:]
                    for taxon in taxa:
                        taxon_to_clusterid[taxon] = clusterid
        """

        return clusterid_to_taxa, taxon_to_clusterid


class output_utils(object):
    """
    Print output cluster and tree files
    """

    def __init__(self, ori_tree_string=None, leaf_node_id_to_leafname=None, taxon_to_clusterid=None, clusterid_to_taxa=None, outfname=None, sensitivity_subsumed_taxa_to_clusterid=False, nosub_taxa_to_clusterid=False):

        self.ori_tree_string = ori_tree_string

        # list of leaves by name
        self.leaf_node_id_to_leafname = leaf_node_id_to_leafname
        self.taxon_list = leaf_node_id_to_leafname.values()

        # convert leaf node_id to leaf name
        self.taxon_to_clusterid = {leaf_node_id_to_leafname[taxon]:clusterid for taxon, clusterid in taxon_to_clusterid.items()}
        self.clusterid_to_taxa = {clusterid:[leaf_node_id_to_leafname[taxon] for taxon in taxa] for clusterid, taxa in clusterid_to_taxa.items()}

        self.outfname = outfname

        self.sensitivity_subsumed_taxa_to_clusterid = sensitivity_subsumed_taxa_to_clusterid
        self.nosub_taxa_to_clusterid = nosub_taxa_to_clusterid

    def cluster_output(self):

        curr_tree_string = self.ori_tree_string

        # get node-ids annotated for internal nodes first
        edited_tree_string = []
        prev_end = 0
        for expr in re.finditer('\)(\d+):', curr_tree_string):
            edited_tree_string.append(curr_tree_string[prev_end:expr.start()+1])
            edited_tree_string.append('[&NODE_ID=\'{}\']:'.format(expr.group(1)))
            prev_end = expr.end()
        edited_tree_string.append(curr_tree_string[prev_end:])
        curr_tree_string = ''.join(edited_tree_string)

        with open('cluster_{}.txt'.format(self.outfname), 'w') as output:

            output.write('CLUSTER\tTAXA\tSensitivity-induced (original cluster-node-id)\tSub-cluster subsumed into parent (original cluster-node-id)\r\n')
            for taxon, clusterid in self.taxon_to_clusterid.items():

                output.write('{}\t{}'.format(clusterid, re.sub("(^'|'$)", '', taxon)))
                if self.sensitivity_subsumed_taxa_to_clusterid != False:
                    output.write('\t{}'.format(self.sensitivity_subsumed_taxa_to_clusterid[taxon] if taxon in self.sensitivity_subsumed_taxa_to_clusterid else ''))
                if self.nosub_taxa_to_clusterid != False:
                    output.write('\t{}'.format(self.nosub_taxa_to_clusterid[taxon] if taxon in self.nosub_taxa_to_clusterid else ''))
                output.write('\r\n')

                taxon_start_index = curr_tree_string.find(taxon)
                if taxon_start_index < 0:
                    raise SystemExit('\r\nERROR: Problem parsing taxon ({}) in tree string.\r\n'.format(taxon))
                else:
                    curr_tree_string = curr_tree_string.replace(taxon, "'{}'[&CLUSTER='{}']".format(re.sub("(^'|'$)", "", taxon), clusterid))

        return curr_tree_string

    def figtree_output(self, modified_tree_string):

        with open('tree_{}.tre'.format(self.outfname), 'w') as output:
            output.write('#NEXUS\r\nBegin taxon;\r\n\tDimensions ntax={};\r\n\t\tTaxlabels\r\n'.format(len(self.taxon_list)))
            for taxon in self.taxon_list:
                if taxon in self.taxon_to_clusterid:
                    output.write("\t\t\t'{}_cluster{}'\r\n".format(re.sub("(^'|'$)", "", taxon), self.taxon_to_clusterid[taxon]))
                else:
                    output.write('\t\t\t{}\r\n'.format(taxon))

            output.write('\t\t\t;\r\nEnd;\r\nBegin trees;\r\n\tTranslate\r\n')

            for i, taxon in enumerate(self.taxon_list):
                if taxon in self.taxon_to_clusterid:
                    output.write("\t\t{:>4} '{}_cluster{}'{}\r\n".format(i+1, re.sub("(^'|'$)", '', taxon), self.taxon_to_clusterid[taxon], '' if i+1 == len(self.taxon_list) else ','))
                else:
                    output.write("\t\t{:>4} '{}'{}\r\n".format(i+1, re.sub("(^'|'$)", '', taxon), '' if i+1 == len(self.taxon_list) else ','))

                taxon_start_index = modified_tree_string.find("'{}'".format(re.sub("(^'|'$)", '', taxon)))
                if taxon_start_index < 0:
                    raise SystemExit('\r\nERROR: Problem parsing taxon ({}) in tree string.\r\n'.format(taxon))
                else:
                    modified_tree_string = modified_tree_string.replace("'{}'".format(re.sub("(^'|'$)", "", taxon)), str(i+1))

            output.write(';\r\ntree TREE1 = {}\r\nEnd;\r\n'.format(modified_tree_string))

    def generate_heatmap(self, ete3_tree, *args):

        for node in ete3_tree:
            if node.is_leaf():
                taxon = node.name

                tc_index = 0
                for args_index in range(0, len(args), 2):
                    taxon_to_classification = args[args_index]
                    class_to_color = args[args_index+1]

                    try:
                        classification = taxon_to_classification[taxon]
                        color = class_to_color[classification]
                    except:
                        classification = ''
                        color = 'white'

                    class_face = ete3.TextFace(classification, ftype='Arial', fsize=2, fgcolor='white')
                    class_face.background.color = color
                    class_face.hz_align = 1
                    class_face.vt_align = 1
                    class_face.margin_left = 5
                    class_face.margin_right = 5
                    node.add_face(class_face, tc_index, "aligned")
                    tc_index += 1

        return ete3_tree

    def prior_output(self, pc_input):
        with open('prior_{}.txt'.format(self.outfname), 'w') as output:
            output.write('Prior_cluster_id\tTaxon\tCurrent_cluster_id\r\n')
            for p, leaves_in_pc in pc_input.items():
                for leaf in leaves_in_pc:
                    leaf = self.leaf_node_id_to_leafname[leaf]
                    output.write('{}\t{}\t'.format(p, leaf))
                    try:
                        output.write('{}\r\n'.format(self.taxon_to_clusterid[leaf]))
                    except:
                        output.write('\r\n')

    def generate_color_scheme(self, class_list):
        from colorsys import hls_to_rgb
        import random
        random.seed(666) # maintain consistent h_list shuffle

        def hls2hex(h, l, s):
            ''' Converts from HLS to RGB '''
            return '#%02x%02x%02x' %tuple(map(lambda x: int(x*255), hls_to_rgb(h, l, s)))
        def get_color(h, s=0.5, l=0.5):
            ''' Returns a RGB color code with the specific Hue, Saturation and
            lightness.'''
            return hls2hex(h, l, s)

        # get widest increment of color hues
        h_increment = 1/(len(class_list))
        h_list = [0+(cl*h_increment) for cl in range(len(class_list))]
        random.shuffle(h_list)
        return {classification:get_color(h_list[cl]) for cl, classification in enumerate(class_list)}

    def ete3_pdf_tree_output(self, pc_input):
        output_tree = ete3.Tree(self.ori_tree_string)
        output_tree.ladderize()

        # treestyle
        ts = ete3.TreeStyle()
        ts.show_leaf_name = False

        # generate color scheme for taxon_to_clusterid
        if len(self.clusterid_to_taxa.keys()) > 1:
            clusterid_to_color = self.generate_color_scheme(self.clusterid_to_taxa.keys())
        else:
            clusterid_to_color = {self.clusterid_to_taxa.keys()[0]:'#ff2929'}

        for n, node in enumerate(output_tree.traverse(strategy='levelorder')):
            if n == 0:
                try:
                    ts.scale_length = float('{:.3f}'.format(node.get_farthest_leaf()[-1]/10))
                except:
                    pass

            if node.is_leaf():
                # color branches
                ns = ete3.NodeStyle()
                ns["size"] = 0 # no node shape

                taxon = node.name
                if taxon in self.taxon_to_clusterid:
                    clusterid = self.taxon_to_clusterid[taxon]
                    ns["hz_line_color"] = clusterid_to_color[clusterid]
                    # write taxon names aligned to the right
                    taxon_name = ete3.TextFace(taxon, ftype='Arial', fsize=2, bold=True, fgcolor=clusterid_to_color[clusterid])
                else:
                    # write taxon names aligned to the right
                    taxon_name = ete3.TextFace(taxon, ftype='Arial', fsize=2, fstyle="italic")

                node.set_style(ns)
                taxon_name.margin_left = 2
                node.add_face(taxon_name, column=0, position='branch-right')

            else:
                ns = ete3.NodeStyle()
                ns["size"] = 0 # no node shape

                # set node style
                node.set_style(ns)

        heatmap_headers = ['Cluster-ID']

        if pc_input != False:
            # prior heat-map
            heatmap_headers.append('Prior')
            taxon_to_pc = {self.leaf_node_id_to_leafname[taxon]:p for p, leaves_of_pc in pc_input.items() for taxon in leaves_of_pc}

            pc_to_color = self.generate_color_scheme(list(set(taxon_to_pc.values())))
            output_tree = self.generate_heatmap(output_tree, self.taxon_to_clusterid, clusterid_to_color, taxon_to_pc, pc_to_color)
        else:
            # normal
            output_tree = self.generate_heatmap(output_tree, self.taxon_to_clusterid, clusterid_to_color)

        # heatmap header
        for lh_index, legend_header in enumerate(heatmap_headers):
            header_face = ete3.TextFace(legend_header, ftype='Arial', fsize=2)
            header_face.hz_align = 1
            header_face.vt_align = 1
            header_face.margin_left = 5
            header_face.margin_right = 5
            ts.aligned_header.add_face(header_face, lh_index)

        # render as pdf
        output_tree.render('pdftree_{}.pdf'.format(self.outfname), tree_style=ts)