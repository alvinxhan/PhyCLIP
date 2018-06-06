#!/usr/bin/env python

from __future__ import division, print_function
from copy import deepcopy as dc
import ete3
import re
import itertools
import numpy as np

from phyclip_modules.stats_utils import qn, multiple_testing_correction, summary_stats
from phyclip_modules.tree_utils import dissociate_subtrees, dissociate_outlier_leaves, post_solver_cleanup, clean_up_sensitivity_induced_clusters
from phyclip_modules.gurobi_solver import gurobi_solver
from phyclip_modules import phyclip_output

if __name__ == '__main__':

    # parse parameters
    import argparse
    parser = argparse.ArgumentParser(description='PhyCLIP v0.1')
    parser.add_argument('-i', '--input_file', required=True, type=str, help='Input file.')
    parser.add_argument('--treeinfo', type=str, help='Master tree information file.')
    parser.add_argument('--collapse_zero_branch_lengths', default=0, choices=[0, 1], type=int, help='Collapse nodes with zero branch lengths of tree prior to running PhyCLIP (default = %(default)s).')
    parser.add_argument('--equivalent_zero_length', default=0.0000001, type=float, help='Maximum branch length to be rounded to zero if the --collapse_zero_branch_lengths flag is passed (advanced option, default = %(default)s).')

    parser.add_argument('--gam_method', choices=['MAD', 'Qn'], default='Qn',  help='Method to estimate robust dispersion measure (default = %(default)s).')
    parser.add_argument('--hypo_test', choices=['Kuiper', 'KS'], default='Kuiper', help='Hypothesis test to use for statistical differentiation of distance distributions (default = %(default)s).')
    parser.add_argument('--preQ', default=0, choices=[0, 1], type=int, help='Perform Benjamini-Hochberg corrections of p-values BEFORE filtering nodes that are < minimum cluster size (advanced option, default = %(default)s).')
    parser.add_argument('--subsume_sensitivity_induced_clusters', default=1, choices=[0, 1], type=int, help='Subsume cluster-size sensitivity-induced clusters into parent cluster (advanced option, default = %(default)s).')
    parser.add_argument('--sensitivity_percentile', default=25, type=int, help='Percentile of cluster size distribution under which a cluster is considered to be sensitivity-induced (advanced option, default = %(default)s percent).')
    parser.add_argument('--ilp_verbose', default=0, choices=[0, 1], type=int, help='ILP solver verbose (default: %(default)s)')
    params = parser.parse_args()

    # check input file format
    infhandle = filter(None, [_.strip() for _ in open(params.input_file, 'rU')])
    # first line of input file should be path to tree file
    treepath = infhandle.pop(0)

    # filenames
    treefname = re.sub('([^/]+/|\.[^.]+$)', '', treepath)
    inputfname = re.sub('([^/]+/|\.[^.]+$)', '', params.input_file)

    try:
        tree = ete3.Tree(treepath, format=5)
        print ('\ntree file...OK')
    except:
        raise TypeError('Invalid tree file given. Check that the correct path to the NEWICK tree file is given in the first line of the input file.')

    # subsequent lines are parameters (cs, fdr, gam)
    parameters = []
    for _, line in enumerate(infhandle):
        try:
            cs, fdr, gam = line.strip().split(',')
            parameters.append((int(cs), float(fdr), float(gam)))
        except:
            raise TypeError('Invalid parameter set format in line {} of input file.\n'.format(_+2))
    print('parameters...OK')

    # ladderize tree
    tree.ladderize()

    # collapse zero branch lengths if required
    if params.collapse_zero_branch_lengths == 1:
        from phyclip_modules.tree_utils import collapse_zero_branch_lengths
        tree = collapse_zero_branch_lengths(tree, params.equivalent_zero_length)

        # write zero branch length collapsed tree as a newick file
        treefname = 'zero-branch-length-collapsed_{}.nwk'.format(treefname)
        tree.write(format=5, outfile=treefname)
        print ('\nnewick tree file with zero branch legnth collapsed...{}'.format(treefname))

    # --- GLOBAL TREE INFORMATION --- #
    print ('\ngetting tree information...')

    # parse treeinfo file if given
    if params.treeinfo:
        from phyclip_modules.tree_utils import parse_treeinfo_file
        try:
            global_leaf_dist_to_node, global_leafpair_to_distance, global_nodepair_to_pval = parse_treeinfo_file(params.treeinfo, params.hypo_test)
            print('master treeinfo file...OK.')
        except:
            raise TypeError('Invalid treeinfo file given.')
    else:
        params.treeinfo = '{}_treeinfo.txt'.format(treefname)
        global_leaf_dist_to_node = {}
        global_leafpair_to_distance = {}
        global_nodepair_to_pval = {}

    from phyclip_modules import get_global_tree_info
    global_tree_info_obj = get_global_tree_info(tree, global_leaf_dist_to_node, global_leafpair_to_distance, global_nodepair_to_pval, params.treeinfo)

    global_tree_string, taxon_list, global_node_to_leaves, global_nindex_to_node, global_node_to_nindex, global_leaf_dist_to_node, global_nodepair_to_dist = global_tree_info_obj.node_indexing()

    global_leafpair_to_distance, global_node_to_pwdist, global_node_to_mean_pwdist, global_node_to_ancestral_nodes, global_node_to_descendant_nodes = global_tree_info_obj.pwdist_dist_and_ancestral_trace(len(taxon_list), global_node_to_leaves, global_nindex_to_node, global_node_to_nindex)

    global_nodepair_to_pval = global_tree_info_obj.get_global_pval(len(global_node_to_nindex), params.hypo_test, global_node_to_leaves, global_node_to_ancestral_nodes, global_node_to_pwdist, global_leafpair_to_distance)

    ### --- GLOBAL TREE INFORMATION --- ###

    # Header of final summary stats output
    statsfname = 'summary-stats_{}.txt'.format(inputfname)
    with open(statsfname, 'w') as output:
        output.write('CS\tFDR\tMAD\tKuiper/KS\tQn/MAD\tq-values\tWithin_cluster_limit\t'
                     '#_of_clustered_sequences\tTotal_no_of_sequences\t%\t#_of_clusters\t'
                     'Mean_cluster_size\tS.D.\tMedian_cluster_size\tMAD\tMin_cluster_size\tMax_cluster_size\t'
                     'Grand_mean_of_mean_pwd\tS.D.\tGrand_mean_of_median_pwd\tS.D.\tMin_mean_pwd\tMax_mean_pwd\t'
                     'Mean_of_inter-cluster_dist\tS.D.\tMedian_of_inter-cluster_dist\tMAD\tMin_inter-cluster_dist\tMax_inter-cluster_dist\n')

    # Start runs
    print ('\nstarting runs...')
    for p_index, (cs, fdr, gam) in enumerate(parameters):
        ### --- PARAMETER-SPECIFIC TREE INFORMATION --- ###
        print ('\ncs{}, fdr{}, gam{}...'.format(cs, fdr, gam))

        # calculate within-cluster limit for given gamma
        med_x = np.median(global_node_to_mean_pwdist.values())
        if params.gam_method == 'Qn':
            mad_x = qn(global_node_to_mean_pwdist.values())
        else:
            mad_x = np.median([x-med_x for x in global_node_to_mean_pwdist.values() if x >= med_x])
        curr_wcl = med_x + (gam * mad_x)

        # remove nodes with leaf counts < cs from consideration
        curr_node_to_leaves = dc({node:leaves for node, leaves in global_node_to_leaves.items() if len(leaves) >= cs})

        # list of nodes
        curr_list_of_ancestral_node = curr_node_to_leaves.keys()[:]

        # multiple-testing correction using BH procedure could be done without filtering for cs (pre) or post filtering for cs (post)
        if params.preQ == 1:
            # if pre, correction needs to be done only once
            if p_index == 0:
                print ('multiple-testing correction (pre)...')
                curr_nodepair_to_qval = multiple_testing_correction(global_nodepair_to_pval)
        else:
            print ('multiple-testing correction (post)...')
            curr_nodepair_to_pval = {(i,j):global_nodepair_to_pval[(i,j)] if (i,j) in global_nodepair_to_pval else global_nodepair_to_pval[(j,i)] for i,j in itertools.combinations(curr_list_of_ancestral_node, 2)}
            curr_nodepair_to_qval = multiple_testing_correction(curr_nodepair_to_pval)

        # update pairwise distance dictionaries
        curr_node_to_pwdist = dc({node:pwdist for node, pwdist in global_node_to_pwdist.items() if node in curr_list_of_ancestral_node})
        curr_node_to_mean_pwdist = dc({node:mean_pwdist for node, mean_pwdist in global_node_to_mean_pwdist.items() if node in curr_list_of_ancestral_node})

        # update nodes' ancestry relations
        curr_node_to_ancestral_nodes = dc({k:list(set(v)&set(curr_list_of_ancestral_node)) for k, v in global_node_to_ancestral_nodes.items() if k in curr_list_of_ancestral_node})
        curr_node_to_descendant_nodes = dc({k:list(set(v)&set(curr_list_of_ancestral_node)) for k, v in global_node_to_descendant_nodes.items() if k in curr_list_of_ancestral_node})

        # dissociate subtrees
        print ('dissociating potential clusters...')
        curr_node_to_leaves, curr_node_to_ancestral_nodes, curr_node_to_descendant_nodes, subtree_nodes = dissociate_subtrees(curr_list_of_ancestral_node, curr_node_to_leaves, curr_node_to_mean_pwdist, curr_node_to_ancestral_nodes, curr_node_to_descendant_nodes, curr_nodepair_to_qval, curr_wcl, global_node_to_ancestral_nodes, global_node_to_descendant_nodes, global_node_to_leaves)

        # update pairwise distance dictionaries/nodes' ancestry relations
        print ('updating tree info...')
        curr_node_to_pwdist = {}
        for node, leaves in curr_node_to_leaves.items():
            # filter by cs
            if len(leaves) < cs:
                del curr_node_to_leaves[node]
                continue
            # if leaf membership of node has changed, recalculate the pairwise distance distribution
            if set(leaves) != set(global_node_to_leaves[node]):
                curr_node_to_pwdist[node] = sorted([global_leafpair_to_distance[(x,y)] for (x,y) in itertools.combinations(curr_node_to_leaves[node], 2)])
            else:
                curr_node_to_pwdist[node] = global_node_to_pwdist[node][:]

        curr_list_of_ancestral_node = curr_node_to_leaves.keys()[:]
        curr_node_to_mean_pwdist = {node:np.mean(pwdist) for node, pwdist in curr_node_to_pwdist.items()}
        curr_node_to_ancestral_nodes = {k:list(set(v)&set(curr_list_of_ancestral_node)) for k,v in curr_node_to_ancestral_nodes.items() if k in curr_list_of_ancestral_node and len(v) > 0}
        curr_node_to_descendant_nodes = {k:list(set(v)&set(curr_list_of_ancestral_node)) for k,v in curr_node_to_descendant_nodes.items() if k in curr_list_of_ancestral_node and len(v) > 0}

        # get leaf to ancestor nodes dictionary
        curr_leaf_to_ancestors = {}
        for node, leaves in curr_node_to_leaves.items():
            for leaf in leaves:
                try:
                    curr_leaf_to_ancestors[leaf].append(node)
                except:
                    curr_leaf_to_ancestors[leaf] = [node]

        # dissociating leaf outliers
        print ('dissociating leaf outliers...')
        curr_node_to_leaves, curr_node_to_ancestral_nodes, curr_node_to_descendant_nodes, curr_nodes_with_leaf_membership_changes = dissociate_outlier_leaves(curr_node_to_leaves, curr_node_to_mean_pwdist, curr_wcl, global_leaf_dist_to_node, curr_node_to_descendant_nodes, curr_node_to_ancestral_nodes, global_nodepair_to_dist, global_node_to_descendant_nodes, global_node_to_ancestral_nodes, curr_list_of_ancestral_node, curr_leaf_to_ancestors)

        # update pairwise distance dictionaries/nodes' ancestry relations again
        print ('updating tree info...')
        curr_node_to_leaves = {k:v for k,v in curr_node_to_leaves.items() if len(v) >= cs}
        curr_list_of_ancestral_node = curr_node_to_leaves.keys()

        # recalculate pairwise distance distribution for nodes with leaf membership changes
        for node in list(set(curr_nodes_with_leaf_membership_changes)&set(curr_list_of_ancestral_node)):
            # no need to sort anymore
            curr_node_to_pwdist[node] = sorted([global_leafpair_to_distance[(x,y)] for (x,y) in itertools.combinations(curr_node_to_leaves[node], 2)])
            curr_node_to_mean_pwdist[node] = np.mean(curr_node_to_pwdist[node])

        curr_node_to_pwdist = {k:v for k,v in curr_node_to_pwdist.items() if k in curr_list_of_ancestral_node}
        curr_node_to_mean_pwdist = {k:v for k,v in curr_node_to_mean_pwdist.items() if k in curr_list_of_ancestral_node}

        curr_node_to_ancestral_nodes = {k:list(set(v)&set(curr_list_of_ancestral_node)) for k,v in curr_node_to_ancestral_nodes.items() if k in curr_list_of_ancestral_node and len(v) > 0}
        curr_node_to_descendant_nodes = {k:list(set(v)&set(curr_list_of_ancestral_node)) for k,v in curr_node_to_descendant_nodes.items() if k in curr_list_of_ancestral_node and len(v) > 0}
        curr_leaf_to_ancestors = {}
        for node, leaves in curr_node_to_leaves.items():
            for leaf in leaves:
                try:
                    curr_leaf_to_ancestors[leaf].append(node)
                except:
                    curr_leaf_to_ancestors[leaf] = [node]

        ### --- PARAMETER-SPECIFIC TREE INFORMATION --- ###

        # Build ILP model and solve
        try:
            # try with Gurobi first
            curr_taxon_to_clusterid = gurobi_solver(curr_node_to_leaves, curr_leaf_to_ancestors, curr_list_of_ancestral_node, curr_nodepair_to_qval, curr_node_to_mean_pwdist, curr_wcl, cs, fdr, params.ilp_verbose)
        except:
            # try GLPK
            raise TypeError('No ILP solvers available.')

        curr_outfname = '{}_{}_cs{}_fdr{}_gam{}_{}'.format(params.gam_method.lower(), params.hypo_test.lower(), str(cs), str(fdr), str(gam), inputfname)
        if curr_taxon_to_clusterid == 'na':
            # continue to next parameter set if no solution
            with open('cluster_{}.txt'.format(curr_outfname), 'w') as output:
                output.write('NO SOLUTION FOUND!\n')
            continue

        curr_clusterid_to_taxa = {}
        for taxon, clusterid in curr_taxon_to_clusterid.items():
            try:
                curr_clusterid_to_taxa[clusterid].append(taxon)
            except:
                curr_clusterid_to_taxa[clusterid] = [taxon]

        # post-ILP solver clean up
        print('cleaning up clusters...')
        curr_clusterid_to_taxa, curr_taxon_to_clusterid, curr_taxon_to_supercluster = post_solver_cleanup(curr_clusterid_to_taxa, curr_taxon_to_clusterid, global_node_to_leaves, global_leaf_dist_to_node, curr_node_to_leaves, global_node_to_descendant_nodes)

        # clean up sensitivity-induced clusters
        if params.subsume_sensitivity_induced_clusters == 0:
            curr_sensitivity_subsumed_taxa_to_clusterid = {}
        else:
            curr_clusterid_to_taxa, curr_taxon_to_clusterid, curr_sensitivity_subsumed_taxa_to_clusterid = clean_up_sensitivity_induced_clusters(curr_clusterid_to_taxa, curr_taxon_to_clusterid, params.sensitivity_percentile, global_node_to_descendant_nodes, global_leafpair_to_distance, curr_wcl)

        # remove any clusters < min cluster size and determine distribution of clusters
        curr_clusterlen_to_frequency = {}
        for clusterid, taxa in curr_clusterid_to_taxa.items():
            clusterlen = len(taxa)
            if clusterlen < cs:
                for taxon in taxa:
                    del curr_taxon_to_clusterid[taxon]
                del curr_clusterid_to_taxa[clusterid]
            else:
                try:
                    curr_clusterlen_to_frequency[clusterlen] += 1
                except:
                    curr_clusterlen_to_frequency[clusterlen] = 1

        curr_clusterlen_distribution = [i for j in [[clusterlen]*frequency for clusterlen, frequency in curr_clusterlen_to_frequency.items()] for i in j]

        # print outputs
        print ('writing outputs...')

        output_obj = phyclip_output(global_tree_string, curr_taxon_to_clusterid, taxon_list, curr_taxon_to_supercluster, curr_sensitivity_subsumed_taxa_to_clusterid, curr_outfname)
        # cluster file
        curr_modified_tree_string = output_obj.cluster_output()
        # figtree annotated tree file
        output_obj.figtree_output(curr_modified_tree_string)

        # append to summary stats output file
        summary_stats(curr_clusterid_to_taxa, global_leafpair_to_distance, global_nodepair_to_dist, curr_clusterlen_distribution, statsfname, len(curr_taxon_to_clusterid), len(taxon_list), cs, fdr, gam, params.hypo_test, params.gam_method, 'pre' if params.preQ == 1 else 'post', curr_wcl)


    print ('\n...all parameters sets analysed.\n')

