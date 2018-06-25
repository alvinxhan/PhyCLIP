#!/usr/bin/env python

from __future__ import division, print_function
import ete3
import re
import itertools
import numpy as np

from phyclip_modules.stats_utils import qn, multiple_testing_correction, summary_stats
from phyclip_modules.tree_utils import parse_newick_tree
from phyclip_modules import get_global_tree_info, phyclip_output, node_leaves_reassociation, clean_up_modules

if __name__ == '__main__':

    # parse parameters
    import argparse
    parser = argparse.ArgumentParser(description='Phylogenetic Clustering by Linear Integer Programming (PhyCLIP) v0.1')
    parser.add_argument('-i', '--input_file', required=True, type=str, help='Input file.')
    parser.add_argument('--treeinfo', type=str, help='Tree information file generated from previous PhyCLIP run.')
    parser.add_argument('--collapse_zero_branch_lengths', default=0, choices=[0, 1], type=int, help='Collapse nodes with zero branch lengths of tree prior to running PhyCLIP (default = %(default)s).')
    parser.add_argument('--equivalent_zero_length', default=0.0000001, type=float, help='Maximum branch length to be rounded to zero if the --collapse_zero_branch_lengths flag is passed (advanced option, default = %(default)s).')
    parser.add_argument('--gam_method', choices=['MAD', 'Qn'], default='Qn',  help='Method to estimate robust dispersion measure (default = %(default)s).')
    parser.add_argument('--hypo_test', choices=['Kuiper', 'KolSmi'], default='Kuiper', help='Hypothesis test to use for statistical differentiation of distance distributions (default = %(default)s).')
    parser.add_argument('--preQ', default=0, choices=[0, 1], type=int, help='Perform Benjamini-Hochberg corrections of p-values BEFORE filtering nodes that are < minimum cluster size (advanced option, default = %(default)s).')
    parser.add_argument('--subsume_sensitivity_induced_clusters', default=1, choices=[0, 1], type=int, help='Subsume cluster-size sensitivity-induced clusters into parent cluster (advanced option, default = %(default)s).')
    parser.add_argument('--sensitivity_percentile', default=25, type=int, help='Percentile of cluster size distribution under which a cluster is considered to be sensitivity-induced (advanced option, default = %(default)s%%).')
    parser.add_argument('--subsume_subclusters', default=0, choices=[0, 1], type=int, help='Subsume sub-clusters into their respective parent clusters (advanced option, default = %(default)s).')
    parser.add_argument('--solver', default='gurobi', choices=['glpk', 'gurobi'], type=str, help='Preferred ILP solver IF more than one solvers are available (default: %(default)s).')
    parser.add_argument('--solver_verbose', default=0, choices=[0, 1], type=int, help='ILP solver verbose (default: %(default)s)')
    parser.add_argument('--solver_check', action='store_true', help='Check available ILP solver(s) installed.')
    params = parser.parse_args()

    print ('{}\n\n{:^72}\n{:^72}\n\n{}'.format(''.join(['-']*72), 'Phylogenetic Clustering by Linear Integer Programming (PhyCLIP)', 'Version 0.1', ''.join(['-']*72)))

    # check solver availability
    available_solvers = []
    try:
        # check for glpsol
        import subprocess
        cmd = ['glpsol', '--help']
        subprocess.call(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        available_solvers.append('glpk')
    except:
        pass

    try:
        import gurobipy
        available_solvers.append('gurobi')
    except:
        pass

    if len(available_solvers) > 0:
        # exit if only to check available ILP solver(s) installed
        if params.solver_check:
            raise SystemExit('\nAvailable solvers...{}\n'.format(', '.join(available_solvers)))
    else:
        raise SystemExit('\nERROR: No supported solvers installed. See manual for details.\n')

    # check input file format
    infhandle = filter(None, [_.strip() for _ in open(params.input_file, 'rU')])
    # first line of input file should be path to tree file
    treepath = infhandle.pop(0)

    # filenames
    treefname = re.sub('([^/]+/|\.[^.]+$)', '', treepath)
    inputfname = re.sub('([^/]+/|\.[^.]+$)', '', params.input_file)

    try:
        newick_tree_string = parse_newick_tree(treepath)
        tree = ete3.Tree(newick_tree_string, format=5)
        print ('\nTree file...OK')
    except:
        raise SystemExit('\nERROR: Invalid tree file. Check that the correct path to the NEWICK tree file is given in the first line of the input file.\n')

    # subsequent lines are parameters (cs, fdr, gam)
    parameters = []
    for _, line in enumerate(infhandle):
        try:
            cs, fdr, gam = line.strip().split(',')
            parameters.append((int(cs), float(fdr), float(gam)))
        except:
            raise SystemExit('\nERROR: Invalid parameter set format in line {} of input file.\n'.format(_+2))
    print('Parameter sets...OK')

    # preferred solver
    if params.solver not in available_solvers:
        print ('\nWARNING: {} is not installed.'.format(params.solver))
        params.solver = available_solvers[0]

    print('\nILP solver...{}'.format(params.solver))

    # ladderize tree
    tree.ladderize()

    # collapse zero branch length
    if params.collapse_zero_branch_lengths == 1:
        from phyclip_modules.tree_utils import collapse_zero_branch_lengths
        tree = collapse_zero_branch_lengths(tree, params.equivalent_zero_length)

        # write zero branch length collapsed tree as a newick file
        treefname = 'zero-branch-length-collapsed_{}.nwk'.format(treefname)
        tree.write(format=5, outfile=treefname)
        print ('\nCollapsed newick tree file...{}'.format(treefname))

    # --- GLOBAL TREE INFORMATION --- #
    print ('\nGetting tree information...')

    # parse treeinfo file if given
    global_leaf_dist_to_node = {}
    global_leafpair_to_distance = {}
    global_nodepair_to_pval = {}

    if params.treeinfo:
        from phyclip_modules.tree_utils import parse_treeinfo_file
        try:
            global_leaf_dist_to_node, global_leafpair_to_distance, global_nodepair_to_pval = parse_treeinfo_file(params.treeinfo, params.hypo_test)
            print('Treeinfo file...OK.')
        except:
            raise SystemExit('\nERROR: Invalid treeinfo file.\n')
    else:
        params.treeinfo = '{}_treeinfo.txt'.format(treefname)

        import os
        if params.treeinfo in os.listdir(os.getcwd()):
            try:
                overwrite_treeinfo = re.search('(y|n)', raw_input('\nWARNING: {} exists in current working directory. Overwrite? (y/n): '.format(params.treeinfo))).group()
            except:
                raise SystemExit('\nERROR: Invalid input.\n')

            if overwrite_treeinfo == 'n':
                from phyclip_modules.tree_utils import parse_treeinfo_file
                try:
                    global_leaf_dist_to_node, global_leafpair_to_distance, global_nodepair_to_pval = parse_treeinfo_file(params.treeinfo, params.hypo_test)
                    print('Treeinfo file...OK.')
                except:
                    raise SystemExit('\nERROR: Invalid treeinfo file.\n')

    global_tree_info_obj = get_global_tree_info(tree, global_leaf_dist_to_node, global_leafpair_to_distance, global_nodepair_to_pval, params.treeinfo, params.hypo_test)

    global_tree_string, taxon_list, global_node_to_leaves, global_nindex_to_node, global_node_to_nindex, global_leaf_dist_to_node, global_nodepair_to_dist, global_node_to_parent_node, global_node_to_mean_child_dist2root = global_tree_info_obj.node_indexing()

    global_leafpair_to_distance, global_node_to_pwdist, global_node_to_mean_pwdist, global_node_to_ancestral_nodes, global_node_to_descendant_nodes, global_leaf_to_ancestors, global_node_to_mean_child_dist2anc = global_tree_info_obj.pwdist_dist_and_ancestral_trace(len(taxon_list), global_node_to_leaves, global_nindex_to_node, global_node_to_nindex, global_node_to_mean_child_dist2root, global_nodepair_to_dist)

    global_nodepair_to_pval = global_tree_info_obj.get_global_pval(len(global_node_to_nindex), params.hypo_test, global_node_to_leaves, global_node_to_ancestral_nodes, global_node_to_pwdist, global_leafpair_to_distance)

    # pre-calculation for within-cluster limit
    med_x = np.median(global_node_to_mean_pwdist.values())
    if params.gam_method == 'Qn':
        mad_x = qn(global_node_to_mean_pwdist.values())
    else:
        mad_x = np.median([x-med_x for x in global_node_to_mean_pwdist.values() if x >= med_x])

    ### --- GLOBAL TREE INFORMATION --- ###

    # Header of final summary stats output
    statsfname = 'summary-stats_{}.txt'.format(inputfname)
    with open(statsfname, 'w') as output:
        output.write('CS\tFDR\tMAD\tKuiper/KS\tQn/MAD\tq-values\tWithin_cluster_limit\tSolution_Index\t'
                     '#_of_clustered_sequences\tTotal_no_of_sequences\t%\t#_of_clusters\t'
                     'Mean_cluster_size\tS.D.\tMedian_cluster_size\tMAD\tMin_cluster_size\tMax_cluster_size\t'
                     'Grand_mean_of_mean_pwd\tS.D.\tGrand_mean_of_median_pwd\tS.D.\tMin_mean_pwd\tMax_mean_pwd\t'
                     'Mean_of_inter-cluster_dist\tS.D.\tMedian_of_inter-cluster_dist\tMAD\tMin_inter-cluster_dist\tMax_inter-cluster_dist\n')

    csgam_to_reassociation_memory = {} # memory of reassociated nodes/leaves for repeated (cs, gam) set

    # Start runs
    print ('\nStarting runs...')
    for p_index, (cs, fdr, gam) in enumerate(parameters):
        ### --- PARAMETER-SPECIFIC TREE INFORMATION --- ###
        print ('\n{:<26}{:^20}{:>26}'.format(''.join(['-']*25), ", ".join(['cs{}'.format(cs), 'fdr{}'.format(fdr), 'gam{}'.format(gam)]), ''.join(['-']*25)))

        # current within-cluster limit
        curr_wcl = med_x + (gam * mad_x)

        # level-order sorted list of nodes with leaves >= cs
        curr_list_of_ancestral_node = sorted([node for node, leaves in global_node_to_leaves.items() if len(leaves) >= cs])

        # reassociate subtrees and leaves based on curr_wcl
        print ('\nReassociating subtrees/leaves...')
        curr_node_to_leaves, curr_node_to_descendant_nodes, curr_node_to_mean_pwdist = [], [], [] # clear memory

        # if (cs, gam) set is repeated, check if we have performed reassociation analyses before to speed up runs
        try:
            curr_node_to_leaves, curr_node_to_descendant_nodes, curr_node_to_mean_pwdist = csgam_to_reassociation_memory[(cs, gam)]
        except:
            nla_object = node_leaves_reassociation(cs, curr_wcl, params.gam_method, curr_list_of_ancestral_node, global_node_to_leaves, global_node_to_descendant_nodes, global_node_to_mean_pwdist, global_node_to_mean_child_dist2anc, global_node_to_parent_node, global_nodepair_to_dist, global_leaf_dist_to_node, global_leaf_to_ancestors)

            curr_node_to_leaves, curr_node_to_descendant_nodes, curr_node_to_mean_pwdist = nla_object.nla_main()

            # save to memory to speed up subsequent runs with the same cs/gam parameter set
            csgam_to_reassociation_memory[(cs, gam)] = [curr_node_to_leaves, curr_node_to_descendant_nodes, curr_node_to_mean_pwdist]

        # update pairwise distance dictionaries/nodes' ancestry relations
        print ('Updating tree info...')
        curr_leaves = []
        for node, leaves in curr_node_to_leaves.items():
            # filter by cs
            if len(leaves) < cs:
                del curr_node_to_leaves[node]
                continue

            curr_leaves = list(set(curr_leaves)|set(leaves))

        curr_list_of_ancestral_node = curr_node_to_leaves.keys()[:]
        curr_node_to_descendant_nodes = {k:list(set(v)&set(curr_list_of_ancestral_node)) for k,v in curr_node_to_descendant_nodes.items() if k in curr_list_of_ancestral_node and len(v) > 0}

        # multiple-testing correction using BH procedure could be done without filtering for cs (pre) or post filtering for cs (post)
        if params.preQ == 1:
            # if pre, correction needs to be done only once
            if p_index == 0:
                print ('Multiple-testing correction (pre)...')
                curr_nodepair_to_qval = multiple_testing_correction(global_nodepair_to_pval)
        else:
            print ('Multiple-testing correction (post)...')
            curr_nodepair_to_pval = {(i,j):global_nodepair_to_pval[(i,j)] if (i,j) in global_nodepair_to_pval else global_nodepair_to_pval[(j,i)] for i,j in itertools.combinations(curr_list_of_ancestral_node, 2)}
            curr_nodepair_to_qval = multiple_testing_correction(curr_nodepair_to_pval)

        ### --- PARAMETER-SPECIFIC TREE INFORMATION --- ###

        # Build ILP model and solve
        if params.solver == 'gurobi':
            from phyclip_modules.gurobi_solver import gurobi_solver
            try:
                all_solutions = gurobi_solver(curr_node_to_leaves, curr_leaves, curr_list_of_ancestral_node, curr_nodepair_to_qval, curr_node_to_mean_pwdist, curr_wcl, cs, fdr, params.solver_verbose)
            except:
                raise SystemExit('\nERROR: Unable to solve ILP model using gurobi. You may try to use other available solvers by the --solver flag.\n')

        else:
            from phyclip_modules.glpk_solver import glpk_solver
            try:
                all_solutions = glpk_solver(curr_node_to_leaves, curr_leaves, curr_list_of_ancestral_node, curr_nodepair_to_qval, curr_node_to_mean_pwdist, curr_wcl, cs, fdr, params.solver_verbose)
            except:
                raise SystemExit('\nERROR: Unable to solve ILP model using glpk. You may try to use other available solvers by the --solver flag.\n')

        if all_solutions == 'na':
            # continue to next parameter set if no solution
            with open('cluster_{}_{}_cs{}_fdr{}_gam{}_sol0_{}'.format(params.gam_method.lower(), params.hypo_test.lower(), str(cs), str(fdr), str(gam), inputfname), 'w') as output:
                output.write('NO OPTIMAL SOLUTION FOUND.')
                print ('\nNO OPTIMAL SOLUTION FOUND.')
            continue # continue to next parameter set

        elif len(all_solutions) > 1:
            print ('\nMultiple ({}) solutions found..'.format(len(all_solutions)))

        # analyse solution and print outputs
        for sol_index, curr_taxon_to_clusterid in enumerate(all_solutions):
            curr_outfname = '{}_{}_cs{}_fdr{}_gam{}_sol{}_{}'.format(params.gam_method.lower(), params.hypo_test.lower(), str(cs), str(fdr), str(gam), sol_index, inputfname)
            curr_clusterid_to_taxa = {}
            for taxon, clusterid in curr_taxon_to_clusterid.items():
                try:
                    curr_clusterid_to_taxa[clusterid].append(taxon)
                except:
                    curr_clusterid_to_taxa[clusterid] = [taxon]

            print('\nCleaning up clusters...')
            cleanup_object = clean_up_modules(curr_node_to_descendant_nodes, global_node_to_parent_node, global_node_to_leaves, global_leafpair_to_distance, curr_node_to_leaves, curr_wcl, cs)

            # ensure that the most descendant-possible node-id is subtending each cluster
            curr_clusterid_to_taxa, curr_taxon_to_clusterid = cleanup_object.most_desc_nodeid_for_cluster(curr_clusterid_to_taxa, curr_taxon_to_clusterid)

            # remove cluster-size sensitivity-induced clusters
            if params.subsume_sensitivity_induced_clusters:
                print ('Removing cluster-size sensitivity-induced clusters...')
                curr_clusterlen_distribution = cleanup_object.get_cluster_size_distribution(curr_clusterid_to_taxa) # determine distribution of clusters
                curr_clusterid_to_taxa, curr_taxon_to_clusterid, curr_sensitivity_subsumed_taxa_to_clusterid = cleanup_object.subsume_subclusters_under_x_percentile(curr_clusterid_to_taxa, curr_taxon_to_clusterid, curr_clusterlen_distribution, params.sensitivity_percentile)

            # de-cluster outliers of descendent cluster that got clustered in ancestor cluster
            curr_clusterid_to_taxa, curr_taxon_to_clusterid = cleanup_object.decluster_outlying_taxa_clustered_to_anc_clusters(curr_clusterid_to_taxa, curr_taxon_to_clusterid)

            # delete any clusters < cs
            curr_clusterid_to_taxa, curr_taxon_to_clusterid = cleanup_object.remove_clusters_below_cs(curr_clusterid_to_taxa, curr_taxon_to_clusterid)

            # further subsume any sub-clusters into parent clusters which taxa set include that of sub-clusters
            if params.subsume_subclusters:
                print ('Subsume sub-clusters within parent...')
                curr_clusterlen_distribution = cleanup_object.get_cluster_size_distribution(curr_clusterid_to_taxa) # get cluster size distribution
                curr_clusterid_to_taxa, curr_taxon_to_clusterid, curr_nosub_taxa_to_clusterid = cleanup_object.subsume_subclusters_under_x_percentile(curr_clusterid_to_taxa, curr_taxon_to_clusterid, curr_clusterlen_distribution, 100)

            # get cluster size distribution
            curr_clusterlen_distribution = cleanup_object.get_cluster_size_distribution(curr_clusterid_to_taxa)

            # print outputs
            print ('Writing outputs...')
            output_obj = phyclip_output(global_tree_string, curr_taxon_to_clusterid, taxon_list, curr_sensitivity_subsumed_taxa_to_clusterid if params.subsume_sensitivity_induced_clusters else False, curr_nosub_taxa_to_clusterid if params.subsume_subclusters else False, curr_outfname)
            # cluster file
            curr_modified_tree_string = output_obj.cluster_output()
            # figtree annotated tree file
            output_obj.figtree_output(curr_modified_tree_string)

            # append to summary stats output file
            summary_stats(curr_clusterid_to_taxa, global_leafpair_to_distance, global_nodepair_to_dist, curr_clusterlen_distribution, statsfname, len(curr_taxon_to_clusterid), len(taxon_list), cs, fdr, gam, params.hypo_test, params.gam_method, 'pre' if params.preQ == 1 else 'post', curr_wcl, sol_index)

    print ('\n...All parameters sets analysed.\n')

