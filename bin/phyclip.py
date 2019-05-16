#!/usr/bin/env python

# Phylogenetic Clustering by Linear Integer Programming (PhyCLIP)
# Authors: Alvin X. Han and Edyth Parker

from __future__ import division
from copy import deepcopy as dc
import ete3
import re
import os
import itertools
import subprocess
import argparse
from multiprocessing import Pool
import numpy as np
import time

import pyximport; pyximport.install()

from phyclip_modulex.pyutilx import parse_newick_tree, parse_input_parameters, clean_up_modules, get_cluster_size_distribution, output_utils
from phyilpx_phyclip import phyilpx_treeinfo, distal_dissociation
from phyilpx_stats import qn, summary_stats

if __name__ == '__main__':
    # parse parameters
    version = 2.0
    parser = argparse.ArgumentParser(description='Phylogenetic Clustering by Linear Integer Programming (PhyCLIP) v{}'.format(version))

    required_args = parser.add_argument_group('Required')
    required_args.add_argument('-i', '--input_file', type=str, help='Input file. See manual for format details.')

    analyses_aid = parser.add_argument_group('Analysis aids')
    analyses_aid.add_argument('--treeinfo', type=str, default=False, help='*_treeinfo.txt file containing between-nodes p-values generated from PREVIOUS PhyCLIP run of the SAME phylogenetic tree.')
    analyses_aid.add_argument('--no_treeinfo', action='store_true', default=False, help='Stop PhyCLIP from generating *_treeinfo.txt file.')
    analyses_aid.add_argument('--pdf_tree', action='store_true', help='PDF tree output annotated with cluster results.')
    analyses_aid.add_argument('--optimise', choices=['intermediate', 'high'], help='PhyCLIP automatically searches for the clustering result from the OPTIMAL input parameter set. See documentation for details. Must have >1 input parameter set in input file.')

    analyses_options = parser.add_argument_group('Analysis options')
    #analyses_options.add_argument('--mode', choices=['genetic', 'geneticX'], default='genetic', help='PhyCLIP clustering mode. \'genetic\' - aims to identify genetic clusters; \'geneticX\' - genetic clustering with heuristics shortcuts for faster analysis; default = %(default)s')
    analyses_options.add_argument('--tree_outgroup', type=str, default=False, help='Taxon (name as appeared in tree) to be set as outgroup for rooting OR \'midpoint\' for mid-point rooting.')
    analyses_options.add_argument('--collapse_zero_branch_length', default=0, choices=[0, 1], type=int, help='Collapse internal nodes with zero branch length of tree before running PhyCLIP (default = %(default)s).')
    analyses_options.add_argument('--equivalent_zero_length', default=0.000001, type=np.float32, help='Maximum branch length to be rounded to zero if the --collapse_zero_branch_length flag is passed (advanced option, default = %(default)s).')

    analyses_options.add_argument('--gam_method', choices=['mad', 'qn'], default='qn',  help='Method to estimate robust dispersion measure (default = %(default)s).')
    analyses_options.add_argument('--hypo_test', choices=['kuiper', 'kolsmi'], default='kuiper', help='Hypothesis test to use for statistical differentiation of distance distributions (default = %(default)s).')
    analyses_options.add_argument('--preQ', default=0, choices=[0, 1], type=int, help='Perform Benjamini-Hochberg corrections of p-values BEFORE filtering nodes that are < minimum cluster size (advanced option, default = %(default)s).')

    cleanup_options = parser.add_argument_group('Clean-up options')
    cleanup_options.add_argument('--subsume_sensitivity_induced_clusters', default=1, choices=[0, 1], type=int, help='Subsume cluster-size sensitivity-induced clusters into parent cluster (default = %(default)s).')
    cleanup_options.add_argument('--sensitivity_percentile', default=25, type=int, help='Percentile of cluster size distribution under which a cluster is considered to be sensitivity-induced (advanced option, default = %(default)s%%).')
    cleanup_options.add_argument('--subsume_subclusters', default=0, choices=[0, 1], type=int, help='Subsume sub-clusters into their respective parent clusters (default = %(default)s).')
    cleanup_options.add_argument('--force', default=0, choices=[0,1], type=int, help='Force-cluster low information putative clusters.')

    prior_options = parser.add_argument_group('Prior options')
    prior_options.add_argument('--prior', type=str, help='Prior information on clustered taxa. File format same as PhyCLIP cluster text file output (i.e. CLUSTER-ID{tab}SEQUENCE-NAME). Only works with gurobi solver, no support for glpk.')
    prior_options.add_argument('--prior_weights', type=str, help='Weights on importance/confidence between prior clusters to remain clustered together. File format: CLUSTER-ID{tab}INT/FLOAT_WEIGHT{\\newline}. Equal weights will be assumed if none is given.')

    solver_options = parser.add_argument_group('Solver options')
    solver_options.add_argument('--solver', default='gurobi', choices=['glpk', 'gurobi'], type=str, help='Preferred ILP solver IF more than one solvers are available (default: %(default)s).')
    solver_options.add_argument('--solver_verbose', default=0, choices=[0, 1], type=int, help='ILP solver verbose (default: %(default)s)')
    solver_options.add_argument('--solver_check', action='store_true', help='Check available ILP solver(s) installed.')
    solver_options.add_argument('--threads', type=int, help='Number of threads to use for Gurobi solver (default = all available nodes).')

    params = parser.parse_args()

    print ('{}\n\n{:^72}\n{:^72}\n\n{}'.format(''.join(['-']*72), 'Phylogenetic Clustering by Linear Integer Programming (PhyCLIP)', 'v{}-dev'.format(version), ''.join(['-']*72)))

    # check solver availability
    available_solvers = {}
    try:
        # check for glpsol
        cmd = ['glpsol', '--version']
        solver_version = 'glpk_{}'.format(re.search('v\d+\.\d+', subprocess.check_output(cmd)).group())
        available_solvers['glpk'] = solver_version
    except:
        pass

    try:
        from gurobipy import gurobi
        solver_version = 'gurobi_v{}'.format('.'.join(map(str, gurobi.version())))
        available_solvers['gurobi'] = solver_version
    except:
        pass

    if len(available_solvers) > 0:
        # exit if only to check available ILP solver(s) installed
        if params.solver_check:
            print ('\nAvailable solvers...{}\n'.format(', '.join(available_solvers.values())))
            exit(0)
    else:
        raise Exception('\nNo supported solvers installed. See manual for details on how to download and install supported ILP solvers.\n')

    # preferred solver
    if params.solver not in available_solvers:
        print ('\nWARNING: {} is not installed.'.format(params.solver))
        params.solver = available_solvers.keys()[0]

    print('\nILP solver...{}'.format(available_solvers[params.solver]))

    # limit number of threads for gurobi solver
    try:
        ncpu = int(params.threads)
    except:
        import multiprocessing as mp
        ncpu = mp.cpu_count()

    # check input file format
    try:
        infhandle = filter(None, [_.strip() for _ in open(params.input_file, 'rU')])
    except:
        raise Exception('\nMissing input file.\n')
    # first line of input file should be path to tree file
    treepath = infhandle.pop(0)

    # filenames
    treefname = re.sub('([^/]+/|\.[^.]+$)', '', treepath)
    inputfname = re.sub('([^/]+/|\.[^.]+$)', '', params.input_file)

    # parse and check newick tree file
    try:
        newick_tree_string = parse_newick_tree(treepath)
    except:
        raise Exception('\nInvalid tree file.\n')

    # subsequent lines are parameters (cs, fdr, gam)
    parameters, max_gam = parse_input_parameters(infhandle)

    # read prior file (if given)
    if params.prior:
        prfhandle = filter(None, [_.strip() for _ in open(params.prior, 'rU')])

        prior_input = {}
        for line in prfhandle:
            try:
                clusterid, taxon = line.split('\t')[:2]
            except:
                raise Exception('\nInvalid prior file given. Check manual or type --help for details on correct file format.\n')

            if clusterid == 'CLUSTER' and taxon == 'TAXA':
                continue
            else:
                taxon = "'{}'".format(re.sub("(^'|'$)", "", taxon))
                try:
                    prior_input[clusterid].append(taxon)
                except:
                    prior_input[clusterid] = [taxon]

        # prior weights
        if params.prior_weights:
            prfhandle = filter(None, [_.strip() for _ in open(params.prior, 'rU')])

            prior_weights = {}
            for line in prfhandle:
                try:
                    clusterid, weight = line.split('\t')
                except:
                    raise Exception('\nInvalid prior weights file given. Check manual or type --help for details on correct file format.\n')

                if clusterid in prior_input:
                    prior_weights[clusterid] = float(weight)
                else:
                    raise Exception('\nCLUSTER-ID {} in prior weights file was not found in prior file.\n'.format(clusterid))

            prior_weights = {clusterid:weight/(sum(prior_weights.values())) for clusterid, weight in prior_weights.items()}
        else:
            prior_weights = {clusterid:1/len(taxa) for clusterid, taxa in prior_input.items()}

        print('Prior file...OK')
    else:
        prior_input = False
        prior_weights = False

    # print mode
    #print ('\nMode...{}'.format(params.mode))

    statsfname = 'summary-stats_{}.txt'.format(inputfname)
    if statsfname in os.listdir(os.getcwd()):
        try:
            overwrite_statsfile = re.search('(y|n)', raw_input('\nWARNING: SUMMARY STATS FILE {} exists in current working directory. Overwrite? (y/n): '.format(statsfname))).group()
        except:
            raise Exception('\nInvalid input.\n')
    else:
        overwrite_statsfile = 'y'

    # write header of summary stats output
    if overwrite_statsfile == 'y':
        with open(statsfname, 'w') as output:
            output.write('Treefile\tCS\tFDR\tMAD\tKuiper/KS\tQn/MAD\tq-values\tForce\tWithin_cluster_limit\tSolution_Index\tPrior_taxa_clustered(%)\t'
                         '#_of_clustered_sequences\tTotal_no_of_sequences\t%_clustered\t#_of_clusters\t'
                         'Mean_cluster_size\tS.D.\tMedian_cluster_size\tMAD\tMin_cluster_size\tMax_cluster_size\t'
                         'Grand_mean_of_mean_pwd\tS.D.\tGrand_mean_of_median_pwd\tS.D.\tMin_mean_pwd\tMax_mean_pwd\t'
                         'Mean_of_inter-cluster_dist\tS.D.\tMedian_of_inter-cluster_dist\tMAD\tMin_inter-cluster_dist\tMax_inter-cluster_dist\t'
                         'solver(version)\n')

    # --- GLOBAL TREE INFORMATION --- #
    print ('\nGetting tree information...')

    phyilpx_obj = phyilpx_treeinfo(newick_tree_string, treefname, params.tree_outgroup, params.collapse_zero_branch_length, params.equivalent_zero_length)
    global_leaf_node_id_to_leafname, global_internal_nodes, original_tree_string, treefname = phyilpx_obj.properties()

    # if cluster prior is given, check that all stated taxa in prior file are present in tree
    if params.prior:
        global_leafname_to_leaf_node_id = {v:k for k, v in global_leaf_node_id_to_leafname.items()}
        prior_input = {clusterid:[global_leafname_to_leaf_node_id[taxon] for taxon in taxa] for clusterid, taxa in prior_input.items()}

        taxa_in_prior_not_present_in_tree = list(set([i for j in prior_input.values() for i in j]) - set(global_leaf_node_id_to_leafname.keys()))
        if len(taxa_in_prior_not_present_in_tree) > 0:
            raise Exception('\nSome taxa in prior file are not present in tree - {}.\n'.format(', '.join([global_leaf_node_id_to_leafname[leaf] for leaf in taxa_in_prior_not_present_in_tree])))

    # check treeinfo file
    if params.no_treeinfo == False: # to print treeinfo file
        if params.treeinfo == False: # no treeinfo file given
            params.treeinfo = '{}_treeinfo.npz'.format(treefname)
            # check if a copy of treeinfo file is already in current working directiory
            if params.treeinfo in os.listdir(os.getcwd()):
                raise Exception('\nA file with the same name as the expected treeinfo file ({}) exists in current working directory. Either remove the file or input it as argument using --treeinfo flag.\n'.format(params.treeinfo))

    # get pairwise distance array (nodepair = all nodes including leaves, leafpair = structured, just leaves)
    global_nodepair_to_dist = phyilpx_obj.get_nodepair_distance()

    # get structured array of leaf dist to node / array of node to leaves (reverse sorted by distance to node)
    global_leaf_dist_to_node, global_node_to_leaves = phyilpx_obj.get_leaf_dist_to_node()

    # structured array of node to parent node
    global_node_to_parent_node = phyilpx_obj.get_node_to_parent_node()

    # ancestry relations/global_node_to_mean_pwdist are python dictionaries
    # global_node_to_mean_child_dist2anc = np.array(N,N)
    global_node_to_descendant_nodes, global_leaf_to_ancestors, global_node_to_mean_child_dist2anc, global_node_to_mean_pwdist = phyilpx_obj.get_ancestral_relations()

    # get global p-values
    global_nodepair_to_pval = phyilpx_obj.get_global_pval({'kuiper':1, 'kolsmi':0}[params.hypo_test], params.treeinfo, params.no_treeinfo)

    # pre-calculation for within-cluster limit
    med_x = np.median(global_node_to_mean_pwdist.values())
    if params.gam_method == 'qn':
        mad_x = qn([x for x in global_node_to_mean_pwdist.values() if x >= med_x])
    else:
        mad_x = np.median([x-med_x for x in global_node_to_mean_pwdist.values() if x >= med_x])

    # Start runs
    print ('\nStarting runs...')
    optimise_dict = {}
    multi_cs_gam_tuples = [(cs, gam) for (cs, fdr, gam) in parameters]
    if len(set(multi_cs_gam_tuples)) < len(multi_cs_gam_tuples): # if there are non-unique
        csgam_to_reassociation_memory = {} # memory of reassociated nodes/leaves for repeated (cs, gam) set
        multi_cs_gam_tuples = [_ for _ in set(multi_cs_gam_tuples) if multi_cs_gam_tuples.count(_) > 1]
    else:
        multi_cs_gam_tuples = []

    for p_index, (cs, fdr, gam) in enumerate(parameters):
        ### --- PARAMETER-SPECIFIC TREE INFORMATION --- ###
        print ('\n{:<26}{:^20}{:>26}'.format(''.join(['-']*25), ", ".join(['cs{:d}'.format(cs), 'fdr{}'.format(str(fdr)), 'gam{}'.format(str(gam))]), ''.join(['-']*25)))

        # current within-cluster limit
        curr_wcl = med_x + (gam * mad_x)

        # level-order sorted list of nodes with leaves >= cs (only reassociate such nodes)
        curr_list_of_ancestral_node = np.array([node for node, leaves in global_node_to_leaves.items() if len(leaves) >= cs], dtype=np.int64)

        # reassociate subtrees and leaves based on curr_wcl
        # if (cs, gam) set is repeated, check if we have performed reassociation analyses before to speed up runs
        try:
            curr_node_to_leaves, curr_node_to_descendant_nodes, curr_node_to_mean_pwdist = csgam_to_reassociation_memory[(cs, gam)]
        except:
            dd_obj = distal_dissociation(cs, curr_wcl, params.gam_method, curr_list_of_ancestral_node, global_node_to_leaves, global_node_to_descendant_nodes, global_node_to_mean_pwdist, global_node_to_mean_child_dist2anc, global_node_to_parent_node, global_nodepair_to_dist, global_leaf_dist_to_node, global_leaf_to_ancestors, global_leaf_node_id_to_leafname)

            curr_node_to_leaves, curr_node_to_descendant_nodes, curr_node_to_mean_pwdist = dd_obj.dd_main()

            for node, leaves in curr_node_to_leaves.items():
                # filter by cs
                if len(leaves) < cs or node not in curr_node_to_mean_pwdist:
                    del curr_node_to_leaves[node]
                    continue

            """
            # ! --- addition --- #
            # get nodes <= curr_wcl
            nodes_below_wcl = [node for node in curr_node_to_leaves.keys() if curr_node_to_mean_pwdist[node] <= curr_wcl]
            upper_tol = np.float32(-1.)

            for node in nodes_below_wcl:
                leaf_pairs = list(itertools.combinations(curr_node_to_leaves[node], 2))
                _pwdist = np.zeros(len(leaf_pairs), dtype = np.float32)

                for _, (i,j) in enumerate(leaf_pairs):
                    _pwdist[_] = global_nodepair_to_dist[(i,j)]

                med_pwdist = np.median(_pwdist)
                rmad_pwdist = np.median([abs(_ - med_pwdist) for _ in _pwdist if _ >= med_pwdist])

                max_pwdist = med_pwdist + 2*rmad_pwdist
                if max_pwdist > upper_tol:
                    upper_tol = max_pwdist

            # homogeneity check
            for node, leaves in curr_node_to_leaves.items():
                try:
                    curr_nodes_and_leaves = set(curr_node_to_descendant_nodes[node])|set(leaves)
                except:
                    curr_nodes_and_leaves = set(leaves)

                children_of_node = list(set(phyilpx_obj.get_children(node))&(curr_nodes_and_leaves))
                if len(children_of_node) < 2:
                    continue

                nodes_to_separate = []
                for child_i, child_j in itertools.combinations(children_of_node, 2):
                    if child_i in global_leaf_to_ancestors:
                        dist_i = global_nodepair_to_dist[(child_i, node)]
                    else:
                        dist_i = global_node_to_mean_child_dist2anc[(child_i, node)]

                    if child_j in global_leaf_to_ancestors:
                        dist_j = global_nodepair_to_dist[(child_j, node)]
                    else:
                        dist_j = global_node_to_mean_child_dist2anc[(child_j, node)]

                    if dist_i + dist_j > upper_tol:
                        nodes_to_separate += [child_i, child_j]

                if len(nodes_to_separate):
                    for child in set(nodes_to_separate):
                        if nodes_to_separate.count(child) == len(children_of_node) - 1:
                            #print node, child

                            if child in global_leaf_to_ancestors: # child is leaf
                                curr_node_to_leaves[node] = np.array(list(set(curr_node_to_leaves[node]) - set([child])), dtype = np.int32)
                            else:
                                try:
                                    curr_node_to_leaves[node] = np.array(list(set(curr_node_to_leaves[node]) - set(curr_node_to_leaves[child])), dtype = np.int32)
                                except:
                                    pass

                                try:
                                    curr_node_to_descendant_nodes[node] = list(set(curr_node_to_descendant_nodes[node]) - set(curr_node_to_descendant_nodes[child]))
                                except:
                                    pass

                    if len(curr_node_to_leaves[node]) < cs:
                        del curr_node_to_leaves[node]
                    else:
                        curr_node_to_mean_pwdist[node] = np.mean([global_nodepair_to_dist[(i,j)] for i, j in itertools.combinations(curr_node_to_leaves[node], 2)])
            # ! --- addition --- #
            """
            # save to memory to speed up subsequent runs with the same cs/gam parameter set
            if (cs, gam) in multi_cs_gam_tuples:
                csgam_to_reassociation_memory[(cs, gam)] = [curr_node_to_leaves, curr_node_to_descendant_nodes, curr_node_to_mean_pwdist]

        # update pairwise distance dictionaries/nodes' ancestry relations
        print ('Updating tree info...')
        curr_leaves = [x for y in curr_node_to_leaves.values() for x in y]
        # level-order sorted list of nodes with leaves >= cs (only reassociate such nodes)
        curr_list_of_ancestral_node = np.sort(np.array(curr_node_to_leaves.keys(), dtype=np.int64))
        curr_node_to_descendant_nodes = {k:list(set(v)&set(curr_list_of_ancestral_node)) for k,v in curr_node_to_descendant_nodes.items() if k in curr_list_of_ancestral_node and len(v) > 0}

        # multiple-testing correction using BH procedure could be done without filtering for cs (pre) or post filtering for cs (post)
        if params.preQ == 0:
            print ('Multiple-testing correction (post)...')
            curr_nodepair_to_qval = phyilpx_obj.multiple_testing_correction(curr_list_of_ancestral_node)
        else:
            # if pre, correction needs to be done only once
            if p_index == 0:
                print ('Multiple-testing correction (pre)...')
                curr_nodepair_to_qval = phyilpx_obj.multiple_testing_correction()

        ### --- PARAMETER-SPECIFIC TREE INFORMATION --- ###

        # Build ILP model and solve
        if params.solver == 'gurobi':
            from phyclip_modulex.gurobi_solverx import gurobi_solver
            all_solutions = gurobi_solver(curr_node_to_leaves, curr_leaves, curr_list_of_ancestral_node, curr_nodepair_to_qval, curr_node_to_mean_pwdist, curr_wcl, cs, fdr, params.solver_verbose, p_index, ncpu, prior_input, prior_weights, curr_node_to_descendant_nodes)

        else:
            from phyclip_modulex.glpk_solverx import glpk_solver
            all_solutions = glpk_solver(curr_node_to_leaves, curr_leaves, curr_list_of_ancestral_node, curr_nodepair_to_qval, curr_node_to_mean_pwdist, curr_wcl, cs, fdr, params.solver_verbose)

        if all_solutions == 'na':
            print ('\nWARNING: No optimal solution found for paramter set ({:d}, {}, {}).'.format(cs, str(fdr), str(gam)))
            if params.solver == 'glpk':
                print ('{:<9}You might want to use GUROBI solver (if available) instead.'.format(''))
            # continue to next parameter set if no solution
            with open(statsfname, 'a') as output:
                output.write('{}\t{:d}\t{}\t{}\tNO OPTIMAL SOLUTION FOUND.\n'.format(treefname, cs, fdr, gam))
            continue # continue to next parameter set

        elif len(all_solutions) > 1:
            print ('\nMultiple ({}) solutions found..'.format(len(all_solutions)))

        # analyse solution and print outputs
        for sol_index, curr_taxon_to_clusterid in all_solutions.items():
            # failed integrality
            if curr_taxon_to_clusterid == False:
                with open(statsfname, 'a') as output:
                    output.write('{}\t{}\t{}\t{}\tNO OPTIMAL SOLUTION FOUND (SOLUTION-INDEX: {}, FAILED INTEGRALITY).\n'.format(treefname, cs, fdr, gam, sol_index))
                continue

            curr_outfname = '{}_{}_cs{}_fdr{}_gam{}_sol{}{}_f{}_{}'.format(params.gam_method.lower(), params.hypo_test.lower(), str(cs), str(fdr), str(gam), sol_index, '_prior' if params.prior else '', params.force, treefname)

            curr_clusterid_to_taxa = {}
            for taxon, clusterid in curr_taxon_to_clusterid.items():
                try:
                    curr_clusterid_to_taxa[clusterid].append(taxon)
                except:
                    curr_clusterid_to_taxa[clusterid] = [taxon]

            print('\nCleaning up clusters...')
            curr_sensitivity_subsumed_taxa_to_clusterid = {}
            curr_nosub_taxa_to_clusterid = {}

            cleanup_object = clean_up_modules(curr_node_to_descendant_nodes, global_node_to_leaves, curr_node_to_leaves, curr_wcl, cs, global_leaf_dist_to_node, global_leaf_to_ancestors, global_node_to_parent_node, global_nodepair_to_dist, global_leaf_node_id_to_leafname)

            # leave-one-out clean up for clusters violating curr_wcl
            curr_clusterid_to_taxa, curr_taxon_to_clusterid = cleanup_object.loo_wcl_violation(curr_clusterid_to_taxa, curr_taxon_to_clusterid, global_leaf_node_id_to_leafname)

            # ensure that the most descendant-possible node-id is subtending each cluster
            curr_clusterid_to_taxa, curr_taxon_to_clusterid = cleanup_object.most_desc_nodeid_for_cluster(curr_clusterid_to_taxa, curr_taxon_to_clusterid)

            # remove cluster-size sensitivity-induced clusters
            if params.subsume_sensitivity_induced_clusters == 1:
                curr_clusterlen_distribution = get_cluster_size_distribution(curr_clusterid_to_taxa) # determine distribution of clusters
                curr_clusterid_to_taxa, curr_taxon_to_clusterid, curr_sensitivity_subsumed_taxa_to_clusterid = cleanup_object.subsume_subclusters_under_x_percentile(curr_clusterid_to_taxa, curr_taxon_to_clusterid, curr_clusterlen_distribution, params.sensitivity_percentile)

            # further subsume any sub-clusters into parent clusters which taxa set include that of sub-clusters
            if params.subsume_subclusters == 1:
                print ('Subsume statistically-homogeneous sub-clusters within parent...')
                curr_clusterlen_distribution = get_cluster_size_distribution(curr_clusterid_to_taxa) # get cluster size distribution
                curr_clusterid_to_taxa, curr_taxon_to_clusterid, curr_nosub_taxa_to_clusterid = cleanup_object.subsume_subclusters_under_x_percentile(curr_clusterid_to_taxa, curr_taxon_to_clusterid, curr_clusterlen_distribution, 100)

            # de-cluster outliers of descendent cluster that got clustered in ancestor cluster
            curr_clusterid_to_taxa, curr_taxon_to_clusterid = cleanup_object.decluster_outlying_taxa_clustered_to_anc_clusters(curr_clusterid_to_taxa, curr_taxon_to_clusterid)

            # force low-information putative clusters to cluster
            if params.force == 1:
                print ('Force low-information nodes to cluster...')
                curr_clusterid_to_taxa, curr_taxon_to_clusterid = cleanup_object.force_lowinfo_cluster(curr_clusterid_to_taxa, curr_taxon_to_clusterid, curr_node_to_mean_pwdist)

                if params.subsume_subclusters == 1:
                    curr_clusterlen_distribution = get_cluster_size_distribution(curr_clusterid_to_taxa) # get cluster size distribution
                    curr_clusterid_to_taxa, curr_taxon_to_clusterid = cleanup_object.subsume_subclusters_under_x_percentile(curr_clusterid_to_taxa, curr_taxon_to_clusterid, curr_clusterlen_distribution, 100)[:2]

            # final check to delete any clusters < cs
            curr_clusterid_to_taxa, curr_taxon_to_clusterid = cleanup_object.remove_clusters_below_cs(curr_clusterid_to_taxa, curr_taxon_to_clusterid)

            # get cluster size distribution
            curr_clusterlen_distribution = get_cluster_size_distribution(curr_clusterid_to_taxa)

            # print outputs
            #print ('Writing outputs (post-clean-up)...')
            print ('Writing outputs...')
            output_obj = output_utils(original_tree_string, global_leaf_node_id_to_leafname, curr_taxon_to_clusterid, curr_clusterid_to_taxa, curr_outfname, curr_sensitivity_subsumed_taxa_to_clusterid if params.subsume_sensitivity_induced_clusters else False, curr_nosub_taxa_to_clusterid if params.subsume_subclusters else False)
            # cluster file
            curr_modified_tree_string = output_obj.cluster_output()

            # figtree annotated tree file
            output_obj.figtree_output(curr_modified_tree_string)

            # append to summary stats output file
            params.mode = '' # single mode of analysis for now 
            curr_coverage, curr_mu_pwd, curr_mu_icd = summary_stats(curr_clusterid_to_taxa, global_nodepair_to_dist, curr_clusterlen_distribution, statsfname, treefname, len(curr_taxon_to_clusterid), len(global_leaf_node_id_to_leafname), cs, fdr, gam, params.hypo_test, params.gam_method, 'pre' if params.preQ == 1 else 'post', params.force, curr_wcl, sol_index, solver_version, params.mode, prior_input)

            # optimise
            if params.optimise:
                optimise_dict[(cs, fdr, gam, sol_index)] = (curr_coverage, curr_mu_pwd, curr_mu_icd)

            # prior output
            if params.prior:
                output_obj.prior_output(prior_input)

            # output pdf tree
            if params.pdf_tree:
                output_obj.ete3_pdf_tree_output(prior_input)

    # find optimal parameter
    if params.optimise:
        if len(parameters) > 1:
            print ('\nSorting for the optimal parameter set...')
            # reverse curr_mu_pwd
            max_curr_mu_pwd = max([v[1] for v in optimise_dict.values()])

            for k, (curr_coverage, curr_mu_pwd, curr_mu_icd) in optimise_dict.items():
                optimise_dict[k] = (curr_coverage, max_curr_mu_pwd-curr_mu_pwd, curr_mu_icd)

            if params.optimise == 'high':
                opt_cs, opt_fdr, opt_gam, opt_sol_index = sorted(optimise_dict.keys(), key=lambda _: (optimise_dict[_][1], optimise_dict[_][2], optimise_dict[_][0]))[-1]
            elif params.optimise == 'intermediate':
                opt_cs, opt_fdr, opt_gam, opt_sol_index = sorted(optimise_dict.keys(), key=lambda _: (optimise_dict[_][0], optimise_dict[_][1], optimise_dict[_][2]))[-1]

            print ('CS = {:d}, FDR = {}, GAMMA = {}, Solution index = {}.'.format(opt_cs, str(opt_fdr), str(opt_gam), opt_sol_index))

            opt_outfname = '{}_{}_cs{}_fdr{}_gam{}_sol{}{}_f{}_{}'.format(params.gam_method.lower(), params.hypo_test.lower(), str(opt_cs), str(opt_fdr), str(opt_gam), opt_sol_index, '_prior' if params.prior else '', str(params.force), treefname)
            subprocess.call('cp tree_{}.tre tree_optimal_parameter_cs{}_fdr{}_gam{}_sol{}_f{}_{}.tre'.format(opt_outfname, str(opt_cs), str(opt_fdr), str(opt_gam), opt_sol_index, str(params.force), treefname), shell=True)
            subprocess.call('cp cluster_{}.txt cluster_optimal_parameter_cs{}_fdr{}_gam{}_sol{}_f{}_{}.txt'.format(opt_outfname, str(opt_cs), str(opt_fdr), str(opt_gam), opt_sol_index, str(params.force), treefname), shell=True)

            if params.pdf_tree:
                subprocess.call('cp pdftree_{}.pdf pdftree_optimal_parameter_cs{}_fdr{}_gam{}_sol{}_f{}_{}.pdf'.format(opt_outfname, str(opt_cs), str(opt_fdr), str(opt_gam), opt_sol_index, str(params.force), treefname), shell=True)

        else:
            print ('\nWARNING: No optimising of input parameters required since only 1 parameter set was analysed.')

    print ('\n...done.\n')
