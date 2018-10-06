#!/usr/bin/env python

# Phylogenetic Clustering by Linear Integer Programming (PhyCLIP)
# Authors: Alvin X. Han and Edyth Parker

from __future__ import division
import ete3
import re
import os
import itertools
import subprocess
import argparse
import numpy as np

from phyclip_modules.stats_utils import qn, multiple_testing_correction, summary_stats, get_cluster_size_distribution
from phyclip_modules.tree_utils import parse_newick_tree
from phyclip_modules import get_global_tree_info, phyclip_output, node_leaves_reassociation, clean_up_modules

if __name__ == '__main__':
    # parse parameters
    version = 1.1
    parser = argparse.ArgumentParser(description='Phylogenetic Clustering by Linear Integer Programming (PhyCLIP) v{}'.format(version))
    
    required_args = parser.add_argument_group('Required')
    required_args.add_argument('-i', '--input_file', type=str, help='Input file. See manual for format details.')
    
    analyses_aid = parser.add_argument_group('Analysis aids')
    analyses_aid.add_argument('--treeinfo', type=str, help='*_treeinfo.txt file generated from PREVIOUS PhyCLIP run of the SAME phylogenetic tree.')
    analyses_aid.add_argument('--no_treeinfo', action='store_true', help='Stop PhyCLIP from generating *_treeinfo.txt file.')
    analyses_aid.add_argument('--pdf_tree', action='store_true', help='PDF tree output annotated with cluster results.')
    analyses_aid.add_argument('--optimise', choices=['intermediate', 'high'], help='PhyCLIP automatically searches for the clustering result from the OPTIMAL input parameter set. See documentation for details. Must have >1 input parameter set in input file.')

    analyses_options = parser.add_argument_group('Analysis options')
    analyses_options.add_argument('--tree_outgroup', type=str, help='Taxon (name as appeared in tree) to be set as outgroup for rooting.')
    analyses_options.add_argument('--midpoint', action='store_true', help='Root tree by mid-point node.')
    analyses_options.add_argument('--collapse_zero_branch_length', default=0, choices=[0, 1], type=int, help='Collapse internal nodes with zero branch length of tree before running PhyCLIP (default = %(default)s).')
    analyses_options.add_argument('--equivalent_zero_length', default=0.000001, type=float, help='Maximum branch length to be rounded to zero if the --collapse_zero_branch_length flag is passed (advanced option, default = %(default)s).')

    analyses_options.add_argument('--subsume_sensitivity_induced_clusters', default=1, choices=[0, 1], type=int, help='Subsume cluster-size sensitivity-induced clusters into parent cluster (default = %(default)s).')
    analyses_options.add_argument('--sensitivity_percentile', default=25, type=int, help='Percentile of cluster size distribution under which a cluster is considered to be sensitivity-induced (advanced option, default = %(default)s%%).')
    analyses_options.add_argument('--subsume_subclusters', default=0, choices=[0, 1], type=int, help='Subsume sub-clusters into their respective parent clusters (default = %(default)s).')

    analyses_options.add_argument('--gam_method', choices=['mad', 'qn'], default='qn',  help='Method to estimate robust dispersion measure (default = %(default)s).')
    analyses_options.add_argument('--hypo_test', choices=['kuiper', 'kolsmi'], default='kuiper', help='Hypothesis test to use for statistical differentiation of distance distributions (default = %(default)s).')
    analyses_options.add_argument('--preQ', default=0, choices=[0, 1], type=int, help='Perform Benjamini-Hochberg corrections of p-values BEFORE filtering nodes that are < minimum cluster size (advanced option, default = %(default)s).')

    prior_options = parser.add_argument_group('Prior options')
    prior_options.add_argument('--prior', type=str, help='Prior information on clustered taxa. File format same as PhyCLIP cluster text file output (i.e. CLUSTER-ID{tab}SEQUENCE-NAME). Only works with gurobi solver, no support for glpk.')
    prior_options.add_argument('--prior_weights', type=str, help='Weights on importance/confidence between prior clusters to remain clustered together. File format: CLUSTER-ID{tab}INT/FLOAT_WEIGHT{\\newline}. Equal weights will be assumed if none is given.')

    solver_options = parser.add_argument_group('Solver options')
    solver_options.add_argument('--solver', default='gurobi', choices=['glpk', 'gurobi'], type=str, help='Preferred ILP solver IF more than one solvers are available (default: %(default)s).')
    solver_options.add_argument('--solver_verbose', default=0, choices=[0, 1], type=int, help='ILP solver verbose (default: %(default)s)')
    solver_options.add_argument('--solver_check', action='store_true', help='Check available ILP solver(s) installed.')

    solver_options.add_argument('--threads', type=int, help='Number of threads to use (default = all).')

    params = parser.parse_args()

    print ('{}\n\n{:^72}\n{:^72}\n\n{}'.format(''.join(['-']*72), 'Phylogenetic Clustering by Linear Integer Programming (PhyCLIP)', 'v{}'.format(version), ''.join(['-']*72)))

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

    # parse newick tree file
    try:
        newick_tree_string = parse_newick_tree(treepath)
    except:
        raise Exception('\nInvalid tree file. Check that the correct path to the NEWICK tree file is given in the first line of the input file.\n')

    try:
        tree = ete3.Tree(newick_tree_string)
    except:
        tree = ete3.Tree(newick_tree_string, format=1)

    print('\nTree file...OK')

    # changing tree root
    if params.tree_outgroup or params.midpoint:

        root_node = ''

        if params.tree_outgroup:
            try:
                root_node = tree&"'{}'".format(params.tree_outgroup)
            except:
                raise Exception('\nGiven outgroup {} did not match any taxon node in tree.\n'.format(params.tree_outgroup))

        if params.midpoint:
            if root_node != '':
                raise Exception('\nOnly ONE rooting method is accepted.\n')
            else:
                root_node = tree.get_midpoint_outgroup()

        tree.set_outgroup(root_node)
        print ('\nWARNING: Tree outgroup changed to {}.'.format(root_node.name if params.tree_outgroup else 'mid-point'))

    tree.ladderize() # ladderize tree
    taxon_list = tree.get_leaf_names() # get list of all taxa

    # collapse zero branch length
    if params.collapse_zero_branch_length == 1:
        from phyclip_modules.tree_utils import collapse_zero_branch_length
        tree = collapse_zero_branch_length(tree, params.equivalent_zero_length)

        if tree == False:
            raise Exception('No branches were collapsed. Check the upper limit of zero-length branches on your tree and adjust accordingly using --equivalent_zero_length')

        # write zero branch length collapsed tree as a newick file
        treefname = 'zero-branch-length-collapsed_{}'.format(treefname)
        tree.write(format=5, outfile='{}.nwk'.format(treefname))
        print ('Collapsed newick tree file...{}.nwk'.format(treefname))

    original_tree_string = tree.write(format=5)

    # subsequent lines are parameters (cs, fdr, gam)
    parameters = []
    for _, line in enumerate(infhandle):
        try:
            cs, fdr, gam = line.strip().split(',')
            try:
                parameters.append((int(cs), float(fdr), float(gam)))
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

        for x in cs_range:
            for y in fdr_range:
                for z in gam_range:
                    parameters.append((x, y, z))

    print('Parameter sets...OK')

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

        # if cluster prior is given, check that all stated taxa in prior file are present in tree
        taxa_in_prior_not_present_in_tree = list(set([i for j in prior_input.values() for i in j]) - set(taxon_list))
        if len(taxa_in_prior_not_present_in_tree) > 0:
            raise Exception('\nSome taxa in prior file are not present in tree - {}.\n'.format(', '.join(taxa_in_prior_not_present_in_tree)))

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

    # preferred solver
    if params.solver not in available_solvers:
        print ('\nWARNING: {} is not installed.'.format(params.solver))
        params.solver = available_solvers.keys()[0]

    print('\nILP solver...{}'.format(available_solvers[params.solver]))

    # only gurobi has prior support
    if params.prior and params.solver != 'gurobi':
        if 'gurobi' in available_solvers:
            params.solver = 'gurobi'
            print ('WARNING: Prior analyses can only be performed using gurobi. Switching to {}...'.format(available_solvers[params.solver]))
        else:
            raise Exception('\nPrior analyses can only be performed using gurobi solver.\n')

    # limit number of threads for parallelization
    try:
        ncpu = int(params.threads)
    except:
        from pathos.helpers import mp
        ncpu = mp.cpu_count()
    print ('Threads...{}'.format(ncpu))

    # check if summary stats file already exist in current working directory
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
            """output.write('Treefile\tCS\tFDR\tMAD\tKuiper/KS\tQn/MAD\tq-values\tWithin_cluster_limit\tSolution_Index\tClean_up\tPrior_taxa_clustered(%)\t'
             '#_of_clustered_sequences\tTotal_no_of_sequences\t%_clustered\t#_of_clusters\t'
             'Mean_cluster_size\tS.D.\tMedian_cluster_size\tMAD\tMin_cluster_size\tMax_cluster_size\t'
             'Grand_mean_of_mean_pwd\tS.D.\tGrand_mean_of_median_pwd\tS.D.\tMin_mean_pwd\tMax_mean_pwd\t'
             'Mean_of_inter-cluster_dist\tS.D.\tMedian_of_inter-cluster_dist\tMAD\tMin_inter-cluster_dist\tMax_inter-cluster_dist\tsolver(version)\n')"""

            output.write('Treefile\tCS\tFDR\tMAD\tKuiper/KS\tQn/MAD\tq-values\tWithin_cluster_limit\tSolution_Index\tPrior_taxa_clustered(%)\t'
                         '#_of_clustered_sequences\tTotal_no_of_sequences\t%_clustered\t#_of_clusters\t'
                        'Mean_cluster_size\tS.D.\tMedian_cluster_size\tMAD\tMin_cluster_size\tMax_cluster_size\t'
                        'Grand_mean_of_mean_pwd\tS.D.\tGrand_mean_of_median_pwd\tS.D.\tMin_mean_pwd\tMax_mean_pwd\t'
                        'Mean_of_inter-cluster_dist\tS.D.\tMedian_of_inter-cluster_dist\tMAD\tMin_inter-cluster_dist\tMax_inter-cluster_dist\tsolver(version)\n')


    # --- GLOBAL TREE INFORMATION --- #
    print ('\nGetting tree information...')

    # parse treeinfo file if given
    global_leaf_dist_to_node = {}
    global_leafpair_to_distance = {}
    global_nodepair_to_pval = {}

    # read treeinfo file if available
    if params.treeinfo:
        from phyclip_modules.tree_utils import parse_treeinfo_file
        try:
            global_leaf_dist_to_node, global_leafpair_to_distance, global_nodepair_to_pval = parse_treeinfo_file(params.treeinfo)
            print('Treeinfo file...OK.')
        except:
            raise Exception('\nInvalid treeinfo file. Try re-running without treeinfo file.\n')
    else:
        # treeinfo file = TRUE
        if params.no_treeinfo == False:
            params.treeinfo = '{}_treeinfo.txt'.format(treefname)

            # check if a copy of treeinfo file is already in current working directiory
            if params.treeinfo in os.listdir(os.getcwd()):
                try:
                    overwrite_treeinfo = re.search('(y|n)', raw_input('\nWARNING: TREEINFO FILE {} exists in current working directory. Overwrite? (y/n): '.format(params.treeinfo))).group()
                except:
                    raise Exception('\nInvalid input.\n')

                # no overwrite over previous treeinfo file.
                if overwrite_treeinfo == 'n':
                    from phyclip_modules.tree_utils import parse_treeinfo_file
                    try:
                        global_leaf_dist_to_node, global_leafpair_to_distance, global_nodepair_to_pval = parse_treeinfo_file(params.treeinfo)
                        print('Treeinfo file...OK.')
                    except:
                        print ('\nWARNING: Invalid treeinfo file in current working directory. Overwriting...')
                        pass

    global_tree_info_obj = get_global_tree_info(tree, global_leaf_dist_to_node, global_leafpair_to_distance, global_nodepair_to_pval, params.treeinfo, params.hypo_test, params.no_treeinfo, ncpu)

    global_tree_string, global_node_to_leaves, global_nindex_to_node, global_node_to_nindex, global_leaf_dist_to_node, global_nodepair_to_dist, global_node_to_parent_node, global_node_to_mean_child_dist2root = global_tree_info_obj.node_indexing()

    global_leafpair_to_distance, global_node_to_pwdist, global_node_to_mean_pwdist, global_node_to_ancestral_nodes, global_node_to_descendant_nodes, global_leaf_to_ancestors, global_node_to_mean_child_dist2anc = global_tree_info_obj.pwdist_dist_and_ancestral_trace(global_node_to_leaves, global_nindex_to_node, global_node_to_nindex, global_node_to_mean_child_dist2root, global_nodepair_to_dist)

    global_nodepair_to_pval = global_tree_info_obj.get_global_pval(params.hypo_test, global_node_to_leaves, global_node_to_ancestral_nodes, global_node_to_pwdist)

    # pre-calculation for within-cluster limit
    med_x = np.median(global_node_to_mean_pwdist.values())
    if params.gam_method == 'Qn':
        mad_x = qn(global_node_to_mean_pwdist.values())
    else:
        mad_x = np.median([x-med_x for x in global_node_to_mean_pwdist.values() if x >= med_x])

    ### --- GLOBAL TREE INFORMATION --- ###

    csgam_to_reassociation_memory = {} # memory of reassociated nodes/leaves for repeated (cs, gam) set

    # Start runs
    print ('\nStarting runs...')
    optimise_dict = {}
    for p_index, (cs, fdr, gam) in enumerate(parameters):
        ### --- PARAMETER-SPECIFIC TREE INFORMATION --- ###
        print ('\n{:<26}{:^20}{:>26}'.format(''.join(['-']*25), ", ".join(['cs{}'.format(cs), 'fdr{}'.format(fdr), 'gam{}'.format(gam)]), ''.join(['-']*25)))

        # current within-cluster limit
        curr_wcl = med_x + (gam * mad_x)

        # level-order sorted list of nodes with leaves >= cs (only reassociate such nodes)
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

            for node, leaves in curr_node_to_leaves.items():
                # filter by cs
                if len(leaves) < cs:
                    del curr_node_to_leaves[node]
                    continue

            # save to memory to speed up subsequent runs with the same cs/gam parameter set
            csgam_to_reassociation_memory[(cs, gam)] = [curr_node_to_leaves, curr_node_to_descendant_nodes, curr_node_to_mean_pwdist]

        # update pairwise distance dictionaries/nodes' ancestry relations
        print ('Updating tree info...')
        curr_leaves = [x for y in curr_node_to_leaves.values() for x in y]
        curr_list_of_ancestral_node = curr_node_to_leaves.keys()[:]
        curr_node_to_descendant_nodes = {k:list(set(v)&set(curr_list_of_ancestral_node)) for k,v in curr_node_to_descendant_nodes.items() if k in curr_list_of_ancestral_node and len(v) > 0}
        #curr_node_to_mean_pwdist = {k:v for k,v in curr_node_to_mean_pwdist.items() if k in curr_list_of_ancestral_node}

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
            #try:
            all_solutions = gurobi_solver(curr_node_to_leaves, curr_leaves, curr_list_of_ancestral_node, curr_nodepair_to_qval, curr_node_to_mean_pwdist, curr_wcl, cs, fdr, prior_input, prior_weights, params.solver_verbose, p_index, ncpu)
            #except:
            #raise Exception('\nUnable to solve ILP model using gurobi. You may try to use other available solvers by the --solver flag.\n')
        else:
            from phyclip_modules.glpk_solver import glpk_solver
            #try:
            all_solutions = glpk_solver(curr_node_to_leaves, curr_leaves, curr_list_of_ancestral_node, curr_nodepair_to_qval, curr_node_to_mean_pwdist, curr_wcl, cs, fdr, params.solver_verbose)
            #except:
            #raise Exception('\nUnable to solve ILP model using glpk. You may try to use other available solvers by the --solver flag.\n')

        if all_solutions == 'na':
            # continue to next parameter set if no solution
            with open(statsfname, 'a') as output:
                output.write('{}\t{}\t{}\t{}\tNO OPTIMAL SOLUTION FOUND.\n'.format(treefname, cs, fdr, gam))
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

            curr_outfname = '{}_{}_cs{}_fdr{}_gam{}_sol{}{}_{}'.format(params.gam_method.lower(), params.hypo_test.lower(), str(cs), str(fdr), str(gam), sol_index, '_prior' if params.prior else '', treefname)
            curr_clusterid_to_taxa = {}
            for taxon, clusterid in curr_taxon_to_clusterid.items():
                try:
                    curr_clusterid_to_taxa[clusterid].append(taxon)
                except:
                    curr_clusterid_to_taxa[clusterid] = [taxon]

            """
            #! --- print pre-clean up results start --- !#

            print ('Writing outputs (pre-clean-up)...')
            pre_clean_up_output_obj = phyclip_output(original_tree_string, global_tree_string, curr_taxon_to_clusterid, curr_clusterid_to_taxa, taxon_list, curr_outfname, 'pre-clean')

            # cluster file
            pre_clean_up_modified_tree_string = pre_clean_up_output_obj.cluster_output()
            # figtree annotated tree file
            pre_clean_up_output_obj.figtree_output(pre_clean_up_modified_tree_string)

            # append to summary stats output file
            # get clusterlen distribution
            pre_clean_up_clusterlen_distribution = get_cluster_size_distribution(curr_clusterid_to_taxa)

            summary_stats(curr_clusterid_to_taxa, global_leafpair_to_distance, global_nodepair_to_dist, pre_clean_up_clusterlen_distribution, statsfname, treefname, len(curr_taxon_to_clusterid), len(taxon_list), cs, fdr, gam, params.hypo_test, params.gam_method, 'pre' if params.preQ == 1 else 'post', curr_wcl, sol_index, 'pre-clean', prior_input, available_solvers[params.solver])

            # ! --- print pre-clean up results end --- !#
            """

            print('\nCleaning up clusters...')
            cleanup_object = clean_up_modules(curr_node_to_descendant_nodes, global_node_to_leaves, global_leafpair_to_distance, curr_node_to_leaves, curr_wcl, cs, global_leaf_dist_to_node, global_leaf_to_ancestors, global_node_to_parent_node, global_nodepair_to_dist)

            # ensure that the most descendant-possible node-id is subtending each cluster
            curr_clusterid_to_taxa, curr_taxon_to_clusterid = cleanup_object.most_desc_nodeid_for_cluster(curr_clusterid_to_taxa, curr_taxon_to_clusterid)

            # leave-one-out clean up for clusters violating wcl
            curr_clusterid_to_taxa, curr_taxon_to_clusterid = cleanup_object.loo_wcl_violation(curr_clusterid_to_taxa, curr_taxon_to_clusterid)

            # remove cluster-size sensitivity-induced clusters
            if params.subsume_sensitivity_induced_clusters:
                curr_clusterlen_distribution = get_cluster_size_distribution(curr_clusterid_to_taxa) # determine distribution of clusters
                curr_clusterid_to_taxa, curr_taxon_to_clusterid, curr_sensitivity_subsumed_taxa_to_clusterid = cleanup_object.subsume_subclusters_under_x_percentile(curr_clusterid_to_taxa, curr_taxon_to_clusterid, curr_clusterlen_distribution, params.sensitivity_percentile)

            # de-cluster outliers of descendent cluster that got clustered in ancestor cluster
            curr_clusterid_to_taxa, curr_taxon_to_clusterid = cleanup_object.decluster_outlying_taxa_clustered_to_anc_clusters(curr_clusterid_to_taxa, curr_taxon_to_clusterid)

            # delete any clusters < cs
            curr_clusterid_to_taxa, curr_taxon_to_clusterid = cleanup_object.remove_clusters_below_cs(curr_clusterid_to_taxa, curr_taxon_to_clusterid)

            # further subsume any sub-clusters into parent clusters which taxa set include that of sub-clusters
            if params.subsume_subclusters:
                print ('Subsume sub-clusters within parent...')
                curr_clusterlen_distribution = get_cluster_size_distribution(curr_clusterid_to_taxa) # get cluster size distribution
                curr_clusterid_to_taxa, curr_taxon_to_clusterid, curr_nosub_taxa_to_clusterid = cleanup_object.subsume_subclusters_under_x_percentile(curr_clusterid_to_taxa, curr_taxon_to_clusterid, curr_clusterlen_distribution, 100)

            # get cluster size distribution
            curr_clusterlen_distribution = get_cluster_size_distribution(curr_clusterid_to_taxa)

            # print outputs
            #print ('Writing outputs (post-clean-up)...')
            print ('Writing outputs...')
            #output_obj = phyclip_output(original_tree_string, global_tree_string, curr_taxon_to_clusterid, curr_clusterid_to_taxa, taxon_list, curr_outfname, 'post-clean', curr_sensitivity_subsumed_taxa_to_clusterid if params.subsume_sensitivity_induced_clusters else False, curr_nosub_taxa_to_clusterid if params.subsume_subclusters else False)
            output_obj = phyclip_output(original_tree_string, global_tree_string, curr_taxon_to_clusterid, curr_clusterid_to_taxa, taxon_list, curr_outfname, curr_sensitivity_subsumed_taxa_to_clusterid if params.subsume_sensitivity_induced_clusters else False, curr_nosub_taxa_to_clusterid if params.subsume_subclusters else False)
            # cluster file
            curr_modified_tree_string = output_obj.cluster_output()
            # figtree annotated tree file
            output_obj.figtree_output(curr_modified_tree_string)

            # append to summary stats output file
            #curr_coverage, curr_mu_pwd, curr_mu_icd = summary_stats(curr_clusterid_to_taxa, global_leafpair_to_distance, global_nodepair_to_dist, curr_clusterlen_distribution, statsfname, treefname, len(curr_taxon_to_clusterid), len(taxon_list), cs, fdr, gam, params.hypo_test, params.gam_method, 'pre' if params.preQ == 1 else 'post', curr_wcl, sol_index, 'post-clean', prior_input, solver_version)
            curr_coverage, curr_mu_pwd, curr_mu_icd = summary_stats(curr_clusterid_to_taxa, global_leafpair_to_distance, global_nodepair_to_dist, curr_clusterlen_distribution, statsfname, treefname, len(curr_taxon_to_clusterid), len(taxon_list), cs, fdr, gam, params.hypo_test, params.gam_method, 'pre' if params.preQ == 1 else 'post', curr_wcl, sol_index, prior_input, solver_version)

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
            print ('\nDetermining the optimal parameter...')
            # reverse curr_mu_pwd
            max_curr_mu_pwd = max([v[1] for v in optimise_dict.values()])

            for k, (curr_coverage, curr_mu_pwd, curr_mu_icd) in optimise_dict.items():
                optimise_dict[k] = (curr_coverage, max_curr_mu_pwd-curr_mu_pwd, curr_mu_icd)

            if params.optimise == 'high':
                opt_cs, opt_fdr, opt_gam, opt_sol_index = sorted(optimise_dict.keys(), key=lambda _: (optimise_dict[_][1], optimise_dict[_][2], optimise_dict[_][0]))[-1]
            elif params.optimise == 'intermediate':
                opt_cs, opt_fdr, opt_gam, opt_sol_index = sorted(optimise_dict.keys(), key=lambda _: (optimise_dict[_][0], optimise_dict[_][1], optimise_dict[_][2]))[-1]

            print ('CS = {}, FDR = {}, GAMMA = {}, Solution index = {}.'.format(opt_cs, opt_fdr, opt_gam, opt_sol_index))

            opt_outfname = '{}_{}_cs{}_fdr{}_gam{}_sol{}{}_{}'.format(params.gam_method.lower(), params.hypo_test.lower(), str(opt_cs), str(opt_fdr), str(opt_gam), opt_sol_index, '_prior' if params.prior else '', treefname)
            subprocess.call('cp tree_{}.tre tree_optimal_parameter_cs{}_fdr{}_gam{}_sol{}_{}.tre'.format(opt_outfname, str(opt_cs), str(opt_fdr), str(opt_gam), opt_sol_index, treefname), shell=True)
            subprocess.call('cp cluster_{}.txt cluster_optimal_parameter_cs{}_fdr{}_gam{}_sol{}_{}.txt'.format(opt_outfname, str(opt_cs), str(opt_fdr), str(opt_gam), opt_sol_index, treefname), shell=True)

            if params.pdf_tree:
                subprocess.call('cp pdftree_{}.pdf pdftree_optimal_parameter_cs{}_fdr{}_gam{}_sol{}_{}.pdf'.format(opt_outfname, str(opt_cs), str(opt_fdr), str(opt_gam), opt_sol_index, treefname), shell=True)

        else:
            print ('\nWARNING: No optimising of input parameters required since only 1 parameter set was analysed.')

    print ('\n...All parameter sets analysed.\n')

