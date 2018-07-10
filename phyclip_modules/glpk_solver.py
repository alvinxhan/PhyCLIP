import itertools
import subprocess

def glpk_solver(node_to_leaves, all_leaves, list_of_ancestral_node, nodepair_to_qval, node_to_mean_pwdist, within_cluster_limit, min_cluster_size, fdr_cutoff, verbose):
    print ('Solving with glpk...')

    # generate model/data file
    with open('phyclip_temp_glpk.mod', 'w') as output:
        # define sets
        output.write('set LEAF_BINARY_INDICES dimen 2;\r\n'
                     'set CURR_TAXON_LIST;\r\n'
                     'set NODE_INDICES;\r\n\r\n')

        # define parameters
        output.write('param QVAL{node_i in NODE_INDICES, node_j in NODE_INDICES: node_i < node_j};\r\n'
                     'param NODE_MEAN_DIST{node in NODE_INDICES};\r\n\r\n')

        # define variables
        output.write('var leaf_decision{(leaf, node) in LEAF_BINARY_INDICES}, binary;\r\n'
                     'var node_decision{node in NODE_INDICES}, binary;\r\n\r\n')

        # optimization function
        output.write('maximize obj: sum{(leaf, node) in LEAF_BINARY_INDICES} leaf_decision[leaf, node];\r\n\r\n')

        # constraints
        # leaf binary should be related to node binary
        output.write('s.t. leaf_node_rs1{(leaf, node) in LEAF_BINARY_INDICES}: leaf_decision[leaf, node] <= node_decision[node];\r\n\r\n')

        # each leaf can only select one node to cluster to
        output.write('s.t. leaf_node_rs2{leaf in CURR_TAXON_LIST}: sum{node in NODE_INDICES: (leaf, node) in LEAF_BINARY_INDICES} leaf_decision[leaf, node] <= 1;\r\n\r\n')

        # leaf always chooses the youngest cluster-able node
        output.write('s.t. leaf_node_rs3{leaf in CURR_TAXON_LIST, node_i in NODE_INDICES, node_j in NODE_INDICES: node_i < node_j and (leaf, node_i) in LEAF_BINARY_INDICES and (leaf, node_j) in LEAF_BINARY_INDICES}: leaf_decision[leaf, node_i] <= 2 - node_decision[node_i] - node_decision[node_j];\r\n\r\n')

        # cluster size constraint
        output.write('s.t. clust_size_cnst{{node in NODE_INDICES}}: ({}+9999)*(node_decision[node] - 1) + {} <= sum{{leaf in CURR_TAXON_LIST: (leaf, node) in LEAF_BINARY_INDICES}} leaf_decision[leaf, node];\r\n\r\n'.format(len(all_leaves), min_cluster_size))

        # inter-cluster constraint
        output.write('s.t. inter_clust_cnst{{node_i in NODE_INDICES, node_j in NODE_INDICES: node_i < node_j}}: 9999*(2 - node_decision[node_i] - node_decision[node_j]) >= QVAL[node_i, node_j] - {};\r\n\r\n'.format(fdr_cutoff))


        # within-cluster constraint
        output.write('s.t. within_clust_cnst{{node in NODE_INDICES}}: 9999*(node_decision[node] - 1) <= {} - NODE_MEAN_DIST[node];\r\n\r\n'.format(within_cluster_limit))

        output.write('solve;\r\n\r\n')

        # write output
        output.write('printf{(leaf, node) in LEAF_BINARY_INDICES}: "%s@@%s@@%s,", leaf, node, leaf_decision[leaf, node] > "phyclip_temp_glpk.out";\r\n\r\n')

        # write data
        output.write('data;\r\n\r\n')
        # write sets
        # leaf_binary_indices
        output.write('set LEAF_BINARY_INDICES := {};\r\n\r\n'.format(' '.join(['("{}", {})'.format(leaf, n) for n, leaves in node_to_leaves.items() for leaf in leaves])))

        # curr_taxon_list
        output.write('set CURR_TAXON_LIST := {};\r\n\r\n'.format(' '.join(['"{}"'.format(leaf) for leaf in all_leaves])))

        # node_indices
        output.write('set NODE_INDICES := {};\r\n\r\n'.format(' '.join(['{}'.format(n) for n in list_of_ancestral_node])))

        # node pairwise mean dist
        output.write('param NODE_MEAN_DIST := {};\r\n\r\n'.format(', '.join(['[{}] {}'.format(n, node_to_mean_pwdist[n]) for n in list_of_ancestral_node])))

        # node pairs to qval
        output.write('param QVAL := {};\r\n\r\n'.format(', '.join(['[{}, {}] {}'.format(n, m, qval) for (n, m), qval in nodepair_to_qval.items() if n < m])))

        output.write('end;\r\n')

    # solve using glpsol
    cmd = ['glpsol', '--math', 'phyclip_temp_glpk.mod']
    if verbose == 0:
        subprocess.call(cmd, stdout=subprocess.PIPE)
    else:
        subprocess.call(cmd)

    try:
        # read output file and return
        fhandle = open('phyclip_temp_glpk.out', 'rU').readlines()
        glpk_output = fhandle.pop(0).strip()
        glpk_output = filter(None, glpk_output.split(','))

        # remove temporary files
        subprocess.call('rm phyclip_temp*', shell=True)

        taxon_to_clusterid = {}
        for _ in glpk_output:
            leaf, node, output_bin = _.split('@@')
            if int(output_bin) == 1:
                taxon_to_clusterid[leaf] = int(node)

        return [taxon_to_clusterid]

    except:
        # remove temporary files
        subprocess.call('rm phyclip_temp*', shell=True)
        return 'na'
