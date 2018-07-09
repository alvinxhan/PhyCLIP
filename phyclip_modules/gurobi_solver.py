from gurobipy import *
import itertools

def gurobi_solver(node_to_leaves, all_leaves, list_of_ancestral_node, nodepair_to_qval, node_to_mean_pwdist, within_cluster_limit, min_cluster_size, fdr_cutoff, prior, pc_weights, verbose):
    print ('Solving with gurobi...')

    # set up indices
    leaf_binary_indices = []
    leaf_to_ancestral_nodes = {}
    if prior:
        p_to_nodes_to_prior_leaves = {}

    for n, leaves in node_to_leaves.items():
        for leaf in leaves:
            leaf_binary_indices.append((leaf, n))

            try:
                leaf_to_ancestral_nodes[leaf].append(n)
            except:
                leaf_to_ancestral_nodes[leaf] = [n]

        if prior:
            for p, leaves_in_prior_cluster in prior.items():
                prior_leaves_subtended = list(set(leaves_in_prior_cluster)&set(leaves))
                if len(prior_leaves_subtended) > 0:
                    try:
                        p_to_nodes_to_prior_leaves[p][n] = prior_leaves_subtended
                    except:
                        p_to_nodes_to_prior_leaves[p] = {n:prior_leaves_subtended}

    if prior:
        # remove any prior clusters where only a single node subtends them
        p_to_nodes_to_prior_leaves = {p:node_to_prior_leaves for p, node_to_prior_leaves in p_to_nodes_to_prior_leaves.items() if len(node_to_prior_leaves) > 1}

        if len(p_to_nodes_to_prior_leaves) == 0:
            prior = False
        else:
            p_to_nodeij_permutations = {p:[(yi, yj) for (yi, yj) in itertools.permutations(node_to_prior_leaves.keys(), 2)] for p, node_to_prior_leaves in p_to_nodes_to_prior_leaves.items()}

    # model
    model = Model('NewModel')

    # verbose
    model.Params.LogToConsole = verbose

    # variables
    node_decision = model.addVars(list_of_ancestral_node, vtype=GRB.BINARY)
    leaf_decision = model.addVars(leaf_binary_indices, vtype=GRB.BINARY)
    if prior:
        prior_taxa_clustered = model.addVars(p_to_nodes_to_prior_leaves.keys(), lb = 0, vtype=GRB.INTEGER)
        # binary denoting 1 if number of prior leaves clustered in yi < yj; 0 if yi >= yj
        y = model.addVars([(p, yi, yj) for p, nodeij_pairs in p_to_nodeij_permutations.items() for (yi, yj) in nodeij_pairs], vtype=GRB.BINARY)
        # binary denoting 1 if sum_j(y[i,j]) == 0 (i.e. nodes with largest counts of prior taxa clustered)
        z = model.addVars([(p, n) for p, node_to_prior_leaves in p_to_nodes_to_prior_leaves.items() for n in node_to_prior_leaves.keys()], vtype=GRB.BINARY)

    # update model
    model.update()

    # constraints
    # leaf binary should be related to node binary
    model.addConstrs(leaf_decision[(leaf,n)] <= node_decision[n] for (leaf, n) in leaf_binary_indices)

    # each leaf can only select one node to cluster to
    model.addConstrs(quicksum(leaf_decision[(leaf, n)] for n in list_of_ancestral_node if (leaf, n) in leaf_binary_indices) <= 1 for leaf in all_leaves)

    # leaf always chooses the youngest cluster-able node
    for leaf, ancestral_nodes in leaf_to_ancestral_nodes.items():
        for n, m in itertools.combinations(ancestral_nodes, 2):
            if n < m:
                model.addConstr(leaf_decision[(leaf, n)] <= 2 - node_decision[n] - node_decision[m])
            else:
                model.addConstr(leaf_decision[(leaf, m)] <= 2 - node_decision[m] - node_decision[n])

    # inter-cluster constraint
    for (n, m) in itertools.combinations(list_of_ancestral_node, 2):
        model.addConstr((2 - node_decision[n] - node_decision[m])*999 >= nodepair_to_qval[(n,m)]-fdr_cutoff)

    for n in list_of_ancestral_node:
        # cluster size constraint
        model.addConstr(len(all_leaves)*(node_decision[n]-1) + min_cluster_size <= quicksum(leaf_decision[(leaf, n)] for leaf in node_to_leaves[n]))
        # within-cluster constraint
        model.addConstr(999*(node_decision[n]-1) <= within_cluster_limit - node_to_mean_pwdist[n])

    # prior
    if prior:
        for p, nodeij_permutations in p_to_nodeij_permutations.items():

            node_to_prior_leaves = p_to_nodes_to_prior_leaves[p]

            for (yi, yj) in nodeij_permutations:
                model.addConstr(quicksum(leaf_decision[(leaf, yi)] for leaf in node_to_prior_leaves[yi]) >= quicksum(leaf_decision[(leaf, yj)] for leaf in node_to_prior_leaves[yj]) - 9999*y[(p, yi, yj)])

                model.addConstr(quicksum(leaf_decision[(leaf, yj)] for leaf in node_to_prior_leaves[yj]) >= quicksum(leaf_decision[(leaf, yi)] for leaf in node_to_prior_leaves[yi]) + 1 - 9999*(1 - y[(p, yi, yj)]))

            for n, prior_leaves_subtended in node_to_prior_leaves.items():
                model.addConstr(
                    quicksum(y[(p, n, m)] for m in node_to_prior_leaves.keys() if n != m) <= 9999*(1-z[(p, n)]))

                model.addConstr(
                    quicksum(y[(p, n, m)] for m in node_to_prior_leaves.keys() if n != m) >= 1 - 9999*z[(p, n)])

            model.addConstr(prior_taxa_clustered[p] <= 9999*quicksum(z[(p, n)] for n in node_to_prior_leaves.keys()))
            for n, prior_leaves_subtended in node_to_prior_leaves.items():
                model.addConstr(prior_taxa_clustered[p] - quicksum(leaf_decision[(leaf, n)] for leaf in prior_leaves_subtended) >= 9999*(z[(p, n)]-1))
                model.addConstr(prior_taxa_clustered[p] - quicksum(leaf_decision[(leaf, n)] for leaf in prior_leaves_subtended) <= 9999*(1-z[(p, n)]))

    # update model
    model.update()

    # objective function - maximize number of strains clustered
    model.ModelSense = GRB.MAXIMIZE
    if prior:
        model.NumObj = 2
        model.setObjectiveN(quicksum(leaf_decision[(leaf, n)] for (leaf, n) in leaf_binary_indices), index=0, priority=1, name='primary')
        model.setObjectiveN(quicksum(pc_weights[p]*prior_taxa_clustered[p] for p in p_to_nodes_to_prior_leaves.keys()), index=1, priority=0, name='prior')
    else:
        model.setObjective(quicksum(leaf_decision[(leaf, n)] for (leaf, n) in leaf_binary_indices))

    # update model
    model.update()

    # optimize
    model.optimize()

    if model.status == GRB.Status.OPTIMAL:
        optimal_solutions = []

        # query optimal solution objective
        optsol_obj = model.getAttr('PoolObjVal')

        # get solution pool size
        solution_pool_size = model.getAttr('SolCount')

        for sol_index in xrange(solution_pool_size):

            # get solution
            model.params.solutionNumber = sol_index
            curr_sol_obj = model.getAttr('PoolObjVal')

            if curr_sol_obj == optsol_obj:
                primaryobj_solution = model.getAttr('Xn', leaf_decision)

                # prior clusters
                if prior:
                    p_to_leaves_clustered = model.getAttr('Xn', prior_taxa_clustered)
                    curr_prior_obj = sum([pc_weights[p]*p_to_leaves_clustered[p] for p in p_to_leaves_clustered.keys()])
                    if sol_index == 0:
                        prior_opt_sol_obj = curr_prior_obj
                    elif curr_prior_obj < prior_opt_sol_obj:
                        break

                taxon_to_clusterid = {leaf:n for (leaf, n) in leaf_binary_indices if primaryobj_solution[(leaf, n)] > 0}

                optimal_solutions.append(taxon_to_clusterid)

        return optimal_solutions
    else:
        return 'na'