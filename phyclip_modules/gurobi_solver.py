from gurobipy import *
import itertools

def gurobi_solver(node_to_leaves, all_leaves, list_of_ancestral_node, nodepair_to_qval, node_to_mean_pwdist, within_cluster_limit, min_cluster_size, fdr_cutoff, verbose):
    print ('Solving with gurobi...')

    # set up indices
    leaf_binary_indices = []
    leaf_to_ancestral_nodes = {}
    for n, leaves in node_to_leaves.items():
        for leaf in leaves:
            leaf_binary_indices.append((leaf, n))

            try:
                leaf_to_ancestral_nodes[leaf].append(n)
            except:
                leaf_to_ancestral_nodes[leaf] = [n]

    # model
    model = Model('NewModel')

    # verbose
    model.Params.LogToConsole = verbose

    # variables
    node_decision = model.addVars(list_of_ancestral_node, vtype=GRB.BINARY)
    leaf_decision = model.addVars(leaf_binary_indices, vtype=GRB.BINARY)

    # objective function - maximize number of strains clustered
    model.ModelSense = GRB.MAXIMIZE
    model.setObjective(quicksum(leaf_decision[(leaf, n)] for (leaf, n) in leaf_binary_indices))
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

    # cluster size constraint
    model.addConstrs((len(all_leaves)+9999)*(node_decision[n]-1) + min_cluster_size <= quicksum(leaf_decision[(leaf, n)] for leaf in node_to_leaves[n]) for n in list_of_ancestral_node)

    # inter-cluster constraint
    for (n, m) in itertools.combinations(list_of_ancestral_node, 2):
        model.addConstr((2 - node_decision[n] - node_decision[m])*9999 >= nodepair_to_qval[(n,m)]-fdr_cutoff)

    # within-cluster constraint
    for n in list_of_ancestral_node:
        model.addConstr(9999*(node_decision[n]-1) <= within_cluster_limit - node_to_mean_pwdist[n])

    # update model
    model.update()

    # optimize
    model.optimize()

    optimal_solutions = []
    if model.status == GRB.Status.OPTIMAL:

        # get solution pool size
        solution_pool_size = model.getAttr('SolCount')

        for sol_index in xrange(solution_pool_size):

            # get solution
            model.params.solutionNumber = sol_index
            solution = model.getAttr('Xn', leaf_decision)
            taxon_to_clusterid = {leaf:n for (leaf, n) in leaf_binary_indices if solution[(leaf, n)] == 1}

            # 1st solution is the optimal solution
            if sol_index == 0:
                optsol_obj = len(taxon_to_clusterid)
            # check if there are multiple optimal solutions
            else:
                if len(taxon_to_clusterid) < optsol_obj:
                    break
            optimal_solutions.append(taxon_to_clusterid)

        return optimal_solutions
    else:
        return 'na'