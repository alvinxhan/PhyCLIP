from gurobipy import *
import itertools

def gurobi_solver(node_to_leaves, leaf_to_ancestors, list_of_ancestral_node, nodepair_to_qval, node_to_mean_pwdist, within_cluster_limit, min_cluster_size, fdr_cutoff, verbose):

    # set up indices
    leaf_binary_indices = [(leaf, n) for n, leaves in node_to_leaves.items() for leaf in leaves]
    curr_taxon_list = leaf_to_ancestors.keys()[:]

    # model
    model = Model('NewModel')

    # verbose
    print ('solver...gurobi.')
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
    model.addConstrs(quicksum(leaf_decision[(leaf, n)] for n in leaf_to_ancestors[leaf]) <= 1 for leaf in curr_taxon_list)

    # cluster size constraint
    model.addConstrs((len(curr_taxon_list)+9999)*(node_decision[n]-1) <= quicksum(leaf_decision[(leaf, n)] for leaf in node_to_leaves[n]) - min_cluster_size for n in list_of_ancestral_node)

    # inter-cluster constraint
    for (n, m) in itertools.combinations(list_of_ancestral_node, 2):
        model.addConstr((2 - node_decision[n] - node_decision[m])*9999 >= nodepair_to_qval[(n,m)]-fdr_cutoff)

    # within-cluster constraint
    for n in list_of_ancestral_node:
        model.addConstr(9999*(node_decision[n]-1) <= within_cluster_limit - node_to_mean_pwdist[n])

    # leaf always chooses the youngest cluster-able node
    for leaf, ancestors_of_leaf in leaf_to_ancestors.items():
        for n, m in itertools.combinations(ancestors_of_leaf, 2):
            if n < m:
                model.addConstr(leaf_decision[(leaf, n)] <= 2 - node_decision[n] - node_decision[m])
            else:
                model.addConstr(leaf_decision[(leaf, m)] <= 2 - node_decision[m] - node_decision[n])

    # update model
    model.update()

    # optimize
    model.optimize()

    if model.status == GRB.Status.OPTIMAL:
        # get solution
        solution = model.getAttr('x', leaf_decision)
        return {leaf: n for (leaf, n) in leaf_binary_indices if solution[(leaf, n)] == 1}
    else:
        return 'na'