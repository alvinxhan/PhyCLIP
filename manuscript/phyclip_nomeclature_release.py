from __future__ import division

if __name__ == '__main__':
    import argparse
    params = argparse.ArgumentParser(description='PhyCLIP cluster nomenclature.')
    params.add_argument('-t', '--tree', type=str, required=True, help='PhyCLIP output tree file.')
    params.add_argument('-p', '--percentile', type=int, default=25, help='Percentile cut-off of cluster size distinguishing between nested and child clusters (default: %(default)s).')
    params = params.parse_args()

    # read tree file
    import re, sys
    fhandle = filter(None, open(params.tree, 'rU'))
    if not re.search('#NEXUS', fhandle.pop(0)):
        sys.exit('\nERROR: Tree must be in NEXUS format.\n')
    print ('\ntree...yes.')

    id_to_taxon = {}
    # parse NEXUS tree
    for line in fhandle:
        try:
            id, taxon = re.search('^(\d+)\s+([^,;]+),*$', line.strip()).group(1, 2)
            id_to_taxon[int(id)] = taxon
        except:
            pass

        try:
            tree = re.search('\([^;]+;', line.strip()).group()
        except:
            pass

    if len(id_to_taxon) == 0 or not tree:
        sys.exit('\nERROR: Invalid tree file given.\n')
    print ('taxa...{}'.format(len(id_to_taxon)))

    # get cluster-taxa dictionaries
    taxon_to_cluster = {}
    cluster_to_taxa = {}
    for id, taxon in id_to_taxon.items():
        try:
            cluster = int(re.search('cluster(\d+)', taxon).group(1))
        except:
            cluster = -1 # clustered
        finally:
            # remove quotes and phyclip cluster annotations
            taxon = re.sub('(^\'|\'$|_cluster\d+)', '', taxon)
            id_to_taxon[id] = taxon

        if cluster > -1:
            taxon_to_cluster[taxon] = cluster
            try:
                cluster_to_taxa[cluster].append(taxon)
            except:
                cluster_to_taxa[cluster] = [taxon]
    if len(taxon_to_cluster) == 0:
        sys.exit('\nERROR: No cluster information parsable from input tree.\n')
    print ('PhyCLIP clusters...{}'.format(len(cluster_to_taxa)))

    # make tree ete3 parsable
    tree = re.sub('\[[^\]]+\]', '', tree)
    new_tree = []
    prev_end = 0
    for expr in re.finditer('[\(,](\d+):', tree):
        new_tree.append(tree[prev_end:expr.start()+1])
        new_tree.append(id_to_taxon[int(expr.group(1))])
        prev_end = expr.end()-1
    new_tree.append(tree[prev_end:])
    tree = ''.join(new_tree)

    import ete3
    tree = ete3.Tree(tree, format=5)
    all_taxa = tree.get_leaf_names()
    # check that we have correctly parsed nexus tree
    if set(all_taxa) != set(id_to_taxon.values()):
        sys.exit('\nCURIOUSER AND CURIOUSER...could not reconcile taxa information and leaves of tree. ASK ALVIN WTF HE DID!\n')

    import string
    alphabets = list(string.ascii_lowercase)

    # level-order
    # cluster numbering is in level order
    levelorder_n_to_node = {}
    levelorder_node_to_n = {}
    levelorder_leafname_to_node = {}
    levelorder_node_to_anc = {}
    levelorder_node_to_leaves = {}
    for n, node in enumerate(tree.traverse()):
        if node.is_leaf():
            levelorder_leafname_to_node[node.get_leaf_names()[0]] = node
        else:
            levelorder_n_to_node[n] = node
            levelorder_node_to_n[node] = n
            levelorder_node_to_anc[n] = [levelorder_node_to_n[anc] for anc in node.get_ancestors() if not anc.is_leaf()]
            levelorder_node_to_leaves[n] = node.get_leaf_names()

    levelorder_node_to_desc = {}
    for node, anc_list in levelorder_node_to_anc.items():
        for anc in anc_list:
            try:
                levelorder_node_to_desc[anc].append(node)
            except:
                levelorder_node_to_desc[anc] = [node]

    # clean up cluster id - must follow most descendant node
    for cluster in sorted(cluster_to_taxa.keys()):
        taxa = cluster_to_taxa[cluster]
        try:
            for node in sorted([cluster] + levelorder_node_to_desc[cluster], reverse=True):
                if set(taxa) <= set(levelorder_node_to_leaves[node]):
                    break
        except:
            continue

        if node != cluster:
            cluster_to_taxa[node] = taxa[:]
            del cluster_to_taxa[cluster]

            for taxon in taxa:
                taxon_to_cluster[taxon] = node

    # determine distribution of clusters
    clusterlen_to_frequency_count = {}
    for cluster in sorted(cluster_to_taxa.keys(), key=lambda _: len(cluster_to_taxa[_])):
        try:
            clusterlen_to_frequency_count[len(cluster_to_taxa[cluster])] += 1
        except:
            clusterlen_to_frequency_count[len(cluster_to_taxa[cluster])] = 1

    clusterlen_distribution = [i for j in [[cluster]*frequency for cluster, frequency in clusterlen_to_frequency_count.items()] for i in j]
    # regard any clusters of length <= 25 percentile lacking evidence to be a potential trajectory
    import numpy as np
    clusterlen_cutoff = int(np.percentile(clusterlen_distribution, params.percentile))
    print ('nested-cluster size cut-off ({}-percentile)...{}'.format(params.percentile, clusterlen_cutoff))

    med_clustersize = np.median(clusterlen_distribution)
    mad_clustersize = np.median([abs(_-med_clustersize) for _ in clusterlen_distribution])*1.4826
    print ('cluster size stats...median={}, mad={:.2f}, mean ={:.2f}; s.d.={:.2f}\ncluster size range...min={}; max={}'.format(med_clustersize, mad_clustersize, np.mean(clusterlen_distribution), np.std(clusterlen_distribution, ddof=1), min(clusterlen_distribution), max(clusterlen_distribution)))

    # determine ancestry of clusters
    cluster_to_desc_clusters = {}
    for cluster in sorted(cluster_to_taxa.keys()):
        try:
            desc_clusters = list(set(levelorder_node_to_desc[cluster])&set(cluster_to_taxa.keys()))
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

    import itertools
    # pre-order
    cluster_to_desc_indices = {}
    cluster_to_nomenclature = {}
    cluster_to_pwdist = {}
    cluster_to_leafdist = {}
    nomenclature_index = 1

    print ('\nassigning...\n{:<10}{}'.format('cluster', 'nomenclature'))
    for node in tree.traverse(strategy='preorder'):
        # skip if leaf node
        if node.is_leaf():
            continue
        n = levelorder_node_to_n[node]

        # skip if n in not a cluster
        if n not in cluster_to_taxa.keys():
            continue

        try:
            immediate_parent = sorted(cluster_to_anc_clusters[n])[-1]
        except:
            immediate_parent = -1

        # cluster has no ancestor cluster
        if immediate_parent < 0:
            nomenclature = str(nomenclature_index)
            cluster_to_desc_indices[nomenclature] = {'nested':0, 'child':1}
            cluster_to_nomenclature[n] = nomenclature
            nomenclature_index += 1
            print ('{:<10}{}'.format(n, nomenclature))
            continue
        else:
            parent_nomenclature = cluster_to_nomenclature[immediate_parent]
            clusterlen = len(cluster_to_taxa[n])
            if clusterlen <= clusterlen_cutoff and n not in cluster_to_desc_clusters:
                # check if all leaves in n have distance <= farthest leaf of parent node
                try:
                    farthest_dist = max(cluster_to_leafdist[immediate_parent])
                except:
                    cluster_to_leafdist[immediate_parent] = [levelorder_leafname_to_node[taxon].get_distance(levelorder_n_to_node[immediate_parent]) for taxon in cluster_to_taxa[immediate_parent]]
                    farthest_dist = max(cluster_to_leafdist[immediate_parent])

                if all(levelorder_leafname_to_node[taxon].get_distance(levelorder_n_to_node[immediate_parent]) <= farthest_dist for taxon in cluster_to_taxa[n]):
                    # nested
                    nomenclature = '{}({})'.format(parent_nomenclature, alphabets[cluster_to_desc_indices[parent_nomenclature]['nested']])
                    cluster_to_desc_indices[parent_nomenclature]['nested'] += 1
                    cluster_to_nomenclature[n] = nomenclature
                else:
                    try:
                        mean_pwdist = np.mean(cluster_to_pwdist[immediate_parent])
                    except:
                        cluster_to_pwdist[immediate_parent] = [levelorder_leafname_to_node[i].get_distance(levelorder_leafname_to_node[j]) for i, j in itertools.combinations(cluster_to_taxa[immediate_parent], 2)]
                        mean_pwdist = np.mean(cluster_to_pwdist[immediate_parent])

                    additional_dist = [levelorder_leafname_to_node[i].get_distance(levelorder_leafname_to_node[j]) for i, j in itertools.product(cluster_to_taxa[n], cluster_to_taxa[immediate_parent])]
                    #print n, np.mean(additional_dist + cluster_to_pwdist[immediate_parent]), mean_pwdist
                    if np.mean(additional_dist + cluster_to_pwdist[immediate_parent]) <= mean_pwdist:
                        # nested
                        nomenclature = '{}({})'.format(parent_nomenclature, alphabets[cluster_to_desc_indices[parent_nomenclature]['nested']])
                        cluster_to_desc_indices[parent_nomenclature]['nested'] += 1
                        cluster_to_nomenclature[n] = nomenclature
                    else:
                        # not nested
                        nomenclature = '{}.{}'.format(parent_nomenclature, cluster_to_desc_indices[parent_nomenclature]['child'])
                        cluster_to_desc_indices[parent_nomenclature]['child'] += 1
                        cluster_to_desc_indices[nomenclature] = {'nested': 0, 'child': 1}
                        cluster_to_nomenclature[n] = nomenclature
            else:
                # not nested
                nomenclature = '{}.{}'.format(parent_nomenclature, cluster_to_desc_indices[parent_nomenclature]['child'])
                cluster_to_desc_indices[parent_nomenclature]['child'] += 1
                cluster_to_desc_indices[nomenclature] = {'nested': 0, 'child': 1}
                cluster_to_nomenclature[n] = nomenclature
            print ('{:<10}{}'.format(n, nomenclature))

    # write tree output
    print ('\nwriting outputs...')
    with open('nomenclature_{}'.format(params.tree), 'w') as output:
        output.write('#NEXUS\nBegin taxon;\n\tDimensions ntax={};\n\t\tTaxlabels\n'.format(len(all_taxa)))

        tree_string = tree.write(format=5)
        for id in sorted(id_to_taxon.keys()):
            taxon = id_to_taxon[id]

            if taxon in taxon_to_cluster:
                cluster = taxon_to_cluster[taxon]
                nomenclature = cluster_to_nomenclature[cluster]
                tree_string = tree_string.replace(taxon, '{}[&CLUSTER=\'{}\']'.format(id, nomenclature))
                taxon = "'{}_cluster{}*'".format(taxon, nomenclature)
            else:
                tree_string = tree_string.replace(taxon, "'{}'".format(id))

            output.write('\t\t\t{}\n'.format(taxon))
        output.write('\t\t\t;\nEnd;\nBegin trees;\n\tTranslate\n')

        for id in sorted(id_to_taxon.keys()):
            taxon = id_to_taxon[id]

            if taxon in taxon_to_cluster:
                cluster = taxon_to_cluster[taxon]
                nomenclature = cluster_to_nomenclature[cluster]
                # tree_string = tree_string.replace(taxon, '{}[&CLUSTER=\'{}\']'.format(id, nomenclature))
                taxon = "'{}_cluster{}*'".format(taxon, nomenclature)



            output.write('\t\t{:>5} {}{}\n'.format(id, taxon, ',' if id < len(all_taxa) else ''))

        output.write(';\ntree TREE = {}\nEnd;'.format(tree_string))

    # write cluster text output
    unclustered_taxa = list(set(all_taxa)-set([i for j in cluster_to_taxa.values() for i in j]))
    with open('nomenclature-cluster_{}.txt'.format(params.tree), 'w') as output:
        output.write('cluster\tnomenclature\ttaxon\n')
        for cluster, taxa in cluster_to_taxa.items():
            nomenclature = cluster_to_nomenclature[cluster]
            for taxon in taxa:
                output.write('{}\t{}\t{}\n'.format(cluster, nomenclature, taxon))

        for taxon in unclustered_taxa:
            output.write('unclustered\tunclustered\t{}\n'.format(taxon))

    print ('...done.\n')
