import re

class phyclip_output(object):
    """
    Print output cluster and tree files
    """
    def __init__(self, tree_string=None, taxon_to_clusterid=None, taxon_list=None, outfname=None, clean_up_status=None, sensitivity_subsumed_taxa_to_clusterid=False, nosub_taxa_to_clusterid=False):
        self.tree_string = tree_string
        self.taxon_to_clusterid = taxon_to_clusterid
        self.taxon_list = taxon_list
        self.sensitivity_subsumed_taxa_to_clusterid = sensitivity_subsumed_taxa_to_clusterid
        self.nosub_taxa_to_clusterid = nosub_taxa_to_clusterid
        self.outfname = outfname
        self.clean_up_status = clean_up_status

    def cluster_output(self):

        curr_tree_string = self.tree_string

        with open('cluster_{}_{}.txt'.format(self.clean_up_status, self.outfname), 'w') as output:

            output.write('CLUSTER\tTAXA\tSensitivity-induced (original cluster-node-id)\tSub-cluster subsumed into parent (original cluster-node-id)\n')
            for taxon, clusterid in self.taxon_to_clusterid.items():


                output.write('{}\t{}'.format(clusterid, taxon))
                if self.sensitivity_subsumed_taxa_to_clusterid != False:
                    output.write('\t{}'.format(self.sensitivity_subsumed_taxa_to_clusterid[taxon] if taxon in self.sensitivity_subsumed_taxa_to_clusterid else ''))
                if self.nosub_taxa_to_clusterid != False:
                    output.write('\t{}'.format(self.nosub_taxa_to_clusterid[taxon] if taxon in self.nosub_taxa_to_clusterid else ''))
                output.write('\n')

                taxon_start_index = curr_tree_string.find(taxon)
                if taxon_start_index < 0:
                    raise SystemExit('\nERROR: Problem parsing taxon ({}) in tree string.\n'.format(taxon))
                else:
                    curr_tree_string = curr_tree_string.replace(taxon, "'{}'[&CLUSTER='{}']".format(re.sub("(^'|'$)", "", taxon), clusterid))

        return curr_tree_string

    def figtree_output(self, modified_tree_string):

        with open('tree_{}_{}.tre'.format(self.clean_up_status, self.outfname), 'w') as output:
            output.write('#NEXUS\nBegin taxon;\n\tDimensions ntax={};\n\t\tTaxlabels\n'.format(len(self.taxon_list)))
            for taxon in self.taxon_list:
                if taxon in self.taxon_to_clusterid:
                    output.write("\t\t\t'{}_cluster{}'\n".format(re.sub("(^'|'$)", "", taxon), self.taxon_to_clusterid[taxon]))
                else:
                    output.write('\t\t\t{}\n'.format(taxon))

            output.write('\t\t\t;\nEnd;\nBegin trees;\n\tTranslate\n')

            for i, taxon in enumerate(self.taxon_list):
                if taxon in self.taxon_to_clusterid:
                    output.write("\t\t{:>4} '{}_cluster{}'{}\n".format(i+1, re.sub("(^'|'$)", '', taxon), self.taxon_to_clusterid[taxon], '' if i+1 == len(self.taxon_list) else ','))
                else:
                    output.write("\t\t{:>4} '{}'{}\n".format(i+1, re.sub("(^'|'$)", '', taxon), '' if i+1 == len(self.taxon_list) else ','))

                taxon_start_index = modified_tree_string.find("'{}'".format(re.sub("(^'|'$)", '', taxon)))
                if taxon_start_index < 0:
                    raise SystemExit('\nERROR: Problem parsing taxon ({}) in tree string.\n'.format(taxon))
                else:
                    modified_tree_string = modified_tree_string.replace("'{}'".format(re.sub("(^'|'$)", "", taxon)), str(i+1))

            output.write(';\ntree TREE1 = {}\nEnd;\n'.format(modified_tree_string))