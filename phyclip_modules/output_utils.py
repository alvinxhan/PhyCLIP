from __future__ import division
import re
import ete3

class phyclip_output(object):
    """
    Print output cluster and tree files
    """
    def __init__(self, ori_tree_string=None, tree_string=None, taxon_to_clusterid=None, clusterid_to_taxa=None, taxon_list=None, outfname=None, clean_up_status=None, sensitivity_subsumed_taxa_to_clusterid=False, nosub_taxa_to_clusterid=False):
        self.ori_tree_string = ori_tree_string
        self.tree_string = tree_string
        self.taxon_to_clusterid = taxon_to_clusterid
        self.clusterid_to_taxa = clusterid_to_taxa
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


                output.write('{}\t{}'.format(clusterid, re.sub("(^'|'$)", '', taxon)))
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

    def prior_output(self, pc_input):
        with open('prior_{}.txt'.format(self.outfname), 'w') as output:
            output.write('Prior_cluster_id\tTaxon\tCurrent_cluster_id\n')
            for p, leaves_in_pc in pc_input.items():
                for leaf in leaves_in_pc:
                    output.write('{}\t{}\t'.format(p, leaf))
                    try:
                        output.write('{}\n'.format(self.taxon_to_clusterid[leaf]))
                    except:
                        output.write('\n')

    def generate_heatmap(self, ete3_tree, *args):

        for node in ete3_tree:
            if node.is_leaf():
                taxon = node.name

                tc_index = 0
                for args_index in range(0, len(args), 2):
                    taxon_to_classification = args[args_index]
                    class_to_color = args[args_index+1]

                    try:
                        classification = taxon_to_classification[taxon]
                        color = class_to_color[classification]
                    except:
                        classification = ''
                        color = 'white'

                    class_face = ete3.TextFace(classification, ftype='Arial', fsize=2, fgcolor='white')
                    class_face.background.color = color
                    class_face.hz_align = 1
                    class_face.vt_align = 1
                    class_face.margin_left = 5
                    class_face.margin_right = 5
                    node.add_face(class_face, tc_index, "aligned")
                    tc_index += 1

        return ete3_tree

    def generate_color_scheme(self, class_list):
        from colorsys import hls_to_rgb
        import random
        random.seed(666) # maintain consistent h_list shuffle

        def hls2hex(h, l, s):
            ''' Converts from HLS to RGB '''
            return '#%02x%02x%02x' %tuple(map(lambda x: int(x*255), hls_to_rgb(h, l, s)))
        def get_color(h, s=0.5, l=0.5):
            ''' Returns a RGB color code with the specific Hue, Saturation and
            lightness.'''
            return hls2hex(h, l, s)

        # get widest increment of color hues
        h_increment = 1/(len(class_list)-1)
        h_list = [0+(cl*h_increment) for cl in xrange(len(class_list))]
        random.shuffle(h_list)
        return {classification:get_color(h_list[cl]) for cl, classification in enumerate(class_list)}

    def ete3_pdf_tree_output(self, pc_input):
        output_tree = ete3.Tree(self.ori_tree_string)
        output_tree.ladderize()

        # treestyle
        ts = ete3.TreeStyle()
        ts.show_leaf_name = False

        # generate color scheme for taxon_to_clusterid
        clusterid_to_color = self.generate_color_scheme(self.clusterid_to_taxa.keys())

        clusterid_traversed = []
        for n, node in enumerate(output_tree.traverse(strategy='levelorder')):
            if n == 0:
                ts.scale_length = float('{:.3f}'.format(node.get_farthest_leaf()[-1]/10))

            if node.is_leaf():
                # color branches
                ns = ete3.NodeStyle()
                ns["size"] = 0 # no node shape

                taxon = node.name
                if taxon in self.taxon_to_clusterid:
                    clusterid = self.taxon_to_clusterid[taxon]
                    ns["hz_line_color"] = clusterid_to_color[clusterid]
                    # write taxon names aligned to the right
                    taxon_name = ete3.TextFace(taxon, ftype='Arial', fsize=2, bold=True, fgcolor=clusterid_to_color[clusterid])
                else:
                    # write taxon names aligned to the right
                    taxon_name = ete3.TextFace(taxon, ftype='Arial', fsize=2, fstyle="italic")

                node.set_style(ns)
                taxon_name.margin_left = 2
                node.add_face(taxon_name, column=0, position='branch-right')

            else:
                ns = ete3.NodeStyle()
                ns["size"] = 0 # no node shape

                # set node style
                node.set_style(ns)

        heatmap_headers = ['Cluster-ID']
        if pc_input != False:
            # prior heat-map
            heatmap_headers.append('Prior')
            taxon_to_pc = {taxon:p for p, leaves_of_pc in pc_input.items() for taxon in leaves_of_pc}

            pc_to_color = self.generate_color_scheme(list(set(taxon_to_pc.values())))
            output_tree = self.generate_heatmap(output_tree, self.taxon_to_clusterid, clusterid_to_color, taxon_to_pc, pc_to_color)
        else:
            # normal
            output_tree = self.generate_heatmap(output_tree, self.taxon_to_clusterid, clusterid_to_color)

        # heatmap header
        for lh_index, legend_header in enumerate(heatmap_headers):
            header_face = ete3.TextFace(legend_header, ftype='Arial', fsize=2)
            header_face.hz_align = 1
            header_face.vt_align = 1
            header_face.margin_left = 5
            header_face.margin_right = 5
            ts.aligned_header.add_face(header_face, lh_index)

        # render as pdf
        output_tree.render('{}{}.pdf'.format('' if pc_input == False else 'prior_', self.outfname), tree_style=ts)

