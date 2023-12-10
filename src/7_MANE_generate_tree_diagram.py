import matplotlib
matplotlib.use('Agg')
import sys
import glob
import re
from ete3 import NCBITaxa, Tree, TreeStyle, NodeStyle, RectFace, TextFace

def generate_tree(species_names):
    ncbi = NCBITaxa()

    name2taxid = ncbi.get_name_translator(species_names)
    taxid2name = {taxid: name for name, taxid_list in name2taxid.items() for taxid in taxid_list}
    taxids = [taxid for taxid_list in name2taxid.values() for taxid in taxid_list]
    tree = ncbi.get_topology(taxids)

    for node in tree.traverse():
        taxid = int(node.name)
        node.add_features(taxid=taxid)
        if node.is_leaf():
            node.name = taxid2name[taxid]
            # node.name.margin_right = 3  # Added bottom margin to move the label a little further from the branch
        lineage_taxids = ncbi.get_lineage(node.taxid)
        lineage_names = ncbi.translate_to_names(lineage_taxids)
        node.add_features(common_name=lineage_names[-1])

    return tree

def read_species_names(file_path):
    with open(file_path, 'r') as file:
        return [line.strip() for line in file]

def color_tree(t, base_color_annot_file, output_file):
    unsorted_color_annot_files = glob.glob(f"{base_color_annot_file}_*.txt")
    color_annot_files = sorted(unsorted_color_annot_files, key=lambda x: int(re.search(r'intron_(\d+)', x).group(1)))
    print(color_annot_files)
    column = 1

    nstyle = NodeStyle()
    nstyle["hz_line_width"] = 2
    nstyle["vt_line_width"] = 2

    for n in t.traverse():
        n.set_style(nstyle)
        if not n.is_leaf():
            common_name_face = TextFace(n.common_name, fsize=10, fgcolor="black")
            n.add_face(common_name_face, column=0, position="branch-top")

    for color_annot_file in color_annot_files:
        intron_number = re.search(r'intron_(\d+)', color_annot_file).group(1)
        intron_number_face = TextFace(intron_number, fsize=10, fgcolor="black")
        with open(color_annot_file, 'r') as file:
            for _ in range(5):
                next(file)
            for line in file:
                if len(line.rstrip().rsplit(None, 1)) > 1:
                    species, color = line.rstrip().rsplit(None, 1)  # Split by last whitespace only
                    # species, color = line.strip().split()  # Split by whitespace
                    nodes = t.search_nodes(name=species)
                    if nodes:
                        node = nodes[0]
                        if color == '#FF0000':
                            face = RectFace(width=15, height=15, fgcolor="Red", bgcolor="Red")
                        else:
                            face = RectFace(width=15, height=15, fgcolor="Blue", bgcolor="Blue")
                        node.add_face(face, column=column, position="aligned")
                        node.add_face(intron_number_face, column=column, position="aligned")
                    else:
                        print(f"Node not found for species: {species}")
        column += 1

    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.show_branch_length = False
    ts.show_branch_support = True

    t.render(output_file, w=183, units="mm", tree_style=ts)

if len(sys.argv) < 4:
    print('Please provide an input file name for the species names, a base name for the color annotation files, and an output file name for the PDF as command-line arguments.')
else:
    species_names = read_species_names(sys.argv[1])
    t = generate_tree(species_names)
    color_tree(t, sys.argv[2], sys.argv[3])

