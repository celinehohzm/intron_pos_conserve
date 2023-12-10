import sys
import glob
from ete3 import Tree, TreeStyle, NodeStyle

# Get the tree file name and base name for the color annotation files from the command line
tree_file = sys.argv[1]
base_color_annot_file = sys.argv[2]

# Load your tree
t = Tree(tree_file)

# Define node styles for the two groups
style1 = NodeStyle()
style1["bgcolor"] = "Red"
style2 = NodeStyle()
style2["bgcolor"] = "Blue"

# Use glob to find all color annotation files that match the base name
color_annot_files = glob.glob(f"{base_color_annot_file}_*.txt")

# Iterate over color annotation files
for color_annot_file in color_annot_files:
    # Open the color annotation file
    with open(color_annot_file, 'r') as file:
        # Skip the header lines
        for _ in range(5):
            next(file)

        # Process the rest of the lines
        for line in file:
            # Split the line into species and color
            species, color = line.strip().split('\t')
            # Find the node for this species
            node = t&species
            # Assign the correct style based on the color
            if color == '#FF0000':
                node.set_style(style1)
            else:
                node.set_style(style2)

# Render the tree to a PDF file
ts = TreeStyle()
ts.show_leaf_name = True
ts.show_branch_length = True
ts.show_branch_support = True
t.render(f'{tree_file}_colored.pdf', w=183, units="mm", tree_style=ts)

