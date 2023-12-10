import sys
from ete3 import NCBITaxa

def generate_tree(species_names, output_file_newick):
    ncbi = NCBITaxa()

    # Get a dictionary of species names to taxids
    name2taxid = ncbi.get_name_translator(species_names)

    # Create a reverse dictionary of taxids to species names
    taxid2name = {taxid: name for name, taxid_list in name2taxid.items() for taxid in taxid_list}

    # Collect the corresponding taxids
    taxids = [taxid for taxid_list in name2taxid.values() for taxid in taxid_list]

    # Get the corresponding tree
    tree = ncbi.get_topology(taxids)

    # Replace taxid with species names for all leaf nodes
    for leaf in tree:
        taxid = int(leaf.name)  # Convert the string to integer
        if taxid in taxid2name:
            leaf.name = taxid2name[taxid]
        else:
            print("not in dic", leaf.name)

    # Write the Newick format of the tree to a file
    tree.write(outfile=output_file_newick, format=1)

def read_species_names(file_path):
    with open(file_path, 'r') as file:
        return [line.strip() for line in file]

# Check that filenames were given as command-line arguments
if len(sys.argv) < 2:
    print('Please provide an input file name, an output file name for the Newick format, and an output file name for the PDF as command-line arguments.')
else:
    species_names = read_species_names(sys.argv[1])
    generate_tree(species_names, sys.argv[2])


