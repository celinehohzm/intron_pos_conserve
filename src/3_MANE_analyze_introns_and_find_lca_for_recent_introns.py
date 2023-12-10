import sys
import re
from Bio import Entrez
from ete3 import NCBITaxa
from collections import defaultdict

def get_species_name(header):
    match = re.search(r'\[(.*?)\]', header)
    if match:
        species_name = match.group(1)
        return species_name
    else:
        return "Unknown"

input_file = sys.argv[1]
output_file = sys.argv[2]
protein_id = sys.argv[3]
output_file_lca = sys.argv[4]

# Initialize the Entrez module from Biopython and provide your email address
Entrez.email = "celinehohzm@gmail.com"
ncbi = NCBITaxa()

# Read the input file
with open(input_file, "r") as infile:
    lines = infile.readlines()

# Separate headers and sequences
headers = []
sequences = []
seq = ""

for line in lines:
    if line.startswith(">"):
        headers.append(line.strip())
        if seq:
            sequences.append(seq)
            seq = ""
    else:
        seq += line.strip()

if seq:
    sequences.append(seq)

# Find intron positions in the first alignment
intron_positions = [i for i, c in enumerate(sequences[0]) if c == "X"]

# Analyze alignments
results = []
for pos in intron_positions:
    intron_name = f"intron_{len(results) + 1}"
    aligned_count = 0
    unaligned_count = 0
    with_x = []
    without_x = []
    for i, seq in enumerate(sequences[1:]):
        # Take slice from seq[pos - 2] to seq[pos + 2] inclusive, and check if "X" is in it
        slice = seq[max(0, pos - 2): min(len(seq), pos + 3)]
        slice_without_dash = [c for c in slice if c != '-']
        species_name = get_species_name(headers[i + 1])
        if "X" in slice_without_dash:
            aligned_count += 1
            with_x.append(species_name)
        else:
            unaligned_count += 1
            without_x.append(species_name)

    results.append({"name": intron_name, "aligned_count": aligned_count, "unaligned_count": unaligned_count, "with_x": with_x, "without_x": without_x})

# Write results to the output file
with open(output_file, "w") as outfile:
    for result in results:
        outfile.write(f"{result['name']} - X aligned count: {result['aligned_count']} - X not aligned count: {result['unaligned_count']}\n")
        outfile.write("Species with X: " + ', '.join(result['with_x']) + "\n")
        outfile.write("Species without X: " + ', '.join(result['without_x']) + "\n\n")


def get_lca(ncbi, species):
    # Get taxids for species
    taxids = [ncbi.get_name_translator([name])[name][0] for name in species if ncbi.get_name_translator([name])]
    # Use get_topology to get a pruned NCBI taxonomy tree
    tree = ncbi.get_topology(taxids, intermediate_nodes=True)
    # The LCA would be the root of the pruned tree
    lca = tree.get_tree_root()
    print("lca_name", tree.get_common_ancestor(taxids).name)
    lca_name = ncbi.get_taxid_translator([int(lca.name)])[int(lca.name)]
    return lca_name

# Inside the loop
with open(output_file_lca, "w") as outfile_lca:
    for result in results:
        if len(result['without_x']) / (len(result['with_x']) + len(result['without_x'])) > 0.2:
            species_with_x = [get_species_name(header) for header in result['with_x']]
            taxids = [ncbi.get_name_translator([name])[name][0] for name in species_with_x if ncbi.get_name_translator([name])]
            # Find common ancestor           
            tree = ncbi.get_topology(taxids)
            common_ancestor = tree.name
            # If you want to know the rank (e.g. kingdom, phylum, etc.) of the common ancestor
            rank = ncbi.get_rank([common_ancestor])
            print("The common ancestor is at the rank: ", rank[common_ancestor])
    outfile_lca.write(f"{protein_id}\t{result['name']}\t{' '.join(result['with_x'])}\n{rank}\n")

