import re
import sys
from Bio import Entrez
from collections import defaultdict
from ete3 import NCBITaxa

ncbi = NCBITaxa()

species_to_taxgroup = {}  # Dictionary to store results

def get_species_name(header):
    match = re.search(r'\[(.*?)\]', header)
    if match:
        species_name = match.group(1)
        return species_name
    else:
        return "Unknown"

def get_taxonomic_group(species_name):
    if species_name in species_to_taxgroup:
        return species_to_taxgroup[species_name]

    name2taxid = ncbi.get_name_translator([species_name])
    taxid = name2taxid[species_name][0] if species_name in name2taxid else None

    tax_class = "Unknown"
    tax_order = "Unknown"
    tax_genus = "Unknown"

    if taxid:
        lineage = ncbi.get_lineage(taxid)
        names = ncbi.get_taxid_translator(lineage)

        for taxid in lineage:
            rank = ncbi.get_rank([taxid])
            if rank[taxid] == "class":
                tax_class = names[taxid]
            if rank[taxid] == "order":
                tax_order = names[taxid]
            if rank[taxid] == "genus":
                tax_genus = names[taxid]

        if tax_class == "Mammalia":
            if tax_order == "Primates":
                species_to_taxgroup[species_name] = f"Mammalia {tax_order} {tax_genus}"
            else:
                species_to_taxgroup[species_name] = f"Mammalia {tax_order}"
        else:
            species_to_taxgroup[species_name] = tax_class
    else:
        species_to_taxgroup[species_name] = "Unknown"

    return species_to_taxgroup[species_name]

input_file = sys.argv[1]
output_file = sys.argv[2]
protein_id = sys.argv[3]
output_file4 = sys.argv[4]

Entrez.email = "celinehohzm@gmail.com"

with open(input_file, "r") as infile:
    lines = infile.readlines()

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

intron_positions = [i for i, c in enumerate(sequences[0]) if c == "X"]

results = []

for pos in intron_positions:
    intron_name = f"intron_{len(results) + 1}"
    aligned_count = 0
    unaligned_count = 0
    tax_groups_with_x = {}
    tax_groups_without_x = {}

    for i, seq in enumerate(sequences[1:]):
        slice = seq[max(0, pos - 2): min(len(seq), pos + 3)]
        slice_without_dash = [c for c in slice if c != '-']
        species_name = get_species_name(headers[i + 1])
        tax_group = get_taxonomic_group(species_name)
        if "X" in slice_without_dash:
            aligned_count += 1
            tax_groups_with_x[tax_group] = tax_groups_with_x.get(tax_group, 0) + 1
        else:
            unaligned_count += 1
            tax_groups_without_x[tax_group] = tax_groups_without_x.get(tax_group, 0) + 1

    results.append({"name": intron_name, "aligned_count": aligned_count, "unaligned_count": unaligned_count, "tax_groups_with_x": tax_groups_with_x, "tax_groups_without_x": tax_groups_without_x})

group_order = {
    'Class': 0,
    'Order': 1,
    'Genus': 2,
    'Mammalia Primates': 3,
    # Add more groups if needed...
}

def group_sort_key(group):
    return group_order.get(group, 999)  # Unknown groups will be put at the end

with open(output_file, "w") as outfile:
    for result in results:
        species_with_x = sorted(result['tax_groups_with_x'].keys(), key=group_sort_key)
        species_without_x = sorted(result['tax_groups_without_x'].keys(), key=group_sort_key)
        outfile.write(f"{result['name']} - X aligned count: {result['aligned_count']} - X not aligned count: {result['unaligned_count']}\n")
        outfile.write("Taxonomic groups with X: " + ', '.join(species_with_x) + "\n")
        outfile.write("Taxonomic groups without X: " + ', '.join(species_without_x) + "\n\n\n\n")

results_with_LCA = []
total_species = len(sequences) - 1

for result in results:
    if result['unaligned_count'] / total_species > 0.1:
        print("result['unaligned_count']", result['unaligned_count'])
        species_with_x = sorted(list(result['tax_groups_with_x'].keys()), key=group_sort_key)
        count_species_with_x = sum(result['tax_groups_with_x'].values())
        species_without_x = sorted(list(result['tax_groups_without_x'].keys()), key=group_sort_key)
        count_species_without_x = sum(result['tax_groups_without_x'].values())

        with open(output_file4, "a") as outfile4:
            outfile4.write(f"Protein_id: {protein_id}, ")
            outfile4.write(f"Intron_number: {result['name'].split('_')[1]}\n")
            outfile4.write("Taxonomic Groups with X: " + ', '.join(species_with_x) + "\n")
            outfile4.write(f"Number of Species in Taxonomic Groups with X: {count_species_with_x}\n")
            outfile4.write("Taxonomic Groups without X: " + ', '.join(species_without_x) + "\n")
            outfile4.write(f"Number of Species in Taxonomic Groups without X: {count_species_without_x}\n\n\n\n")

