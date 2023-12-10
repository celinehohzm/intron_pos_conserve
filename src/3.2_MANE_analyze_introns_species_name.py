import sys
import re
from Bio import Entrez
from ete3 import NCBITaxa
from collections import defaultdict, Counter

def get_species_name(header):
    match = re.search(r'\[(.*?)\]', header)
    if match:
        return match.group(1)
    else:
        return "Unknown"

input_file = sys.argv[1]
mane_protein_id = sys.argv[2]
output_file = sys.argv[3]
output_file2 = sys.argv[4]
output_file3 = sys.argv[5]

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

total_species = len(sequences)
intron_positions = [i for i, c in enumerate(sequences[0]) if c == "X"]

# Initialize a dictionary to track proteins that have at least one "X"
proteins_with_x = {}

# Analyze alignments
results = []
for pos in intron_positions:
    intron_name = f"intron_{len(results) + 1}"
    aligned_count = 0
    unaligned_count = 0
    with_x = []
    without_x = []
    for i, seq in enumerate(sequences[1:]):
        slice = seq[max(0, pos - 2): min(len(seq), pos + 3)]
        slice_without_dash = [c for c in slice if c != '-']
        if "X" in slice_without_dash:
            aligned_count += 1
            with_x.append(headers[i + 1])
            proteins_with_x[headers[i + 1]] = True  # Update proteins_with_x dictionary
        else:
            unaligned_count += 1
            without_x.append(headers[i + 1])

    if mane_protein_id in with_x:
        results.append({"name": intron_name, "aligned_count": aligned_count, "unaligned_count": unaligned_count, "with_x": with_x, "without_x": without_x})

# Find common lineage for potential recent gained introns
def get_lineages(species_names):
    ncbi = NCBITaxa()
    lineages = defaultdict(list)

    for species in species_names:
        name2taxid = ncbi.get_name_translator([species])
        if name2taxid:
            taxid = name2taxid[species][0]
            lineage = ncbi.get_lineage(taxid)
            names = ncbi.get_taxid_translator(lineage)
            lineages[species] = [names[taxid] for taxid in lineage]
        else:
            print(f"Species {species} not found in the database.")

    return lineages


def best_lineage_for_species_with_x(lineages_with_x, lineages_without_x):
    # Count taxa for species with X and without X
    lineage_counts_with_x = Counter(taxon for lineage in lineages_with_x.values() for taxon in lineage)
    lineage_counts_without_x = Counter(taxon for lineage in lineages_without_x.values() for taxon in lineage)

    # Calculate the score for each taxon as count_with_x - count_without_x
    scores = {taxon: lineage_counts_with_x.get(taxon, 0) - lineage_counts_without_x.get(taxon, 0) for taxon in set(lineage_counts_with_x) | set(lineage_counts_without_x)}

    # Find the taxon with the maximum score
    best_taxon = max(scores, key=scores.get)

    return best_taxon


# Inside the block that writes to output files
with open(output_file, "w") as outfile:
    for result in results:
        # Filter out proteins that have zero "X" aligned
        result['with_x'] = [header for header in result['with_x'] if header in proteins_with_x]
        result['without_x'] = [header for header in result['without_x'] if header in proteins_with_x]
        
        if result['aligned_count'] > 0:
            outfile.write(f"{result['name']} - X aligned count: {result['aligned_count']} - X not aligned count: {result['unaligned_count']}\n")
            outfile.write("Species with X: " + ', '.join([get_species_name(header) for header in result['with_x']]) + "\n")
            outfile.write("Species without X: " + ', '.join([get_species_name(header) for header in result['without_x']]) + "\n\n")


with open(output_file2, "w") as outfile, open(output_file3, "a") as outfile3:
    for result in results:
        if result['aligned_count'] > 0:  # Add this line
            if result['unaligned_count'] / total_species > 0.1:
                outfile.write(f"{result['name']} - X aligned count: {result['aligned_count']} - X not aligned count: {result['unaligned_count']}\n")
                species_with_x.extend([get_species_name(header) for header in result['with_x']])
                species_without_x.extend([get_species_name(header) for header in result['without_x']])
                lineages_with_x = get_lineages(species_with_x)
                lineages_without_x = get_lineages(species_without_x)
                best_lineage = best_lineage_for_species_with_x(lineages_with_x, lineages_without_x)

                outfile.write("Best Lineage for Species with X but not with Species without X: " + best_lineage + "\n")
                outfile.write("Species with X: " + ', '.join([get_species_name(header) for header in result['with_x']]) + "\n")
                outfile.write("Species without X: " + ', '.join([get_species_name(header) for header in result['without_x']]) + "\n\n")

                outfile3.write(f"{mane_protein_id}_{result['name']}\t{best_lineage}\n")
                outfile3.write(f"- X aligned count: {result['aligned_count']} - X not aligned count: {result['unaligned_count']}\n")