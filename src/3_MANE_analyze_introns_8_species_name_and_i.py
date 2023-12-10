import sys
import re
from Bio import Entrez

def get_species_name(header):
    match = re.search(r'\[(.*?)\]', header)
    if match:
        species_name = match.group(1)
        words = species_name.split()
        if len(words) >= 2:
            species_id = words[0][:3] + '_' + words[1][:2]
        else:
            species_id = species_name[:5].replace(" ", "_")
        return species_id
    else:
        return "Unknown"

input_file = sys.argv[1]
output_file = sys.argv[2]

# Initialize the Entrez module from Biopython and provide your email address
Entrez.email = "celinehohzm@gmail.com"

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
        if seq[pos] == "X":
            aligned_count += 1
            with_x.append(headers[i + 1])
        else:
            unaligned_count += 1
            without_x.append(headers[i + 1])

    results.append({"name": intron_name, "aligned_count": aligned_count, "unaligned_count": unaligned_count, "with_x": with_x, "without_x": without_x})

# Write results to the output file
with open(output_file, "w") as outfile:
    for result in results:
        outfile.write(f"{result['name']} - X aligned count: {result['aligned_count']} - X not aligned count: {result['unaligned_count']}\n")
        outfile.write("Species with X: " + ', '.join([get_species_name(header) for header in result['with_x']]) + "\n")
        outfile.write("Species without X: " + ', '.join([get_species_name(header) for header in result['without_x']]) + "\n\n")

