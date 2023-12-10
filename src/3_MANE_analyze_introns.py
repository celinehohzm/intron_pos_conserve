import sys
import re
import time
from Bio import Entrez

# Define input file and output file names
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

    # Query the NCBI Taxonomy database to find taxonomy IDs for the species
    with_x_species = list(set(re.findall(r"\[(.*?)\]", ", ".join(with_x))))
    without_x_species = list(set(re.findall(r"\[(.*?)\]", ", ".join(without_x))))
    with_x_taxids = []
    without_x_taxids = []

    # Batch the queries into groups of 10 species at a time
    with_x_batches = [with_x_species[i:i + 10] for i in range(0, len(with_x_species), 10)]
    without_x_batches = [without_x_species[i:i + 10] for i in range(0, len(without_x_species), 10)]

    for batch in with_x_batches:
        term = " OR ".join(batch)
        handle = Entrez.esearch(db="taxonomy", term=term)
        record = Entrez.read(handle)
        if record["IdList"]:
            with_x_taxids += record["IdList"]
        time.sleep(0.5)  # Sleep for 0.5 seconds between queries

    for batch in without_x_batches:
        term = " OR ".join(batch)
        handle = Entrez.esearch(db="taxonomy", term=term)
        record = Entrez.read(handle)
        if record["IdList"]:
            without_x_taxids += record["IdList"]
        time.sleep(0.5)  # Sleep for 0.5 seconds between queries

    results.append({"name": intron_name, "aligned_count": aligned_count, "unaligned_count": unaligned_count, "with_x": with_x_taxids, "without_x": without_x_taxids})

# Write results to the output file
with open(output_file, "w") as outfile:
    for result in results:
        outfile.write(f"{result['name']} - X aligned count: {result['aligned_count']} - X not aligned count: {result['unaligned_count']}\n")
        outfile.write(f"Taxonomy IDs with X: {', '.join(result['with_x'])}\n")
        outfile.write(f"Taxonomy IDs without X: {', '.join(result['without_x'])}\n\n")

