import sys
import re
import time
import subprocess
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

# Create a Kraken database
# !kraken2-build --standard --db kraken2-db

# Classify sequences using Kraken
for result in results:
    with_x_taxids = result['with_x']
    without_x_taxids = result['without_x']
    all_taxids = with_x_taxids + without_x_taxids

    # Write the sequences to a temporary file
    with open("tmp.fa", "w") as tmpfile:
        for header, sequence in zip(headers[1:], sequences[1:]):
            tmpfile.write(f">{header}\n{sequence}\n")

    # Classify the sequences with Kraken
    # subprocess.run(['/ccb/sw/bin/kraken2', '--db', 'kraken2-db', 'tmp.fa', '>', 'tmp.out'])
    subprocess.run('/ccb/sw/bin/kraken2 --db kraken2-db tmp.fa > tmp.out', shell=True)


    # Parse the Kraken output and create a taxonomy tree
    with_x_tree = {}
    without_x_tree = {}

    with open("tmp.out", "r") as krakenout:
        for line in krakenout:
            fields = line.strip().split("\t")
            taxid = fields[2]
            if taxid in all_taxids:
                if fields[0] in with_x_taxids:
                    if taxid not in with_x_tree:
                        with_x_tree[taxid] = {}
                    taxon = Entrez.read(Entrez.efetch(db="taxonomy", id=taxid))
                    lineage = taxon["LineageEx"]
                    for i in range(len(lineage)):
                        rank = lineage[i]["Rank"]
                        name = lineage[i]["ScientificName"]
                        if rank in ["species", "genus", "family", "order", "class", "phylum", "superkingdom"]:
                            if name not in with_x_tree[taxid]:
                                with_x_tree[taxid][name] = {}
                            if i < len(lineage) - 1:
                                with_x_tree[taxid] = with_x_tree[taxid][name]
                            else:
                                with_x_tree[taxid][name] = {"count": fields[1], "children": {}}
                else:
                    if taxid not in without_x_tree:
                        without_x_tree[taxid] = {}
                    taxon = Entrez.read(Entrez.efetch(db="taxonomy", id=taxid))
                    lineage = taxon["LineageEx"]
                    for i in range(len(lineage)):
                        rank = lineage[i]["Rank"]
                        name = lineage[i]["ScientificName"]
                        if rank in ["species", "genus", "family", "order", "class", "phylum", "superkingdom"]:
                            if name not in without_x_tree[taxid]:
                                without_x_tree[taxid][name] = {}
                            if i < len(lineage) - 1:
                                without_x_tree[taxid] = without_x_tree[taxid][name]
                            else:
                                without_x_tree[taxid][name] = {"count": fields[1], "children": {}}

    # Write the taxonomy tree to a file
    with open(f"{result['name']}_with_x_taxonomy.txt", "w") as withxfile:
        withxfile.write("WITH X:\n\n")
        withxfile.write(str(with_x_tree))

    with open(f"{result['name']}_without_x_taxonomy.txt", "w") as withoutxfile:
        withoutxfile.write("WITHOUT X:\n\n")
        withoutxfile.write(str(without_x_tree))
