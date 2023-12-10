import sys
from Bio import Entrez
import re

# Set email address for Entrez
Entrez.email = "celinehohzm@gmail.com"
Entrez.api_key = "e8a704e313c1befab6c289256a5db4f16608"

# Parse input and output file names
input_file = sys.argv[1]
output_file = sys.argv[2]

# Read the input file
with open(input_file, "r") as infile:
    lines = infile.readlines()

# Process each intron
results = []
for i in range(0, len(lines), 4):
    # Extract information for the current intron
    intron_name = lines[i].split()[0]

    try:
        x_aligned_count = int(lines[i].split()[5])
        x_not_aligned_count = int(lines[i].split()[11])
    except ValueError:
        print(f"Error encountered when processing the line: {lines[i]}")
        print(lines[i].split())
        continue

    with_x = [int(re.search(r'\d+', taxid).group()) for taxid in lines[i+1].split(":")[1].split(",")]
    without_x = [int(re.search(r'\d+', taxid).group()) for taxid in lines[i+2].split(":")[1].split(",")]

    # with_x = [int(taxid.strip()) for taxid in lines[i+1].split(":")[1].split(",")]
    # without_x = [int(taxid.strip()) for taxid in lines[i+2].split(":")[1].split(",")]

    # Retrieve taxonomy information for each taxid
    taxonomy_with_x = {}
    for taxid in with_x:
        handle = Entrez.efetch(db="Taxonomy", id=str(taxid), retmode="xml")
        record = Entrez.read(handle)
        if record:  # Check if the record is not empty
            taxonomy_with_x[taxid] = record[0]["ScientificName"]
        else:
            print(f"No taxonomy record found for taxid: {taxid}")
        # taxonomy_with_x[taxid] = record["ScientificName"]

    taxonomy_without_x = {}
    for taxid in without_x:
        handle = Entrez.efetch(db="Taxonomy", id=str(taxid), retmode="xml")
        record = Entrez.read(handle)
        if record:  # Check if the record is not empty
            taxonomy_without_x[taxid] = record[0]["ScientificName"]
        else:
            print(f"No taxonomy record found for taxid: {taxid}")
        # taxonomy_without_x[taxid] = record["ScientificName"]

    # Add the taxonomy information to the results
    results.append({"intron_name": intron_name, "x_aligned_count": x_aligned_count, "x_not_aligned_count": x_not_aligned_count, "taxonomy_with_x": taxonomy_with_x, "taxonomy_without_x": taxonomy_without_x})

# Write the results to the output file
with open(output_file, "w") as outfile:
    for result in results:
        outfile.write(f"{result['intron_name']} - X aligned count: {result['x_aligned_count']} - X not aligned count: {result['x_not_aligned_count']}\n")
        outfile.write(f"Taxonomy IDs with X: {', '.join([str(taxid) for taxid in result['taxonomy_with_x']])}\n")
        outfile.write("Taxonomy with X:\n")
        for taxid in result['taxonomy_with_x']:
            outfile.write(f"{taxid}\t{result['taxonomy_with_x'][taxid]}\n")
        outfile.write(f"Taxonomy IDs without X: {', '.join([str(taxid) for taxid in result['taxonomy_without_x']])}\n")
        outfile.write("Taxonomy without X:\n")
        for taxid in result['taxonomy_without_x']:
            outfile.write(f"{taxid}\t{result['taxonomy_without_x'][taxid]}\n")
        outfile.write("\n")
