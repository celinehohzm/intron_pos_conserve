#!/usr/bin/env python3
import sys

if len(sys.argv) != 3:
    print("Usage: python script_name.py <MANE_refseq_protein_file> <gene2refseq_file>")
    sys.exit(1)

input_file = sys.argv[1]
gene2refseq_file = sys.argv[2]

# Placeholder: Read the gene2refseq file and map protein accessions to taxids
protein_to_taxid = {}
with open(gene2refseq_file, 'r') as f:
    # Skipping header line
    next(f)
    for line in f:
        parts = line.strip().split("\t")
        tax_id = parts[0]
        #print(parts[5])
        protein_accession = parts[5]  # Get only the base accession without the version
        protein_to_taxid[protein_accession] = tax_id

# Read the input file and process data
with open(input_file, 'r') as f:
    data = {}
    added_taxids_per_key = {}  # Dictionary to keep track of taxids that have been processed for each key
    for line in f:
        line = line.strip()
        key, value = line.split("\t")
        # Extract the accession number from the value
        protein_accession = value.split("|")[1]  # Get only the base accession without the version
        tax_id = protein_to_taxid.get(protein_accession, None)

        if tax_id:
            if key not in added_taxids_per_key:
                added_taxids_per_key[key] = set()

            if tax_id not in added_taxids_per_key[key]:
                added_taxids_per_key[key].add(tax_id)
                if key not in data:
                    data[key] = []
                data[key].append(protein_accession)  # Modified this line

# Write to separate output files
for key in data:
    output_file = f"{key}_orthologs_vertebrates_n_primates_proteinids.txt"
    with open(output_file, 'w') as f:
        f.write(key + "\n")
        f.write("\n".join(data[key]))

print("Done!")

