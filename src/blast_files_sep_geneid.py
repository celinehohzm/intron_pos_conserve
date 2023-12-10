#!/usr/bin/env python3
import sys

if len(sys.argv) != 3:
    print("Usage: python script_name.py <MANE_refseq_protein_file> <gene2refseq_file>")
    sys.exit(1)

input_file = sys.argv[1]
gene2refseq_file = sys.argv[2]

# Read the gene2refseq file and map protein accessions to GeneIDs
protein_to_geneid = {}
with open(gene2refseq_file, 'r') as f:
    # Skipping header line
    next(f)
    for line in f:
        parts = line.strip().split("\t")
        gene_id = parts[1]
        protein_accession = parts[5].split(".")[0]  # Get only the base accession without the version
        protein_to_geneid[protein_accession] = gene_id

# Read the input file and process data
with open(input_file, 'r') as f:
    data = {}
    added_geneids_per_key = {}  # Dictionary to keep track of geneids that have been processed for each key
    key_to_original_protein = {}  # Dictionary to map gene_ids back to their original protein_ids
    
    for line in f:
        line = line.strip()
        key, value = line.split("\t")
        # Extract the accession number from the value
        protein_accession_value = value.split("|")[1].split(".")[0]  # Get only the base accession without the version
        gene_id_value = protein_to_geneid.get(protein_accession_value, None)
        key_gene_id = protein_to_geneid.get(key.split(".")[0], key)  # Convert key to its gene_id
        
        # Store the original protein ID for each gene_id
        key_to_original_protein[key_gene_id] = key
        
        if gene_id_value:
            if key_gene_id not in added_geneids_per_key:
                added_geneids_per_key[key_gene_id] = set()

            if gene_id_value not in added_geneids_per_key[key_gene_id]:
                added_geneids_per_key[key_gene_id].add(gene_id_value)
                if key_gene_id not in data:
                    data[key_gene_id] = []
                data[key_gene_id].append(gene_id_value)  # Add gene_id instead of protein_accession

# Write to separate output files
for key_gene_id, values in data.items():
    # Use the key's original protein ID for the file name
    original_key = key_to_original_protein[key_gene_id]
    output_file = f"blast_files/{original_key}_orthologs_vertebrates_n_primates_geneids.txt"
    with open(output_file, 'w') as f:
        f.write(key_gene_id + "\n")
        f.write("\n".join(values))

print("Done!")

