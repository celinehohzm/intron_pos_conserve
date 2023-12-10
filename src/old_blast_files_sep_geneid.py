#!/usr/bin/env python3

# Define file paths
input_file = "blast_files/MANE_refseq_protein_chunk18_orthologs.txt"
gene2refseq_file = "gene2refseq"

# Create a mapping from protein accession numbers to GeneID from the gene2refseq file
protein_to_gene = {}
with open(gene2refseq_file, 'r') as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) > 6:  # Ensure the line is long enough to avoid out of range errors
            gene_id = parts[1]
            protein_accession = parts[5]
            protein_to_gene[protein_accession] = gene_id

# Read the input file and process data
data = {}
with open(input_file, 'r') as f:
    for line in f:
        line = line.strip()
        protein_key, value = line.split("\t")
        
        # Extract the gene corresponding to the value protein
        gene_key = protein_to_gene.get(protein_key, protein_key)
        gene_value = protein_to_gene.get(value.split("|")[1], value)

        if protein_key in data:
            data[protein_key].append(gene_value)


# Write to separate output files
for protein_key in data:
    output_file = f"test_blast_files/old_{protein_key}_orthologs_geneids.txt"
    with open(output_file, 'w') as f:
        f.write("\n".join(data[protein_key]))

print("Done!")
