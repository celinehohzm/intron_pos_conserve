def extract_descriptions(protein_ids_file, fasta_file, output_file):
    # Read protein IDs into a set for faster lookup
    with open(protein_ids_file, 'r') as f:
        protein_ids = {line.strip() for line in f}

    # Open the output file
    with open(output_file, 'w') as outfile:
        # Go through the FASTA file
        with open(fasta_file, 'r') as fasta:
            for line in fasta:
                # Check if the line starts with '>'
                if line.startswith('>'):
                    # Extract the protein ID (without the '>')
                    protein_id = line[1:].split()[0]
                    if protein_id in protein_ids:
                        # Find the index of the last '[' character, to extract the description
                        index_of_last_bracket = line.rfind('[')
                        # Extract the description part
                        description = line[1:index_of_last_bracket].split(protein_id)[1].strip()
                        # Write the protein ID and description to the output file
                        outfile.write(f"{protein_id}\t{description}\n")

# Replace 'path_to_protein_ids.txt' and 'path_to_fasta.faa' with your actual file paths
extract_descriptions('filtered_ALL_recent_introns_uniq_proteinids.txt', 
                     'MANE.GRCh38.v1.0.refseq_protein.faa', 
                     'ALL_recent_protein_descriptions.txt')

