import os

# Define the root directory and filename
root_dir = "/home/choh1/intron_pos_conserv/blastdb"  # Modify this to your root directory
protein_ids_file = "filtered_ALL_recent_introns_uniq_proteinids.txt"
output_file = "ALL_recent_introns_sequences.fasta"

def extract_fasta_from_file(file_path):
    """Extract fasta sequence from a given file."""
    with open(file_path, 'r') as f:
        return f.read()

def main():
    # Read the protein IDs from the file
    with open(protein_ids_file, 'r') as f:
        protein_ids = [line.strip() for line in f]

    # Extract and combine sequences
    combined_sequences = ""
    for protein_id in protein_ids:
        fasta_file_path = os.path.join(root_dir, protein_id, f"{protein_id}_.fasta")
        if os.path.isfile(fasta_file_path):
            combined_sequences += extract_fasta_from_file(fasta_file_path)
            combined_sequences += "\n"  # Separate sequences with a newline

    # Write combined sequences to the output file
    with open(output_file, 'w') as f:
        f.write(combined_sequences)

if __name__ == "__main__":
    main()

