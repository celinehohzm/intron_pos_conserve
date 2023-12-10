# Usage: python extract_proteins.py last_round_proteinids.txt MANE.GRCh38.v1.0.refseq_protein.faa output_sequences.faa

import sys

def read_protein_ids(protein_id_file):
    """Read protein IDs from a file and return them as a set."""
    with open(protein_id_file, 'r') as f:
        return set(line.strip() for line in f)

def extract_sequences(protein_ids, fasta_file, output_file):
    """Extract sequences with IDs in protein_ids from fasta_file and write to output_file."""
    with open(fasta_file, 'r') as fasta, open(output_file, 'w') as output:
        write_sequence = False
        for line in fasta:
            if line.startswith('>'):
                protein_id = line.split()[0][1:]  # Get the ID part, remove '>'
                write_sequence = protein_id in protein_ids
            if write_sequence:
                output.write(line)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python extract_proteins.py protein_ids_file fasta_file output_file")
        sys.exit(1)

    protein_ids_file, fasta_file, output_file = sys.argv[1], sys.argv[2], sys.argv[3]
    protein_ids = read_protein_ids(protein_ids_file)
    extract_sequences(protein_ids, fasta_file, output_file)
    print(f"Protein sequences extracted to {output_file}")

