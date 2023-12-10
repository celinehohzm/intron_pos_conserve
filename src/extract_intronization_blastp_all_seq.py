import os
import sys
import re
import subprocess
from Bio import SeqIO

# Check if the correct number of arguments is provided
if len(sys.argv) != 3:
    print("Usage: python extract_sequences.py input_file.txt recentintrons_exons_of_intron.fa")
    sys.exit(1)

input_file = sys.argv[1]
sequence_file = sys.argv[2]

# Function to fetch protein sequences using entrezpy
def fetch_protein_sequence(protein_id):
    try:
        result = subprocess.run(['/home/choh1/anaconda3/envs/gtf/lib/python3.10/site-packages/entrezpy', 'efetch', '-db', 'protein', '-id', protein_id, '-format', 'fasta'],
                                capture_output=True, text=True, check=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error fetching sequence for {protein_id}: {e}")
        return None

def main():
    current_header = ""
    current_proteins = []

    with open(input_file, 'r') as infile:
        for line in infile:
            line = line.strip()
            fields = line.split()
            header = fields[0]
            protein_id = fields[1].split('|')[1]

            if header != current_header:
                # New header, process the previous one and reset the current proteins list
                if current_header:
                    # Fetch and save sequences of the current proteins
                    output_file = f"recentintrons_{current_header}_orthologs_aa_seq.fa"
                    with open(output_file, 'a') as outfile:
                        for pid in current_proteins:
                            seq = fetch_protein_sequence(pid)
                            if seq:
                                outfile.write(seq)
                current_header = header
                current_proteins = [protein_id]
            else:
                # Same header, add the protein id to the list
                current_proteins.append(protein_id)

    # Process the last header and its proteins
    if current_header:
        output_file = f"recentintrons_{current_header}_orthologs_aa_seq.fa"
        with open(output_file, 'a') as outfile:
            for pid in current_proteins:
                seq = fetch_protein_sequence(pid)
                if seq:
                    outfile.write(seq)

    # Extract the corresponding sequences from the given sequence file
    output_directory = os.path.dirname(sequence_file)
    with open(input_file, 'r') as infile:
        for line in infile:
            line = line.strip()
            fields = line.split()
            header = fields[0]
            seq_id = header.split(',')[0]
            seq_id = re.sub(r'[^\w\s]', '', seq_id)  # Remove non-alphanumeric characters
            sequence = None

            # Find the corresponding sequence in the sequence file
            with open(sequence_file, 'r') as seqfile:
                for record in SeqIO.parse(seqfile, 'fasta'):
                    if record.id.startswith(seq_id):
                        sequence = record.seq
                        break

            # Save the sequence to the appropriate file
            if sequence:
                output_file = os.path.join(output_directory, f"recentintrons_{header}_orthologs_aa_seq.fa")
                with open(output_file, 'w') as outfile:
                    outfile.write(f">{record.description}\n{sequence}\n")

if __name__ == "__main__":
    main()

