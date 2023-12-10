import sys
from Bio import SeqIO

def extract_sequences(input_file, output_file, start_position, end_position):
    with open(output_file, 'w') as outfile:
        for record in SeqIO.parse(input_file, "fasta"):
            # Adjusting for zero-based indexing in Python
            extracted_seq = record.seq[start_position-1:end_position]
            outfile.write(f">{record.id}\n{extracted_seq}\n")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python script.py <input_fasta> <output_fasta> <start_pos> <end_pos>")
        sys.exit(1)

    input_fasta = sys.argv[1]
    output_fasta = sys.argv[2]
    start_pos = int(sys.argv[3])
    end_pos = int(sys.argv[4])

    extract_sequences(input_fasta, output_fasta, start_pos, end_pos)

