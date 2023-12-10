import os
import sys
import subprocess

def read_blast_output(filename):
    """Reads a BLAST output file and returns a dictionary with the second column as keys."""
    results = {}
    with open(filename, 'r') as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) > 10:
                key = parts[1]  # Second column
                results[key] = parts
    return results

def fetch_sequence(protein_id, start, stop):
    """Fetches a protein sequence using efetch."""
    result = subprocess.run(
        ['efetch', '-db', 'protein', '-id', protein_id, '-format', 'fasta', '-seq_start', str(start), '-seq_stop', str(stop)],
        capture_output=True, text=True
    )
    return result.stdout

def compare_files(file1, file2, protein_id, intron_len):
    """Compares two BLAST output files and prints the desired output."""
    file1_data = read_blast_output(file1)
    file2_data = read_blast_output(file2)

    header = "Ortho_protein\tseq2_start\tseq1_end\tseq1_2_gap\tseq1_start\tseq2_end\tseq1\tseq2\n"

    for key in file1_data:
        if key in file2_data:
            filename = protein_id + '_blasts_outputs.txt'
            # Check if the file already exists
            file_exists = os.path.isfile(filename)

            with open(filename, 'a') as f:
                # Write the header if the file does not exist
                if not file_exists:
                    f.write(header)

                # Write the data
                f.write(f"{key}\t{file2_data[key][8]}\t{file1_data[key][9]}\t{int(file2_data[key][8]) - int(file1_data[key][9])}\t{int(file1_data[key][8])}\t{int(file2_data[key][9])}\t{file1_data[key][-1]}\t{file2_data[key][-1]}\n")
                
                # Fetch and print the sequence
                if int(file2_data[key][8]) - int(file1_data[key][9]) > 0.03*intron_len:
                    sequence = fetch_sequence(key, int(file1_data[key][8]), int(file2_data[key][9]))
                    with open(protein_id + '_exons_of_intron.fa', 'a') as f:
                        f.write(f"{sequence}\n")
                    print(sequence)

# Check if the protein ID is provided as a command-line argument
if len(sys.argv) > 2:
    protein_id = sys.argv[1]
    intron_len = int(sys.argv[2])
    compare_files(protein_id + '_blast1.txt', protein_id + '_blast2.txt', protein_id, intron_len)
else:
    print("Usage: script.py <protein_id>")

