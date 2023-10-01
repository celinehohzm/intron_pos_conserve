import sys

def parse_fasta(fasta_file):
    with open(fasta_file, 'r') as f:
        contents = f.read().splitlines()

    fasta_dict = {}
    current_header = ""

    for line in contents:
        if line.startswith('>'):
            current_header = line
            fasta_dict[current_header] = ""
        else:
            fasta_dict[current_header] += line

    return fasta_dict

def extract_sequence_by_id(protein_id, fasta_dict):
    for header, sequence in fasta_dict.items():
        if protein_id in header:
            return header, sequence
    return None, None

def write_output_file(protein_id, header, sequence):
    if header and sequence:
        with open(f"{protein_id}_.fasta", 'w') as output_file:
            output_file.write(f"{header}\n{sequence}\n")
            print(f"Sequence for protein ID {protein_id} has been written to {protein_id}.fasta")
    else:
        print(f"Protein ID {protein_id} not found in the input file.")

def main():
    if len(sys.argv) != 3:
        print("Usage: python script.py <protein_id>")
        sys.exit(1)

    protein_id = sys.argv[1]
    input_file = sys.argv[2]

    fasta_dict = parse_fasta(input_file)
    header, sequence = extract_sequence_by_id(protein_id, fasta_dict)
    write_output_file(protein_id, header, sequence)

if __name__ == "__main__":
    main()
