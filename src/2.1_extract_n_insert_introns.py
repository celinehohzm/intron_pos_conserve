import sys

def new_extract_aa_intron_positions(table_path):

    with open(table_path, 'r') as file:
        lines = file.readlines()

    coding_lengths = []
    intron_lengths = []
    for line in lines[2:]:
        split_line = line.split()
        if len(split_line) == 7: # ensure there are at least 6 columns
            coding_lengths.append(int(split_line[-2]))
            intron_lengths.append(int(split_line[-1]))
        if len(split_line) == 6:
            print(split_line[-1])
            coding_lengths.append(int(split_line[-1]))

    print("coding lengths", coding_lengths)
    print("intron_lengths", intron_lengths)

    protein_coding_len = [round(x / 3) for x in coding_lengths]
    print("protein_coding_len", protein_coding_len)

    protein_intron_len = [round(x / 3) for x in intron_lengths]
    print("protein_intron_len", protein_intron_len)

    for i in range(1, len(protein_coding_len)):
        protein_coding_len[i] += protein_coding_len[i - 1]
    print("protein_coding_len_sum", protein_coding_len)

    intron_pos = [x + 1 for x in protein_coding_len[:-1]]
    print("intron_positions", intron_pos)

    return intron_pos


def insert_introns_into_sequence(sequence, intron_positions):
    sequence_with_introns = list(sequence)

    for pos in intron_positions:
        sequence_with_introns.insert(pos, 'X')

    return ''.join(sequence_with_introns)

def extract_amino_acid_sequence(file_content):
    lines = file_content.split('\n')
    sequence_lines = [line for line in lines if not line.startswith(">")]
    return "".join(sequence_lines)

def get_fasta_header(file_content):
    lines = file_content.split('\n')
    header_line = [line for line in lines if line.startswith(">")]
    return header_line[0] if header_line else ""

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print('Usage: python insert_introns.py <table.txt> <amino_acid_fasta.txt> <first_protein_id>')
        sys.exit(1)

    table_path = sys.argv[1]
    aa_fasta_path = sys.argv[2]
    first_protein_id = sys.argv[3]

    try:
        aa_intron_positions = new_extract_aa_intron_positions(table_path)

        with open(aa_fasta_path, 'r') as f:
            file_content = f.read()

        fasta_header = get_fasta_header(file_content)
        amino_acid_sequence = extract_amino_acid_sequence(file_content)

        sequence_with_introns = insert_introns_into_sequence(amino_acid_sequence, aa_intron_positions)

        with open(first_protein_id + "_sequence_w_introns.fa", 'a') as f:
            f.write(f"{fasta_header}\n")
            f.write(f"{sequence_with_introns}\n")

        with open(first_protein_id + "_intron_positions.txt", 'a') as f:
            f.write(f"{fasta_header}\n")
            f.write(', '.join(map(str, aa_intron_positions)) + '\n')


    except Exception as e:
        print(f"An error occurred: {e}")
