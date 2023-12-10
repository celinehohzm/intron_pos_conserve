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

    total_protein_len = protein_coding_len[-1]
    print("total_protein_len", total_protein_len)
    intron_pos_score = [ round((x / total_protein_len) * 100, 2) for x in intron_pos]
    print("intron_pos_score", intron_pos_score)

    return intron_pos, intron_lengths, intron_pos_score

def insert_introns_into_sequence(sequence, intron_positions):
    sequence_with_introns = list(sequence)

    for pos in intron_positions:
        sequence_with_introns.insert(pos, 'X')

    return ''.join(sequence_with_introns)

def get_sequence_from_fasta(fasta_content, protein_id):
    # Splitting the content into lines
    lines = fasta_content.strip().split('\n')
    # Iterating through the lines
    for i, line in enumerate(lines):
        if line.startswith(">") and protein_id in line:
            # Return the next line if the protein_id is found
            return line, lines[i + 1]
    return None


if __name__ == '__main__':
    if len(sys.argv) != 5:
        print('Usage: python insert_introns.py <table.txt> <amino_acid_fasta.txt> <protein_id> <first_protein_id>')
        sys.exit(1)

    table_path = sys.argv[1]
    aa_fasta_path = sys.argv[2]
    protein_id = sys.argv[3]
    first_protein_id = sys.argv[4]

    try:
        aa_intron_positions, intron_lengths, intron_pos_score = new_extract_aa_intron_positions(table_path)

        with open(aa_fasta_path, 'r') as f:
            file_content = f.read()

        fasta_header, amino_acid_sequence = get_sequence_from_fasta(file_content, protein_id)

        sequence_with_introns = insert_introns_into_sequence(amino_acid_sequence, aa_intron_positions)


        print("amino_acid_sequence", amino_acid_sequence)
        print("sequence_with_introns", sequence_with_introns)


        with open(first_protein_id + "_sequence_w_introns.fa", 'a') as f:
            f.write(f"{fasta_header}\n")
            f.write(f"{sequence_with_introns}\n")

        with open(first_protein_id + "_intron_positions.txt", 'a') as f:
            f.write(f"{fasta_header}\n")
            f.write(', '.join(map(str, aa_intron_positions)) + '\n')

        with open(first_protein_id + "_intron_lengths.txt", 'a') as f:
            f.write(f"{fasta_header}\n")
            f.write(', '.join(map(str, intron_lengths)) + '\n')

        with open(first_protein_id + "_intron_position_score.txt", 'a') as f:
            f.write(f"{fasta_header}\n")
            f.write(', '.join(map(str, intron_pos_score)) + '\n')

    except Exception as e:
        print(f"An error occurred: {e}")
