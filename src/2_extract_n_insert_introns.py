import sys

def new_extract_aa_intron_positions(table_path):

    with open(table_path, 'r') as file:
        lines = file.readlines()

    coding_lengths = []
    for line in lines[2:]:
        split_line = line.split()
        print(split_line[-2])
        # if len(split_line) >= 6: # ensure there are at least 6 columns
        #if split_line[5].strip() != "": # check if the 6th column (coding lengths) is empty
        coding_lengths.append(int(split_line[-2]))
        #else:
            #raise ValueError("The 'coding lengths' column is missing a value.")
        # else:
        #     raise ValueError("There are not enough columns in the table.")

    print("table coding len", coding_lengths)

    #genomic_intervals = [list(map(int, line.split()[0].split('-'))) for line in lines[2:]]
    #print("genomic interval", genomic_intervals)

    #reversed_order = False
    #if genomic_intervals[0][1] < genomic_intervals[0][0]:
       # reversed_order = True
    #print(reversed_order)

    #if len(genomic_intervals) < 1:
       # raise ValueError("There should be at least two genomic intervals in the input table file.")

    protein_coding_len = [round(x / 3) for x in coding_lengths]
    print("protein_coding_len", protein_coding_len)

    #if reversed_order == True:
     #   protein_coding_len = protein_coding_len[::-1]
    #print("protein_coding_len_after", protein_coding_len)

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
        aa_intron_positions = new_extract_aa_intron_positions(table_path)

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


    except Exception as e:
        print(f"An error occurred: {e}")
