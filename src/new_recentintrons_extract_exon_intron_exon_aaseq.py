import sys

def extract_aa_intron_positions_and_exons(table_path):
    with open(table_path, 'r') as file:
        lines = file.readlines()

    gene_interval_coding = []
    coding_lengths = []
    intron_lengths = []
    for line in lines[2:]:
        split_line = line.split()
        if len(split_line) == 7: # ensure there are at least 6 columns
            coding_lengths.append(split_line[-2])
            gene_interval_coding.append(split_line[-4].split('-'))
            intron_lengths.append(int(split_line[-1])//3)
        if len(split_line) == 6:
            ###print(split_line[-1])
            coding_lengths.append(split_line[-2])
            gene_interval_coding.append(split_line[-4].split('-'))
            #coding_lengths = [list(map(int, line.split()[6].split('-')))]
    ###print("coding lengths", coding_lengths)
    ###print("intron lengths", intron_lengths)

    # calculate exon_positions
    exon_positions = []
    cumulative_sum = 0
    for length in coding_lengths:
        # Since the starting position is cumulative_sum + 1
        aa_length = int(length)//3
        start = cumulative_sum 
        # Ending position is the cumulative_sum + current aa_length
        end = cumulative_sum + aa_length 
        exon_positions.append([start, end])
        # Update cumulative_sum
        cumulative_sum = end
        
    #gene_interval_coding = [list(map(int, line.split()[3].split('-'))) for line in lines[2:]]
    #print(gene_interval_coding)

    # # Convert genomic intervals into exon positions
    interval_coding_aa_positions = [(int(x[0])//3, int(x[1])//3) for x in gene_interval_coding]

    # # Calculate intron positions from exon positions
    intron_positions = []
    for i in range(len(interval_coding_aa_positions)-1):
        intron_start = interval_coding_aa_positions[i][1] + 1
        intron_end = interval_coding_aa_positions[i+1][0] - 1
        intron_positions.append((intron_start, intron_end))

    ###print("exon_positions", exon_positions)
    print("interval_coding_aa_positions", interval_coding_aa_positions)
    ###print("intron postions", intron_positions)
    #print("intron_lengths", intron_lengths)

    return exon_positions, intron_positions, intron_lengths



def extract_exon_sequences(sequence, exon_positions, intron_num, intron_lengths):
    # Subtract one because we are using 0-based indexing
    intron_num = int(intron_num)
    flanking_exon_positions = [exon_positions[intron_num-1], exon_positions[intron_num]]
    flanking_exon_sequences = [sequence[pos[0]:pos[1]] for pos in flanking_exon_positions]
    intron_sequence = "X" * int(intron_lengths[intron_num-1])
    print(int(intron_lengths[intron_num-1]))
    ###print("len seq", len(sequence))
    ###print(flanking_exon_positions)
    ###print(flanking_exon_sequences)

    new_sequence = flanking_exon_sequences[0] + intron_sequence + flanking_exon_sequences[1]
    print("flanking_exon_sequences[0]", flanking_exon_sequences[0])
    print("flanking_exon_sequences[1]", flanking_exon_sequences[1])
    ###print(new_sequence)

    return new_sequence, flanking_exon_sequences

def extract_amino_acid_sequence(file_content):
    lines = file_content.split('\n')
    sequence_lines = [line for line in lines if not line.startswith(">")]
    return "".join(sequence_lines)

if __name__ == '__main__':
    if len(sys.argv) != 5:
        print('Usage: python recentintrons_extract_exon_intron_exon_aaseq.py <table.txt> <protein_id> <intron_num> <amino_acid_seq.fa>')
        sys.exit(1)

    table_path = sys.argv[1]
    protein_id = sys.argv[2]
    intron_num = int(sys.argv[3])
    aa_seq_path = sys.argv[4]

    try:
        exon_positions, intron_positions, intron_lengths = extract_aa_intron_positions_and_exons(table_path)

        with open(aa_seq_path, 'r') as f:
            file_content = f.read()
        print("intron_positions", intron_positions)

        amino_acid_sequence = extract_amino_acid_sequence(file_content)
        new_sequences, flanking_exon_sequences = extract_exon_sequences(amino_acid_sequence, exon_positions, intron_num, intron_lengths)

        with open(protein_id + '_to_be_aligned_for_stevem.fa', 'w') as f:
            f.write(f">{protein_id},{intron_num}\n")
            f.write(f"{new_sequences}\n")

    except Exception as e:
        print(f"An error occurred: {e}")
