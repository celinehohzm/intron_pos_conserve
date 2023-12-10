def extract_sequences(input_file, sequence_file, output_file):
    # Read the headers from the input_file
    with open(input_file, 'r') as f:
        headers = set(line.split()[0] for line in f)

    # Extract sequences corresponding to the headers from sequence_file
    with open(sequence_file, 'r') as f:
        sequences = f.read().split('>')

    # Write the selected sequences to output_file
    with open(output_file, 'w') as f:
        for seq in sequences:
            lines = seq.strip().split('\n')
            header = lines[0]
            sequence = '\n'.join(lines[1:])
            if header in headers:
                f.write(f'>{header}\n{sequence}\n')


if __name__ == '__main__':
    input_file = 'recentintrons_exons_of_intron_orthologs_POTENTIAL_INTRONIZATION.txt'
    sequence_file = 'recentintrons_exons_of_intron.fa'
    output_file = 'output_sequences.fa'
    extract_sequences(input_file, sequence_file, output_file)

