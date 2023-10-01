import sys
import os
import gzip
from Bio import SeqIO

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

    # calculate exon_positions
    exon_aa_positions = []
    cumulative_sum = 0
    for length in protein_coding_len:
        start = cumulative_sum + 1
        end = length
        exon_aa_positions.append([start, end])
        cumulative_sum = end
    print(exon_aa_positions)
    # Transforming the list
    exon_intron_exon_pos = [(exon_aa_positions[i][0], exon_aa_positions[i][1], exon_aa_positions[i+1][1]) for i in range(len(exon_aa_positions)-1)]
    return intron_pos, exon_intron_exon_pos, protein_intron_len


def insert_introns_into_sequence(sequence, intron_positions):
    sequence_with_introns = list(sequence)
    for pos in intron_positions:
        sequence_with_introns.insert(pos, 'X')
    return ''.join(sequence_with_introns)

def get_sequence_from_fasta_db(fasta_file_path, protein_id, idx_filename="fasta.idx"):
    if not os.path.exists(idx_filename):
        _ = SeqIO.index_db(idx_filename, fasta_file_path, "fasta")

    fasta_index = SeqIO.index_db(idx_filename)
    if protein_id in fasta_index:
        record = fasta_index[protein_id]
        return record.description, str(record.seq)
    return None, None


if __name__ == '__main__':
    if len(sys.argv) != 5:
        print('Usage: python insert_introns.py <table.txt> <amino_acid_fasta.txt> <protein_id> <first_protein_id>')
        sys.exit(1)

    table_path = sys.argv[1]
    aa_fasta_path = sys.argv[2]
    protein_id = sys.argv[3]
    first_protein_id = sys.argv[4]

    try:
        aa_intron_positions, exon_intron_exon_pos, protein_intron_len = new_extract_aa_intron_positions(table_path)

        fasta_header, amino_acid_sequence = get_sequence_from_fasta_db(aa_fasta_path, protein_id)
        if not fasta_header or not amino_acid_sequence:
            raise ValueError(f"Protein ID {protein_id} not found in FASTA file.")

        sequence_with_introns = insert_introns_into_sequence(amino_acid_sequence, aa_intron_positions)

        print("amino_acid_sequence", amino_acid_sequence)
        print("sequence_with_introns", sequence_with_introns)

        with open(first_protein_id + "_sequence_w_introns.fa", 'a') as f:
            f.write(f">{fasta_header}\n")
            f.write(f"{sequence_with_introns}\n")

        with open(first_protein_id + "_intron_positions.txt", 'a') as f:
            f.write(f">{fasta_header}\n")
            f.write(', '.join(map(str, aa_intron_positions)) + '\n')

    except Exception as e:
        print(f"An error occurred: {e}")
