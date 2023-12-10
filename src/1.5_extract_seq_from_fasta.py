import sys

def determine_strand_orientation(genomic_intervals):
    strand_orientation = "positive" if genomic_intervals[0][0] < genomic_intervals[1][0] else "negative"
    return strand_orientation

def read_fasta(file):
    """Generator function to read a fasta file and yield header and sequence"""
    with open(file, 'r') as f:
        header = None
        sequence = ''
        for line in f:
            if line.startswith('>'):
                if header is not None:
                    yield header, sequence
                    sequence = ''
                header = line.rstrip()
            else:
                sequence += line.rstrip()
        if header is not None:
            yield header, sequence

def extract_sequence(file_content):
    lines = file_content.split('\n')
    sequence_lines = [line for line in lines if not line.startswith(">")]
    return "".join(sequence_lines)

def get_fasta_header(file_content):
    lines = file_content.split('\n')
    header_line = [line for line in lines if line.startswith(">")]
    return header_line[0] if header_line else ""

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage: python insert_introns.py <table.txt> <amino_acid_fasta.txt>')
        sys.exit(1)

    # dna_fasta_path = sys.argv[1]
    aa_fasta_path = sys.argv[1]

    # Read in the fasta file and extract the sequences with "whole genome shotgun sequence" in the header
    # for header, sequence in read_fasta(dna_fasta_path):
    #     dna_header= header
    #     dna_sequence = sequence
    #     break
    #
    # with open("org_dna_sequence.fasta", "w") as output_file:
    #     output_file.write(f"{dna_header}\n{dna_sequence}\n")

    with open(aa_fasta_path, 'r') as f:
        aa_content = f.read()
    aa_header = get_fasta_header(aa_content)
    aa_sequence = extract_sequence(aa_content)
    print(f"{aa_header}\n{aa_sequence}")

    #with open("org_aa_sequence.fa", "w") as output_file:
    #    output_file.write(f"{aa_header}\n{aa_sequence}\n")
