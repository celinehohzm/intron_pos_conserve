from Bio import SeqIO

def split_fasta(input_file, output_pattern, split_count):
    records = list(SeqIO.parse(input_file, 'fasta'))
    records_per_file = len(records) // split_count
    for i in range(split_count):
        start = i * records_per_file
        if i == split_count - 1:
            # Include the remainder sequences in the last file
            end = len(records)
        else:
            end = start + records_per_file
        SeqIO.write(records[start:end], output_pattern.format(i+1), 'fasta')

input_fasta = "ALL_recent_introns_sequences.fasta"
output_pattern = "ALL_recent_introns_sequences_split_{}.fasta"  # {} will be replaced with file number
split_count = 3

split_fasta(input_fasta, output_pattern, split_count)

