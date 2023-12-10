def split_fasta(filename, sequences_per_file):
    with open(filename, 'r') as f:
        count = 0
        file_num = 1
        output_file = None

        for line in f:
            # If the line is a header and count is a multiple of sequences_per_file, start a new file
            if line.startswith('>') and count % sequences_per_file == 0:
                if output_file:
                    output_file.close()
                output_filename = f"MANE_refseq_protein_chunk{file_num}.fasta"
                output_file = open(output_filename, 'w')
                file_num += 1
            output_file.write(line)
            if line.startswith('>'):
                count += 1
        if output_file:
            output_file.close()

filename = "MANE.GRCh38.v1.0.refseq_protein.faa"
sequences_per_file = 956
split_fasta(filename, sequences_per_file)
