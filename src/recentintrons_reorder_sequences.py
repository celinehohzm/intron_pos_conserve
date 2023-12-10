import re
import sys

def parse_introns_file(introns_file, intron_num):
    with open(introns_file, 'r') as file:
        introns_data = file.read()

    pattern = fr'intron_{intron_num} - X aligned count: \d+ - X not aligned count: \d+\n[\s\S]*?Species with X: ([\s\S]*?)\nSpecies without X: ([\s\S]*?)(?=\n\n|\Z)'
    match = re.search(pattern, introns_data)

    if match:
        species_with_x_list = match.group(1).split(', ')
        species_without_x_list = match.group(2).split(', ')
    else:
        raise ValueError(f"Intron number {intron_num} not found in the file.")

    return species_with_x_list, species_without_x_list

def parse_fasta_file(fasta_file, species_with_x, species_without_x):
    with open(fasta_file, 'r') as file:
        fasta_data = file.read().split('>')[1:]  # Split and remove empty first element

    sequences = {}
    for entry in fasta_data:
        lines = entry.split('\n')
        header = lines[0]
        sequence = ''.join(lines[1:])
        species_name = re.search(r'\[([^\]]+)\]', header).group(1)

        # Insert Y or N right after the protein ID
        protein_id = re.search(r'(\w+\.\d+)', header).group(1)
        suffix = '_Y' if species_name in species_with_x else '_N' if species_name in species_without_x else ''
        modified_header = header.replace(protein_id, protein_id + suffix)

        sequences[species_name] = '>' + modified_header + '\n' + sequence

    return sequences

def reorder_sequences(species_order, sequences):
    ordered_sequences = []
    for species in species_order:
        if species in sequences:
            ordered_sequences.append(sequences[species])

    # Ensure Homo sapiens is first if present
    if 'Homo sapiens' in sequences:
        ordered_sequences.insert(0, sequences['Homo sapiens'])

    return ordered_sequences

def main():
    if len(sys.argv) != 3:
        print("Usage: script.py <protein_id> <intron_num>")
        sys.exit(1)

    protein_id = sys.argv[1]
    intron_num = sys.argv[2]

    introns_file = f'{protein_id}/{protein_id}_recent_introns.txt'
    fasta_file = f'{protein_id}/{protein_id}_sequence_w_introns.fa'
    #fasta_file = f'{protein_id}/{protein_id}_reorder_aln_seq_w_introns.afa'
    output_file = f'{protein_id}/{protein_id}_{intron_num}_reordered_sequences.fa'

    species_with_x, species_without_x = parse_introns_file(introns_file, intron_num)
    sequences = parse_fasta_file(fasta_file, species_with_x, species_without_x)
    ordered_sequences = reorder_sequences(species_with_x + species_without_x, sequences)

    with open(output_file, 'w') as file:
        file.write('\n'.join(ordered_sequences))

if __name__ == "__main__":
    main()

