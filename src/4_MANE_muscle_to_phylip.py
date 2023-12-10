import sys
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import re
from collections import defaultdict

input_file = sys.argv[1]
output_file = sys.argv[2]


def generate_unique_id(species_name, species_count):
    words = species_name.split()
    if len(words) >= 2:
        unique_id = words[0][:min(3, len(words[0]))] + words[1][:min(6, len(words[1]))].capitalize()
    else:
        unique_id = species_name

    species_count[unique_id] += 1
    if species_count[unique_id] > 1:
        unique_id += f"{species_count[unique_id]}"

    return unique_id


# Read the input FASTA file
alignment = AlignIO.read(input_file, "fasta")

# Create a new alignment with unique identifiers
unique_alignment = MultipleSeqAlignment([])

# Create a defaultdict to count occurrences of each species name
species_count = defaultdict(int)

for i, record in enumerate(alignment):
    # Extract the species name from the record description using a regular expression
    species_match = re.search(r'\[(.*?)\]', record.description)

    if species_match:
        species_name = species_match.group(1)
    else:
        species_name = record.id

    # Use the generate_unique_id() function to generate the unique identifier
    unique_id = generate_unique_id(species_name, species_count)

    record.id = unique_id
    unique_alignment.append(record)

# Write the alignment with unique identifiers to the output PHYLIP file
AlignIO.write(unique_alignment, output_file, "phylip")

