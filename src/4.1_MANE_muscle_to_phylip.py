import sys
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import re
from collections import defaultdict

input_file = sys.argv[1]
output_file = sys.argv[2]
id_mapping_file = sys.argv[3]  # Add an additional command-line argument for the ID mapping file

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

# Create a dictionary to store the mapping between shortened and original identifiers
id_mapping = {}

for i, record in enumerate(alignment):
    # Extract the species name from the record description using a regular expression
    species_match = re.search(r'\[(.*?)\]', record.description)

    if species_match:
        species_name = species_match.group(1)
    else:
        species_name = record.id

    # Use the generate_unique_id() function to generate the unique identifier
    unique_id = generate_unique_id(species_name, species_count)

    # Create a shortened identifier
    short_id = f"sp{i}"

    # Store the mapping between the shortened and original identifiers
    id_mapping[short_id] = unique_id

    # Use the shortened identifier in the output alignment
    record.id = short_id
    unique_alignment.append(record)

# Write the alignment with shortened identifiers to the output PHYLIP file
AlignIO.write(unique_alignment, output_file, "phylip")

# Write the ID mapping to a separate file
with open(id_mapping_file, "w") as f:
    for short_id, original_id in id_mapping.items():
        f.write(f"{short_id}\t{original_id}\n")

