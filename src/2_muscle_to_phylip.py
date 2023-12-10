from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import sys
import re

input_file = sys.argv[1]
output_file = sys.argv[2]

# Read the input FASTA file
alignment = AlignIO.read(input_file, "fasta")

# Create a new alignment with unique identifiers
unique_alignment = MultipleSeqAlignment([])

for i, record in enumerate(alignment):
    # Extract the species name from the record description using a regular expression
    species_match = re.search(r'\[(.*?)\]', record.description)
    
    # Use the species name as the identifier or use the original ID if the species name is not found
    if species_match:
        unique_id = species_match.group(1).replace(" ", "_")
    else:
        unique_id = record.id
    
    # Append a unique index to the identifier to avoid duplicates
    unique_id = f"{unique_id}_{i}"
    record.id = unique_id
    unique_alignment.append(record)

# Write the alignment with unique identifiers to the output PHYLIP file
AlignIO.write(unique_alignment, output_file, "phylip")

