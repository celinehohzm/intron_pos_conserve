
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

# Read the input FASTA file
alignment = AlignIO.read(input_file, "fasta")

# Create a new alignment with unique identifiers
unique_alignment = MultipleSeqAlignment([])

for i, record in enumerate(alignment):
    # Truncate the identifier to 7 characters and append a unique index
    unique_id = record.id[:5] 
    record.id = unique_id
    unique_alignment.append(record)

# Write the alignment with unique identifiers to the output PHYLIP file
AlignIO.write(unique_alignment, output_file, "phylip")

