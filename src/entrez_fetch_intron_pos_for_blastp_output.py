from Bio import Entrez
import re

Entrez.email = "celinehohzm@gmail.com"  # Always tell NCBI who you are

def fetch_intron_positions(protein_id):
    handle = Entrez.efetch(db="protein", id=protein_id, rettype="gb", retmode="text")
    record = handle.read()
    intron_positions = []

    for match in re.finditer("intron", record):
        start_pos = record.rfind("     ", 0, match.start())
        position_line = record[start_pos:match.start()].strip()
        
        # Extracting the start and end positions
        start, end = map(int, position_line.split(".."))
        intron_positions.append((start, end))
    
    return intron_positions

# Get the unique protein IDs
protein_ids = set([line.split()[1].split('|')[1] for line in open('blast_output.txt')])

for protein_id in protein_ids:
    positions = fetch_intron_positions(protein_id)
    print(f"{protein_id}: {positions}")

