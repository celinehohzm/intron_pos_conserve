#!/bin/bash

# Path to the gene2refseq file
GENE2REFSEQ="gene2refseq"

# Loop through each chunk in the blast_files directory
for CHUNK_FILE in blast_files/MANE_refseq_protein_chunk*_orthologs.txt; do
    python blast_files_sep_geneid.py "$CHUNK_FILE" $GENE2REFSEQ
done

echo "All chunks processed!"

