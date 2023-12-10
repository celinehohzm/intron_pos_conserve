#!/bin/bash

# Check if a protein ID was provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 PROTEIN_ID"
    exit 1
fi

# Get the protein ID from command line argument
protein_id=$1

# Fetch the amino_acid_fasta.txt
efetch -db protein -id "${protein_id}" -format fasta
