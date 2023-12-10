#!/bin/bash

# Check if a protein ID was provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 PROTEIN_ID"
    exit 1
fi

# Get the protein ID from command line argument
protein_id=$1

# Fetch the whole_table.txt
esearch -db protein -query "${protein_id}" | elink -target gene | efetch -format gene_table > whole_table.txt

# Extract the table from whole_table.txt
python3 1_extract_table.py whole_table.txt ${protein_id}

