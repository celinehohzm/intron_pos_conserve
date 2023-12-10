#!/bin/bash

# File containing your protein IDs
PROTEIN_IDS="filtered_ALL_recent_introns_proteinids.txt"

# Loop through each protein ID
while read -r PROTEIN_ID; do
    # Use elink to link protein IDs to gene IDs, then esearch and efetch to get the gene name
    GENE_NAME=$(efetch -db protein -id "$PROTEIN_ID" -format docsum | 
                xtract -pattern DocumentSummary -element Caption)

    echo "$PROTEIN_ID -> $GENE_NAME"
done < "$PROTEIN_IDS"

