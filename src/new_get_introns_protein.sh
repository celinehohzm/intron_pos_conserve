#!/bin/bash

echo "Starting script..."

# Read each line of the input file
IFS=$'\n' read -d '' -r -a lines < "NP_001013860.1_proteinids.txt"
echo "Read lines from filtered_ALL_recent_introns.txt"

# Loop through each protein ID in the array
for line in "${lines[@]}"
do
    echo "Processing line: $line"

    # Skip empty lines or lines starting with '-'
    if [[ -z "$line" || "$line" == -* ]]; then
            echo "Skipping empty line or line starting with '-'"
            continue
    fi

    # Extract protein_id and intron_number from the line
    if [[ "$line" =~ ^(NP_[0-9]+\.[0-9]+)[[:space:]]+([0-9]+)[[:space:]]+([A-Za-z]+) ]]; then
            protein_id=${BASH_REMATCH[1]}
            intron_number=${BASH_REMATCH[2]}
            echo "Extracted Protein ID: ${protein_id}, Intron number: ${intron_number}"
    else
            echo "Invalid line format: $line"
            continue
    fi

    echo "Running 1_extract_table.py for ${protein_id}"
    python3 1_extract_table.py ${first_protein_id}/${protein_id}_gene_table.txt ${protein_id} > recentintrons_extracted_table.txt

    echo "Running recentintrons_extract_exon_intron_exon_aaseq.py for ${protein_id}"
    python3 new_recentintrons_extract_exon_intron_exon_aaseq.py recentintrons_extracted_table.txt ${protein_id} ${intron_number} ${first_protein_id}/${protein_id}_sequences.fasta

done
