#!/bin/bash

# Read each line of the input file
IFS=$'\n' read -d '' -r -a lines < "test_ALL_recent_introns_new.txt"
export PATH=$PATH
echo $PATH

# Loop through each line
for line in "${lines[@]}"
do
    # Skip empty lines or lines starting with '-'
    if [[ -z "$line" || "$line" == -* ]]; then
        continue
    fi

    # Extract protein_id and intron_number from the line
    if [[ "$line" =~ ([^[:space:]]+)_intron_([0-9]+) ]]; then
        protein_id=${BASH_REMATCH[1]}
        intron_number=${BASH_REMATCH[2]}
    else
        echo "Invalid line format: $line"
        continue
    fi

    echo "Processing Protein ID: $protein_id, Intron number: $intron_number"

    # Fetch the whole_table.txt
    echo "Fetching whole table..."
    cmd="esearch -db protein -query \"${protein_id}\" | elink -target gene | efetch -format gene_table > RECENTINTRONS_whole_table.txt"
    echo "Command: $cmd"
    eval $cmd

    # Extract the table from whole_table.txt
    echo "Extracting table..."
    read genomic_id start_pos end_pos reversed_order <<< $(python3 recentintrons_extract_table.py RECENTINTRONS_whole_table.txt ${protein_id} ${intron_number})
    echo $genomic_id $start_pos $end_pos $reversed_order

    # Extract transcript's DNA sequence
    cmd="efetch -db nucleotide -id ${genomic_id} -format fasta > RECENTINTRONS_genomic_dnaseq.fa"
    echo "Command: $cmd"
    eval $cmd

    echo "Extracting intron DNA sequence..."
    #python3 recentintrons_extract_intron_dna_seq.py ${protein_id} ${intron_number} RECENTINTRONS_extracted_table.txt RECENTINTRONS_genomic_dnaseq.fa >> recentintrons_new_multi_dnaseq.fa

    python3 recentintrons_extract_intron_dnaseq.py RECENTINTRONS_genomic_dnaseq.fa ${protein_id} ${intron_number} ${start_pos} ${end_pos} ${reversed_order} test_recentintrons_dnaseq_new.fa
    echo "Completed Protein ID: $protein_id, Intron number: $intron_number"
    echo "----------------------------------"

    # Move files to a separate directory if needed
    # mv RECENTINTRONS_* RECENTINTRONS/

done
