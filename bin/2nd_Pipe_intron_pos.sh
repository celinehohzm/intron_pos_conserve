#!/bin/bash
set -x
if [ $# -ne 1 ]; then
    echo "Usage: $0 <orthologs_protein_id_list.txt>"
    exit 1
fi

orthologs_protein_id_list=$1

# Define the path to the scripts
SCRIPTS_DIR="$(pwd)/src"

# Read the first protein ID and store it in a variable
first_protein_id=$(head -n 1 "$orthologs_protein_id_list")
echo "First protein ID: $first_protein_id"

# Clear the output file
echo "" > sequence_w_introns.fa

# Loop through each protein ID in the list
for protein_id in $(cat "$orthologs_protein_id_list"); do

    echo "Processing ortholog protein ${protein_id}"

    # Fetch the amino_acid_fasta.txt
    efetch -db protein -id "${protein_id}" -format fasta > amino_acid.fa

    # Fetch the whole_table.txt
    esearch -db protein -query "${protein_id}" | elink -target gene | efetch -format gene_table > whole_table.txt

    # Extract the table from whole_table.txt
    python3 "${SCRIPTS_DIR}/1_extract_table.py" whole_table.txt ${protein_id} > extracted_table.txt

    # Process multi-line fasta files, outputs org_aa_sequence.fasta
    python3 "${SCRIPTS_DIR}/1.5_extract_seq_from_fasta.py" amino_acid.fa > org_aa_sequence.fa

    cat org_aa_sequence.fa >> ${first_protein_id}_org_aa_sequence.fa

    # Outputs ${first_protein_id}_sequence_w_introns.fasta and ${first_protein_id}_intron_positions.txt
    python3 "${SCRIPTS_DIR}/2.1_extract_n_insert_introns.py" extracted_table.txt org_aa_sequence.fa ${first_protein_id}

    echo "Completed processing ${protein_id}"

done

