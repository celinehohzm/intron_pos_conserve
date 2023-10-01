#!/bin/bash
set -e  # Stop script on the first error

if [ $# -ne 2 ]; then
    echo "Usage: $0 <protein_id_list.txt> <first_protein_id>"
    exit 1
fi

orthologs_protein_id_list=$1
first_protein_id=$2

# Debug: Print input arguments
echo "Input arguments:"
echo "Orthologs protein ID list: $orthologs_protein_id_list"
echo "First protein ID: $first_protein_id"

# Read the first protein ID and store it in a variable
echo "First protein ID: $first_protein_id"

# Debug: Print the content of the protein ID list file
echo "Content of protein ID list file:"
cat $orthologs_protein_id_list
echo "-----------------------"

# Fetch the amino_acid_fasta.txt
if efetch -db protein -id "$(cat $orthologs_protein_id_list)" -format fasta | \
   awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);} END {printf("\n");}' > ${first_protein_id}_sequences.fasta; then
    echo "Fetched amino_acid_fasta.txt successfully"
else
    echo "Failed to fetch amino_acid_fasta.txt"
    exit 1
fi

# Debug: Print fetched amino_acid_fasta.txt content
echo "Content of amino_acid_fasta.txt:"
cat ${first_protein_id}_sequences.fasta
echo "-----------------------"

# Fetch the whole_table.txt
if cat $orthologs_protein_id_list | epost -db protein | elink -target gene | efetch -format uid > ${first_protein_id}_gene_ids.txt; then
    echo "Fetched whole_table.txt successfully"
else
    echo "Failed to fetch whole_table.txt"
    exit 1
fi

# Debug: Print fetched whole_table.txt content
echo "Content of whole_table.txt:"
cat ${first_protein_id}_gene_ids.txt
echo "-----------------------"

chunk_size=5
# Split the input file into smaller chunks
split -l $chunk_size ${first_protein_id}_gene_ids.txt ${first_protein_id}_gene_ids_chunk_

# Process each chunk
for chunk in ${first_protein_id}_gene_ids_chunk_*; do
    if efetch -db gene -id "$(cat $chunk)" -format gene_table >> ${first_protein_id}_gene_table.txt; then
        echo "Processed chunk $chunk successfully"
    else
        echo "Failed to process chunk $chunk"
        exit 1
    fi
    sleep 5  # Wait for 5 seconds before the next request
done

# Cleanup the chunks
rm ${first_protein_id}_gene_ids_chunk_*

# Loop through each protein ID in the list
for protein_id in $(cat "$orthologs_protein_id_list"); do
    echo "Processing ortholog protein ${protein_id}"

    if python3 1_extract_table.py ${first_protein_id}_gene_table.txt ${protein_id} > full_extracted_table.txt; then
        echo "Extracted table for ${protein_id} successfully"
    else
        echo "Failed to extract table for ${protein_id}"
        exit 1
    fi

    if python3 2.1_extract_n_insert_introns.py full_extracted_table.txt ${first_protein_id}_sequences.fasta ${protein_id} ${first_protein_id}; then
        echo "Inserted introns for ${protein_id} successfully"
    else
        echo "Failed to insert introns for ${protein_id}"
        exit 1
    fi

    echo "Completed processing ${protein_id}"
done

echo "Script completed successfully"

