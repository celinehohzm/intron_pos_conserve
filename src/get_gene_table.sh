#!/bin/bash

# Set your email and API key (Already set in bash_profile)
# export NCBI_API_KEY="ffc920f74d224ccc6e9e89822ed3ee355b08"
# export NCBI_EMAIL="celinehohzm@gmail.com"
echo "API Key: $NCBI_API_KEY"


chunk_size=5

for gene_id_file in blast_files/*_orthologs_vertebrates_n_primates_geneids.txt; do
    # Remove path and extract protein ID from filename
    file_name=$(basename "$gene_id_file")
    first_protein_id="${file_name/_orthologs_vertebrates_n_primates_geneids.txt/}"

    #first_protein_id=$(basename "$gene_id_file" "_orthologs_vertebrates_n_primates_geneids.txt")

    # Check if the gene table file already exists
    if [[ -f "${first_protein_id}_gene_table.txt" ]]; then
        echo "Gene table for ${first_protein_id} already exists. Skipping..."
        continue
    fi

    echo "Processing $first_protein_id ..."

    # If you want to reset the gene table for each protein ID, uncomment the next line:
    # > ${first_protein_id}_gene_table.txt

    # Split the input file into smaller chunks
    split -l $chunk_size "$gene_id_file" "${first_protein_id}_gene_ids_chunk_"

    # Process each chunk
    for chunk in ${first_protein_id}_gene_ids_chunk_*; do
        efetch -db gene -id "$(cat $chunk)" -format gene_table >> ${first_protein_id}_gene_table.txt
        sleep 3  # Wait for 3 seconds before the next request
    done

    # Cleanup the chunks
    rm ${first_protein_id}_gene_ids_chunk_*
done

mv "${first_protein_id}_gene_table.txt" "new_gene_tables/"

echo "All files processed."

