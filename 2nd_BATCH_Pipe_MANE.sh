#!/bin/bash
# set -x
# if [ $# -ne 2 ]; then
#     echo "Usage: $0 <protein_id_list.txt> <first_protein_id>"
#     exit 1
# fi

orthologs_protein_id_list=$1
first_protein_id=$2

# Read the first protein ID and store it in a variable
# first_protein_id=$(head -n 1 "$orthologs_protein_id_list")
echo "First protein ID: $first_protein_id"

# Fetch the amino_acid_fasta.txt
efetch -db protein -id "$(cat $orthologs_protein_id_list)" -format fasta | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);} END {printf("\n");}' > ${first_protein_id}_sequences.fasta
# efetch -db protein -id "$(cat NP_000005.3/NP_000005.3_orthologs_protein_ids.txt)" -format fasta | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);} END {printf("\n");}' > NP_000005.3_sequences.fasta

# Fetch the whole_table.txt
cat $orthologs_protein_id_list | epost -db protein | elink -target gene | efetch -format uid > ${first_protein_id}_gene_ids.txt
# Fetching all gene tables is too much for ncbi, so let's split it up
## efetch -db gene -id "$(cat ${first_protein_id}_gene_ids.txt)" -format gene_table > ${first_protein_id}_gene_table.txt

chunk_size=5
# Split the input file into smaller chunks
split -l $chunk_size ${first_protein_id}_gene_ids.txt ${first_protein_id}_gene_ids_chunk_

# Process each chunk
for chunk in ${first_protein_id}_gene_ids_chunk_*; do
	efetch -db gene -id "$(cat $chunk)" -format gene_table >> ${first_protein_id}_gene_table.txt
	sleep 5  # Wait for 5 seconds before the next request
done

# Cleanup the chunks
rm ${first_protein_id}_gene_ids_chunk_*

# Loop through each protein ID in the list
for protein_id in $(cat "$orthologs_protein_id_list"); do

    echo "Processing ortholog protein ${protein_id}"

    # Extract the table from whole_table.txt
    python3 1_extract_table.py ${first_protein_id}_gene_table.txt ${protein_id} > full_extracted_table.txt
    # python3 1_extract_table.py NP_000005.3_gene_table.txt NP_000005.3 > extracted_table.txt

    # Outputs ${first_protein_id}_sequence_w_introns.fasta and ${first_protein_id}_intron_positions.txt
    python3 2_extract_n_insert_introns.py full_extracted_table.txt ${first_protein_id}_sequences.fasta ${protein_id} ${first_protein_id}
    # python3 2.1_extract_n_insert_introns.py extracted_table.txt NP_000005.3_sequences.fasta NP_000005.3

    echo "Completed processing ${protein_id}"

done
