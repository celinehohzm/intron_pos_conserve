#!/bin/bash

input_file="test_gene_ids.txt"
chunk_size=10
output_file="test_gene_table.txt"

# Split the input file into smaller chunks
split -l $chunk_size $input_file ${input_file}_chunk_

# Process each chunk
for chunk in ${input_file}_chunk_*; do
    efetch -db gene -id "$(cat $chunk)" -format gene_table >> $output_file
    sleep 5  # Wait for 5 seconds before the next request
done

# Cleanup the chunks
rm ${input_file}_chunk_*

