#!/bin/bash

echo "Starting script..."

# Read each line of the input file
IFS=$'\n' read -d '' -r -a lines < "filtered_ALL_recent_introns.txt"
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
            tax_group=${BASH_REMATCH[3]}
            echo "Extracted Protein ID: ${protein_id}, Intron number: ${intron_number}, Tax group: ${tax_group}"
    else
            echo "Invalid line format: $line"
            continue
    fi

    echo "Running 1_extract_table.py for ${protein_id}"
    python3 1_extract_table.py ${protein_id}/${protein_id}_gene_table.txt ${protein_id} > recentintrons_extracted_table.txt

    echo "Running recentintrons_extract_exon_intron_exon_aaseq.py for ${protein_id}"
    intron_len=$(python3 recentintrons_extract_exon_intron_exon_aaseq.py recentintrons_extracted_table.txt ${protein_id} ${intron_number} ${protein_id}/${protein_id}_sequences.fasta)

    echo "Intron length: ${intron_len}"

    output=$(/home/choh1/ncbi-blast-2.14.0+/bin/get_species_taxids.sh -n ${tax_group})
    # Extract the taxid using grep and awk
    taxid=$(echo "$output" | grep -oP 'Taxid : \K\d+')

    # Check if taxid is extracted
    if [ -z "$taxid" ]; then
    	echo "Taxid not found."
    	exit 1
    fi
    
    echo "Extracted Taxid: $taxid"

    # Use the taxid in the next command
    /home/choh1/ncbi-blast-2.14.0+/bin/get_species_taxids.sh -t "$taxid" > "${taxid}.txids"
    echo "Taxid list saved to ${taxid}.txids"

    echo "Running BLAST for blast1.fa excluding ${tax_group}"
    /home/choh1/ncbi-blast-2.14.0+/bin/blastp -db /home/choh1/ncbi-blast-2.14.0+/refseq_protein_db/refseq_protein -query blast1.fa -negative_taxidlist ${taxid}.txids -outfmt "6 std staxids scomname ssciname sseq" -out ${protein_id}_blast1.txt
    #/ccb/sw/bin/blastp -query blast1.fa -db refseq_protein -out blast1.txt -outfmt "6 std staxids scomname ssciname sseq" -entrez_query "NOT ${tax_group}[Organism]" -remote

    echo "Running BLAST for blast2.fa excluding ${tax_group}"
    /home/choh1/ncbi-blast-2.14.0+/bin/blastp -db /home/choh1/ncbi-blast-2.14.0+/refseq_protein_db/refseq_protein -query blast2.fa -negative_taxidlist ${taxid}.txids -outfmt "6 std staxids scomname ssciname sseq" -out ${protein_id}_blast2.txt 
    #/ccb/sw/bin/blastp -query blast2.fa -db refseq_protein -out blast2.txt -outfmt "6 std staxids scomname ssciname sseq" -entrez_query "NOT ${tax_group}[Organism]" -remote

    echo "Running compare_recentintrons_blast_outputs_for_intronizations.py for ${protein_id}"
    python compare_recentintrons_blast_outputs_for_intronizations.py ${protein_id} ${intron_len}

done

echo "Script completed."

