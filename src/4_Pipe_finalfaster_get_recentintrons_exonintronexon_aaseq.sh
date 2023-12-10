#!/bin/bash


# Read each line of the input file
IFS=$'\n' read -d '' -r -a lines < "filtered_ALL_recent_introns.txt"
export PATH=$PATH
echo $PATH

# Loop through each protein ID in the array
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

    	echo "Processing Protein ID: ${protein_id}, Intron number: ${intron_number}"

		python3 1_extract_table.py ${protein_id}/${protein_id}_gene_table.txt ${protein_id} > recentintrons_extracted_table.txt

        # BLASTP this recentintrons WITHOUT XXX, Outputs amino acid sequences of exon+intron+exon  intron_num in recentintrons_exons_of_intron_0911.fa
        python3 recentintrons_extract_exon_intron_exon_aaseq.py recentintrons_extracted_table.txt ${protein_id} ${intron_number} ${protein_id}/${protein_id}_.fasta

		## Show Steven this to verify, outputs recentintrons WITH XXX, outputs into recentintrons_exons_of_intronXXXX_0911.fa
		python3 recentintrons_extract_exon_intronXXX_exon_aaseq.py recentintrons_extracted_table.txt ${protein_id} ${intron_number} ${protein_id}/${protein_id}_.fasta
done


#/home/choh1/ncbi-blast-2.14.0+/bin/blastp -query recentintrons_exons_of_intron_0911.fa -db /home/choh1/ncbi-blast-2.14.0+/refseq_protein_db/refseq_protein -outfmt "6 qseqid sseqid pident length mismatch gapopen gaps qstart qend sstart send evalue bitscore qcovs qlen slen" -out recentintrons_potentialintronization_0911.txt -num_threads 4
#awk '$6 < 1 && $15 < $16' recentintrons_potentialintronization_0911.txt > filtered_recentintrons_potentialintronization_0911.txt
# 1 gap open ($6 > 0), query sequence < subject sequence ($15 < $16)
# dont want 0 gap open, 
 
