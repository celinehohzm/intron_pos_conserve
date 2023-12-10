#!/bin/bash

echo "Inside blast_run.sh script..."

# Path configurations
BLASTP_PATH="/home/choh1/ncbi-blast-2.14.0+/bin/blastp"
DB_PATH="/home/choh1/ncbi-blast-2.14.0+/refseq_protein_db/refseq_protein"
TAXIDLIST_PATH="/home/choh1/ncbi-blast-2.14.0+/bin/vertebrates_taxids.txt"

input_file=$1
filename=$(basename $input_file)
output_file="/home/choh1/intron_pos_conserv/blastdb/blast_files/${filename%.fasta}_orthologs.txt"

# Extracting the number of threads from the second argument.
num_threads=$2

echo "    Input file: $input_file"
echo "    Number of threads: $num_threads"
echo "    Output file will be: $output_file"

# Debug messages to verify the paths and file contents
echo "Checking paths and file contents..."
echo "    BLASTP Path: $BLASTP_PATH"
echo "    DB Path: $DB_PATH"
echo "    TAXIDLIST Path: $TAXIDLIST_PATH"

# Check if paths are accessible
if [[ ! -x "$BLASTP_PATH" ]]; then
    echo "    Error: BLASTP path is not executable or doesn't exist!"
fi

if [[ ! -f "$DB_PATH.pto" ]]; then
    echo "    Error: BLAST Database files not found!"
fi

if [[ ! -f "$TAXIDLIST_PATH" ]]; then
    echo "    Error: Taxid list file doesn't exist!"
fi

echo "Starting BLAST..."
$BLASTP_PATH -query $input_file \
             -db $DB_PATH \
             -out $output_file \
             -taxidlist $TAXIDLIST_PATH \
             -outfmt "6 qseqid sseqid" \
             -num_threads $num_threads 2>&1 | tee -a "/home/choh1/intron_pos_conserv/blastdb/blast_files/${filename%.fasta}_BLAST_debug.log"

echo "BLAST finished for $input_file!"

# Check the output file contents
if [[ ! -s "$output_file" ]]; then
    echo "Warning: BLAST output file is empty!"
fi

echo "Debug messages finished. BLAST finished for $input_file!"

