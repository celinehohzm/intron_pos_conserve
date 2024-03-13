#!/bin/bash
set -x
if [ $# -ne 1 ]; then
    echo "Usage: $0 <protein_id_list.txt>"
    exit 1
fi

# Activate environment, adjust as necessary if using Conda or another environment manager
# source path/to/environment/bin/activate

MANE_protein_id_list=$1

# Clear the output files
echo "" > sequence_w_introns.fasta
echo "" > aln_seq_w_introns.afa

# Define the path to the scripts and data
SCRIPTS_DIR="$(pwd)/src"
DATA_DIR="$(pwd)/data" # Adjust this if your data is located elsewhere
BIN_DIR="$(pwd)/bin"

# Loop through each protein ID in the list
for MANE_protein_id in $(cat "$MANE_protein_id_list"); do

    echo "Processing MANE protein ${MANE_protein_id}"

    # Create a directory for the protein_id
    mkdir -p "${MANE_protein_id}"

    # Use the Python scripts from the src directory
    python3 "${SCRIPTS_DIR}/1_MANE_pipe_extract_protein_fasta_from_multifasta.py" ${MANE_protein_id} "${DATA_DIR}/MANE.GRCh38.v1.0.refseq_protein.faa"

    # Assuming BLAST is installed globally or available in the user's environment
    blastp -query "${MANE_protein_id}.fasta" -db "${DATA_DIR}/refseq_protein_db/refseq_protein" -out "${MANE_protein_id}_orthologs.txt" -taxidlist "${DATA_DIR}/vertebrates_except_primates_taxids.txt" -outfmt "6 std staxids scomname ssciname"

    # Continue using the scripts for processing
    python3 "${SCRIPTS_DIR}/1.5_MANE_filter_alignment.py" ${MANE_protein_id}_orthologs.txt ${MANE_protein_id}_orthologs_protein_ids.txt

    # Execute the second pipeline script
    "${BIN_DIR}/2nd_Pipe_intron_pos.sh" ${MANE_protein_id}_orthologs_protein_ids.txt

    # Further processing
    grep -v '^$' "${MANE_protein_id}_sequence_w_introns.fa" | grep -v '^XX' > "${MANE_protein_id}_cleaned_sequence_w_introns.fa"

    # Assuming MUSCLE is installed globally or available in the user's environment
    muscle -align "${MANE_protein_id}_cleaned_sequence_w_introns.fa" -output aln_seq_w_introns.afa

    python3 "${SCRIPTS_DIR}/2_MANE_reorder_alignment.py" aln_seq_w_introns.afa ${MANE_protein_id}_reorder_aln_seq_w_introns.afa ${MANE_protein_id}

    python3 "${SCRIPTS_DIR}/3.1_MANE_analyze_introns_species_name.py" ${MANE_protein_id}_reorder_aln_seq_w_introns.afa ${MANE_protein_id} ${MANE_protein_id}_intron_analysis.txt ${MANE_protein_id}_recent_introns.txt ALL_recent_introns.txt

    python3 "${SCRIPTS_DIR}/4_MANE_extract_species_names_from_muscle_aln.py" ${MANE_protein_id}_reorder_aln_seq_w_introns.afa ${MANE_protein_id}_species_names.txt

    python3 "${SCRIPTS_DIR}/5_MANE_make_color_annot.py" ${MANE_protein_id}_intron_analysis.txt ${MANE_protein_id}

    python3 "${SCRIPTS_DIR}/7_generate_tree_and_diagram.py" "${MANE_protein_id}/${MANE_protein_id}_species_names.txt" "${MANE_protein_id}/${MANE_protein_id}_color_annot_intron" "${MANE_protein_id}/${MANE_protein_id}_tree.pdf"

    mv ${MANE_protein_id}_* "${MANE_protein_id}/"

    echo "Completed processing MANE protein ${MANE_protein_id}"

done


