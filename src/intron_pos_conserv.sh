#!/bin/bash
set -x
if [ $# -ne 1 ]; then
    echo "Usage: $0 <protein_id_list.txt>"
    exit 1
fi

MANE_protein_id_list=$1

echo "" > full_aln_seq_w_introns.afa

# Loop through each protein ID in the list
for MANE_protein_id in $(cat "$MANE_protein_id_list"); do

    echo "Processing MANE protein ${MANE_protein_id}"

    # Create a directory for the protein_id
    mkdir -p "${MANE_protein_id}"

    # outputs ${MANE_protein_id}_.fasta
    python3 1_MANE_pipe_extract_protein_fasta_from_multifasta.py ${MANE_protein_id} ../data/MANE.GRCh38.v1.0.refseq_protein.faa

    # outputs blastp orthologs
    blastp -query ${MANE_protein_id}_.fasta -db refseq_protein -entrez_query "txid7742[ORGN] NOT txid9443[ORGN]" -out ${MANE_protein_id}_orthologs.txt -outfmt "6 std staxids scomname ssciname" -remote

    # filter alignment
    python3 1.5_MANE_filter_alignment.py ${MANE_protein_id}_orthologs.txt ${MANE_protein_id}_orthologs_protein_ids.txt

    # insert MANE protein id to the top of the orthologs protein file
    echo "${MANE_protein_id}" > full_temp.txt && cat ${MANE_protein_id}_orthologs_protein_ids.txt >> full_temp.txt && mv full_temp.txt ${MANE_protein_id}_orthologs_protein_ids.txt

     # outputs sequence_w_introns.fasta
    ./2nd_BATCH_Pipe_MANE.sh ${MANE_protein_id}_orthologs_protein_ids.txt ${MANE_protein_id}

    # clean the sequence_w_introns.fa file
    grep -v '^$'  ${MANE_protein_id}_sequence_w_introns.fa | grep -v '^XX' > ${MANE_protein_id}_cleaned_sequence_w_introns.fa

    # MUSCLE alignment
    muscle -align ${MANE_protein_id}_cleaned_sequence_w_introns.fa -output full_aln_seq_w_introns.afa

    # reorder MUSCLE alignment to set human protein to the top
    python3 2_MANE_reorder_alignment.py full_aln_seq_w_introns.afa ${MANE_protein_id}_reorder_aln_seq_w_introns.afa ${MANE_protein_id}

    # analyze intron species names and outputs the recent introns and their best lineage
    python3 3.1_MANE_analyze_introns_species_name.py ${MANE_protein_id}_reorder_aln_seq_w_introns.afa ${MANE_protein_id} ${MANE_protein_id}_intron_analysis.txt ${MANE_protein_id}_recent_introns.txt ALL_recent_introns_new.txt

    # extracts list of species names froom reordered alignment
    python3 4_MANE_extract_species_names_from_muscle_aln.py ${MANE_protein_id}_reorder_aln_seq_w_introns.afa ${MANE_protein_id}_species_names.txt

    python3 5_MANE_make_color_annot.py ${MANE_protein_id}_intron_analysis.txt ${MANE_protein_id}
    
    python3 7_MANE_generate_tree_and_diagram.py ${MANE_protein_id}/${MANE_protein_id}_species_names.txt ${MANE_protein_id}/${MANE_protein_id}_color_annot_intron_ ${MANE_protein_id}

    mv ${MANE_protein_id}_* "${MANE_protein_id}/"

    echo "Completed processing MANE protein ${MANE_protein_id}"

    rm full_aln_seq_w_introns.afa
    rm full_temp.txt
done
