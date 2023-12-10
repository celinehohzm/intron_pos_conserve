#!/bin/bash

# Go through all directories
for protein_id in */; do
    # Remove trailing "/"
    protein_id=${protein_id%/}
    echo "Running script for ${protein_id}"
    python 7_generate_tree_and_diagram.py "${protein_id}/${protein_id}_species_names.txt" "${protein_id}/${protein_id}_color_annot_intron" "new_tree_diagrams/${protein_id}_tree.pdf"
done
