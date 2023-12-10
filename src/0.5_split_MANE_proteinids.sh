#!/bin/bash

# count the total number of lines
total_lines=$(wc -l < MANE.GRCh38.v1.0.refseq_protein_ids.txt)

# calculate lines per file
((lines_per_file = (total_lines + 19) / 20))  # this ensures that the division rounds up

# split the file
split -l $lines_per_file --numeric-suffixes=1 --additional-suffix=.txt MANE.GRCh38.v1.0.refseq_protein_ids.txt MANE_protein_ids_chunk

