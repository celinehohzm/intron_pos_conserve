# Tracking Recent Intron Gains in Primates

## Goal: 
To keep track of the number of recent intron gains in primates

## Why is this important??
Can improve annotation - 
If introns are conserved in a certain faimily, we can infer that other species who belong to that family have this intron too!

## Input: 
A file with any number of NCBI protein accession numbers (eg. NP_000005.3)

## Output:
- A phylogenetic tree diagram that shows which species have this intron conserved (red = not conserved, blue = conserved)
- An analysis file that shows, for each intron in that protein:
  - species list that has that intron aligned,
  - species list that has does not have that intron aligned,
  - the subgroup in the tree that best fits that intron gain

## What counts as an intron gain event?



## Methods:
1. Get list of similar/orthologous proteins via BlastP-ing MANE protein
2. Keep only most similar protein for each species
3. Get each orthologous proteins’ genomic interval coding coordinates from NCBI via Entrez direct
4. Add “X” at the intron positions into the amino acid sequence
5. Multi-sequence align MANE + all of its orthologous proteins via MUSCLE 
6. For each intron, count how many species have introns aligned/not unaligned at that position
7. Build tree that shows intron conservation via ETEToolkit


