# Tracking Recent Intron Gains in Primates

## Goal: 
To keep track of the number of recent intron gains in primates

## Why is this important??
Can improve annotation - \
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
- Non-conserved in a minimum of 10% of all examined orthologous species.
- We allotted "intron conservation scores" to each topology clusters following the criteria of maximized intron preservation and minimized non-conservation
- Subgroup with highest intron conservation score is chosen as the timeline where the intron was first gained.


## Methods:
1. Get list of similar/orthologous proteins via BlastP-ing MANE protein
2. Keep only most similar protein for each species
3. Get each orthologous proteins’ genomic interval coding coordinates from NCBI via Entrez direct
4. Add “X” at the intron positions into the amino acid sequence
5. Multi-sequence align MANE + all of its orthologous proteins via MUSCLE 
6. For each intron, count how many species have introns aligned/not unaligned at that position
7. Build tree that shows intron conservation via ETEToolkit


## Steps to run
### 1. Creating environment
```
conda env create -f environment.yml
```
```
conda activate intron_conserv_env
```

### 2. Prepare protein_ids.txt input file
Include all your NCBI protein ids in the file, eg. ```test/test_proteinid.txt```

### 3. Run intron_pos_conserv
``` ./intron_pos_conserv.sh ../test/test_proteinid.txt```

### 4. See results
A directory will be created for each separate protein_id.  \
Each directory contains: 
1. A NCBI Taxa Tree with color annotation that shows if each intron in the protein is aligned or not aligned to the human protein introns (red = aligned, blue = not aligned)
2. An intron analysis file for each protein that shows for each intron, what species have that intron aligned, and which species do not
3. A ALL_recent_introns.txt file that tracks all the recent introns (useful if you have more than one protein id in the input file)

