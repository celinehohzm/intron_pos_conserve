import os
import matplotlib.pyplot as plt
import numpy as np
from ete3 import NCBITaxa
import collections

# Initialize NCBITaxa object
ncbi = NCBITaxa()

# Name of the file containing intron information
input_filename = 'ALL_recent_introns.txt'

def is_float(value):
    try:
        float(value)
        return True
    except ValueError:
        return False

# Determine the "age" of Vertebrata
vertebrata_taxid = ncbi.get_name_translator(["Vertebrata"])["Vertebrata"][0]
vertebrata_age = len(ncbi.get_lineage(vertebrata_taxid))

# We will store the frequency of introns for each taxa age in this dictionary
taxa_age_intron_frequency = collections.defaultdict(int)

# Open the input file and process each intron entry
with open(input_filename, 'r') as file:
    for intron_info in file:
        intron_info = intron_info.strip()
        count_info = next(file).strip()

        parts = intron_info.split()
        protein_id = parts[0].split('_intron_')[0]
        intron_number = int(parts[0].split('_intron_')[1])
        taxonomic_group = parts[1]

        aligned_count = int(count_info.split('X aligned count: ')[1].split(' ')[0])
        not_aligned_count = int(count_info.split('X not aligned count: ')[1])

        total_count = aligned_count + not_aligned_count
        confidence_score = (not_aligned_count / total_count) * 100

        if confidence_score > 50:
            name2taxid = ncbi.get_name_translator([taxonomic_group])
            if taxonomic_group in name2taxid and name2taxid[taxonomic_group]:
                taxid = name2taxid[taxonomic_group][0]
                lineage = ncbi.get_lineage(taxid)
                taxa_age = len(lineage)

                # Skip taxa that are older than or equal to Vertebrata
                if taxa_age <= vertebrata_age:
                    continue

                # Increment the count of introns for this taxa age
                taxa_age_intron_frequency[taxa_age] += 1

                print(f'Processed Intron: {protein_id} {intron_number} {taxonomic_group} {confidence_score:.2f}')

# Preparing data for the plot
taxa_ages = list(taxa_age_intron_frequency.keys())
intron_frequencies = list(taxa_age_intron_frequency.values())

# Creating the scatter plot without labels
plt.figure(figsize=(10,10))
plt.scatter(taxa_ages, intron_frequencies, alpha=0.5)

plt.title('Intron Frequency vs Taxa Age')
plt.xlabel('Taxa Age (Number of Ancestors in Lineage)')
plt.ylabel('Intron Frequency')

plt.ylim(0, max(intron_frequencies) + 10)

plt.savefig('intron_frequency_vs_taxa_age_no_labels.pdf')
plt.show()

