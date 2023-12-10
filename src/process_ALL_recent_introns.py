import os
import matplotlib.pyplot as plt
import numpy as np
from ete3 import NCBITaxa
import requests

# Create an NCBITaxa instance
ncbi = NCBITaxa()


# Name of the file containing intron information
input_filename = 'ALL_recent_introns.txt'

def is_earlier_than_vertebrates(taxonomic_group):
    """
    Check if the taxonomic group appears before Vertebrates in the taxonomy tree
    """
    # Convert taxonomic group name to taxid
    name2taxid = ncbi.get_name_translator([taxonomic_group])
    if not name2taxid or taxonomic_group not in name2taxid:
        return False  # Cannot find taxid for this taxonomic group, so it's safer to return False
    taxid = name2taxid[taxonomic_group][0]

    # Retrieve the lineage of the given taxid
    lineage = ncbi.get_lineage(taxid)
    names = ncbi.get_taxid_translator(lineage)
    lineage_names = [names[taxid] for taxid in lineage]

    # return True if earlier than Vertberates
    return "Vertebrata" not in lineage_names 

def is_float(value):
    try:
        float(value)
        return True
    except ValueError:
        return False

intron_scores = []
intron_lengths = []
divergence_times = []

proteins_to_skip = set()

# First pass to identify proteins based on our two criteria
with open(input_filename, 'r') as file:
    last_protein_id = None
    consecutive_introns_count = 0
    introns_seen = []

    for intron_info in file:
        intron_info = intron_info.strip()
        count_info = next(file).strip()

        parts = intron_info.split()
        protein_id = parts[0].split('_intron_')[0]
        #print(parts)
        intron_number = int(parts[0].split('_intron_')[1])

        # Logic for consecutive introns
        if last_protein_id == protein_id:
            if intron_number == last_intron_number + 1:
                consecutive_introns_count += 1
            else:
                consecutive_introns_count = 1
        else:
            consecutive_introns_count = 1

        if consecutive_introns_count >= 3:
            proteins_to_skip.add(protein_id)

        last_protein_id = protein_id
        last_intron_number = intron_number

        # Logic for proteins with only 1 intron or introns 1 and 2
        if protein_id != last_protein_id and last_protein_id is not None:
            if (1 in introns_seen) or (len(introns_seen) >= 2 and set(introns_seen) == {1, 2}):
                proteins_to_skip.add(last_protein_id)
            introns_seen = []

        introns_seen.append(intron_number)
        last_protein_id = protein_id

    if (1 in introns_seen) or (len(introns_seen) >= 2 and set(introns_seen) == {1, 2}):
        proteins_to_skip.add(last_protein_id)

# Now process the introns, skipping those from proteins identified above
with open(input_filename, 'r') as file:
    for intron_info in file:
        intron_info = intron_info.strip()
        count_info = next(file).strip()

        parts = intron_info.split()
        protein_id = parts[0].split('_intron_')[0]
        intron_number = int(parts[0].split('_intron_')[1])
        taxonomic_group = parts[1]

        #print("NP_000005.3" in proteins_to_skip)
            #print(taxonomic_group.lower() == "root") 
            #print(protein_id in proteins_to_skip) 
            #print(is_earlier_than_vertebrates(taxonomic_group))
        
        # reset proteins to skip to take into account intron 1 and intron 1,2
        #proteins_to_skip = []
        if taxonomic_group.lower() == "root" or protein_id in proteins_to_skip or is_earlier_than_vertebrates(taxonomic_group):
            continue

        aligned_count = int(count_info.split('X aligned count: ')[1].split(' ')[0])
        not_aligned_count = int(count_info.split('X not aligned count: ')[1])

        total_count = aligned_count + not_aligned_count
        confidence_score = (not_aligned_count / total_count) * 100
        #print("aligned_count", aligned_count, "not_aligned_count", not_aligned_count, "confidence_score", confidence_score)

        if confidence_score > 50:
            score_filename = f'{protein_id}/{protein_id}_intron_position_score.txt'
            if os.path.exists(score_filename):
                with open(score_filename, 'r') as score_file:
                    lines = score_file.readlines()
                    for i, line in enumerate(lines):
                        if protein_id in line:
                            scores = lines[i+1].strip().split(',')
                            scores = [score for score in scores if is_float(score.strip())]
                            if intron_number <= len(scores):
                                intron_score = float(scores[intron_number - 1])
                                intron_scores.append(intron_score)
                                print(f'{protein_id}\t{intron_number}\t{taxonomic_group}\t{confidence_score:.2f}\t{intron_score}')
                            break

            length_filename = f'{protein_id}/{protein_id}_intron_lengths.txt'
            if os.path.exists(length_filename):
                with open(length_filename, 'r') as length_file:
                    lines = length_file.readlines()
                    for i, line in enumerate(lines):
                        if protein_id in line:
                            lengths = lines[i+1].strip().split(',')
                            lengths = [length for length in lengths if is_float(length.strip())]
                            if intron_number <= len(lengths):
                                intron_length = float(lengths[intron_number - 1])
                                intron_lengths.append(intron_length)
                            break

        # At the end of the loop, if you've switched to a different protein, reset the skip_protein flag
        if last_protein_id != protein_id:
            skip_protein = False

# Creating the position score histogram plot
plt.figure()
plt.hist(intron_scores, bins=range(0, 105, 5), edgecolor='black')
plt.title('Recently-Gained Intron Count vs Intron Position Score')
plt.xlabel('Intron Position Score = Recent Intron Position / Full Protein Length')
plt.ylabel('Intron Count')
plt.savefig('recent_introns_position_score_histogram.pdf')

# Determine the maximum length
max_length = max(intron_lengths)

# Creating the lengths histogram plot with logarithmic scale x-axis
plt.figure()
bins = np.logspace(np.log10(min(intron_lengths)), np.log10(max(intron_lengths)), 50)  # Create log-spaced bins
plt.hist(intron_lengths, bins=bins, edgecolor='black')
plt.xscale('log')  # Set the x-axis to a logarithmic scale
plt.gca().set_xlim(left=min(intron_lengths))  # Adjust the lower bound of the x-axis
plt.title('Recently-Gained Intron Count vs Intron Lengths (Log Scale)')
plt.xlabel('Intron Lengths (log scale)')
plt.ylabel('Intron Count')
plt.savefig('recent_introns_lengths_histogram_log_scale_adjusted.pdf')


# Creating the lengths histogram plot with logarithmic scale x-axis
plt.figure()
bins = np.logspace(np.log10(min(intron_lengths)), np.log10(max(intron_lengths)), 50)  # Create log-spaced bins
plt.hist(intron_lengths, bins=bins, edgecolor='black')
plt.gca().set_xlim(left=min(intron_lengths))  # Adjust the lower bound of the x-axis
plt.title('Recently-Gained Intron Count vs Intron Lengths')
plt.xlabel('Intron Lengths')
plt.ylabel('Intron Count')
plt.savefig('recent_introns_lengths_histogram.pdf')

# Creating the lengths histogram plot NO log scale x-axis
plt.figure()
bins = np.logspace(np.log10(min(intron_lengths)), np.log10(max(intron_lengths)), 50)  # Cre
