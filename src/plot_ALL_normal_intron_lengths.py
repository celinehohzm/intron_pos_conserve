import os
import re
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter

# Root directory where the protein directories begin
root_dir = "/home/choh1/intron_pos_conserv/blastdb"  # Modify this to point to your root directory

def extract_first_numbers_from_file(file_path):
    with open(file_path, 'r') as f:
        content = f.readlines()
        for line in content:
            if line.startswith(">"):
                continue
            number_strings = [n for n in re.split(r',\s*', line.strip()) if n]
            numbers = [int(n) for n in number_strings]
            return numbers
    return []

def main():
    all_first_numbers = []

    for subdir in [d for d in os.listdir(root_dir) if d.startswith("NP_0")]:
        file_path = os.path.join(root_dir, subdir, f"{subdir}_intron_lengths.txt")
        if os.path.isfile(file_path):
            first_numbers = extract_first_numbers_from_file(file_path)
            if first_numbers:
                all_first_numbers.extend(first_numbers)

    # Count occurrences of each intron length
    intron_counts = Counter(all_first_numbers)

    # Extract intron lengths and their counts for plotting
    intron_lengths = list(intron_counts.keys())
    intron_counts_list = list(intron_counts.values())

    # Create the histogram plot
    plt.figure()
    bins = np.logspace(np.log10(min(intron_lengths)), np.log10(max(intron_lengths)), 50)
    plt.hist(intron_lengths, weights=intron_counts_list, bins=bins, edgecolor='black')
    plt.xscale('log')
    plt.gca().set_xlim(left=min(intron_lengths))
    plt.title('Normal Intron Count vs Intron Lengths')
    plt.xlabel('Intron Lengths (log scale)')
    plt.ylabel('Intron Count')
    plt.savefig('normal_introns_lengths_histogram_log_scale_adjusted.pdf')

if __name__ == "__main__":
    main()

