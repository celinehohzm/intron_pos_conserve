import matplotlib.pyplot as plt
from collections import defaultdict

# Given divergence_times dictionary
divergence_times = {
    "Actinopterygii": 450,
    "Amniota": 321,
    "Artiodactyla": 62.5,
    "Boreoeutheria": 90,
    "Callithrix": 22.5,
    "Caniformia": 48.5,
    "Carnivora": 48.5,
    "Catarrhini": 27.5,
    "Cebidae": 25,
    "Cercopithecidae": 22.5,
    "Euarchontoglires": 82.5,
    "Euteleostomi": 450,
    "Eutheria": 105,
    "Felidae": 60,
    "Feliformia": 57.5,
    "Glires": 75,
    "Gnathostomata": 500,
    "Gorilla": 8,
    "Haplorrhini": 42.5,
    "Hominidae": 17,
    "Homininae": 7.5,
    "Hominoidea": 17.5,
    "Hylobatidae": 17.5,
    "Laurasiatheria": 90,
    "Mammalia": 187.5,
    "Metatheria": 170,
    "Muroidea": 25,
    "Primates": 60,
    "Pteropus": 55,
    "Rodentia": 75,
    "Sarcopterygii": 400,
    "Sauropsida": 321,
    "Simiiformes": 35,
    "Teleostomi": 500,
    "Testudinoidea": 260,
    "Tetrapoda": 365,
    "Theria": 170,
    "Vertebrata": 500,
    "Whippomorpha": 57.5,
    "Xenopus": 260
}

# Read the taxon names from the file and count their frequency
taxon_counts = defaultdict(int)

with open("unique_filtered_ALL_recent_introns_taxon.txt", "r") as f:  # replace with your file name if different
    for line in f:
        taxon_name = line.strip() 
        if taxon_name in divergence_times:
            taxon_counts[taxon_name] += 1

# Extract taxons and their counts, sorted by their divergence times
sorted_taxons = sorted(taxon_counts.keys(), key=lambda x: divergence_times[x])
sorted_counts = [taxon_counts[taxon] for taxon in sorted_taxons]
sorted_times = [divergence_times[taxon] for taxon in sorted_taxons]

# Plotting
plt.figure(figsize=(12, 6))
plt.bar(sorted_times, sorted_counts, align='center', tick_label=sorted_taxons)
plt.xticks(rotation=90)
plt.xlabel("Divergence Time")
plt.ylabel("Taxon Frequency")
plt.title("Taxon Frequency vs Divergence Time")
plt.tight_layout()
plt.show()

