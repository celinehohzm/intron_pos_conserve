
#!/usr/bin/env python3

input_file = "blast_files/MANE_refseq_protein_chunk18_orthologs.txt"

# Read the input file and process data
with open(input_file, 'r') as f:
    data = {}
    for line in f:
        line = line.strip()
        key, value = line.split("\t")
        # Extract the accession number from the value
        value = value.split("|")[1]
        if key not in data:
            data[key] = []
        data[key].append(value)

# Write to separate output files
for key in data:
    output_file = f"test_blast_files/old_{key}_orthologs_proteinids.txt"
    with open(output_file, 'w') as f:
        f.write(key + "\n")
        f.write("\n".join(data[key]))

print("Done!")
