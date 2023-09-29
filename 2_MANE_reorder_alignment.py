import sys

input_file = sys.argv[1]
output_file = sys.argv[2]
protein_id = sys.argv[3]

with open(input_file, "r") as infile:
    lines = infile.readlines()

with open(output_file, "w") as outfile:
    protein_id_found = False
    i = 0
    while i < len(lines):
        if lines[i].startswith(">") and protein_id in lines[i]:
            outfile.write(lines[i])
            i += 1
            while i < len(lines) and not lines[i].startswith(">"):
                outfile.write(lines[i])
                i += 1
            protein_id_found = True
            break
        i += 1

    if not protein_id_found:
        print(f"Protein ID {protein_id} not found in the input file.")
        exit(1)

    for line in lines:
        if not (line.startswith(">") and protein_id in line):
            outfile.write(line)
