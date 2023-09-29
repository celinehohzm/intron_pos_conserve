import sys
import re
from collections import defaultdict

def generate_unique_id(species_name, species_count):
    words = species_name.split()
    if len(words) >= 2:
        unique_id = words[0][:min(3, len(words[0]))] + words[1][:min(6, len(words[1]))].capitalize()
    else:
        unique_id = species_name

    species_count[unique_id] += 1
    if species_count[unique_id] > 1:
        unique_id += f"{species_count[unique_id]}"

    return unique_id

def process_input_file(input_file, protein_id):
    with open(input_file, 'r') as file:
        content = file.readlines()

    intron_data = {}
    intron_name = None

    for line in content:
        if line.startswith("intron"):
            intron_name = line.split()[0]
            intron_data[intron_name] = {"x": [], "no_x": []}
        elif "Species with X" in line:
            intron_data[intron_name]["x"] = line.split(': ')[1].split(', ')
        elif "Species without X" in line:
            intron_data[intron_name]["no_x"] = line.split(': ')[1].split(', ')

    for intron_name, data in intron_data.items():
        output_file = f"{protein_id}_color_annot_{intron_name}.txt"
        species_count = defaultdict(int)

        with open(output_file, 'w') as output:
            output.write('DATASET_COLORSTRIP\n')
            output.write('SEPARATOR TAB\n')
            output.write('DATASET_LABEL\tMy_Colors\n')
            output.write('COLOR\t#ff0000\n')
            output.write('DATA\n')

            for species in data["x"]:
                # unique_id = generate_unique_id(species.strip(), species_count)
                output.write(f'{species.strip()}\t#FF0000\n')

            for species in data["no_x"]:
                # unique_id = generate_unique_id(species.strip(), species_count)
                output.write(f'{species.strip()}\t#0000FF\n')

if __name__ == '__main__':
    input_file = sys.argv[1]
    protein_id = sys.argv[2]
    process_input_file(input_file, protein_id)

