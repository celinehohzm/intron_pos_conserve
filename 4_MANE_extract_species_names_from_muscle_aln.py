import sys
import re

def extract_species_names(input_file, output_file):
    with open(input_file, "r") as in_file, open(output_file, "w") as out_file:
        for line in in_file:
            if line.startswith('>'):
                species = re.search('\[(.*?)\]', line)
                if species:
                    out_file.write(species.group(1) + '\n')

if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    extract_species_names(input_file, output_file)

