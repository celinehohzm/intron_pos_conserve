import sys

def filter_alignment_file(input_file, output_proteinid):
    seen_organisms = set()
    filtered_speciesnames = []
    filtered_proteinid = []

    with open(input_file, 'r') as infile:
        for line in infile:
            # Split the line into columns
            columns = line.strip().split('\t')
            
            # Extract the scientific name (column 15)
            organism = columns[14]
            
            # If the organism has not been seen before, add it to the set and keep the row
            if organism not in seen_organisms:
                seen_organisms.add(organism)
                # filtered_rows.append(line)
                filtered_speciesnames.append(columns[14] + '\n')
                filtered_proteinid.append(columns[1] + '\n')


    # Write the filtered rows to the output file
    with open(output_proteinid, 'w') as outfile:
        outfile.writelines(filtered_proteinid)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python filter_alignment.py input_alignment_file.txt output_filtered_alignment_file.txt")
        sys.exit(1)

    input_file = sys.argv[1]
    output_proteinid = sys.argv[2]
    #output_speciesnames = sys.argv[3]
    filter_alignment_file(input_file, output_proteinid)

