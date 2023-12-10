import sys
import re
from Bio import Entrez

def read_alignment(file_path):
    print(f"Reading alignment file: {file_path}")
    with open(file_path, 'r') as file:
        lines = file.readlines()

    headers = []
    sequences = []
    current_seq = ''

    for line in lines:
        if line.startswith('>'):
            headers.append(line.strip())
            if current_seq:
                sequences.append(current_seq)
                current_seq = ''
        else:
            current_seq += line.strip()

    sequences.append(current_seq)  # Add the last sequence

    print(f"Alignment file read successfully: {len(headers)} sequences found")
    return headers, sequences

def format_alignment(headers, sequences):
    print(f"Formatting alignment...")
    formatted_alignment = []
    for header, sequence in zip(headers, sequences):
        species_name = re.findall(r"\[(.*?)\]", header)[0]
        tax_id = get_taxonomy_id(species_name)
        formatted_header = f"{header} (TaxID: {tax_id})"
        formatted_alignment.append(f"{formatted_header:<20} {sequence}")

    print(f"Alignment formatted successfully")
    return formatted_alignment

def write_output(formatted_alignment, output_file):
    print(f"Writing output to file: {output_file}")
    with open(output_file, 'w') as file:
        for line in formatted_alignment:
            file.write(line + '\n')
    print(f"Output written successfully")

def get_taxonomy_id(species_name):
    Entrez.email = "celinehohzm@gmail.com" # Change this to your email address
    print(f"Querying NCBI Taxonomy database for: {species_name}")
    handle = Entrez.esearch(db="taxonomy", term=species_name)
    record = Entrez.read(handle)
    if record["IdList"]:
        print(f"Taxonomy ID found: {record['IdList'][0]}")
        return record["IdList"][0]
    else:
        print(f"No taxonomy ID found")
        return "not found"
    time.sleep(3)

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python alignment_formatter.py input_file output_file")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    headers, sequences = read_alignment(input_file)
    formatted_alignment = format_alignment(headers, sequences)
    write_output(formatted_alignment, output_file)

