import sys

def extract_table(input_data, protein_id):
    # Find the transcript variant with the specified protein ID
    target_line = f"protein {protein_id}"
    start_index = input_data.find(target_line)

    # Extract the corresponding table
    start_table = input_data.find("Genomic", start_index)
    end_table = input_data.find("\n\n", start_table)
    table = input_data[start_table:end_table]

    return table

if __name__ == "__main__":
    input_file = sys.argv[1]
    protein_id = sys.argv[2]

    with open(input_file, 'r') as f:
        input_data = f.read()

    extracted_table = extract_table(input_data, protein_id)
    print(extracted_table)
