import sys

def read_alignment(file_path):
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

    return headers, sequences

def format_alignment(headers, sequences):
    formatted_alignment = []
    for header, sequence in zip(headers, sequences):
        formatted_alignment.append(f"{header:<20} {sequence}")
    return formatted_alignment

def write_output(formatted_alignment, output_file):
    with open(output_file, 'w') as file:
        for line in formatted_alignment:
            file.write(line + '\n')

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python alignment_formatter.py input_file output_file")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    headers, sequences = read_alignment(input_file)
    formatted_alignment = format_alignment(headers, sequences)
    write_output(formatted_alignment, output_file)

