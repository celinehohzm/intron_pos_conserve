from Bio.Align import substitution_matrices

# Get the original BLOSUM62 matrix
original_matrix = substitution_matrices.load("BLOSUM62")

# Create a dictionary to store the modified matrix
custom_matrix = {}

# Define a heavy penalty for gap openings
gap_penalty = -10000

# Modify the matrix to penalize gap openings
for pair, score in original_matrix.items():
    if '-' in pair:
        custom_matrix[pair] = gap_penalty
    else:
        custom_matrix[pair] = score

# Modify the score for gap-to-gap alignment
custom_matrix[('-', '-')] = 1

# Write the modified matrix to a file
amino_acids = sorted(list("ARNDCQEGHILKMFPSTWYVBZX-"))

with open("custom_matrix.txt", "w") as f:
    # Write the header row
    f.write("   " + "  ".join(amino_acids) + "\n")

    # Write the matrix rows
    for aa1 in amino_acids:
        row = [aa1]
        for aa2 in amino_acids:
            pair = (aa1, aa2)
            if pair not in custom_matrix:
                pair = (aa2, aa1)
            if pair in custom_matrix:
                row.append(str(custom_matrix[pair]).rjust(4))
            else:
                row.append("   .")
        f.write("  ".join(row) + "\n")

