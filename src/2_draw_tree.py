import sys
import re
import matplotlib.pyplot as plt
from collections import defaultdict

if len(sys.argv) < 2:
    print("Usage: python script.py introns_data.txt")
    sys.exit()

input_file = sys.argv[1]

def parse_names_dmp(names_dmp_file):
    names = {}
    with open(names_dmp_file, 'r') as f:
        for line in f:
            columns = [col.strip() for col in line.split('|')]
            tax_id, name, _, name_class, _ = columns
            if name_class == "scientific name":
                names[int(tax_id)] = name
    return names

def parse_nodes_dmp(nodes_dmp_file):
    parents = {}
    with open(nodes_dmp_file, 'r') as f:
        for line in f:
            columns = [col.strip() for col in line.strip().split('|')]
            tax_id, parent_tax_id, *_ = columns
            parents[tax_id] = parent_tax_id
    return parents

def find_lineage(tax_id, parents):
    lineage = [tax_id]
    while tax_id != parents.get(tax_id, None):
        try:
            tax_id = parents[tax_id]
        except KeyError:
            print(f"Error processing tax_id: {tax_id}, tax_id not found in parents dictionary")
            return []
        lineage.append(tax_id)
    return lineage


def parse_input_file(filename):
    introns_taxon_ids = {}
    intron_name = None
    aligned_ids = []
    not_aligned_ids = []

    with open(filename, 'r') as file:
        for line in file:
            if line.startswith("intron_"):
                if intron_name is not None:
                    introns_taxon_ids[intron_name] = (aligned_ids, not_aligned_ids)
                intron_name = line.split(" ")[0]
                aligned_ids = []
                not_aligned_ids = []
            elif "Taxonomy IDs with X:" in line:
                aligned_ids = [int(x.strip()) for x in line.replace("Taxonomy IDs with X:", "").split(",")]
            elif "Taxonomy IDs without X:" in line:
                not_aligned_ids = [int(x.strip()) for x in line.replace("Taxonomy IDs without X:", "").split(",")]

        if intron_name is not None:
            introns_taxon_ids[intron_name] = (aligned_ids, not_aligned_ids)

    return introns_taxon_ids

introns_taxon_ids = parse_input_file(input_file)

names_dmp_file = 'names.dmp'
nodes_dmp_file = 'nodes.dmp'

names = parse_names_dmp(names_dmp_file)
parents = parse_nodes_dmp(nodes_dmp_file)

def draw_tree(lineages, names, output_file):
    tree = defaultdict(list)
    for lineage in lineages:
        for parent, child in zip(lineage[:-1], lineage[1:]):
            tree[parent].append(child)

    def visit(node, depth=0):
        children = tree[node]
        yield node, depth
        for child in children:
            yield from visit(child, depth + 1)

    visited = list(visit('1'))
    nodes, depths = zip(*visited)

    fig, ax = plt.subplots(figsize=(10, len(nodes) / 2))
    for node, depth in visited:
        ax.plot([depth, depth + 1], [nodes.index(node)] * 2, 'k')
        for child in tree[node]:
            ax.plot([depth + 1] * 2, [nodes.index(node), nodes.index(child)], 'k')
        ax.text(depth, nodes.index(node), names[int(node)], ha='right' if depth == 0 else 'center', va='center')

    ax.axis('off')
    plt.savefig(output_file)
    plt.close(fig)

for intron_name, (x_aligned_ids, x_not_aligned_ids) in introns_taxon_ids.items():
    lineages_aligned = []
    for tax_id in x_aligned_ids:
        try:
            lineage = find_lineage(tax_id, parents)
            lineages_aligned.append(lineage)
        except KeyError as e:
            print(f"Error processing tax_id: {tax_id}, error: {e}")

    lineages_not_aligned = [find_lineage(tax_id, parents) for tax_id in x_not_aligned_ids]

    draw_tree(lineages_aligned, names, f"{intron_name}_aligned_taxonomy_tree.png")
    draw_tree(lineages_not_aligned, names, f"{intron_name}_not_aligned_taxonomy_tree.png")


