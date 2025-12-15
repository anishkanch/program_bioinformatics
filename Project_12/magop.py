#!/usr/bin/env python3

import sys
import argparse
from magnumopus.ispcr import create_dict
from magnumopus.mapping import map_reads_to_ref
from magnumopus.nw import needleman_wunsch
import numpy as np
from biotite.sequence.phylo import neighbor_joining

def parse_args():
    parse = argparse.ArgumentParser()
    parse.add_argument("-a", "--assembly", nargs="*")
    parse.add_argument("-p", "--primer", required=True)
    parse.add_argument("-r", "--reads", nargs="*")
    parse.add_argument("-s", "--ref", nargs="*")
    parse.add_argument("-m", "--mas", type=int, default=500)
    return parse.parse_args()

args = parse_args()

if args.reads and not args.ref:
    print("--reads requires --ref", file=sys.stderr)
    sys.exit(1)

primer_file = args.primer
max_amplicon_size = args.mas

# assembly to ispcr
amplicons = {}

if args.assembly:
    for assembly in args.assembly:
        amplicons.update(create_dict(primer_file, assembly, max_amplicon_size))

if args.ref:
    for ref_file in args.ref:
        amplicons.update(create_dict(primer_file, ref_file, max_amplicon_size))

# reads to ispcr
def clean_reads(read_list):
    pairs = {}
    for r in read_list:
        sample = r.replace("_1.fastq", "").replace("_2.fastq", "")
        if sample not in pairs:
            pairs[sample] = []
        pairs[sample].append(r)
    return pairs

if args.reads:
    pairs = clean_reads(args.reads)

    for s, files in pairs.items():
        if len(files) != 2:
            continue
        
        r1 = [f for f in files if "_1.fastq" in f][0]
        r2 = [f for f in files if "_2.fastq" in f][0]

        sam_object = map_reads_to_ref(args.ref[0], r1, r2)

        consensus = sam_object.best_consensus(fasta=False)
        
        tmp = f"{s}_consensus_tmp.fasta"
        # FIX: write to the filename in tmp, not a literal "tmp"
        with open(tmp, "w") as file:
            file.write(f">{s}\n{consensus.upper()}\n")

        amplicons.update(create_dict(primer_file, tmp, max_amplicon_size))


if not amplicons:
    print("ERROR: No amplicons found. Exiting.", file=sys.stderr)
    sys.exit(1)

# orientation
def reverse_complement(seq):
    comp = ''
    for base in seq:
        base = base.upper()
        if base == "A":
            comp += "T"
        elif base == "T":
            comp += "A"
        elif base == "C":
            comp += "G"
        elif base == "G":
            comp += "C"
    return comp[::-1]


spot = next(iter(amplicons.values()))
orientation = {}

for name, seq in amplicons.items():
    # Needleman–Wunsch forward
    (_, _), score_forward = needleman_wunsch(spot, seq, 1, -1, -1)
    reverse_seq = reverse_complement(seq)
    # Needleman–Wunsch reverse
    (_, _), score_reverse = needleman_wunsch(spot, reverse_seq, 1, -1, -1)
    
    if score_reverse > score_forward:
        orientation[name] = reverse_seq
    else:
        orientation[name] = seq

#clean list 
def clean_label(name):
    name = name.split("/")[-1] #directory
    name = name.split(":")[0] #before first colon
    name = name.replace("_consensus_tmp", "") #temporary
    if name.endswith(".fastq"):
        name = name.replace(".fastq", "") #fasta extensions (only fna and fastq in data folder given)

    for c in [" ", ",", "(", ")", ";"]: #illegal characters
        name = name.replace(c, "_")
    return name

cleaned = {}
for old, seq in orientation.items():
    cleaned[clean_label(old)] = seq

orientation = cleaned

# distance
distances = {}

for name1 in orientation:
    for name2 in orientation:
        if name1 == name2:
            continue
        if (name1, name2) in distances or (name2, name1) in distances:
            continue

        seq1 = orientation[name1]
        seq2 = orientation[name2]

        (alignment1, alignment2), score = needleman_wunsch(seq1, seq2, 1, -1, -1)

        mismatches = 0
        for i in range(len(alignment1)):
            if alignment1[i] != alignment2[i]:
                mismatches += 1
        
        distance = mismatches / len(alignment1)

        distances[(name1, name2)] = distance

# neighbor joining Biotite Neighbor
names = []
for name in orientation:
    names.append(name)

n = len(names)

matrix = []
for i in range(n):
    matrix.append([0] * n)

for pair, dist in distances.items():
    a, b = pair

    i = names.index(a)
    j = names.index(b)

    matrix[i][j] = dist
    matrix[j][i] = dist

arr = np.array(matrix)

tree = neighbor_joining(arr)

print(tree.to_newick(labels=names))