#!/usr/bin/env python3

#read a file that is in the command argument (DNA in fasta format)
#print the sequences to the terminal without the header
#add a pipe symbol between identical bases and spaces between differing

#usage will be pretty_align.py <FASTA file>

import sys
input_file = sys.argv[1]

seq1 = ""
seq2 = ""
match = ""
count = 0

file = open(input_file)
lines = file.readlines()

for line in lines:
    line = line.strip()
    if line.startswith(">"):
        count += 1
        continue
    if count == 1:
        seq1 += line
    elif count == 2:
        seq2 += line

for i in range(len(seq1)):
    if seq1[i] == seq2[i]:
        match += "|"
    else:
        match += " "

print(seq1)
print(match)
print(seq2)


