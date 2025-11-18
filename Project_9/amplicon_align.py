#!/usr/bin/env python3

import argparse
from magnumopus import ispcr, needleman_wunsch

#Create command line interface directly from Exercise9 PDF
pars = argparse.ArgumentParser(description="Perform in-silico PCR on two assemblies and align the amplicons")
pars.add_argument('-1', dest = 'ASSEMBLY1', required=True, help='Path to the first assembly file')
pars.add_argument('-2', dest = 'ASSEMBLY2', required=True, help='Path to the second assembly file')
pars.add_argument('-p', dest = 'PRIMERS', required=True, help='Path to the primer file')
pars.add_argument('-m', dest = 'MAX_AMPLICON_SIZE', type=int, required=True, help='maximum amplicon size for isPCR')
pars.add_argument('--match', dest = 'MATCH',type=int, required=True, help='match score to use in alignment')
pars.add_argument('--mismatch',dest = 'MISMATCH', type=int, required=True, help='mismatch penalty to use in alignment')
pars.add_argument('--gap', dest = 'GAP', type=int, required=True, help='gap penalty to use in alignment')
args = pars.parse_args()

#isPCR using ispcr function we built, dest will be the names we refer to to pull from
#isPCR needs primer_file: str, assembly_file: str, max_amplicon_size: int
amp1 = ispcr(args.PRIMERS, args.ASSEMBLY1, args.MAX_AMPLICON_SIZE)
amp2 = ispcr(args.PRIMERS, args.ASSEMBLY2, args.MAX_AMPLICON_SIZE)
#split, remove the headerline, merge the two sequences onto one line
seq1 = ''.join(amp1.strip().split('>')[-1].splitlines()[1:])
seq2 = ''.join(amp2.strip().split('>')[-1].splitlines()[1:])

#align sequences using needleman_wunsch
#needleman_wunsch needs seq_a: str, seq_b: str, match: int, mismatch: int, gap: int
#align it first with both forward
(aligna1, alignb1), score1 = needleman_wunsch(seq1, seq2, args.MATCH, args.MISMATCH, args.GAP)
#then align it with seq2 in reverse
compliment = ''
for b in seq2:
    if b == "G": 
        compliment = compliment + "C" 
    elif b == "C": 
        compliment = compliment + "G" 
    elif b == "A": 
        compliment = compliment + "T" 
    elif b == "T": 
        compliment = compliment + "A"   
seq2_rev = compliment[::-1]  
(aligna2, alignb2), score2 = needleman_wunsch(seq1, seq2_rev, args.MATCH, args.MISMATCH, args.GAP)

#check for the highest score between score 1 and score 2. that is the correct alignment
if score2 > score1:
    aligna1 = aligna2
    alignb1 = alignb2
    score1 = score2

print(aligna1)
print(alignb1)
print(score1)