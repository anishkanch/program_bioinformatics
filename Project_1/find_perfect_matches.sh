#!/bin/bash

#cd "/Users/ark/Documents/GT/BINF/Fall/Program Bioinformatics/Project 2"

#basic usage of blastn: 
#blastn -query <query fasta file> -subject <subject fasta file> -task <task> - outfmt <desired format> [-out <outfile>]

#usage of your script should be: 
#./find_perfect_matches.sh <query file> <subject file> <output file>

#HINT
#Think about what information you need to identify a perfect sequence match. In order for two sequences to be considered identical, they must have 100% sequence identity and equal length. Which output fields do you need to add (if any) to the BLAST output to determine which matches meet those criteria? Don't hard code a specic value for the length. Get it from the BLAST output.

#information below came from https://www.metagenomics.wiki/tools/blast/blastn-output-format-6
#100% identity depends on pident (should be 100)
#equal length depends on length = qlen

#load in file from usage of script
query_file="$1"
subject_file="$2"
output_file="$3"

#blastn, $1 will be pident, $2 will be length, $3 will be qlen
#we want $1 (pident) to be 100 and for $2(length) to equal $3(qlen), did some research on how to use two conditions with awk via stackoverflow: Awk AND operator
blastn -query "$query_file" -subject "$subject_file" -task blastn-short -outfmt "6 pident length qlen" | awk '$1 == 100 && $2==$3' > "$output_file"

#print number of perfect matches (will be length of output_file (wc -l for line count, not word count))
wc -l < "$output_file"
