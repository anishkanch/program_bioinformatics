#!/bin/bash

#.fna means nucleic acid (DNA/RNA) fasta file whereas .faa means amino acid fasta file
#we want header to contain aditional information, specficially, the accession number to the beginning of the header

#usage of your script should be change_headers.sh <input file> <output file>

#cd "/Users/ark/Documents/GT/BINF/Fall/Program Bioinformatics/Project 2"

#obtain file names
input_file="$1"
output_file="$2"

#value of accession number
#(because the script and the file are in different folders, we are first removing the entire file path up until the accession number), used stackoverflow to learn "Remove a fixed prefix/suffix from a string in Bash"
temporary="${input_file##*_data/}"
#then we are removing everything after the accession number (.fna or .faa)
accession_number="${temporary%.*}"

#add accession number to beginning of each FASTA file, used stackoverflow for this to learn more about the sed function "Add words at beginning and end of a FASTA header line with sed"
sed "s/^>/>$accession_number/" "$input_file" > "$output_file"